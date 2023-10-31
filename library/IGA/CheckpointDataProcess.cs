using CenterSpace.NMath.Core;
using DEMSoft.Common;
using System;
using System.Collections.Generic;
using System.IO;
using System.Threading.Tasks;

namespace DEMSoft.IGA
{
  class CheckpointDataProcess
  {
    public double LastPrevFi0;
    public double LastPrevFi1;
    public int CurrentStepLoad;
    public int CurrentSubStep;
    public int CurrentTotalStep;
    public int CurrentTotalStepSave;
    public DoubleVector un0;
    public AbstractModel model;
    public bool IsFinished;
    //public List<DoubleVector> DisplacementTime;
    public List<double> Time;
    public List<double[]> MonitorDataStore;
    //public List<List<DataGaussPoint[,]>> patch;
    public double[][][,] currentPhase2d;
    public double[][][,] xiEpsilon2d;
    public double[][][,,] currentPhase3d;
    public double[][][,,] xiEpsilon3d;
    //public DoubleVector[][][,] currentStress;
    //public DoubleVector[][][,] currentStrain;
    //formatter.Serialize(stream, new DoubleVector[] { u, v});
    private string path;
    public CheckpointDataProcess(string path)
    {
      this.path = path;
    }
    public void SaveData()
    {
      string filenameCheckpointParameter = path + "restore_parameter.msgpack";
      //string filenameCheckpointTime = path + "restore_time.msgpack";
      string filenameCheckpointUn0 = path + "restore_un0.msgpack";
      string fileNameTime = path + "time.msgpack";
      string filenameCheckpointMonitorDataStore = path + "restore_monitor.msgpack";
      string filenameCheckpointGausspoint = path + "restore_gausspoint.msgpack";
      string directory = Path.GetDirectoryName(filenameCheckpointParameter);
      if (!Directory.Exists(directory))
        Directory.CreateDirectory(directory);
      if (File.Exists(filenameCheckpointParameter))
        File.Delete(filenameCheckpointParameter);
      //if (File.Exists(filenameCheckpointTime))
      //  File.Delete(filenameCheckpointTime);
      if (File.Exists(filenameCheckpointUn0))
        File.Delete(filenameCheckpointUn0);
      if (File.Exists(filenameCheckpointMonitorDataStore))
        File.Delete(filenameCheckpointMonitorDataStore);
      if (File.Exists(filenameCheckpointGausspoint))
        File.Delete(filenameCheckpointGausspoint);

      int numPatch = model.CountPatch();
      if (model.StructureDimension == Dimension.Plane || model.StructureDimension == Dimension.Plate || model.StructureDimension == Dimension.Shell)
      {
        currentPhase2d = new double[numPatch][][,];
        xiEpsilon2d = new double[numPatch][][,];
      }
      else
      {
        currentPhase3d = new double[numPatch][][,,];
        xiEpsilon3d = new double[numPatch][][,,];
      }
      for (int i = 0; i < numPatch; i++)
      {
        int numElement = model.GetPatch(i).CountElements();
        if (model.StructureDimension == Dimension.Plane || model.StructureDimension == Dimension.Plate || model.StructureDimension == Dimension.Shell)
        {
          currentPhase2d[i] = new double[numElement][,];
          xiEpsilon2d[i] = new double[numElement][,];
        }
        else
        {
          currentPhase3d[i] = new double[numElement][,,];
          xiEpsilon3d[i] = new double[numElement][,,];
        }
        if (!AbstractModel.IsParallelProcesing)
        {
          for (int j = 0; j < numElement; j++)
          {
            StoreGausspointValue(i, j);
          }
        }
        else
        {
          var degreeOfParallelism = AbstractModel.NumberOfCPUs;
          var tasks = new Task[degreeOfParallelism];
          //object monitor = new object();
          for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
          {
            // capturing taskNumber in lambda wouldn't work correctly
            int taskNumberCopy = taskNumber;

            tasks[taskNumber] = Task.Factory.StartNew(
                 () =>
                 {
                   var max = numElement * (taskNumberCopy + 1) / degreeOfParallelism;
                   for (int j = numElement * taskNumberCopy / degreeOfParallelism; j < max; j++)
                   {
                     StoreGausspointValue(i, j);
                   }
                 });
          }
          Task.WaitAll(tasks);
        }
      }

      // Presist to file
      IOIGA.WriteMessagePackData(filenameCheckpointParameter, CurrentStepLoad, CurrentSubStep, CurrentTotalStep, CurrentTotalStepSave, false, LastPrevFi0, LastPrevFi1);
      IOIGA.WriteMessagePackDoubleVector(filenameCheckpointUn0, un0);
      IOIGA.WriteMessagePackData(fileNameTime, model.Time);
      //List<double[]> disTime = new List<double[]>();
      //for (int i = 0; i < model.DisplacementTime.Count; i++)
      //{
      //  disTime.Add(model.DisplacementTime[i].ToArray());
      //}
      if (model.GetMonitorDataHistorism().MonitorDataStore != null)
      {
        IOIGA.WriteMessagePackData(filenameCheckpointMonitorDataStore, model.GetMonitorDataHistorism().MonitorDataStore);
      }
      if (model.StructureDimension == Dimension.Plane || model.StructureDimension == Dimension.Plate || model.StructureDimension == Dimension.Shell)
      {
        IOIGA.WriteMessagePackData(filenameCheckpointGausspoint, currentPhase2d, xiEpsilon2d);
      }
      else
      {
        IOIGA.WriteMessagePackData(filenameCheckpointGausspoint, currentPhase3d, xiEpsilon3d);
      }
    }

    public void RestoreLastCheckpointData(string path)
    {
      string filenameCheckpointParameter = path + "restore_parameter.msgpack";
      string filenameCheckpointUn0 = path + "restore_un0.msgpack";
      string fileNameTime = path + "time.msgpack";
      string filenameCheckpointMonitorDataStore = path + "restore_monitor.msgpack";
      string filenameCheckpointGausspoint = path + "restore_gausspoint.msgpack";

      object[] cpParameter = IOIGA.ReadMessagePackData<object>(filenameCheckpointParameter);
      CurrentStepLoad = Convert.ToInt16(cpParameter[0]);
      CurrentSubStep = Convert.ToInt16(cpParameter[1]);
      CurrentTotalStep = Convert.ToInt16(cpParameter[2]);
      CurrentTotalStepSave = Convert.ToInt16(cpParameter[3]);
      IsFinished = Convert.ToBoolean(cpParameter[4]);
      LastPrevFi0 = Convert.ToDouble(cpParameter[5]);
      LastPrevFi1 = Convert.ToDouble(cpParameter[6]);

      un0 = IOIGA.ReadMessagePackDoubleVector(filenameCheckpointUn0);
      Time = IOIGA.ReadMessagePackData<List<double>>(fileNameTime)[0];
      //DisplacementTime = new List<DoubleVector>();
      List<double[]>[] cpMonitor = IOIGA.ReadMessagePackData<List<double[]>>(filenameCheckpointMonitorDataStore);
      MonitorDataStore = cpMonitor[0];
      //List<double[]> disTimeTemp = cpMonitor[1];
      //for (int i = 0; i < disTimeTemp.Count; i++)
      //{
      //  DisplacementTime.Add(new DoubleVector(disTimeTemp[i]));
      //}
      if (model.StructureDimension == Dimension.Plane || model.StructureDimension == Dimension.Plate || model.StructureDimension == Dimension.Shell)
      {
        double[][][][,] cpGausspoint2d = IOIGA.ReadMessagePackData<double[][][,]>(filenameCheckpointGausspoint);
        currentPhase2d = cpGausspoint2d[0];
        xiEpsilon2d = cpGausspoint2d[1];
      }
      else
      {
        double[][][][,,] cpGausspoint3d = IOIGA.ReadMessagePackData<double[][][,,]>(filenameCheckpointGausspoint);
        currentPhase3d = cpGausspoint3d[0];
        xiEpsilon3d = cpGausspoint3d[1];
      }
      int numPatch = model.CountPatch();
      for (int i = 0; i < numPatch; i++)
      {
        int numElement = model.GetPatch(i).CountElements();
        if (!AbstractModel.IsParallelProcesing)
        {
          for (int j = 0; j < numElement; j++)
          {
            RestoreGausspointValue(i, j);
          }
        }
        else
        {
          var degreeOfParallelism = AbstractModel.NumberOfCPUs;
          var tasks = new Task[degreeOfParallelism];
          //object monitor = new object();
          for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
          {
            // capturing taskNumber in lambda wouldn't work correctly
            int taskNumberCopy = taskNumber;

            tasks[taskNumber] = Task.Factory.StartNew(
                 () =>
                             {
                               var max = numElement * (taskNumberCopy + 1) / degreeOfParallelism;
                               for (int j = numElement * taskNumberCopy / degreeOfParallelism; j < max; j++)
                               {
                                 RestoreGausspointValue(i, j);
                               }
                             });
          }
          Task.WaitAll(tasks);
        }
      }
    }
    private void StoreGausspointValue(int i, int j)
    {
      if (model.StructureDimension == Dimension.Plane || model.StructureDimension == Dimension.Plate || model.StructureDimension == Dimension.Shell)
      {
        AbstractElement2D elem = (AbstractElement2D)model.GetPatch(i).GetElement(j);
        int numGauss = elem.GetNumberOfGaussPointOnEachDirection();
        currentPhase2d[i][j] = new double[numGauss, numGauss];
        xiEpsilon2d[i][j] = new double[numGauss, numGauss];
        //currentStress[i][j] = new DoubleVector[numGauss, numGauss];
        //currentStrain[i][j] = new DoubleVector[numGauss, numGauss];
        for (int ii = 0; ii < elem.GetNumberOfGaussPointOnEachDirection(); ii++)
        {
          for (int jj = 0; jj < elem.GetNumberOfGaussPointOnEachDirection(); jj++)
          {
            GaussPoints gp = elem.GetGaussPoint(ii, jj);
            //currentStrain[i][j][ii, jj] = gp.lastStrain;
            //currentStress[i][j][ii, jj] = gp.lastStress;
            object v = gp.GetValue(DataInGausspoint.lastPhase);
            if (v != null)
              currentPhase2d[i][j][ii, jj] = (double)v;

            object v1 = gp.GetValue(DataInGausspoint.xiEpsilon);
            if (v1 != null)
              xiEpsilon2d[i][j][ii, jj] = (double)v1;
          }
        }
      }
      else
      {
        AbstractElement3D elem = (AbstractElement3D)model.GetPatch(i).GetElement(j);
        int numGauss = elem.GetNumberOfGaussPointOnEachDirection();
        currentPhase3d[i][j] = new double[numGauss, numGauss, numGauss];
        xiEpsilon3d[i][j] = new double[numGauss, numGauss, numGauss];
        for (int ii = 0; ii < elem.GetNumberOfGaussPointOnEachDirection(); ii++)
        {
          for (int jj = 0; jj < elem.GetNumberOfGaussPointOnEachDirection(); jj++)
          {
            for (int kk = 0; kk < elem.GetNumberOfGaussPointOnEachDirection(); kk++)
            {
              GaussPoints gp = elem.GetGaussPoint(ii, jj, kk);
              object v = gp.GetValue(DataInGausspoint.lastPhase);
              if (v != null)
                currentPhase3d[i][j][ii, jj, kk] = (double)v;

              object v1 = gp.GetValue(DataInGausspoint.xiEpsilon);
              if (v1 != null)
                xiEpsilon3d[i][j][ii, jj, kk] = (double)v1;
            }
          }
        }
      }
    }
    private void RestoreGausspointValue(int i, int j)
    {
      if (model.StructureDimension == Dimension.Plane || model.StructureDimension == Dimension.Plate || model.StructureDimension == Dimension.Shell)
      {
        AbstractElement2D elem = (AbstractElement2D)model.GetPatch(i).GetElement(j);
        int numGauss = elem.GetNumberOfGaussPointOnEachDirection();

        for (int ii = 0; ii < elem.GetNumberOfGaussPointOnEachDirection(); ii++)
        {
          for (int jj = 0; jj < elem.GetNumberOfGaussPointOnEachDirection(); jj++)
          {
            GaussPoints gp = elem.GetGaussPoint(ii, jj);
            //gp.currentStrain = currentStrain[i][j][ii, jj];
            //gp.currentStress = currentStress[i][j][ii, jj];
            gp.SetValue(DataInGausspoint.currentPhase, currentPhase2d[i][j][ii, jj]);
            gp.SetValue(DataInGausspoint.lastPhase, currentPhase2d[i][j][ii, jj]);
            gp.SetValue(DataInGausspoint.xiEpsilon, xiEpsilon2d[i][j][ii, jj]);
            //gp.lastPhase = gp.currentPhase = currentPhase[i][j][ii, jj];
            //gp.xiEpsilon = xiEpsilon[i][j][ii, jj];
          }
        }
      }
      else
      {
        AbstractElement3D elem = (AbstractElement3D)model.GetPatch(i).GetElement(j);
        int numGauss = elem.GetNumberOfGaussPointOnEachDirection();

        for (int ii = 0; ii < elem.GetNumberOfGaussPointOnEachDirection(); ii++)
        {
          for (int jj = 0; jj < elem.GetNumberOfGaussPointOnEachDirection(); jj++)
          {
            for (int kk = 0; kk < elem.GetNumberOfGaussPointOnEachDirection(); kk++)
            {
              GaussPoints gp = elem.GetGaussPoint(ii, jj, kk);
              gp.SetValue(DataInGausspoint.currentPhase, currentPhase3d[i][j][ii, jj, kk]);
              gp.SetValue(DataInGausspoint.lastPhase, currentPhase3d[i][j][ii, jj, kk]);
              gp.SetValue(DataInGausspoint.xiEpsilon, xiEpsilon3d[i][j][ii, jj, kk]);
            }
          }
        }
      }
    }
  }
}
