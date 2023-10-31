using System.Collections.Generic;
using CenterSpace.NMath.Core;
using System.Threading.Tasks;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Abstract model of piezoelectric
  /// </summary>
  public abstract class AbstractModelPiezoelectric : AbstractModelStructure
  {
    protected List<int> listActuatorPatch;
    protected List<int> listSensorPatch;
    /// <summary>
    /// Constructor class
    /// </summary>
    public AbstractModelPiezoelectric(TypeAnalysisModel type, Dimension dimension, string pathProject = "temp", string nameProject = "project")
        : base(type, dimension, pathProject, nameProject)
    {
      if (StructureDimension == Dimension.Plane)
      {
        listComputeResult.Add(Result.VOLT);
      }
      else if (StructureDimension == Dimension.Solid)
      {
        listComputeResult.Add(Result.VOLT);
      }
      TypeModel = TypeModelProblem.Piezoelectric;
      typeOfFieldsInMultifield.Add(TypeFields.Electric);
    }
    public void AddActuatorPatch(params int[] idPatch)
    {
      if (listActuatorPatch == null)
        listActuatorPatch = new List<int>();
      foreach (int i in idPatch)
        listActuatorPatch.Add(i);
    }
    public void AddSensorPatch(params int[] idPatch)
    {
      if (listSensorPatch == null)
        listSensorPatch = new List<int>();
      foreach (int i in idPatch)
        listSensorPatch.Add(i);
    }
    public void PreProcessing()
    {
      base.PreProcessing();
      tArray = new int[typeOfFieldsInMultifield.Count][];
      tArrayConstraint = new int[typeOfFieldsInMultifield.Count][];
      tArrayUnconstraint = new int[typeOfFieldsInMultifield.Count][];
      if (!IsParallelProcesing)
      {
        tArray[0] = GetTArrayGlobalMechanic();
        tArray[1] = GetTArrayGlobalElectric();
        if (bcMatrix != null)
        {
          if (listActuatorPatch == null && listSensorPatch == null)
          {
            tArrayUnconstraint[0] = GetTArrayGlobalMechanicUnConstraint();
            tArrayConstraint[0] = GetTArrayGlobalMechanicConstraint();
            tArrayUnconstraint[1] = GetTArrayGlobalElectricUnConstraint();
            tArrayConstraint[1] = GetTArrayGlobalElectricConstraint();
          }
        }
        else
        {
          tArrayUnconstraint[0] = tArray[0];
          tArrayUnconstraint[1] = tArray[1];
        }
      }
      else
      {
        if (bcMatrix != null)
        {
          if (listActuatorPatch == null && listSensorPatch == null)
          {
            Task<int[]>[] tasks ={
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalMechanic()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalElectric()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalMechanicUnConstraint()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalMechanicConstraint()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalElectricUnConstraint()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalElectricConstraint())};
            //Block until all tasks complete.
            Task.WaitAll(tasks);
            tArray[0] = tasks[0].Result;
            tArray[1] = tasks[1].Result;
            tArrayUnconstraint[0] = tasks[2].Result;
            tArrayConstraint[0] = tasks[3].Result;
            tArrayUnconstraint[1] = tasks[4].Result;
            tArrayConstraint[1] = tasks[5].Result;
          }
        }
        else
        {
          Task<int[]>[] tasks ={
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalMechanic()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalElectric()) };
          //Block until all tasks complete.
          Task.WaitAll(tasks);
          tArray[0] = tasks[0].Result;
          tArray[1] = tasks[1].Result;
          tArrayUnconstraint[0] = tArray[0];
          tArrayUnconstraint[1] = tArray[1];
        }
      }
    }
    public Structure2DState StressState
    {
      get; set;
    }
    protected DoubleVector GetConstraint(int[] tArrayConstraint)
    {
      if (tArrayConstraint.Length != 0)
      {
        DoubleVector u0 = new DoubleVector(tArrayConstraint.Length);
        for (int i = 0; i < u0.Length; i++)
        {
          int idx = MatrixTool.IndexOfCol(bcMatrix, 0, tArrayConstraint[i], AbstractModel.NumberOfCPUs);
          u0[i] = bcMatrix[idx, 1];
        }
        return u0;
      }
      else
        return null;
    }

    private int[] GetTArrayGlobalMechanic()
    {
      List<int> listTArrayGlobalMechanic = new List<int>();
      for (int i = 0; i < listPatch.Count; i++)
      {
        AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
        int nen = pa.GetCountGlobalBasisFunctions();
        for (int k = 0; k < nen; k++)
        {
          for (int j = 0; j < (int)StructureDimension; j++)
          {
            int temp = pa.GetIDInGlobal(j, k);
            if (temp != -1)
            {
              if (listTArrayGlobalMechanic.IndexOf(temp) == -1)
              {
                listTArrayGlobalMechanic.Add(temp);
              }
            }
          }
        }
      }
      return listTArrayGlobalMechanic.ToArray();
    }

    private int[] GetTArrayGlobalMechanicConstraint()
    {
      List<int> listTArrayGlobalMechanic = new List<int>();
      for (int i = 0; i < listPatch.Count; i++)
      {
        AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
        int nen = pa.GetCountGlobalBasisFunctions();
        for (int k = 0; k < nen; k++)
        {
          for (int j = 0; j < (int)StructureDimension; j++)
          {
            int temp = pa.GetIDInGlobal(j, k);
            if (temp != -1)
            {
              int idxConstraint = MatrixTool.IndexOfCol(bcMatrix, 0, temp, AbstractModel.NumberOfCPUs);
              if (idxConstraint != -1)
                if ((listTArrayGlobalMechanic.IndexOf(temp) == -1) && (bcMatrix[idxConstraint, 2] != (int)StructureDimension/* khong phai la bac tu do volt*/))
                {
                  listTArrayGlobalMechanic.Add(temp);
                }
            }
          }
        }
      }
      return listTArrayGlobalMechanic.ToArray();
    }

    private int[] GetTArrayGlobalMechanicUnConstraint()
    {
      List<int> listTArrayGlobalMechanic = new List<int>();
      for (int i = 0; i < listPatch.Count; i++)
      {
        AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
        int nen = pa.GetCountGlobalBasisFunctions();
        for (int k = 0; k < nen; k++)
        {
          for (int j = 0; j < (int)StructureDimension; j++)
          {
            int temp = pa.GetIDInGlobal(j, k);
            if (temp != -1)
            {
              int idxConstraint = MatrixTool.IndexOfCol(bcMatrix, 0, temp, AbstractModel.NumberOfCPUs);
              if ((listTArrayGlobalMechanic.IndexOf(temp) == -1) && (idxConstraint == -1))
              {
                listTArrayGlobalMechanic.Add(temp);
              }
            }
          }
        }
      }
      return listTArrayGlobalMechanic.ToArray();
    }

    private int[] GetTArrayGlobalElectric()
    {
      List<int> listTArrayGlobalElectric = new List<int>();
      for (int i = 0; i < listPatch.Count; i++)
      {
        AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
        if (TypeModel == TypeModelProblem.Piezoelectric)
        {
          int nen = pa.GetCountGlobalBasisFunctions();
          for (int k = 0; k < nen; k++)
          {
            int temp = pa.GetIDInGlobal((int)StructureDimension, k);
            if (temp != -1)
            {
              if (listTArrayGlobalElectric.IndexOf(temp) == -1)
              {
                listTArrayGlobalElectric.Add(temp);
              }
            }
          }
        }
      }
      return listTArrayGlobalElectric.ToArray();
    }
    private int[] GetTArrayGlobalElectricConstraint()
    {
      List<int> listTArrayGlobalElectric = new List<int>();
      for (int i = 0; i < listPatch.Count; i++)
      {
        AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
        if (TypeModel == TypeModelProblem.Piezoelectric)
        {
          int nen = pa.GetCountGlobalBasisFunctions();
          for (int k = 0; k < nen; k++)
          {
            int temp = pa.GetIDInGlobal((int)StructureDimension, k);
            if (temp != -1)
            {
              int idxConstraint = MatrixTool.IndexOfCol(bcMatrix, 0, temp, AbstractModel.NumberOfCPUs);
              if (idxConstraint != -1)
                if ((listTArrayGlobalElectric.IndexOf(temp) == -1) && (bcMatrix[idxConstraint, 2] == (int)StructureDimension/*bac tu do volt*/))
                {
                  listTArrayGlobalElectric.Add(temp);
                }
            }
          }
        }
      }
      return listTArrayGlobalElectric.ToArray();
    }

    private int[] GetTArrayGlobalElectricUnConstraint()
    {
      List<int> listTArrayGlobalElectric = new List<int>();
      for (int i = 0; i < listPatch.Count; i++)
      {
        AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
        if (TypeModel == TypeModelProblem.Piezoelectric)
        {
          int nen = pa.GetCountGlobalBasisFunctions();
          for (int k = 0; k < nen; k++)
          {
            int temp = pa.GetIDInGlobal((int)StructureDimension, k);
            if (temp != -1)
            {
              int idxConstraint = MatrixTool.IndexOfCol(bcMatrix, 0, temp, AbstractModel.NumberOfCPUs);
              if ((listTArrayGlobalElectric.IndexOf(temp) == -1) && (idxConstraint == -1))
              {
                listTArrayGlobalElectric.Add(temp);
              }
            }
          }
        }
      }
      return listTArrayGlobalElectric.ToArray();
    }

    //protected void SetUGlobalUnconstraint(DoubleVector uGlobal)
    //{
    //  for (int i = 0; i < listPatch.Count; i++)
    //  {
    //    AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
    //    ControlPoint[] cps = pa.GetAllControlPoints();
    //    //if (!IsParallelProcesing)
    //    //{
    //    for (int j = 0; j < cps.Length; j++)
    //    {
    //      ControlPoint currentCps = cps[j];
    //      for (int k = 0; k < (int)StructureDimension; k++)
    //      {
    //        int tArray = currentCps.GetTArrayGlobal()[k];
    //        if (tArray != -1)
    //        {
    //          int tArrayUnconstraint = Array.IndexOf(tArrayMechanicUnconstraint, tArray);
    //          if (tArrayUnconstraint != -1)
    //          {
    //            currentCps.SetULocal(k, uGlobal[tArrayUnconstraint]);
    //            currentCps.AddResult(listComputeResult[k], uGlobal[tArrayUnconstraint]);
    //          }
    //        }
    //      }
    //    }
    //    //}
    //    //else
    //    //{
    //    //    Parallel.For(0, cps.Length - 1, j =>
    //    //    {
    //    //        ControlPoint currentCps = cps[j];
    //    //        for (int k = 0; k < pa.GetCountField(); k++)
    //    //        {
    //    //            int tArray = currentCps.GetTArrayGlobal()[k];
    //    //            if (tArray != -1)
    //    //                currentCps.SetULocal(k, uGlobal[tArray]);
    //    //        }
    //    //    });
    //    //}
    //  }
    //}

    //protected void SetPGlobalUnconstraint(DoubleVector uGlobal)
    //{
    //  for (int i = 0; i < listPatch.Count; i++)
    //  {
    //    AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
    //    ControlPoint[] cps = pa.GetAllControlPoints();
    //    //if (!IsParallelProcesing)
    //    //{
    //    for (int j = 0; j < cps.Length; j++)
    //    {
    //      ControlPoint currentCps = cps[j];
    //      int tArray = currentCps.GetTArrayGlobal()[(int)StructureDimension];
    //      if (tArray != -1)
    //      {
    //        int tArrayUnconstraint = Array.IndexOf(tArrayElectricUnconstraint, tArray);
    //        if (tArrayUnconstraint != -1)
    //        {
    //          currentCps.SetULocal((int)StructureDimension, uGlobal[tArrayUnconstraint]);
    //          currentCps.AddResult(listComputeResult[(int)StructureDimension], uGlobal[tArrayUnconstraint]);
    //        }
    //      }
    //    }
    //    //}
    //    //else
    //    //{
    //    //    Parallel.For(0, cps.Length - 1, j =>
    //    //    {
    //    //        ControlPoint currentCps = cps[j];
    //    //        for (int k = 0; k < pa.GetCountField(); k++)
    //    //        {
    //    //            int tArray = currentCps.GetTArrayGlobal()[k];
    //    //            if (tArray != -1)
    //    //                currentCps.SetULocal(k, uGlobal[tArray]);
    //    //        }
    //    //    });
    //    //}
    //  }
    //}
  }
}
