//using System;
//using System.Linq;
//using DEMSoft.Drawing;
//using CenterSpace.NMath.Core;
//using System.Collections.Generic;
//using DEMSoft.NURBS;
//using DEMSoft.Function;
//using System.IO;
//using System.Threading.Tasks;
//using CenterSpace.NMath.Matrix;
//using System.Diagnostics;

//namespace DEMSoft.IGA
//{
  
//  /// <summary>
//  /// Implement the load controlled Newton-Raphson solver for non-linear problems
//  /// </summary>
//  class SolverLoadControlled
//  {
//    private bool bConverged = true;
//    private int type; // 0 for Newton-Raphson, 1 for Modified Newton-Raphson
//    private AbstractModelStructure model;
//    private CheckpointDataProcess checkpointData;
//    internal static bool IsSolveSecondStaggered;
//    private int[][] tArrayConstraintNoneZero;
//    private bool conditionRestoreCheckpoint;
//    internal static double deltaT;
//    internal DoubleCsrSparseMatrix TMatrixSparse;
//    internal double prevF, lastPrevF0, lastPrevF1;
//    /// <summary>
//    /// Constructor class
//    /// </summary>
//    /// <param name="model"></param>
//    /// <param name="type"></param>
//    public SolverLoadControlled(AbstractModelStructure model, int type)
//    {
//      this.type = type;
//      this.model = model;
//      checkpointData = new CheckpointDataProcess();
//      checkpointData.model = model;
//      conditionRestoreCheckpoint = !model.IsRunFromInitial && File.Exists(model.filenameCheckpoint[0]);
//      if (conditionRestoreCheckpoint)
//      {
//        checkpointData.RestoreLastCheckpointData(model.filenameCheckpoint);
//      }
//      Solve();
//    }

//    /// <summary>
//    /// Solve
//    /// </summary>
//    private void Solve()
//    {
//      int totalStep = ((model.IsRunFromInitial) || (checkpointData.IsFinished) || !conditionRestoreCheckpoint) ? 1 : checkpointData.CurrentTotalStep;//1;
//      int step = ((model.IsRunFromInitial) || (checkpointData.IsFinished) || !conditionRestoreCheckpoint) ? 1 : checkpointData.CurrentStepLoad;//1;
//      int numberOfTotalStep = 0;
//      int numberOfStepLoad = model.numberOfSubstepLoad.Length;

//      model.deltaT = new double[numberOfStepLoad];
//      for (int i = 0; i < numberOfStepLoad; i++)
//      {
//        numberOfTotalStep += model.numberOfSubstepLoad[i];
//        model.deltaT[i] = 1.0/*[s]*/ / model.numberOfSubstepLoad[i];
//      }

//      int countData = model.GetMonitorDataHistorism().CountData();
//      model.DisplacementTime = new DoubleMatrix(model.CountDOF(), numberOfTotalStep);
//      model.Time = new DoubleVector(numberOfTotalStep);
//      if (countData > 0)
//        model.MonitorDataStore = new DoubleMatrix(numberOfTotalStep, countData);
//      DoubleVector un0 = new DoubleVector(model.CountDOF());
//      DoubleVector un1 = null;
//      if ((!model.IsRunFromInitial) && (!checkpointData.IsFinished) && conditionRestoreCheckpoint)
//      {
//        totalStep = checkpointData.CurrentTotalStep;
//        step = checkpointData.CurrentStepLoad;
//        model.DisplacementTime = checkpointData.DisplacementTime;
//        model.Time = checkpointData.Time;
//        if (countData > 0)
//          model.MonitorDataStore = checkpointData.MonitorDataStore;
//        un0 = checkpointData.un0;
//        lastPrevF0 = checkpointData.LastPrevFi0;
//        lastPrevF1 = checkpointData.LastPrevFi1;
//      }

//      string[] filenameCheckpoint = model.filenameCheckpoint;
//      string fileNameDataTime = model.FileNameDataTime;
//      string fileNameTime = model.FileNameTime;
//      while (step <= numberOfStepLoad)
//      {
//        if (model.IsWriteLogFile)
//        {
//          IOIGA.WriteLogFileAndConsole(model.logFile, true, "-----------------------------------------------------------");
//          IOIGA.WriteLogFileAndConsole(model.logFile, true, "Load step " + step + " : ");
//        }
//        int numberOfSubstepLoad = model.numberOfSubstepLoad[step - 1];
//        //Compute increment vector of field
//        DoubleVector[] DNoneZeroVector = new DoubleVector[model.typeOfFieldsInMultifield.Count];
//        tArrayConstraintNoneZero = new int[model.typeOfFieldsInMultifield.Count][];

//        for (int i = 0; i < model.typeOfFieldsInMultifield.Count; i++)
//        {
//          ComputeDisplacementVector(step - 1, model.tArrayConstraint[i], out DNoneZeroVector[i], out tArrayConstraintNoneZero[i]);

//          /////////////////////////////////////////////////////////////
//          ///////Store only increment none zero///////
//          /////////////////////////////////////////////////////////////
//          if (DNoneZeroVector[i] != null)
//          {
//            if (DNoneZeroVector[i].Length > 0 && i == 0/*Structure*/)
//              DNoneZeroVector[i] = DNoneZeroVector[i] / numberOfSubstepLoad;
//          }
//        }

//        int substep = ((model.IsRunFromInitial) || (checkpointData.IsFinished) || !conditionRestoreCheckpoint) ? 1 : checkpointData.CurrentSubStep;//1;

//        while (substep <= numberOfSubstepLoad)
//        {
//          if (model.IsWriteLogFile)
//          {
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "-----------------------------------------------------------");
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "SubLoad step " + substep + " : ");
//          }

//          // Apply inhomogenous boundary condition
//          for (int i = 0; i < model.typeOfFieldsInMultifield.Count; i++)
//          {
//            if (tArrayConstraintNoneZero[i] != null)
//              for (int ii = 0; ii < tArrayConstraintNoneZero[i].Length; ii++)
//              {
//                //Chua dung cho dkb phase=1
//                un0[tArrayConstraintNoneZero[i][ii]] += DNoneZeroVector[i][ii];
//              }
//          }

//          model.SetUGlobal(un0);

//          //////////////////////////////////////////////////////////////////
//          //double loadFactor = (double)step / (double)numLoadStep;
//          //DoubleVector fStep = loadFactor * f;
//          //////////////////////////////////////////////////////////////////
//          double time = step - 1 + substep * model.deltaT[step - 1];
//          deltaT = model.deltaT[step - 1];
//          model.Time[totalStep - 1] = time;
//          DoubleVector fStep = null;
//          model.AssemblyTractionVector(out fStep, time);
//          DoubleVector FR = null;
//          // Solve
//          if (AbstractModel.IsMonolithic)
//          {
//            un1 = SolveOneSubstepLoad(un0, fStep, out FR, time);
//            lastPrevF0 = prevF;
//          }
//          else
//          {
//            // Solve phase field
//            IsSolveSecondStaggered = true;
//            un1 = SolveOneSubstepLoad(un0, null, out DoubleVector FR1, time);
//            lastPrevF0 = prevF;
//            // Solve displacement
//            IsSolveSecondStaggered = false;
//            un1 = SolveOneSubstepLoad(un1, fStep, out FR, time);
//            lastPrevF1 = prevF;
//          }

//          //////////////////////////////
//          ////// update gausspoint variables ////
//          ///////////////////////////////
//          if (model.GetPatch(0).Material.TypeMaterialStructure == EngineeringData.TypeMaterialStructure.Plasticity || AbstractModel.TypeModel == TypeModelProblem.PhaseField)
//            switch (model.StructureDimension)
//            {
//              case Dimension.Plane:
//                for (int i = 0; i < model.CountPatch(); i++)
//                {
//                  AbstractPatch patch = model.GetPatch(i);
//                  int numElement = patch.CountElements();
//                  if (!AbstractModel.IsParallelProcesing)
//                  {
//                    for (int j = 0; j < numElement; j++)
//                    {
//                      AbstractElement2D elem = (AbstractElement2D)patch.GetElement(j);
//                      int numberGps = elem.GetNumberOfGaussPointOnEachDirection();
//                      for (int jj = 0; jj < numberGps; jj++)
//                        for (int ii = 0; ii < numberGps; ii++)
//                        {
//                          elem.GetGaussPoint(ii, jj).UpdateVariables();
//                        }
//                    }
//                  }
//                  else
//                  {
//                    var degreeOfParallelism = AbstractModel.NumberOfCPUs;
//                    var tasks = new Task[degreeOfParallelism];
//                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
//                    {
//                      // capturing taskNumber in lambda wouldn't work correctly
//                      int taskNumberCopy = taskNumber;

//                      tasks[taskNumber] = Task.Factory.StartNew(
//                           () =>
//                           {
//                             var max = numElement * (taskNumberCopy + 1) / degreeOfParallelism;
//                             for (int j = numElement * taskNumberCopy / degreeOfParallelism; j < max; j++)
//                             {
//                               AbstractElement2D elem = (AbstractElement2D)patch.GetElement(j);
//                               int numberGps = elem.GetNumberOfGaussPointOnEachDirection();
//                               for (int jj = 0; jj < numberGps; jj++)
//                                 for (int ii = 0; ii < numberGps; ii++)
//                                 {
//                                   elem.GetGaussPoint(ii, jj).UpdateVariables();
//                                 }
//                             }
//                           });
//                    }
//                    Task.WaitAll(tasks);
//                  }
//                }
//                break;
//              case Dimension.Solid:
//                for (int i = 0; i < model.CountPatch(); i++)
//                {
//                  AbstractPatch patch = model.GetPatch(i);
//                  int numElement = patch.CountElements();
//                  if (!AbstractModel.IsParallelProcesing)
//                  {
//                    for (int j = 0; j < numElement; j++)
//                    {
//                      AbstractElement3D elem = (AbstractElement3D)patch.GetElement(j);
//                      int numberGps = elem.GetNumberOfGaussPointOnEachDirection();
//                      for (int kk = 0; kk < numberGps; kk++)
//                        for (int jj = 0; jj < numberGps; jj++)
//                          for (int ii = 0; ii < numberGps; ii++)
//                          {
//                            elem.GetGaussPoint(ii, jj, kk).UpdateVariables();
//                          }
//                    }
//                  }
//                  else
//                  {
//                    var degreeOfParallelism = AbstractModel.NumberOfCPUs;
//                    var tasks = new Task[degreeOfParallelism];
//                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
//                    {
//                      // capturing taskNumber in lambda wouldn't work correctly
//                      int taskNumberCopy = taskNumber;

//                      tasks[taskNumber] = Task.Factory.StartNew(
//                           () =>
//                           {
//                             var max = numElement * (taskNumberCopy + 1) / degreeOfParallelism;
//                             for (int j = numElement * taskNumberCopy / degreeOfParallelism; j < max; j++)
//                             {
//                               AbstractElement3D elem = (AbstractElement3D)patch.GetElement(j);
//                               int numberGps = elem.GetNumberOfGaussPointOnEachDirection();
//                               for (int kk = 0; kk < numberGps; kk++)
//                                 for (int jj = 0; jj < numberGps; jj++)
//                                   for (int ii = 0; ii < numberGps; ii++)
//                                   {
//                                     elem.GetGaussPoint(ii, jj, kk).UpdateVariables();
//                                   }
//                             }
//                           });
//                    }
//                    Task.WaitAll(tasks);
//                  }
//                }
//                break;
//            }
//          ///////////////////////////////////////////////////////
//          if (!bConverged)
//          {

//            //Console.WriteLine("Not converged !!");
//            if (model.MonitorDataStore != null)
//            {
//              IOIGA.WriteExcelResult(model.PathProject + "monitor_" + model.NameProject + ".xls", model.MonitorDataStore);
//            }
//            model.stopWatchWholeModel.Stop();
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "-----------------------------------------------------------");
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "Not converged !!");
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "Elapsed time = " + model.stopWatchWholeModel.Elapsed);
//            IOIGA.WriteResult(model.logFile, model.PathProject + model.NameProject + ".log");
//            // Give a message: not converged
//            //throw new NotImplementedException("Not converged !!");
//            // Stop the solver
//            break;
//          }
//          if (model.IsWriteLogFile)
//          {
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "umax = " + un1.Max());
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "FR = " + FR.Sum());
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "Elapsed time = " + model.stopWatchWholeModel.Elapsed);
//          }
//          MatrixTool.SetColomn(model.DisplacementTime, un1, totalStep - 1);

//          /////////////////////////////////////////////////////////////////////////////
//          if (countData > 0)
//          {
//            MonitorDataHistorism monitorDataHistorism = model.GetMonitorDataHistorism();
//            DoubleMatrix data = model.MonitorDataStore;
//            for (int j = 0; j < countData; j++)
//            {
//              MonitorDataNode nodeData = monitorDataHistorism.GetMonitorDataNode(j);
//              double[] xi = nodeData.Xi;
//              Result result = nodeData.ResultType;
//              AbstractPatch patch = nodeData.Patch;
//              switch (result)
//              {
//                case Result.UX:
//                case Result.UY:
//                case Result.UZ:
//                case Result.PHASEFIELD:
//                case Result.VOLT:
//                case Result.TEMP:
//                  if (patch is PatchStructure2D)
//                  {
//                    data[totalStep - 1, j] = ((PatchStructure2D)patch).GetApproximateAt(result, xi[0], xi[1]);
//                  }
//                  else if (patch is PatchStructure3D)
//                  {
//                    data[totalStep - 1, j] = ((PatchStructure3D)patch).GetApproximateAt(result, xi[0], xi[1], xi[2]);
//                  }
//                  break;
//                case Result.USUM:
//                  if (patch is PatchStructure2D)
//                  {
//                    double[] u = new double[2];
//                    u[0] = ((PatchStructure2D)patch).GetApproximateAt(Result.UX, xi[0], xi[1]);
//                    u[1] = ((PatchStructure2D)patch).GetApproximateAt(Result.UY, xi[0], xi[1]);
//                    data[totalStep - 1, j] = Math.Sqrt(u[0] * u[0] + u[1] + u[1]);
//                  }
//                  else if (patch is PatchStructure3D)
//                  {
//                    double[] u = new double[3];
//                    u[0] = ((PatchStructure3D)patch).GetApproximateAt(Result.UX, xi[0], xi[1], xi[2]);
//                    u[1] = ((PatchStructure3D)patch).GetApproximateAt(Result.UY, xi[0], xi[1], xi[2]);
//                    u[2] = ((PatchStructure3D)patch).GetApproximateAt(Result.UZ, xi[0], xi[1], xi[2]);
//                    data[totalStep - 1, j] = Math.Sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
//                  }
//                  break;
//                case Result.SIGMAXX:
//                case Result.SIGMAYY:
//                case Result.SIGMAXY:
//                case Result.SIGMAZZ:
//                case Result.SIGMAYZ:
//                case Result.SIGMAXZ:
//                case Result.SIGMAEQV:
//                  if (result != Result.SIGMAEQV)
//                  {
//                    model.Extrapolation(result);
//                  }
//                  else
//                  {
//                    model.Extrapolation(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAZZ, Result.SIGMAXY, Result.SIGMAYZ, Result.SIGMAXZ);
//                  }
//                  if ((patch is PatchStructure3D) || ((patch is PatchStructure2D) && (result != Result.SIGMAYZ && result != Result.SIGMAXZ)))
//                  {
//                    data[totalStep - 1, j] = patch.GetApproximateAt(result, xi);
//                  }
//                  break;
//                case Result.EPSILONXX:
//                case Result.EPSILONYY:
//                case Result.EPSILONXY:
//                case Result.EPSILONZZ:
//                case Result.EPSILONYZ:
//                case Result.EPSILONXZ:
//                case Result.EPSILONEQV:
//                  if (result != Result.EPSILONEQV)
//                  {
//                    model.Extrapolation(result);
//                  }
//                  else
//                  {
//                    model.Extrapolation(Result.EPSILONXX, Result.EPSILONYY, Result.EPSILONZZ, Result.EPSILONXY, Result.EPSILONYZ, Result.EPSILONXZ);
//                  }
//                  if ((patch is PatchStructure3D) || ((patch is PatchStructure2D) && (result != Result.EPSILONYZ && result != Result.EPSILONXZ)))
//                  {
//                    data[totalStep - 1, j] = ((PatchStructure2D)patch).GetApproximateAt(result, xi);
//                  }
//                  break;
//                //case Result.EPSILONPXX:
//                //case Result.EPSILONPYY:
//                //case Result.EPSILONPXY:
//                //case Result.EPSILONPZZ:
//                //case Result.EPSILONPYZ:
//                //case Result.EPSILONPXZ:
//                //	ExtrapolationPlasticStrain();
//                //	if (patch is Structure2DPatch)
//                //	{
//                //		if (result != Result.EPSILONPYZ && result != Result.EPSILONPXZ)
//                //			data[totalStep - 1, j] = ((Structure2DPatch)patch).GetStrainAt(xi[0], xi[1], result);
//                //	}
//                //	else if (patch is Structure3DPatch)
//                //	{
//                //		data[totalStep - 1, j] = ((Structure3DPatch)patch).GetStrainAt(xi[0], xi[1], xi[2], result);
//                //	}
//                //	break;
//                //case Result.EPSILONPEQV:
//                //	ExtrapolationEquivalentPlasticStrain();
//                //	if (patch is Structure2DPatch)
//                //	{
//                //		data[totalStep - 1, j] = ((Structure2DPatch)patch).GetStrainAt(xi[0], xi[1], result);
//                //	}
//                //	else if (patch is Structure3DPatch)
//                //	{
//                //		data[totalStep - 1, j] = ((Structure3DPatch)patch).GetStrainAt(xi[0], xi[1], xi[2], result);
//                //	}
//                //	break;
//                case Result.REACFORCESUM:
//                  data[totalStep - 1, j] = FR.Sum();
//                  break;
//              }
//            }
//          }
//          ////////////////////////////////////////////////////////////////////////////


//          if (model.IsSaveStepByStep)
//          {
//            IOIGA.WriteSoapFormatter(fileNameDataTime, model.DisplacementTime);
//            IOIGA.WriteSoapFormatter(fileNameTime, model.Time);
//          }

//          un0 = un1;
//          substep++;
//          totalStep++;
//          if ((substep - 1) < numberOfSubstepLoad && ((substep - 1) % model.NumberOfStepStorageNonlinear == 0))
//          {
//            checkpointData.CurrentStepLoad = step;
//            checkpointData.CurrentSubStep = substep;
//            checkpointData.CurrentTotalStep = totalStep;
//            checkpointData.un0 = un0;
//            checkpointData.LastPrevFi0 = lastPrevF0;
//            checkpointData.LastPrevFi1 = lastPrevF1;
//            //store.IsFinished = IsFinished;
//            checkpointData.SaveData(filenameCheckpoint);
//            if (model.MonitorDataStore != null)
//            {
//              IOIGA.WriteExcelResult(model.PathProject + "monitor_" + model.NameProject + ".xls", model.MonitorDataStore);
//            }
//          }
//        }
//        if (!bConverged)
//          break;
//        step++;
//        if ((step - 1) < numberOfStepLoad)
//        {
//          checkpointData.CurrentStepLoad = step;
//          checkpointData.CurrentSubStep = 1;
//          checkpointData.CurrentTotalStep = totalStep;
//          checkpointData.un0 = un0;
//          checkpointData.LastPrevFi0 = lastPrevF0;
//          checkpointData.LastPrevFi1 = lastPrevF1;
//          //store.IsFinished = IsFinished;
//          checkpointData.SaveData(filenameCheckpoint);
//          if (model.MonitorDataStore != null)
//          {
//            IOIGA.WriteExcelResult(model.PathProject + "monitor_" + model.NameProject + ".xls", model.MonitorDataStore);
//          }
//        }
//      }
//      if (!model.IsSaveStepByStep)
//      {
//        IOIGA.WriteSoapFormatter(fileNameDataTime, model.DisplacementTime);
//        IOIGA.WriteSoapFormatter(fileNameTime, model.Time);
//      }
//      if (model.MonitorDataStore != null)
//      {
//        IOIGA.WriteExcelResult(model.PathProject + "monitor_" + model.NameProject + ".xls", model.MonitorDataStore);
//      }
//      for (int i = 0; i < filenameCheckpoint.Length; i++)
//      {
//        if (File.Exists(filenameCheckpoint[i]))
//          File.Delete(filenameCheckpoint[i]);
//      }
//      model.IsFinished = true;
//      //checkpointData.SaveData(model.filenameCheckpoint);
//    }

//    /// <summary>
//    /// Solve for one load step. External load: FStep = lambda*F (F is the global load vector)
//    /// </summary>
//    /// <param name="un"></param>
//    /// <param name="Fstep"></param>
//    /// <returns></returns>
//    private DoubleVector SolveOneSubstepLoad(DoubleVector un, DoubleVector Fstep, out DoubleVector FR, double time)
//    {
//      FR = null;
//      // Last converged values: un
//      double err1 = 1;
//      int nIter = 0;
//      bool bUpdate = true;
//      DoubleVector ui = un;
//      //DoubleVector unReduce = null;
//      int[] tArray = null;
//      if (!IsSolveSecondStaggered)
//      {
//        tArray = model.tArray[0];
//      }
//      else
//      {
//        tArray = model.tArray[1];
//      }
//      double wnorm = 0;
//      prevF = 0;
//      //if (AbstractModel.IsMonolithic)
//      //{
//      //  unReduce = un;
//      //}
//      //else
//      //{
//      //  unReduce = MatrixTool.GetSubVector(un, tArray);
//      //}

//      //switch (model.criteria)
//      //{
//      //  case CriterionConvergence.ForceCriterion:
//      //    wnorm = 1e-4;// Math.Sqrt(MatrixFunctions.Dot(rGlobal, rGlobal));Structural (<1e-2), thermal (<1e-6), Emg(<1e-12)
//      //    break;
//      //  case CriterionConvergence.EnergyCriterion:
//      //    wnorm = 1.0;// Math.Sqrt(MatrixFunctions.Dot(duGlobal, rGlobal));
//      //    break;
//      //}

//      while ((err1 > model.convergenceOption.TOL) && (nIter < model.convergenceOption.maximumIteration))
//      {
//        nIter++;
//        // Solve for each iteration
//        if (this.type == 0)
//        {
//          // Newton-Raphson: update tangent 						// stiffness at every iteration
//          bUpdate = true;
//        }
//        else if (this.type == 1)
//        {
//          // Modified Newton-Raphson: only 						// update tangent stiffness at the 						// first iteration (iter = 0)
//          if (nIter == 1) bUpdate = true;
//          else bUpdate = false;
//        }

//        DoubleVector ui1 = null;
//        DoubleVector fi = null;
//        DoubleVector f = null;
//        if (!model.IsSparseData)
//          ui1 = SolveOneIterationDense(Fstep, ui, ref err1, out FR, out f, out fi, bUpdate);
//        else
//          ui1 = SolveOneIterationSparse(Fstep, ui, ref err1, out FR, out f, out fi, bUpdate);

//        if (nIter == 1)
//        {
//          var du = MatrixTool.GetSubVector(ui1 - ui, tArray);
//          switch (model.convergenceOption.criterionConvergence)
//          {
//            case CriterionConvergence.DisplacementCriterion:
//              wnorm = Math.Sqrt(MatrixFunctions.Dot(du, du));
//              break;
//            case CriterionConvergence.ForceCriterion:
//              var g0 = fi;
//              wnorm = Math.Sqrt(MatrixFunctions.Dot(g0, g0));
//              break;
//            case CriterionConvergence.EnergyCriterion:
//              wnorm = MatrixFunctions.Dot(du, fi);// Math.Sqrt(MatrixFunctions.Dot(duGlobal, rGlobal));
//              break;
//            case CriterionConvergence.ResidualCriterion:
//              var fn = (AbstractModel.IsMonolithic || IsSolveSecondStaggered) ? lastPrevF0 : lastPrevF1;
//              wnorm = Math.Abs(Math.Sqrt(MatrixFunctions.Dot(fi, fi)) - fn);//Math.Abs(Math.Sqrt(MatrixFunctions.Dot(f, f)) - fn);
//              break;
//          }
//        }

//        if (wnorm != 0 && err1 != 0)
//          err1 = err1 / wnorm;
//        if (model.IsWriteLogFile)
//        {
//          IOIGA.WriteLogFileAndConsole(model.logFile, true, "\terror = " + err1);
//        }
//        ui = ui1;

//      }

//      if (model.IsWriteLogFile)
//      {
//        IOIGA.WriteLogFileAndConsole(model.logFile, true, "\tNumber of iteration = " + nIter);
//      }
//      if (nIter == model.convergenceOption.maximumIteration && !model.convergenceOption.IsSkipUnconvergenceStep)
//      {
//        this.bConverged = false;
//        //return null;
//      }
//      //else
//      //{
//      //  if (model.GetPatch(0).Material.TypeMaterialStructure == EngineeringData.TypeMaterialStructure.Plasticity || AbstractModel.TypeModel == TypeModelProblem.PhaseField)
//      //    switch (model.StructureDimension)
//      //    {
//      //      case Dimension.Plane:
//      //        for (int i = 0; i < model.CountPatch(); i++)
//      //        {
//      //          AbstractPatch patch = model.GetPatch(i);
//      //          int numElement = patch.CountElements();
//      //          if (!AbstractModel.IsParallelProcesing)
//      //          {
//      //            for (int j = 0; j < numElement; j++)
//      //            {
//      //              AbstractElement2D elem = (AbstractElement2D)patch.GetElement(j);
//      //              int numberGps = elem.GetNumberOfGaussPointOnEachDirection();
//      //              for (int jj = 0; jj < numberGps; jj++)
//      //                for (int ii = 0; ii < numberGps; ii++)
//      //                {
//      //                  elem.GetGaussPoint(ii, jj).UpdateVariables();
//      //                }
//      //            }
//      //          }
//      //          else
//      //          {
//      //            var degreeOfParallelism = AbstractModel.NumberOfCPUs;
//      //            var tasks = new Task[degreeOfParallelism];
//      //            for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
//      //            {
//      //              // capturing taskNumber in lambda wouldn't work correctly
//      //              int taskNumberCopy = taskNumber;

//      //              tasks[taskNumber] = Task.Factory.StartNew(
//      //                   () =>
//      //                   {
//      //                     var max = numElement * (taskNumberCopy + 1) / degreeOfParallelism;
//      //                     for (int j = numElement * taskNumberCopy / degreeOfParallelism; j < max; j++)
//      //                     {
//      //                       AbstractElement2D elem = (AbstractElement2D)patch.GetElement(j);
//      //                       int numberGps = elem.GetNumberOfGaussPointOnEachDirection();
//      //                       for (int jj = 0; jj < numberGps; jj++)
//      //                         for (int ii = 0; ii < numberGps; ii++)
//      //                         {
//      //                           elem.GetGaussPoint(ii, jj).UpdateVariables();
//      //                         }
//      //                     }
//      //                   });
//      //            }
//      //            Task.WaitAll(tasks);
//      //          }
//      //        }
//      //        break;
//      //      case Dimension.Solid:
//      //        for (int i = 0; i < model.CountPatch(); i++)
//      //        {
//      //          AbstractPatch patch = model.GetPatch(i);
//      //          int numElement = patch.CountElements();
//      //          if (!AbstractModel.IsParallelProcesing)
//      //          {
//      //            for (int j = 0; j < numElement; j++)
//      //            {
//      //              AbstractElement3D elem = (AbstractElement3D)patch.GetElement(j);
//      //              int numberGps = elem.GetNumberOfGaussPointOnEachDirection();
//      //              for (int kk = 0; kk < numberGps; kk++)
//      //                for (int jj = 0; jj < numberGps; jj++)
//      //                  for (int ii = 0; ii < numberGps; ii++)
//      //                  {
//      //                    elem.GetGaussPoint(ii, jj, kk).UpdateVariables();
//      //                  }
//      //            }
//      //          }
//      //          else
//      //          {
//      //            var degreeOfParallelism = AbstractModel.NumberOfCPUs;
//      //            var tasks = new Task[degreeOfParallelism];
//      //            for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
//      //            {
//      //              // capturing taskNumber in lambda wouldn't work correctly
//      //              int taskNumberCopy = taskNumber;

//      //              tasks[taskNumber] = Task.Factory.StartNew(
//      //                   () =>
//      //                   {
//      //                     var max = numElement * (taskNumberCopy + 1) / degreeOfParallelism;
//      //                     for (int j = numElement * taskNumberCopy / degreeOfParallelism; j < max; j++)
//      //                     {
//      //                       AbstractElement3D elem = (AbstractElement3D)patch.GetElement(j);
//      //                       int numberGps = elem.GetNumberOfGaussPointOnEachDirection();
//      //                       for (int kk = 0; kk < numberGps; kk++)
//      //                         for (int jj = 0; jj < numberGps; jj++)
//      //                           for (int ii = 0; ii < numberGps; ii++)
//      //                           {
//      //                             elem.GetGaussPoint(ii, jj, kk).UpdateVariables();
//      //                           }
//      //                     }
//      //                   });
//      //            }
//      //            Task.WaitAll(tasks);
//      //          }
//      //        }
//      //        break;
//      //    }
//      //}
//      return ui;
//    }

//    /// <summary>
//    /// Solve for one iteration.
//    /// </summary>
//    /// <param name="Fstep"></param>
//    /// <param name="ui"></param>
//    /// <param name="err1"></param>
//    /// <param name="bUpdate"></param>
//    /// <returns></returns>
//    private DoubleVector SolveOneIterationDense(DoubleVector Fstep, DoubleVector ui, ref double err1, out DoubleVector FR, out DoubleVector f, out DoubleVector fi, bool bUpdate)
//    {
//      DoubleMatrix kGlobal = null;
//      DoubleVector internalForceGlobal = null;
//      FR = null;

//      // Update Stiffness matrix and Internal force
//      if (bUpdate == true)
//      {
//        //model.AssemblyStiffnessMatrixAndInternalForce(out kGlobal, out residualGlobal);
//        model.AssemblyStiffnessMatrix(out kGlobal);
//        model.AssemblyInternalForce(out internalForceGlobal);
//      }
//      else
//      {
//        model.AssemblyInternalForce(out internalForceGlobal);
//      }
//      fi = internalForceGlobal;
//      DoubleVector rGlobal = -internalForceGlobal;
//      if (Fstep != null)
//      {
//        rGlobal += Fstep;
//      }
//      f = rGlobal;
//      DoubleVector rGlobal1 = internalForceGlobal;
//      //int numDOFUnConstraint = 0;
//      // Compute length of Unconstraint DOF
//      int[] tArray = null;
//      if (AbstractModel.IsMonolithic)
//      {
//        //if (model.tArrayUnconstraint[0].Length > 0)
//        //{
//        //  for (int i = 0; i < model.tArrayUnconstraint.Length; i++)
//        //    tArray = tArray.Union(model.tArrayUnconstraint[i]).ToArray();
//        //}
//        //else
//        //{
//        //  for (int i = 0; i < model.tArray.Length; i++)
//        //    tArray = tArray.Union(model.tArray[i]).ToArray();
//        //}
//        tArray = model.tArray[0];
//        for (int i = 1; i < model.tArray.Length; i++)
//          tArray = tArray.Union(model.tArray[i]).ToArray();
//      }
//      else
//      {
//        if (!IsSolveSecondStaggered)
//        {
//          tArray = model.tArray[0];
//        }
//        else
//        {
//          tArray = model.tArray[1];
//        }
//      }
//      f = MatrixTool.GetSubVector(f, tArray);
//      fi = MatrixTool.GetSubVector(fi, tArray);

//      int[] tArrayIndependent = tArray.Intersect(model.tArrayIndependent).ToArray();
//      int[] tArrayDependent = tArray.Intersect(model.tArrayDependent).ToArray();
//      int[] tArrayUnCoupleAndUnconstraint = tArray.Intersect(model.tArrayUncouple).ToArray();

//      DoubleMatrix Kuu = null;
//      DoubleMatrix Kus = null;
//      DoubleMatrix Kum = null;
//      DoubleMatrix Kmm = null;
//      DoubleVector Ru = null;
//      DoubleVector Rm = null;
//      if (tArrayUnCoupleAndUnconstraint.Length != 0)
//      {
//        Kuu = MatrixTool.GetSubMatrix(kGlobal, tArrayUnCoupleAndUnconstraint, tArrayUnCoupleAndUnconstraint);
//        Ru = MatrixTool.GetSubVector(rGlobal, tArrayUnCoupleAndUnconstraint);
//        if (tArrayDependent.Length != 0)
//        {
//          Kus = MatrixTool.GetSubMatrix(kGlobal, tArrayUnCoupleAndUnconstraint, tArrayDependent);
//        }
//      }

//      DoubleMatrix Tsm = null;
//      if (tArrayDependent.Length != 0 && tArrayIndependent.Length != 0)
//      {
//        Tsm = MatrixTool.GetSubMatrix(model.TDense, tArrayDependent, tArrayIndependent);
//        DoubleMatrix TsmTranspose = Tsm.Transpose();
//        DoubleMatrix Kms = MatrixTool.GetSubMatrix(kGlobal, tArrayIndependent, tArrayDependent);
//        DoubleMatrix Kss = MatrixTool.GetSubMatrix(kGlobal, tArrayDependent, tArrayDependent);
//        DoubleVector Rs = MatrixTool.GetSubVector(rGlobal, tArrayDependent);
//        Kmm = MatrixTool.GetSubMatrix(kGlobal, tArrayIndependent, tArrayIndependent) + MatrixFunctions.Product(TsmTranspose, Kms.Transpose() + MatrixFunctions.Product(Kss, Tsm)) + MatrixFunctions.Product(Kms, Tsm);
//        Rm = MatrixTool.GetSubVector(rGlobal, tArrayIndependent) + MatrixFunctions.Product(TsmTranspose, Rs);
//        if (tArrayUnCoupleAndUnconstraint.Length != 0)
//          Kum = MatrixTool.GetSubMatrix(kGlobal, tArrayUnCoupleAndUnconstraint, tArrayIndependent) + MatrixFunctions.Product(Kus, Tsm);
//      }
//      kGlobal = new DoubleMatrix(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length, tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
//      rGlobal = new DoubleVector(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
//      if (!AbstractModel.IsParallelProcesing)
//      {
//        for (int i = 0; i < tArrayUnCoupleAndUnconstraint.Length; i++)
//        {
//          for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
//          {
//            kGlobal[i, j] = Kuu[i, j];
//          }
//          if (Kum != null)
//          {
//            for (int j = 0; j < tArrayIndependent.Length; j++)
//            {
//              kGlobal[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
//            }
//          }
//          rGlobal[i] = Ru[i];
//        }
//        for (int i = 0; i < tArrayIndependent.Length; i++)
//        {
//          for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
//          {
//            kGlobal[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
//          }

//          for (int j = 0; j < tArrayIndependent.Length; j++)
//          {
//            kGlobal[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
//          }
//          rGlobal[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
//        }
//      }
//      else
//      {
//        var degreeOfParallelism = AbstractModel.NumberOfCPUs;
//        var tasks = new Task[degreeOfParallelism];
//        for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
//        {
//          // capturing taskNumber in lambda wouldn't work correctly
//          int taskNumberCopy = taskNumber;

//          tasks[taskNumber] = Task.Factory.StartNew(
//               () =>
//               {
//                 var max = tArrayUnCoupleAndUnconstraint.Length * (taskNumberCopy + 1) / degreeOfParallelism;
//                 for (int i = tArrayUnCoupleAndUnconstraint.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
//                 {
//                   for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
//                   {
//                     kGlobal[i, j] = Kuu[i, j];
//                   }
//                   if (Kum != null)
//                   {
//                     for (int j = 0; j < tArrayIndependent.Length; j++)
//                     {
//                       kGlobal[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
//                     }
//                   }
//                   rGlobal[i] = Ru[i];
//                 }
//               });
//        }
//        Task.WaitAll(tasks);

//        for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
//        {
//          // capturing taskNumber in lambda wouldn't work correctly
//          int taskNumberCopy = taskNumber;

//          tasks[taskNumber] = Task.Factory.StartNew(
//               () =>
//               {
//                 var max = tArrayIndependent.Length * (taskNumberCopy + 1) / degreeOfParallelism;
//                 for (int i = tArrayIndependent.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
//                 {
//                   for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
//                   {
//                     kGlobal[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
//                   }

//                   for (int j = 0; j < tArrayIndependent.Length; j++)
//                   {
//                     kGlobal[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
//                   }
//                   rGlobal[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
//                 }
//               });
//        }
//        Task.WaitAll(tasks);
//      }
//      int[] tArrayIndependentInNew = new int[tArrayIndependent.Length];
//      int c = 0;
//      for (int i = tArrayUnCoupleAndUnconstraint.Length; i < tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length; i++)
//      {
//        tArrayIndependentInNew[c++] = i;
//      }

//      int[] tArrayUncoupleUnconstraintInNew = new int[tArrayUnCoupleAndUnconstraint.Length];
//      for (int i = 0; i < tArrayUnCoupleAndUnconstraint.Length; i++)
//      {
//        tArrayUncoupleUnconstraintInNew[i] = i;
//      }

//      DoubleVector duLocalUM = null;
//      DoubleVector duLocalU = null;
//      DoubleVector duLocalM = null;
//      DoubleVector duLocalS = null;
//      if (rGlobal.Length != 0)
//      {
//        duLocalUM = MatrixFunctions.Solve(kGlobal, rGlobal);
//        duLocalU = MatrixTool.GetSubVector(duLocalUM, tArrayUncoupleUnconstraintInNew);
//        duLocalM = MatrixTool.GetSubVector(duLocalUM, tArrayIndependentInNew);
//        if (duLocalM != null)
//          duLocalS = MatrixFunctions.Product(Tsm, duLocalM);
//      }

//      DoubleVector duGlobal = new DoubleVector(model.CountDOF());
//      if (!AbstractModel.IsParallelProcesing)
//      {
//        for (int i = 0; i < tArrayUncoupleUnconstraintInNew.Length; i++)
//        {
//          duGlobal[tArrayUnCoupleAndUnconstraint[i]] = duLocalU[i];
//        }
//        for (int i = 0; i < tArrayIndependentInNew.Length; i++)
//        {
//          duGlobal[tArrayIndependent[i]] = duLocalM[i];
//        }
//        for (int i = 0; i < tArrayDependent.Length; i++)
//        {
//          duGlobal[tArrayDependent[i]] = duLocalS[i];
//        }
//      }
//      else
//      {
//        var degreeOfParallelism = AbstractModel.NumberOfCPUs;
//        var tasks = new Task[degreeOfParallelism];
//        for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
//        {
//          // capturing taskNumber in lambda wouldn't work correctly
//          int taskNumberCopy = taskNumber;

//          tasks[taskNumber] = Task.Factory.StartNew(
//               () =>
//               {
//                 var max = tArrayUncoupleUnconstraintInNew.Length * (taskNumberCopy + 1) / degreeOfParallelism;
//                 for (int i = tArrayUncoupleUnconstraintInNew.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
//                 {
//                   duGlobal[tArrayUnCoupleAndUnconstraint[i]] = duLocalU[i];
//                 }
//               });
//        }
//        Task.WaitAll(tasks);
//        for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
//        {
//          // capturing taskNumber in lambda wouldn't work correctly
//          int taskNumberCopy = taskNumber;

//          tasks[taskNumber] = Task.Factory.StartNew(
//               () =>
//               {
//                 var max = tArrayIndependentInNew.Length * (taskNumberCopy + 1) / degreeOfParallelism;
//                 for (int i = tArrayIndependentInNew.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
//                 {
//                   duGlobal[tArrayIndependent[i]] = duLocalM[i];
//                 }
//               });
//        }
//        Task.WaitAll(tasks);

//        for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
//        {
//          // capturing taskNumber in lambda wouldn't work correctly
//          int taskNumberCopy = taskNumber;

//          tasks[taskNumber] = Task.Factory.StartNew(
//               () =>
//               {
//                 var max = tArrayDependent.Length * (taskNumberCopy + 1) / degreeOfParallelism;
//                 for (int i = tArrayDependent.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
//                 {
//                   duGlobal[tArrayDependent[i]] = duLocalS[i];
//                 }
//               });
//        }
//        Task.WaitAll(tasks);
//      }

//      err1 = ErrorCalculationIteration(rGlobal, duGlobal, f, fi);
//      if (err1 == double.NaN)
//      {
//        this.bConverged = false;
//        //return null;
//      }
//      if (AbstractModel.IsMonolithic)
//      {
//        FR = MatrixTool.GetSubVector(rGlobal1, tArrayConstraintNoneZero[0]);
//      }
//      else
//      {
//        if (!IsSolveSecondStaggered)
//        {
//          FR = MatrixTool.GetSubVector(rGlobal1, tArrayConstraintNoneZero[0]);
//        }
//      }
//      DoubleVector uGlobal = model.GetUGlobal() + duGlobal;
//      model.SetUGlobal(uGlobal);
//      return uGlobal;
//    }

//    private double ErrorCalculationIteration(DoubleVector rGlobal, DoubleVector duGlobal, DoubleVector f, DoubleVector fi)
//    {
//      double err1 = 0;
//      switch (model.convergenceOption.criterionConvergence)
//      {
//        case CriterionConvergence.DisplacementCriterion:
//          err1 = Math.Sqrt(MatrixFunctions.Dot(duGlobal, duGlobal));
//          break;
//        case CriterionConvergence.ForceCriterion:
//          err1 = Math.Sqrt(MatrixFunctions.Dot(f, f));
//          break;
//        case CriterionConvergence.EnergyCriterion:
//          err1 = MatrixFunctions.Dot(duGlobal, rGlobal);
//          break;
//        case CriterionConvergence.ResidualCriterion:
//          var curF = Math.Sqrt(MatrixFunctions.Dot(f, f));
//          err1 = Math.Abs(curF - prevF);
//          prevF = curF;
//          break;
//      }

//      return err1;
//    }

//    private DoubleVector SolveOneIterationSparse(DoubleVector Fstep, DoubleVector ui, ref double err1, out DoubleVector FR, out DoubleVector f, out DoubleVector fi, bool bUpdate)
//    {
//      DoubleCsrSparseMatrix kGlobalSparse = (DoubleCsrSparseMatrix)model.kGlobalSparse.Clone();
//      DoubleVector internalForceGlobal = null;
//      FR = null;

//      // Update Stiffness matrix and Internal force
//      if (bUpdate == true)
//      {
//        //model.AssemblyStiffnessMatrixAndInternalForce(out kGlobal, out residualGlobal);
//        model.AssemblyStiffnessMatrix(ref kGlobalSparse);
//        model.AssemblyInternalForce(out internalForceGlobal);
//      }
//      else
//      {
//        model.AssemblyInternalForce(out internalForceGlobal);
//      }
//      fi = internalForceGlobal;
//      DoubleVector rGlobal = -internalForceGlobal;
//      if (Fstep != null)
//      {
//        rGlobal += Fstep;
//      }
//      f = (DoubleVector)rGlobal.Clone();
//      DoubleVector rGlobal1 = internalForceGlobal;

//      int[] tArray = null;
//      if (AbstractModel.IsMonolithic)
//      {
//        tArray = model.tArray[0];
//        for (int i = 1; i < model.tArray.Length; i++)
//          tArray = tArray.Union(model.tArray[i]).ToArray();
//      }
//      else
//      {
//        if (!IsSolveSecondStaggered)
//        {
//          tArray = model.tArray[0];
//        }
//        else
//        {
//          tArray = model.tArray[1];
//        }
//      }
//      f = MatrixTool.GetSubVector(f, tArray);
//      fi = MatrixTool.GetSubVector(fi, tArray);

//      model.ApplyTAndg(ref kGlobalSparse, ref rGlobal, null);
//      DoubleVector duGlobal = MatrixFunctions.Solve(kGlobalSparse, rGlobal);

//      duGlobal = MatrixFunctions.Product(model.TSparse, duGlobal) + model.g;
//      err1 = ErrorCalculationIteration(rGlobal, duGlobal, f, fi);
//      if (err1 == double.NaN)
//      {
//        this.bConverged = false;
//        //return null;
//      }
//      if (AbstractModel.IsMonolithic)
//      {
//        FR = MatrixTool.GetSubVector(rGlobal1, tArrayConstraintNoneZero[0]);
//      }
//      else
//      {
//        if (!IsSolveSecondStaggered)
//        {
//          FR = MatrixTool.GetSubVector(rGlobal1, tArrayConstraintNoneZero[0]);
//        }
//      }
//      DoubleVector uGlobal = model.GetUGlobal() + duGlobal;
//      model.SetUGlobal(uGlobal);
//      return uGlobal;
//    }
//    private void ComputeDisplacementVector(int stepLoad, int[] tArrayContraint, out DoubleVector DUNoneZeroVector, out int[] tArrayContraintNoneZero)
//    {
//      List<double> DUNone = new List<double>();
//      List<int> tArrayNone = new List<int>();
//      DUNoneZeroVector = null;
//      tArrayContraintNoneZero = null;
//      foreach (AbstractConstraintValue c in model.listConstraint)
//      {
//        double valueConstraint = 0;
//        if ((c.GetValueConstraint() is double) || (c.GetValueConstraint() is ConstantFunctionRToR))
//          valueConstraint = c.GetValueConstraint(stepLoad + 1) - c.GetValueConstraint(stepLoad);
//        ControlPoint[] cpsel = c.GetControlPointsConstrainted();
//        double[] valueFunctionConstraint = null;
//        if (c is ConstraintValueEdge2D && (c.GetValueConstraint() is FunctionRToR) && (!(c.GetValueConstraint() is ConstantFunctionRToR)))
//        {
//          valueFunctionConstraint = (new DoubleVector(c.ComputeValueConstraintOnControlPoints(stepLoad + 1)) - new DoubleVector(c.ComputeValueConstraintOnControlPoints(stepLoad))).ToArray();
//        }
//        else if (c is ConstraintValueFace3D)
//        {

//        }

//        for (int k = 0; k < cpsel.Length; k++)
//        {
//          if ((c.GetValueConstraint() is FunctionRToR) && (!(c.GetValueConstraint() is ConstantFunctionRToR)))
//          {
//            valueConstraint = valueFunctionConstraint[k];
//          }
//          int[] tArrayCp = cpsel[k].GetTArrayGlobal();
//          int tArray = tArrayCp[c.GetIDOfField()];
//          int indexConstraint = Array.IndexOf(tArrayContraint, tArray);
//          if (indexConstraint != -1 && valueConstraint != 0 && tArrayNone.IndexOf(tArray) == -1/*not exist in list*/)
//          {
//            DUNone.Add(valueConstraint);
//            tArrayNone.Add(tArray);
//          }
//        }
//      }
//      if (DUNone.Count > 0)
//      {
//        DUNoneZeroVector = new DoubleVector(DUNone.ToArray());
//        tArrayContraintNoneZero = tArrayNone.ToArray();
//      }
//    }
//  }
//}
