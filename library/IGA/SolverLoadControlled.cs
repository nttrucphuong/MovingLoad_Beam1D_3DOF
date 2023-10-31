//using System;
//using System.Linq;
//using CenterSpace.NMath.Core;
//using System.Collections.Generic;
//using System.Threading.Tasks;
//using CenterSpace.NMath.Matrix;

//namespace DEMSoft.IGA
//{

//  /// <summary>
//  /// Implement the load controlled Newton-Raphson solver for non-linear problems
//  /// </summary>
//  internal class SolverLoadControlled : AbstractNonlinearSolver
//  {
//    /// <summary>
//    /// Constructor class
//    /// </summary>
//    /// <param name="model"></param>
//    /// <param name="type"></param>
//    public SolverLoadControlled(AbstractModelStructure model, TypeNonlinearSolver type)
//      : base(model, type)
//    {
//    }

//    /// <summary>
//    /// Solve
//    /// </summary>
//    protected override void Solve(TypeNonlinearSolver type)
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
//      //model.DisplacementTime = new List<DoubleVector>();// (model.CountDOF(), numberOfTotalStep);
//      model.Time = new List<double>();//new DoubleVector(numberOfTotalStep);
//      if (countData > 0)
//        model.MonitorDataStore = new List<double[]>();
//      DoubleVector un0 = new DoubleVector(model.CountDOF());
//      DoubleVector un1 = null;
//      if ((!model.IsRunFromInitial) && (!checkpointData.IsFinished) && conditionRestoreCheckpoint)
//      {
//        totalStep = checkpointData.CurrentTotalStep;
//        step = checkpointData.CurrentStepLoad;
//        //model.DisplacementTime = checkpointData.DisplacementTime;
//        model.Time = checkpointData.Time;
//        if (countData > 0)
//          model.MonitorDataStore = checkpointData.MonitorDataStore;
//        un0 = checkpointData.un0;
//      }

//      //string[] filenameCheckpoint = model.filenameCheckpoint;
//      string fileNameTime = model.FileNameTime;
//      while (step <= numberOfStepLoad)
//      {
//        if (model.IsWriteLogFile)
//        {
//          IOIGA.WriteLogFileAndConsole(model.logFile, true, "-----------------------------------------------------------");
//          IOIGA.WriteLogFileAndConsole(model.logFile, true, "Load step " + step + " : ");
//        }
//        int numberOfSubstepLoad = model.numberOfSubstepLoad[step - 1];
//        int substep = ((model.IsRunFromInitial) || (checkpointData.IsFinished) || !conditionRestoreCheckpoint) ? 1 : checkpointData.CurrentSubStep;//1;

//        while (substep <= numberOfSubstepLoad)
//        {
//          if (model.IsWriteLogFile)
//          {
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "-----------------------------------------------------------");
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "SubLoad step " + substep + " : ");
//          }
//          //////////////////////////////////////////////////////////////////
//          //double loadFactor = (double)step / (double)numLoadStep;
//          //DoubleVector fStep = loadFactor * f;
//          //////////////////////////////////////////////////////////////////
//          double time = step - 1 + substep * model.deltaT[step - 1];
//          deltaT = model.deltaT[step - 1];
//          model.Time.Add(time);
//          DoubleVector fStep = null;
//          model.AssemblyTractionVector(out fStep, time);
//          DoubleVector FR = null;


//          // Solve
//          //un1 = SolveOneSubstepLoad(un0, fStep, out FR, time);

//          switch (AbstractModel.TypeSchemeNonlinearSolver)
//          {
//            case TypeSchemeNonlinearSolver.Monolithic:
//            case TypeSchemeNonlinearSolver.SingleStaggered:
//              un1 = SolveOneSubstepLoad(un0, fStep, out FR, time);
//              break;
//            case TypeSchemeNonlinearSolver.Staggered:
//              // Solve phase field
//              IsSolvePhasefieldStaggered = true;
//              un1 = SolveOneSubstepLoad(un0, null, out DoubleVector FR1, time);
//              // Solve displacement
//              IsSolvePhasefieldStaggered = false;
//              un1 = SolveOneSubstepLoad(un1, fStep, out FR, time);
//              break;
//          }

//          UpdateConvergedGaussPoints();
//          ///////////////////////////////////////////////////////
//          if (!bConverged)
//          {
//            if (model.MonitorDataStore != null)
//            {
//              IOIGA.WriteExcelResult(model.PathProject + "monitor_" + model.NameProject + ".xls", model.MonitorDataStore);
//            }
//            model.stopWatchWholeModel.Stop();
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "-----------------------------------------------------------");
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "Not converged !!");
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "Elapsed time = " + model.stopWatchWholeModel.Elapsed);
//            IOIGA.WriteResult(model.logFile, model.PathProject + model.NameProject + ".log");
//            break;
//          }
//          if (model.IsWriteLogFile)
//          {
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "umax = " + un1.Max());
//            IOIGA.WriteLogFileAndConsole(model.logFile, true, "Elapsed time = " + model.stopWatchWholeModel.Elapsed);
//          }
//          //model.DisplacementTime.Add(un1);
//          ComputeMonitorResults(totalStep, countData, out double valueReactionForce);
//          ////////////////////////////////////////////////////////////////////////////


//          un0 = un1;
//          substep++;
//          totalStep++;
//          if ((substep - 1) < numberOfSubstepLoad && ((substep - 1) % model.NumberOfStepStorageNonlinear == 0))
//          {
//            checkpointData.CurrentStepLoad = step;
//            checkpointData.CurrentSubStep = substep;
//            checkpointData.CurrentTotalStep = totalStep;
//            checkpointData.un0 = un0;
//            //store.IsFinished = IsFinished;
//            checkpointData.SaveData();
//            if (model.MonitorDataStore != null)
//            {
//              IOIGA.WriteExcelResult(model.PathProject + "monitor_" + model.NameProject + ".xls", model.MonitorDataStore);
//            }
//            if (model.IsSaveStepByStep)
//            {
//              //IOIGA.WriteListDoubleVectorMessagePackFormatter(fileNameDataTime, model.DisplacementTime);
//              string fileNameDataTime = model.PathProject + "data_" + (totalStep - 1) + ".msgpack";
//              IOIGA.WriteMessagePackDoubleVector(fileNameDataTime, un0);
//              IOIGA.WriteMessagePackData(fileNameTime, model.Time);
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
//          //store.IsFinished = IsFinished;
//          checkpointData.SaveData();
//          if (model.MonitorDataStore != null)
//          {
//            IOIGA.WriteExcelResult(model.PathProject + "monitor_" + model.NameProject + ".xls", model.MonitorDataStore);
//          }
//          if (model.IsSaveStepByStep)
//          {
//            //IOIGA.WriteListDoubleVectorMessagePackFormatter(fileNameDataTime, model.DisplacementTime);
//            string fileNameDataTime = model.PathProject + "data_" + (totalStep - 1) + ".msgpack";
//            IOIGA.WriteMessagePackDoubleVector(fileNameDataTime, un0);
//            IOIGA.WriteMessagePackData(fileNameTime, model.Time);
//          }
//        }
//      }
//      if (!model.IsSaveStepByStep)
//      {
//        //IOIGA.WriteListDoubleVectorMessagePackFormatter(fileNameDataTime, model.DisplacementTime);
//        string fileNameDataTime = model.PathProject + "data_" + (totalStep - 1) + ".msgpack";
//        IOIGA.WriteMessagePackDoubleVector(fileNameDataTime, un0);
//        IOIGA.WriteMessagePackData(fileNameTime, model.Time);
//      }
//      if (model.MonitorDataStore != null)
//      {
//        IOIGA.WriteExcelResult(model.PathProject + "monitor_" + model.NameProject + ".xls", model.MonitorDataStore);
//      }
//      //model.IsFinished = true;
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

//      int nIter = 0;

//      DoubleVector ui = un;

//      int[] tArray = null;

//      switch (AbstractModel.TypeSchemeNonlinearSolver)
//      {
//        case TypeSchemeNonlinearSolver.Monolithic:
//        case TypeSchemeNonlinearSolver.SingleStaggered:
//          tArray = model.tArray[0];
//          for (int i = 1; i < model.tArray.Length; i++)
//            tArray = tArray.Union(model.tArray[i]).ToArray();
//          break;
//        case TypeSchemeNonlinearSolver.Staggered:
//          if (!IsSolvePhasefieldStaggered)
//          {
//            tArray = model.tArray[0];
//          }
//          else
//          {
//            tArray = model.tArray[1];
//          }
//          break;
//      }

//      DoubleVector internalForceGlobal = null;
//      model.AssemblyInternalForce(out internalForceGlobal);
//      DoubleVector residualForce = Fstep - internalForceGlobal;
//      DoubleCsrSparseMatrix kGlobalSparse = null;
//      DoubleSymmetricMatrix kGlobal = null;
//      if (!model.IsSparseData)
//      {
//        // Update Stiffness matrix and Internal force
//        if (bUpdate == false)
//        {
//          model.AssemblyStiffnessMatrix(out kGlobal);
//        }
//      }
//      else
//      {
//        kGlobalSparse = (DoubleCsrSparseMatrix)model.kGlobalSparse.Clone();
//        if (bUpdate == false)
//        {
//          model.AssemblyStiffnessMatrix(ref kGlobalSparse);
//        }
//      }

//      double[] err1 = new double[model.tArray.Length];
//      for (int i = 0; i < err1.Length; i++)
//      {
//        err1[i] = 1;
//      }
//      double err = 1;
//      double wnorm = 0;
//      switch (model.convergenceOption.criterionConvergence)
//      {
//        case CriterionConvergence.ForceCriterion:
//          wnorm = 1 + Math.Sqrt(MatrixFunctions.Dot(Fstep, Fstep));
//          err1 = ErrorCalculationIteration(null, residualForce, internalForceGlobal, internalForceGlobal, tArray);
//          if (err1 == double.NaN)
//          {
//            this.bConverged = false;
//            //return null;
//          }
//          if (wnorm != 0 && err1 != 0)
//            err1 = err1 / wnorm;
//          break;
//      }

//      while ((err1 > model.convergenceOption.TOL) && (nIter < model.convergenceOption.maximumIteration))
//      {
//        nIter++;
//        // Solve for each iteration

//        DoubleVector ui1 = null;

//        if (!model.IsSparseData)
//        {
//          ui1 = SolveOneIterationDense(kGlobal, Fstep, ui, ref err1, ref residualForce);
//        }
//        else
//        {
//          ui1 = SolveOneIterationSparse(kGlobalSparse, Fstep, ui, ref err1, ref residualForce);
//        }
//        if (nIter == 1)
//        {
//          var du = ui1 - ui;
//          switch (model.convergenceOption.criterionConvergence)
//          {
//            case CriterionConvergence.DisplacementCriterion:
//              wnorm = 1 + Math.Sqrt(MatrixFunctions.Dot(du, du));
//              break;
//            case CriterionConvergence.EnergyCriterion:
//              wnorm = 1 + MatrixFunctions.Dot(du, Fstep);// Math.Sqrt(MatrixFunctions.Dot(duGlobal, rGlobal));
//              break;
//            case CriterionConvergence.ResidualCriterion:
//              //var fn = (AbstractModel.IsMonolithic || IsSolveSecondStaggered) ? lastPrevF0 : lastPrevF1;
//              //wnorm = Math.Abs(Math.Sqrt(MatrixFunctions.Dot(fi, fi)) - fn);//Math.Abs(Math.Sqrt(MatrixFunctions.Dot(f, f)) - fn);
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
//      }
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
//    protected override DoubleVector SolveOneIterationDense(DoubleSymmetricMatrix kGlobal, DoubleVector Fstep, DoubleVector ui, ref double[] err1, ref DoubleVector residualForce)
//    {
//      if (bUpdate == true)
//      {
//        model.AssemblyStiffnessMatrix(out kGlobal);
//      }
//      //DoubleVector residualForceClone = (DoubleVector)residualForce.Clone();
//      // Compute length of Unconstraint DOF
//      int[] tArray = null;

//      switch (AbstractModel.TypeSchemeNonlinearSolver)
//      {
//        case TypeSchemeNonlinearSolver.Monolithic:
//        case TypeSchemeNonlinearSolver.SingleStaggered:
//          tArray = model.tArray[0];
//          for (int i = 1; i < model.tArray.Length; i++)
//            tArray = tArray.Union(model.tArray[i]).ToArray();
//          break;
//        case TypeSchemeNonlinearSolver.Staggered:
//          if (!IsSolvePhasefieldStaggered)
//          {
//            tArray = model.tArray[0];
//          }
//          else
//          {
//            tArray = model.tArray[1];
//          }
//          break;
//      }

//      int[] tArrayIndependent = tArray.Intersect(model.tArrayIndependent).ToArray();
//      int[] tArrayDependent = tArray.Intersect(model.tArrayDependent).ToArray();
//      int[] tArrayUnCoupleAndUnconstraint = tArray.Intersect(model.tArrayUncouple).ToArray();

//      DoubleSymmetricMatrix Kuu = null;
//      DoubleMatrix Kus = null;
//      DoubleMatrix Kum = null;
//      DoubleMatrix Kmm = null;
//      DoubleVector Ru = null;
//      DoubleVector Rm = null;

//      DoubleVector gS = null;
//      if (tArrayDependent.Length != 0)
//      {
//        gS = MatrixTool.GetSubVector(model.g, tArrayDependent);
//      }
//      if (tArrayUnCoupleAndUnconstraint.Length != 0)
//      {
//        Kuu = MatrixTool.GetSubMatrix(kGlobal, tArrayUnCoupleAndUnconstraint);
//        Ru = MatrixTool.GetSubVector(residualForce, tArrayUnCoupleAndUnconstraint);
//        if (tArrayDependent.Length != 0)
//        {
//          Kus = MatrixTool.GetSubMatrix(kGlobal, tArrayUnCoupleAndUnconstraint, tArrayDependent);
//          Ru += -(MatrixFunctions.Product(Kus, gS));
//        }
//      }

//      DoubleMatrix Tsm = null;
//      if (tArrayDependent.Length != 0 && tArrayIndependent.Length != 0)
//      {
//        Tsm = MatrixTool.GetSubMatrix(model.TDense, tArrayDependent, tArrayIndependent);
//        DoubleMatrix TsmTranspose = Tsm.Transpose();
//        DoubleMatrix Kms = MatrixTool.GetSubMatrix(kGlobal, tArrayIndependent, tArrayDependent);
//        DoubleMatrix Kss = MatrixTool.GetSubMatrix(kGlobal, tArrayDependent, tArrayDependent);
//        DoubleVector Rs = MatrixTool.GetSubVector(residualForce, tArrayDependent);
//        Kmm = MatrixTool.GetSubMatrix(kGlobal, tArrayIndependent, tArrayIndependent) + MatrixFunctions.Product(TsmTranspose, Kms.Transpose() + MatrixFunctions.Product(Kss, Tsm)) + MatrixFunctions.Product(Kms, Tsm);
//        Rm = MatrixTool.GetSubVector(residualForce, tArrayIndependent) + MatrixFunctions.Product(TsmTranspose, Rs);
//        if (tArrayUnCoupleAndUnconstraint.Length != 0)
//          Kum = MatrixTool.GetSubMatrix(kGlobal, tArrayUnCoupleAndUnconstraint, tArrayIndependent) + MatrixFunctions.Product(Kus, Tsm);
//        Rm += -(NMathFunctions.Product(Kms, gS) + NMathFunctions.Product(TsmTranspose, NMathFunctions.Product(Kss, gS)));
//      }
//      DoubleMatrix lhs = new DoubleMatrix(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length, tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
//      DoubleVector rhs = new DoubleVector(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
//      if (!AbstractModel.IsParallelProcesing)
//      {
//        for (int i = 0; i < tArrayUnCoupleAndUnconstraint.Length; i++)
//        {
//          for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
//          {
//            lhs[i, j] = Kuu[i, j];
//          }
//          if (Kum != null)
//          {
//            for (int j = 0; j < tArrayIndependent.Length; j++)
//            {
//              lhs[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
//            }
//          }
//          rhs[i] = Ru[i];
//        }
//        for (int i = 0; i < tArrayIndependent.Length; i++)
//        {
//          for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
//          {
//            lhs[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
//          }

//          for (int j = 0; j < tArrayIndependent.Length; j++)
//          {
//            lhs[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
//          }
//          rhs[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
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
//                     lhs[i, j] = Kuu[i, j];
//                   }
//                   if (Kum != null)
//                   {
//                     for (int j = 0; j < tArrayIndependent.Length; j++)
//                     {
//                       lhs[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
//                     }
//                   }
//                   rhs[i] = Ru[i];
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
//                     lhs[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
//                   }

//                   for (int j = 0; j < tArrayIndependent.Length; j++)
//                   {
//                     lhs[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
//                   }
//                   rhs[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
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
//      DoubleVector duLocalS = gS;
//      if (rhs.Length != 0)
//      {
//        duLocalUM = MatrixFunctions.Solve(lhs, rhs);
//        duLocalU = MatrixTool.GetSubVector(duLocalUM, tArrayUncoupleUnconstraintInNew);
//        duLocalM = MatrixTool.GetSubVector(duLocalUM, tArrayIndependentInNew);
//        if (duLocalM != null)
//          duLocalS += MatrixFunctions.Product(Tsm, duLocalM);
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
//      DoubleVector uGlobal = model.GetUGlobal() + duGlobal;
//      model.SetUGlobal(uGlobal);

//      model.AssemblyInternalForce(out DoubleVector internalForceGlobal);
//      residualForce = Fstep - internalForceGlobal;

//      err1 = ErrorCalculationIteration(duGlobal, residualForce, internalForceGlobal, internalForceGlobal, tArray);
//      if (err1 == double.NaN)
//      {
//        this.bConverged = false;
//        //return null;
//      }
//      return uGlobal;
//    }

//    protected override DoubleVector SolveOneIterationSparse(DoubleCsrSparseMatrix kGlobalSparse, DoubleVector Fstep, DoubleVector ui, ref double err1, ref DoubleVector residualForce)
//    {
//      // Update Stiffness matrix and Internal force
//      if (bUpdate == true)
//      {
//        model.AssemblyStiffnessMatrix(ref kGlobalSparse);
//      }

//      DoubleVector residualForceClone = (DoubleVector)residualForce.Clone();

//      model.ApplyTAndg(ref kGlobalSparse, ref residualForceClone, null);
//      DoubleVector duGlobal = MatrixFunctions.Solve(kGlobalSparse, residualForceClone);

//      duGlobal = MatrixFunctions.Product(model.TSparse, duGlobal) + model.g;
//      DoubleVector uGlobal = model.GetUGlobal() + duGlobal;
//      model.SetUGlobal(uGlobal);

//      model.AssemblyInternalForce(out DoubleVector internalForceGlobal);
//      residualForce = Fstep - internalForceGlobal;

//      int[] tArray = null;
//      //switch(AbstractModel.TypeSchemeNonlinearSolver)
//      //{
//      //  case TypeSchemeNonlinearSolver.Monolithic:
//      //  case TypeSchemeNonlinearSolver.SingleStaggered:
//      //    break;
//      //  case TypeSchemeNonlinearSolver.Staggered:
//      //    break;
//      //}
//      switch (AbstractModel.TypeSchemeNonlinearSolver)
//      {
//        case TypeSchemeNonlinearSolver.Monolithic:
//        case TypeSchemeNonlinearSolver.SingleStaggered:
//          tArray = model.tArray[0];
//          for (int i = 1; i < model.tArray.Length; i++)
//            tArray = tArray.Union(model.tArray[i]).ToArray();
//          break;
//        case TypeSchemeNonlinearSolver.Staggered:
//          if (!IsSolvePhasefieldStaggered)
//          {
//            tArray = model.tArray[0];
//          }
//          else
//          {
//            tArray = model.tArray[1];
//          }
//          break;
//      }

//      err1 = ErrorCalculationIteration(duGlobal, residualForce, internalForceGlobal, internalForceGlobal, tArray);
//      if (err1 == double.NaN)
//      {
//        this.bConverged = false;
//        //return null;
//      }
//      return uGlobal;
//    }
//  }
//}

using System;
using System.Linq;
using CenterSpace.NMath.Core;
using System.Collections.Generic;
using System.Threading.Tasks;
using DEMSoft.Common;

namespace DEMSoft.IGA
{

  /// <summary>
  /// Implement the load controlled Newton-Raphson solver for non-linear problems
  /// </summary>
  internal class SolverLoadControlled : AbstractNonlinearSolver
  {
    /// <summary>
    /// Constructor class
    /// </summary>
    /// <param name="model"></param>
    /// <param name="type"></param>
    public SolverLoadControlled(AbstractModelStructure model, TypeNonlinearSolver type)
      : base(model, type)
    {
    }

    /// <summary>
    /// Solve
    /// </summary>
    protected override void Solve(TypeNonlinearSolver type)
    {
      int totalStep = ((model.IsRunFromInitial) || (checkpointData.IsFinished) || !conditionRestoreCheckpoint) ? 1 : checkpointData.CurrentTotalStep;//1;
      int step = ((model.IsRunFromInitial) || (checkpointData.IsFinished) || !conditionRestoreCheckpoint) ? 1 : checkpointData.CurrentStepLoad;//1;
      int totalStepSaved = ((model.IsRunFromInitial) || (checkpointData.IsFinished) || !conditionRestoreCheckpoint) ? 0 : checkpointData.CurrentTotalStepSave;//1;
      int numberOfTotalStep = 0;
      int numberOfStepLoad = model.numberOfSubstepLoad.Length;

      model.deltaT = new double[numberOfStepLoad];
      for (int i = 0; i < numberOfStepLoad; i++)
      {
        numberOfTotalStep += model.numberOfSubstepLoad[i];
        model.deltaT[i] = 1.0/*[s]*/ / model.numberOfSubstepLoad[i];
      }
      MonitorDataHistorism monitorDataHistorism = model.GetMonitorDataHistorism();
      int countData = monitorDataHistorism.CountData();
      model.Time = new List<double>();
      if (countData > 0)
        monitorDataHistorism.MonitorDataStore = new List<double[]>();
      DoubleVector un0 = new DoubleVector(model.CountDOF());
      DoubleVector un1 = null;
      if ((!model.IsRunFromInitial) && (!checkpointData.IsFinished) && conditionRestoreCheckpoint)
      {
        totalStep = checkpointData.CurrentTotalStep;
        step = checkpointData.CurrentStepLoad;
        totalStepSaved = checkpointData.CurrentTotalStepSave;
        model.Time = checkpointData.Time;
        if (countData > 0)
          monitorDataHistorism.MonitorDataStore = checkpointData.MonitorDataStore;
        un0 = checkpointData.un0;
        lastPrevF0 = checkpointData.LastPrevFi0;
        lastPrevF1 = checkpointData.LastPrevFi1;
        if (model.IsWriteLogFile)
        {
          model.logFile = IOIGA.ReadLogFile(model.PathProject + model.NameProject + ".log");
        }
      }
      if (monitorDataHistorism.listOfCurveInfo != null)
        if (monitorDataHistorism.listOfCurveInfo.Count > 0)
        {
          monitorDataHistorism.RunGraph(conditionRestoreCheckpoint);
        }

      bool isRunIntoLoopSubstep = false;
      while (step <= numberOfStepLoad)
      {
        if (model.IsWriteLogFile)
        {
          IOIGA.WriteLogFileAndConsole(model.logFile, true, "-----------------------------------------------------------");
          IOIGA.WriteLogFileAndConsole(model.logFile, true, "Load step " + step + " : ");
        }
        int numberOfSubstepLoad = model.numberOfSubstepLoad[step - 1];
        //Compute increment vector of field
        DoubleVector[] DNoneZeroVector = new DoubleVector[model.typeOfFieldsInMultifield.Count];
        tArrayConstraintNoneZero = new int[model.typeOfFieldsInMultifield.Count][];

        for (int i = 0; i < model.typeOfFieldsInMultifield.Count; i++)
        {
          ComputeDisplacementVector(step - 1, model.tArrayConstraint[i], out DNoneZeroVector[i], out tArrayConstraintNoneZero[i]);

          /////////////////////////////////////////////////////////////
          ///////Store only increment none zero///////
          /////////////////////////////////////////////////////////////
          if (DNoneZeroVector[i] != null)
          {
            if (DNoneZeroVector[i].Length > 0 && i == 0/*Structure*/)
              DNoneZeroVector[i] = DNoneZeroVector[i] / numberOfSubstepLoad;
          }
        }
        int substep = ((model.IsRunFromInitial) || (checkpointData.IsFinished) || !conditionRestoreCheckpoint) ? 1 : checkpointData.CurrentSubStep;//1;
        bool isSave = false;
        while (substep <= numberOfSubstepLoad)
        {
          isRunIntoLoopSubstep = true;
          isSave = false;
          if (model.IsWriteLogFile)
          {
            IOIGA.WriteLogFileAndConsole(model.logFile, true, "-----------------------------------------------------------");
            IOIGA.WriteLogFileAndConsole(model.logFile, true, "SubLoad step " + substep + " : ");
          }
          //////////////////////////////////////////////////////////////////
          // Apply inhomogenous boundary condition
          for (int i = 0; i < model.typeOfFieldsInMultifield.Count; i++)
          {
            if (tArrayConstraintNoneZero[i] != null)
              for (int ii = 0; ii < tArrayConstraintNoneZero[i].Length; ii++)
              {
                //Chua dung cho dkb phase=1
                un0[tArrayConstraintNoneZero[i][ii]] += DNoneZeroVector[i][ii];
              }
          }
          //if (model is ModelStructurePhaseFieldStatic)
          //  ((ModelStructurePhaseFieldStatic)model).ExtrapolationInitialPhasefieldGausspoint(ref un0);
          model.SetUGlobal(un0);

          //////////////////////////////////////////////////////////////////
          double time = step - 1 + substep * model.deltaT[step - 1];
          deltaT = model.deltaT[step - 1];
          model.Time.Add(time);
          DoubleVector fStep = null;
          model.AssemblyTractionVector(out fStep, time);

          DoubleVector residualForce = null;
          // Solve
          switch (AbstractModel.TypeSchemeNonlinearSolver)
          {
            case TypeSchemeNonlinearSolver.Monolithic:
            case TypeSchemeNonlinearSolver.SingleStaggered:
              un1 = SolveOneSubstepLoad(un0, fStep, out residualForce, time);
              if (!bConverged)
              {
                model.stopWatchWholeModel.Stop();
                IOIGA.WriteLogFileAndConsole(model.logFile, true, "-----------------------------------------------------------");
                IOIGA.WriteLogFileAndConsole(model.logFile, true, "Not converged !!");
                IOIGA.WriteLogFileAndConsole(model.logFile, true, "Elapsed time = " + model.stopWatchWholeModel.Elapsed);
                IOIGA.WriteLogFile(model.logFile, model.PathProject + model.NameProject + ".log");
                break;
              }
              UpdateConveredGaussPointsEachField();
              lastPrevF0 = prevF;
              break;
            case TypeSchemeNonlinearSolver.Staggered:

              // Solve phase field
              IsSolvePhasefieldStaggered = true;
              un1 = SolveOneSubstepLoad(un0, fStep, out DoubleVector residualForce1, time);
              if (!bConverged)
              {
                model.stopWatchWholeModel.Stop();
                IOIGA.WriteLogFileAndConsole(model.logFile, true, "-----------------------------------------------------------");
                IOIGA.WriteLogFileAndConsole(model.logFile, true, "Not converged !!");
                IOIGA.WriteLogFileAndConsole(model.logFile, true, "Elapsed time = " + model.stopWatchWholeModel.Elapsed);
                IOIGA.WriteLogFile(model.logFile, model.PathProject + model.NameProject + ".log");
                break;
              }
              UpdateConveredGaussPointsEachField();
              lastPrevF0 = prevF;

              // Solve displacement
              IsSolvePhasefieldStaggered = false;
              un1 = SolveOneSubstepLoad(un1, fStep, out residualForce, time);
              if (!bConverged)
              {
                model.stopWatchWholeModel.Stop();
                IOIGA.WriteLogFileAndConsole(model.logFile, true, "-----------------------------------------------------------");
                IOIGA.WriteLogFileAndConsole(model.logFile, true, "Not converged !!");
                IOIGA.WriteLogFileAndConsole(model.logFile, true, "Elapsed time = " + model.stopWatchWholeModel.Elapsed);
                IOIGA.WriteLogFile(model.logFile, model.PathProject + model.NameProject + ".log");
                break;
              }
              UpdateConveredGaussPointsEachField();
              lastPrevF1 = prevF;
              break;
          }

          ///////////////////////////////////////////////////////
          UpdateConvergedGaussPoints();
          ///////////////////////////////////////////////////////
          model.residualForce = residualForce;
          ComputeMonitorResults();
          if (monitorDataHistorism.IsGraphRealTime)
          {
            monitorDataHistorism.UpdateGraph();
          }
          if (model.IsWriteLogFile)
          {
            IOIGA.WriteLogFileAndConsole(model.logFile, true, "umax = " + un1.Max());
            IOIGA.WriteLogFileAndConsole(model.logFile, true, "F = " + fStep.TwoNorm());
            WriteOutputMonitorData(monitorDataHistorism);
            IOIGA.WriteLogFileAndConsole(model.logFile, true, "Elapsed time = " + model.stopWatchWholeModel.Elapsed);
          }

          ////////////////////////////////////////////////////////////////////////////

          un0 = un1;
          substep++;
          totalStep++;
          if ((substep - 1) < numberOfSubstepLoad && ((substep - 1) % model.NumberOfStepStorageNonlinear == 0))
          {
            StoreData(totalStep, step, ref totalStepSaved, un0, substep, ref isSave);
          }
        }
        if (!bConverged)
          break;
        step++;
        if ((step - 1) < numberOfStepLoad)
        {
          StoreData(totalStep, step, ref totalStepSaved, un0, 1, ref isSave);
        }
        if ((step - 1) == numberOfStepLoad && isRunIntoLoopSubstep == true)
        {
          StoreData(totalStep, step - 1, ref totalStepSaved, un0, substep, ref isSave);
        }
      }
      if (!model.IsSaveStepByStep)
      {
        IOIGA.WriteMessagePackDoubleVector(model.FileNameDataTime, un0);
        IOIGA.WriteMessagePackData(model.FileNameTime, model.Time);
        if (monitorDataHistorism.MonitorDataStore != null)
        {
          IOIGA.WriteXMLResult(monitorDataHistorism.MonitorDataStore, model.PathProject + "monitor_" + model.NameProject + ".xml");
        }
      }
      //model.IsFinished = true;
    }

    private void WriteOutputMonitorData(MonitorDataHistorism monitorDataHistorism)
    {
      if (monitorDataHistorism.IndexColOutputMonitor != null)
      {
        int numberOutputData = monitorDataHistorism.IndexColOutputMonitor.Length;
        for (int ii = 0; ii < numberOutputData; ii++)
        {
          int index = monitorDataHistorism.IndexColOutputMonitor[ii];
          MonitorDataNode node = monitorDataHistorism.GetMonitorDataNode(index);
          IOIGA.WriteLogFileAndConsole(model.logFile, true, node.ResultType.ToString() + " = " + monitorDataHistorism.MonitorDataStore[monitorDataHistorism.GetLastIndexMonitorDataStore()][index]);
        }
      }
    }

    private void StoreData(int totalStep, int step, ref int totalStepSaved, DoubleVector un0, int substep, ref bool isSave)
    {

      //store.IsFinished = IsFinished;
      if (model.GetMonitorDataHistorism().MonitorDataStore != null)
      {
        IOIGA.WriteXMLResult(model.GetMonitorDataHistorism().MonitorDataStore, model.PathProject + "monitor_" + model.NameProject + ".xml");
      }
      if (model.IsSaveStepByStep && !isSave)
      {
        totalStepSaved++;
        string fileNameDataTime = model.PathProject + "data_" + (totalStepSaved) + ".msgpack";
        //IOIGA.WriteListDoubleVectorMessagePackFormatter(fileNameDataTime, model.DisplacementTime);
        IOIGA.WriteMessagePackDoubleVector(fileNameDataTime, un0);
        IOIGA.WriteMessagePackData(model.FileNameTime, model.Time);
        isSave = true;
        IOIGA.WriteLogFile(model.logFile, model.PathProject + model.NameProject + ".log");
      }
      checkpointData.CurrentStepLoad = step;
      checkpointData.CurrentSubStep = substep;
      checkpointData.CurrentTotalStep = totalStep;
      checkpointData.CurrentTotalStepSave = totalStepSaved;
      checkpointData.un0 = un0;
      checkpointData.LastPrevFi0 = lastPrevF0;
      checkpointData.LastPrevFi1 = lastPrevF1;
      checkpointData.SaveData();
    }

    /// <summary>
    /// Solve for one load step. External load: FStep = lambda*F (F is the global load vector)
    /// </summary>
    /// <param name="un"></param>
    /// <param name="Fstep"></param>
    /// <returns></returns>
    private DoubleVector SolveOneSubstepLoad(DoubleVector un, DoubleVector Fstep, out DoubleVector residualForce, double time)
    {
      // Last converged values: un
      int nIter = 0;

      DoubleVector ui = un;

      int[] tArray = null;
      switch (AbstractModel.TypeSchemeNonlinearSolver)
      {
        case TypeSchemeNonlinearSolver.Monolithic:
        case TypeSchemeNonlinearSolver.SingleStaggered:
          tArray = model.tArray[0];
          for (int i = 1; i < model.tArray.Length; i++)
            tArray = tArray.Union(model.tArray[i]).ToArray();
          break;
        case TypeSchemeNonlinearSolver.Staggered:
          if (!IsSolvePhasefieldStaggered)
          {
            tArray = model.tArray[0];
          }
          else
          {
            tArray = model.tArray[1];
          }
          break;
      }


      DoubleVector internalForceGlobal = null;
      model.AssemblyInternalForce(out internalForceGlobal);
      residualForce = Fstep - internalForceGlobal;
      DoubleCsrSparseMatrix kGlobalSparse = null;
      DoubleMatrix kGlobal = null;
      if (!model.IsSparseData)
      {
        // Update Stiffness matrix and Internal force
        if (bUpdate == false)
        {
          model.AssemblyStiffnessMatrix(out kGlobal);
        }
      }
      else
      {
        kGlobalSparse = (DoubleCsrSparseMatrix)model.kGlobalSparse.Clone();
        if (bUpdate == false)
        {
          model.AssemblyStiffnessMatrix(ref kGlobalSparse);
        }
      }

      double[] err1 = new double[model.tArray.Length];
      for (int i = 0; i < err1.Length; i++)
      {
        err1[i] = 1;
      }
      double err = 1;
      double wnorm = 1;
      double wnorm1 = 1;
      double wnorm2 = 1;
      prevF = 0;
      //err1 = ErrorCalculationIteration(new DoubleVector(2), residualForce, internalForceGlobal, tArray);
      while ((err > model.convergenceOption.TOL) && (nIter < model.convergenceOption.maximumIteration))
      {
        nIter++;
        // Solve for each iteration
        DoubleVector ui1 = null;

        if (!model.IsSparseData)
        {
          ui1 = SolveOneIterationDense(kGlobal, Fstep, ui, ref err1, ref residualForce);
        }
        else
        {
          ui1 = SolveOneIterationSparse(kGlobalSparse, Fstep, ui, ref err1, ref residualForce);
        }
        if (nIter == 1)
        {
          var du = ui1 - ui;
          DoubleVector reduceDu = null;
          DoubleVector reduceResidualForce = null;
          switch (model.convergenceOption.criterionConvergence)
          {
            case CriterionConvergence.DisplacementCriterion:
              reduceDu = MatrixTool.GetSubVector(du, tArray, AbstractModel.NumberOfCPUs);
              wnorm = Math.Sqrt(MatrixFunctions.Dot(reduceDu, reduceDu));
              break;
            case CriterionConvergence.ForceCriterion:
              //DoubleVector reduceResidualForce1 = null;
              //  reduceResidualForce1 = MatrixTool.GetSubVector(residualForce, tArray);
              //double aa = reduceResidualForce1.Sum();


              //DoubleVector reduceResidualForce = null;
              //if (AbstractModel.IsMonolithic)
              //{
              //  int[] noConstraintTArray = tArray.Except(model.tArrayConstraint[0]).Except(model.tArrayConstraint[1]).ToArray();
              //  reduceResidualForce = MatrixTool.GetSubVector(residualForce, noConstraintTArray);
              //}
              //else
              //{
              //  if (!IsSolveSecondStaggered)
              //  {
              //    int[] noConstraintTArray = tArray.Except(model.tArrayConstraint[0]).ToArray();
              //    reduceResidualForce = MatrixTool.GetSubVector(residualForce, noConstraintTArray);
              //  }
              //  else
              //  {
              //    int[] noConstraintTArray = tArray.Except(model.tArrayConstraint[1]).ToArray();
              //    reduceResidualForce = MatrixTool.GetSubVector(residualForce, noConstraintTArray);
              //  }
              //}
              //wnorm = Math.Sqrt(MatrixFunctions.Dot(reduceResidualForce, reduceResidualForce));
              wnorm = 1;
              break;
            case CriterionConvergence.EnergyCriterion:
              reduceDu = MatrixTool.GetSubVector(du, tArray, AbstractModel.NumberOfCPUs);
              wnorm = Math.Abs(MatrixFunctions.Dot(reduceDu, MatrixTool.GetSubVector(internalForceGlobal - residualForce, tArray, AbstractModel.NumberOfCPUs)));
              break;
            case CriterionConvergence.ResidualCriterion:

              double fn;
              if ((AbstractModel.TypeSchemeNonlinearSolver == TypeSchemeNonlinearSolver.Monolithic) || (AbstractModel.TypeSchemeNonlinearSolver == TypeSchemeNonlinearSolver.SingleStaggered) || (IsSolvePhasefieldStaggered))
              {
                fn = lastPrevF0;
              }
              else
              {
                fn = lastPrevF1;
              }
              //var fn = (AbstractModel.IsMonolithic || IsSolvePhasefieldStaggered) ? lastPrevF0 : lastPrevF1;
              reduceResidualForce = MatrixTool.GetSubVector(residualForce, tArray, AbstractModel.NumberOfCPUs);
              wnorm = Math.Abs(Math.Sqrt(MatrixFunctions.Dot(reduceResidualForce, reduceResidualForce)) - fn);//Math.Abs(Math.Sqrt(MatrixFunctions.Dot(f, f)) - fn);
              break;
            case CriterionConvergence.DisplacementCriterion2:
              reduceDu = MatrixTool.GetSubVector(du, tArray, AbstractModel.NumberOfCPUs);
              wnorm = 1 + MatrixFunctions.Dot(reduceDu, reduceDu);
              break;
            case CriterionConvergence.DisplacementCriterion22:
              reduceDu = MatrixTool.GetSubVector(du, tArray, AbstractModel.NumberOfCPUs);
              wnorm = 1 + Math.Sqrt(MatrixFunctions.Dot(reduceDu, reduceDu));
              break;
            case CriterionConvergence.ResidualCriterion2:
              reduceResidualForce = MatrixTool.GetSubVector(residualForce, tArray, AbstractModel.NumberOfCPUs);
              wnorm = 1 + MatrixFunctions.Dot(reduceResidualForce, reduceResidualForce);
              break;
            case CriterionConvergence.DisplacementCriterion3:
              DoubleVector reduceUn = MatrixTool.GetSubVector(un, tArray, AbstractModel.NumberOfCPUs);
              wnorm = Math.Sqrt(MatrixFunctions.Dot(reduceUn, reduceUn));
              if (wnorm == 0)
                wnorm = 1;
              break;
            case CriterionConvergence.DisplacementCriterion5:
              wnorm = 1;
              break;
            case CriterionConvergence.DisplacementCriterion4:
              int[] noConstraintTArray = null;
              switch (AbstractModel.TypeSchemeNonlinearSolver)
              {
                case TypeSchemeNonlinearSolver.Monolithic:
                case TypeSchemeNonlinearSolver.SingleStaggered:
                  int[] noConstraintTArray0 = model.tArray[0].Except(model.tArrayConstraint[0]).ToArray();
                  int[] noConstraintTArray1 = model.tArray[1].Except(model.tArrayConstraint[1]).ToArray();
                  DoubleVector reducedu0 = MatrixTool.GetSubVector(un, noConstraintTArray0, AbstractModel.NumberOfCPUs);
                  DoubleVector reducedu1 = MatrixTool.GetSubVector(un, noConstraintTArray1, AbstractModel.NumberOfCPUs);
                  wnorm1 = Math.Sqrt(MatrixFunctions.Dot(reducedu0, reducedu0));
                  wnorm2 = Math.Sqrt(MatrixFunctions.Dot(reducedu1, reducedu1));
                  //wnorm = Math.Max(Math.Sqrt(MatrixFunctions.Dot(reducedu0, reducedu0)), Math.Sqrt(MatrixFunctions.Dot(reducedu1, reducedu1)));
                  break;
                case TypeSchemeNonlinearSolver.Staggered:
                  if (!IsSolvePhasefieldStaggered)
                  {
                    noConstraintTArray = tArray.Except(model.tArrayConstraint[0]).ToArray();
                    reduceDu = MatrixTool.GetSubVector(un, noConstraintTArray, AbstractModel.NumberOfCPUs);
                  }
                  else
                  {
                    noConstraintTArray = tArray.Except(model.tArrayConstraint[1]).ToArray();
                    reduceDu = MatrixTool.GetSubVector(un, noConstraintTArray, AbstractModel.NumberOfCPUs);
                  }
                  wnorm = Math.Sqrt(MatrixFunctions.Dot(reduceDu, reduceDu));
                  break;
              }
              if (wnorm == 0)
                wnorm = 1;
              if (wnorm1 == 0)
                wnorm1 = 1;
              if (wnorm2 == 0)
                wnorm2 = 1;
              break;
            case CriterionConvergence.ResidualCriterion3:
              ////int[] noConstraintTArray = null;
              //switch (AbstractModel.TypeSchemeNonlinearSolver)
              //{
              //  case TypeSchemeNonlinearSolver.Monolithic:
              //  case TypeSchemeNonlinearSolver.SingleStaggered:
              //    //noConstraintTArray = tArray.Except(model.tArrayConstraint[0]).Except(model.tArrayConstraint[1]).ToArray();
              //    //noConstraintTArray = tArray.Except(model.tArrayConstraint[0]).ToArray();
              //    int[] noConstraintTArray0 = model.tArray[0].Except(model.tArrayConstraint[0]).ToArray();
              //    int[] noConstraintTArray1 = model.tArray[1].Except(model.tArrayConstraint[1]).ToArray();
              //    DoubleVector reduceResidualForce0 = MatrixTool.GetSubVector(residualForce, noConstraintTArray0);
              //    DoubleVector reduceResidualForce1 = MatrixTool.GetSubVector(residualForce, noConstraintTArray1);
              //    wnorm1 = Math.Sqrt(MatrixFunctions.Dot(reduceResidualForce0, reduceResidualForce0));
              //    wnorm2 = Math.Sqrt(MatrixFunctions.Dot(reduceResidualForce1, reduceResidualForce1));
              //    break;
              //  case TypeSchemeNonlinearSolver.Staggered:
              //    if (!IsSolvePhasefieldStaggered)
              //    {
              //      noConstraintTArray = tArray.Except(model.tArrayConstraint[0]).ToArray();
              //      reduceResidualForce = MatrixTool.GetSubVector(residualForce, noConstraintTArray);
              //    }
              //    else
              //    {
              //      noConstraintTArray = tArray.Except(model.tArrayConstraint[1]).ToArray();
              //      reduceResidualForce = MatrixTool.GetSubVector(residualForce, noConstraintTArray);
              //    }
              //    wnorm = Math.Sqrt(MatrixFunctions.Dot(reduceResidualForce, reduceResidualForce));
              //    break;
              //}
              wnorm1 = wnorm2 = wnorm = 1;
              break;
          }
        }

        if (model.convergenceOption.criterionConvergence != CriterionConvergence.DisplacementCriterion4
          || model.convergenceOption.criterionConvergence != CriterionConvergence.ResidualCriterion3)
        {
          err = err1[0];
          if (wnorm != 0 && err1[0] != 0)
          {
            err = err / wnorm;
          }
        }
        else
        {
          err1[0] = err1[0] / wnorm1;
          err1[1] = err1[1] / wnorm2;
          err = err1.Max();
        }
        if (model.IsWriteLogFile)
        {
          IOIGA.WriteLogFileAndConsole(model.logFile, true, "\terror = " + err);
        }
        ui = ui1;
      }

      if (model.IsWriteLogFile)
      {
        IOIGA.WriteLogFileAndConsole(model.logFile, true, "\tNumber of iteration = " + nIter);
      }
      if (nIter == model.convergenceOption.maximumIteration && !model.convergenceOption.IsSkipUnconvergenceStep)
      {
        this.bConverged = false;
      }
      return ui;
    }

    /// <summary>
    /// Solve for one iteration.
    /// </summary>
    /// <param name="Fstep"></param>
    /// <param name="ui"></param>
    /// <param name="err1"></param>
    /// <param name="bUpdate"></param>
    /// <returns></returns>
    protected override DoubleVector SolveOneIterationDense(DoubleMatrix kGlobal, DoubleVector Fstep, DoubleVector ui, ref double[] err1, ref DoubleVector residualForce)
    {
      if (bUpdate == true)
      {
        model.AssemblyStiffnessMatrix(out kGlobal);
      }
      DoubleVector internalForceGlobali = Fstep - residualForce;
      // Compute length of Unconstraint DOF
      int[] tArray = null;
      switch (AbstractModel.TypeSchemeNonlinearSolver)
      {
        case TypeSchemeNonlinearSolver.Monolithic:
        case TypeSchemeNonlinearSolver.SingleStaggered:
          tArray = model.tArray[0];
          for (int i = 1; i < model.tArray.Length; i++)
            tArray = tArray.Union(model.tArray[i]).ToArray();
          break;
        case TypeSchemeNonlinearSolver.Staggered:
          if (!IsSolvePhasefieldStaggered)
          {
            tArray = model.tArray[0];
          }
          else
          {
            tArray = model.tArray[1];
          }
          break;
      }

      int[] tArrayIndependent = tArray.Intersect(model.tArrayIndependent).ToArray();
      int[] tArrayDependent = tArray.Intersect(model.tArrayDependent).ToArray();
      int[] tArrayUnCoupleAndUnconstraint = tArray.Intersect(model.tArrayUncouple).ToArray();
      int[] tArrayIndependentInNew = new int[tArrayIndependent.Length];
      int c = 0;
      for (int i = tArrayUnCoupleAndUnconstraint.Length; i < tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length; i++)
      {
        tArrayIndependentInNew[c++] = i;
      }

      int[] tArrayUncoupleUnconstraintInNew = new int[tArrayUnCoupleAndUnconstraint.Length];
      for (int i = 0; i < tArrayUnCoupleAndUnconstraint.Length; i++)
      {
        tArrayUncoupleUnconstraintInNew[i] = i;
      }
      DoubleVector gS = null;
      if (tArrayDependent.Length != 0)
      {
        gS = MatrixTool.GetSubVector(model.g, tArrayDependent, AbstractModel.NumberOfCPUs);
      }
      DoubleVector Ru = null;
      DoubleVector Rm = null;

      DoubleMatrix Kuu = null;
      DoubleMatrix Kus = null;
      DoubleMatrix Kum = null;
      DoubleMatrix Kmm = null;
      if (tArrayUnCoupleAndUnconstraint.Length != 0)
      {
        Kuu = MatrixTool.GetSubMatrix(kGlobal, tArrayUnCoupleAndUnconstraint, tArrayUnCoupleAndUnconstraint, AbstractModel.NumberOfCPUs);
        Ru = MatrixTool.GetSubVector(residualForce, tArrayUnCoupleAndUnconstraint, AbstractModel.NumberOfCPUs);
        if (tArrayDependent.Length != 0)
        {
          Kus = MatrixTool.GetSubMatrix(kGlobal, tArrayUnCoupleAndUnconstraint, tArrayDependent, AbstractModel.NumberOfCPUs);
          Ru += -(NMathFunctions.Product(Kus, gS));
        }
      }

      DoubleMatrix Tsm = null;
      if (tArrayDependent.Length != 0 && tArrayIndependent.Length != 0)
      {
        Tsm = MatrixTool.GetSubMatrix(model.TDense, tArrayDependent, tArrayIndependent, AbstractModel.NumberOfCPUs);
        DoubleMatrix TsmTranspose = Tsm.Transpose();
        DoubleMatrix Kms = MatrixTool.GetSubMatrix(kGlobal, tArrayIndependent, tArrayDependent, AbstractModel.NumberOfCPUs);
        DoubleMatrix Kss = MatrixTool.GetSubMatrix(kGlobal, tArrayDependent, tArrayDependent, AbstractModel.NumberOfCPUs);
        DoubleVector Rs = MatrixTool.GetSubVector(residualForce, tArrayDependent, AbstractModel.NumberOfCPUs);
        Kmm = MatrixTool.GetSubMatrix(kGlobal, tArrayIndependent, tArrayIndependent, AbstractModel.NumberOfCPUs) + NMathFunctions.Product(TsmTranspose, Kms.Transpose() + NMathFunctions.Product(Kss, Tsm)) + NMathFunctions.Product(Kms, Tsm);
        Rm = MatrixTool.GetSubVector(residualForce, tArrayIndependent, AbstractModel.NumberOfCPUs) + NMathFunctions.Product(TsmTranspose, Rs);
        if (tArrayUnCoupleAndUnconstraint.Length != 0)
          Kum = MatrixTool.GetSubMatrix(kGlobal, tArrayUnCoupleAndUnconstraint, tArrayIndependent, AbstractModel.NumberOfCPUs) + NMathFunctions.Product(Kus, Tsm);
        Rm += -(NMathFunctions.Product(Kms, gS) + NMathFunctions.Product(TsmTranspose, NMathFunctions.Product(Kss, gS)));
      }
      DoubleSymmetricMatrix lhs = new DoubleSymmetricMatrix(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
      DoubleVector rhs = new DoubleVector(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int i = 0; i < tArrayUnCoupleAndUnconstraint.Length; i++)
        {
          for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
          {
            lhs[i, j] = Kuu[i, j];
          }
          if (Kum != null)
          {
            for (int j = 0; j < tArrayIndependent.Length; j++)
            {
              lhs[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
            }
          }
          rhs[i] = Ru[i];
        }
        for (int i = 0; i < tArrayIndependent.Length; i++)
        {
          for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
          {
            lhs[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
          }

          for (int j = 0; j < tArrayIndependent.Length; j++)
          {
            lhs[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
          }
          rhs[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
        }
      }
      else
      {
        var degreeOfParallelism = AbstractModel.NumberOfCPUs;
        var tasks = new Task[degreeOfParallelism];
        for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
        {
          // capturing taskNumber in lambda wouldn't work correctly
          int taskNumberCopy = taskNumber;

          tasks[taskNumber] = Task.Factory.StartNew(
               () =>
               {
                 var max = tArrayUnCoupleAndUnconstraint.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                 for (int i = tArrayUnCoupleAndUnconstraint.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                 {
                   for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
                   {
                     lhs[i, j] = Kuu[i, j];
                   }
                   if (Kum != null)
                   {
                     for (int j = 0; j < tArrayIndependent.Length; j++)
                     {
                       lhs[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
                     }
                   }
                   rhs[i] = Ru[i];
                 }
               });
        }
        Task.WaitAll(tasks);

        for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
        {
          // capturing taskNumber in lambda wouldn't work correctly
          int taskNumberCopy = taskNumber;

          tasks[taskNumber] = Task.Factory.StartNew(
               () =>
               {
                 var max = tArrayIndependent.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                 for (int i = tArrayIndependent.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                 {
                   for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
                   {
                     lhs[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
                   }

                   for (int j = 0; j < tArrayIndependent.Length; j++)
                   {
                     lhs[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
                   }
                   rhs[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
                 }
               });
        }
        Task.WaitAll(tasks);
      }

      DoubleVector duLocalUM = null;
      DoubleVector duLocalU = null;
      DoubleVector duLocalM = null;
      DoubleVector duLocalS = gS;
      if (rhs.Length != 0)
      {
        duLocalUM = MatrixFunctions.Solve(lhs, rhs);
        duLocalU = MatrixTool.GetSubVector(duLocalUM, tArrayUncoupleUnconstraintInNew, AbstractModel.NumberOfCPUs);
        duLocalM = MatrixTool.GetSubVector(duLocalUM, tArrayIndependentInNew, AbstractModel.NumberOfCPUs);
        if (duLocalM != null)
          duLocalS += MatrixFunctions.Product(Tsm, duLocalM);
      }

      DoubleVector duGlobal = new DoubleVector(model.CountDOF());
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int i = 0; i < tArrayUncoupleUnconstraintInNew.Length; i++)
        {
          duGlobal[tArrayUnCoupleAndUnconstraint[i]] = duLocalU[i];
        }
        for (int i = 0; i < tArrayIndependentInNew.Length; i++)
        {
          duGlobal[tArrayIndependent[i]] = duLocalM[i];
        }
        for (int i = 0; i < tArrayDependent.Length; i++)
        {
          duGlobal[tArrayDependent[i]] = duLocalS[i];
        }
      }
      else
      {
        var degreeOfParallelism = AbstractModel.NumberOfCPUs;
        var tasks = new Task[degreeOfParallelism];
        for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
        {
          // capturing taskNumber in lambda wouldn't work correctly
          int taskNumberCopy = taskNumber;

          tasks[taskNumber] = Task.Factory.StartNew(
               () =>
               {
                 var max = tArrayUncoupleUnconstraintInNew.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                 for (int i = tArrayUncoupleUnconstraintInNew.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                 {
                   duGlobal[tArrayUnCoupleAndUnconstraint[i]] = duLocalU[i];
                 }
               });
        }
        Task.WaitAll(tasks);
        for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
        {
          // capturing taskNumber in lambda wouldn't work correctly
          int taskNumberCopy = taskNumber;

          tasks[taskNumber] = Task.Factory.StartNew(
               () =>
               {
                 var max = tArrayIndependentInNew.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                 for (int i = tArrayIndependentInNew.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                 {
                   duGlobal[tArrayIndependent[i]] = duLocalM[i];
                 }
               });
        }
        Task.WaitAll(tasks);

        for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
        {
          // capturing taskNumber in lambda wouldn't work correctly
          int taskNumberCopy = taskNumber;

          tasks[taskNumber] = Task.Factory.StartNew(
               () =>
               {
                 var max = tArrayDependent.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                 for (int i = tArrayDependent.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                 {
                   duGlobal[tArrayDependent[i]] = duLocalS[i];
                 }
               });
        }
        Task.WaitAll(tasks);
      }
      DoubleVector uGlobal = model.GetUGlobal() + duGlobal;
      model.SetUGlobal(uGlobal);

      model.AssemblyInternalForce(out DoubleVector internalForceGlobal);
      residualForce = Fstep - internalForceGlobal;

      err1 = ErrorCalculationIteration(duGlobal, residualForce, internalForceGlobal, internalForceGlobali, tArray);
      for (int i = 0; i < err1.Length; i++)
      {
        if (err1[i] == double.NaN)
        {
          this.bConverged = false;
          //return null;
        }
      }
      return uGlobal;
    }

    protected override DoubleVector SolveOneIterationSparse(DoubleCsrSparseMatrix kGlobalSparse, DoubleVector Fstep, DoubleVector ui, ref double[] err1, ref DoubleVector residualForce)
    {
      // Update Stiffness matrix and Internal force
      if (bUpdate == true)
      {
        model.AssemblyStiffnessMatrix(ref kGlobalSparse);
      }
      DoubleVector internalForceGlobali = Fstep - residualForce;
      DoubleVector rhs = residualForce;//(DoubleVector)residualForce.Clone();
      DoubleCsrSparseMatrix lhs = (DoubleCsrSparseMatrix)kGlobalSparse.Clone();
      model.ApplyTAndg(ref lhs, ref rhs, null);
      DoubleVector duGlobal = MatrixFunctions.Solve(lhs, rhs);

      duGlobal = MatrixFunctions.Product(model.TSparse, duGlobal) + model.g;
      DoubleVector uGlobal = model.GetUGlobal() + duGlobal;
      model.SetUGlobal(uGlobal);

      model.AssemblyInternalForce(out DoubleVector internalForceGlobal);
      residualForce = Fstep - internalForceGlobal;

      int[] tArray = null;
      switch (AbstractModel.TypeSchemeNonlinearSolver)
      {
        case TypeSchemeNonlinearSolver.Monolithic:
        case TypeSchemeNonlinearSolver.SingleStaggered:
          tArray = model.tArray[0];
          for (int i = 1; i < model.tArray.Length; i++)
            tArray = tArray.Union(model.tArray[i]).ToArray();
          break;
        case TypeSchemeNonlinearSolver.Staggered:
          if (!IsSolvePhasefieldStaggered)
          {
            tArray = model.tArray[0];
          }
          else
          {
            tArray = model.tArray[1];
          }
          break;
      }

      err1 = ErrorCalculationIteration(duGlobal, residualForce, internalForceGlobali, internalForceGlobal, tArray);
      for (int i = 0; i < err1.Length; i++)
      {
        if (err1[i] == double.NaN)
        {
          this.bConverged = false;
          //return null;
        }
      }
      return uGlobal;
    }
  }
}
