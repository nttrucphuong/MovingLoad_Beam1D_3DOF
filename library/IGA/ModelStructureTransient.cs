using System;
using CenterSpace.NMath.Core;
using System.Collections.Generic;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Structure 2D problem (plane stress, plane strain). Default is plane stress
  /// </summary>
  public class ModelStructureTransient : AbstractModelStructure
  {
    /// <summary>
    /// Constructor class
    /// </summary>
    public ModelStructureTransient(Dimension structureDimension, string pathProject = "temp", string nameProject = "project")
                 : base(TypeAnalysisModel.Transient, structureDimension, pathProject, nameProject)
    {
      //listInitialConstraint = new System.Collections.Generic.List<AbstractConstraintValue>();
    }

    /// <summary>
    /// Solve problem
    /// </summary>
    public override void Solve()
    {
      if (listPatch.Count == 0)
        throw new NullReferenceException("Model must be assigned Patch");

      WriteInformationProblem("Transient structural module");
      //DoubleMatrix kGlobal = null;
      //DoubleMatrix mGlobal = null;
      DoubleVector uGlobal = null;
      DoubleVector rGlobal = null;

      //////////////////////////////////////////////////////
      //double delta = 0.5;
      //double alpha = 0.25;

      //AssemblyStiffnessMatrix(ref kGlobal);
      //AssemblyMassMatrix(ref mGlobal);

      //// kNew = K + b1*M + b4* C
      //DoubleMatrix kNew = mGlobal + alpha * kGlobal;//+ a1 * cGlobal));
      //double[,] uuu = new double[countDOF, 3];// Store u, du, ddu
      //displacementTime = DoubleMatrix(countDOF, numberOfStep);
      //DoubleVector time = DoubleVector(numberOfStepSave);
      //DoubleVector q1 = null;
      //int stepSave = 0;
      //for (int step = 1; step <= numberOfStep; step++)
      //{
      //    AssemblyTractionVector(ref rGlobal, step * deltaT);
      //    DoubleVector q0 = DoubleVector(countDOF);
      //    DoubleVector dq0 = DoubleVector(countDOF);
      //    DoubleVector ddq0 = DoubleVector(countDOF);
      //    DoubleVector ddq1 = DoubleVector(countDOF);
      //    //DoubleVector cRight = new DoubleVector(countDOF);
      //    for (int i = 0; i < countDOF; i++)
      //    {
      //        q0[i] = uuu[i, 0];
      //        dq0[i] = uuu[i, 1];
      //        ddq0[i] = uuu[i, 2];
      //    }
      //    rGlobal -= kGlobal * q0 + kGlobal * deltaT * dq0 + (1.0 / 2.0 - alpha) * kGlobal * ddq0 * deltaT * deltaT;
      //    ApplyBoundaryCondition(ref kNew, ref rGlobal, step * deltaT);
      //    solver = new Solver(kNew, rGlobal);
      //    if (IsSparseData)
      //        ddq1 = solver.Solve(SolverType.TFQMR, TOL, maxIter);
      //    else
      //        ddq1 = solver.Solve(SolverType.QR);
      //    // Compute velocity, acceleration
      //    for (int i = 0; i < countDOF; i++)
      //    {
      //        double du = dq0[i] + (1 - delta) * dq0[i] * deltaT + delta * ddq1[i] * deltaT;
      //        double u = q0[i] + dq0[i] * deltaT + (1.0 / 2.0 - alpha) * ddq0[i] * deltaT * deltaT + alpha * ddq1[i] * deltaT * deltaT;
      //        uuu[i, 0] = u;
      //        uuu[i, 1] = du;
      //        uuu[i, 2] = ddq1[i];
      //    }



      //double gamma = 0.5;
      //double beta = 0.25;
      //double b1 = 1.0 / (beta * deltaT * deltaT);
      //double b2 = 1.0 / (beta * deltaT);
      //double b3 = 1.0 / (2.0 * beta) - 1;// beta - 1.0 / 2.0;
      //double b4 = gamma * deltaT * b1;
      //double b5 = 1 + gamma * deltaT * b2;
      //double b6 = deltaT * (1.0 - gamma);
      //double a7 = deltaT * gamma;
      //AssemblyStiffnessMatrix(ref kGlobal);
      //AssemblyMassMatrix(ref mGlobal);

      //// kNew = K + b1*M + b4* C
      //DoubleMatrix kNew = (kGlobal + b1 * mGlobal);//+ a1 * cGlobal));
      //DoubleMatrix uuu = ApplyInitialBoundaryCondition();// new double[countDOF, 3];// Store u, du, ddu
      //displacementTime = DoubleMatrix(countDOF, numberOfStep);
      //DoubleVector time = DoubleVector(numberOfStepSave);
      //int stepSave = 0;
      //for (int step = 1; step <= numberOfStep; step++)
      //{
      //    AssemblyTractionVector(ref rGlobal, step * deltaT);
      //    DoubleVector mRight = DoubleVector(countDOF);
      //    //DoubleVector cRight = new DoubleVector(countDOF);
      //    for (int i = 0; i < countDOF; i++)
      //    {
      //        mRight[i] = b1 * uuu[i, 0] + b2 * uuu[i, 1] + b3 * uuu[i, 2];
      //        //cRight[i] = a1 * uuu[i, 0] + a4 * uuu[i, 1] + a5 * uuu[i, 2];
      //    }
      //    rGlobal += mGlobal * mRight;// +((DoubleVector)(cGlobal * cRight));
      //    DoubleMatrix kNewClone = 1.0 * kNew;
      //    ApplyBoundaryCondition(ref kNewClone, ref rGlobal, step * deltaT);
      //    solver = new Solver(kNewClone, rGlobal);
      //    if (IsSparseData)
      //        uGlobal = solver.Solve(SolverType.TFQMR, TOL, maxIter);
      //    else
      //        uGlobal = solver.Solve(SolverType.QR);
      //    // Compute velocity, acceleration
      //    for (int i = 0; i < countDOF; i++)
      //    {
      //        double ddu = b1 * (uGlobal[i] - uuu[i, 0]) - b2 * uuu[i, 1] - b3 * uuu[i, 2];
      //        double du = uuu[i, 1] + b6 * uuu[i, 2] + a7 * ddu;//uuu[i, 1] + a6 * uuu[i, 1] + a7 * uuu[i, 2];
      //        uuu[i, 0] = uGlobal[i];
      //        uuu[i, 1] = du;
      //        uuu[i, 2] = ddu;
      //    }

      ////ANSYS
      double delta = 0.5;
      double alpha = 0.25;

      DoubleMatrix kGlobalDense = null; // stiffness
      DoubleMatrix mGlobalDense = null; // "mass"
      // Assembly
      if (!IsSparseData)
      {
        AssemblyStiffnessMatrix(out kGlobalDense);
        AssemblyMassMatrix(out mGlobalDense);
      }
      else
      {
        AssemblyStiffnessMatrix(ref kGlobalSparse);
        AssemblyMassMatrix(ref mGlobalSparse);
      }

      int stepSave = 0;
      int numberOfStepSaveLocal = 0;
      for (int step = 1; step <= deltaT.Length; step++)
      {
        numberOfStepSaveLocal += numberOfStepSave[step - 1];
      }
      //DisplacementTime = new List<DoubleVector>();// (countDOF, numberOfStepSaveLocal);
      Time = new List<double>();//new DoubleVector(numberOfStepSaveLocal);
      int countData = GetMonitorDataHistorism().CountData();
      if (countData > 0)
        GetMonitorDataHistorism().MonitorDataStore = new List<double[]>();//new DoubleMatrix(numberOfStepSaveLocal, countData);
      double time = 0;
      DoubleMatrix uuu = ApplyInitialBoundaryCondition();//DoubleMatrix(countDOF, 3);// Store u, du, ddu

      for (int step = 1; step <= deltaT.Length; step++)
      {
        double deltaTLocal = deltaT[step - 1];

        if (IsWriteLogFile)
        {
          IOIGA.WriteLogFileAndConsole(logFile, true, "-----------------------------------------------------------");
          IOIGA.WriteLogFileAndConsole(logFile, true, "Load step " + step + " : delta=" + deltaTLocal);
        }
        double a0 = 1.0 / (alpha * deltaTLocal * deltaTLocal);
        double a1 = delta / (alpha * deltaTLocal);
        double a2 = 1.0 / (alpha * deltaTLocal);
        double a3 = 1.0 / (2.0 * alpha) - 1.0;
        double a4 = delta / alpha - 1.0;
        double a5 = deltaTLocal / 2.0 * (delta / alpha - 2.0);
        double a6 = deltaTLocal * (1 - delta);
        double a7 = delta * deltaTLocal;
        // kNew = K + b1*M + b4* C
        DoubleMatrix kNewDense = null; // M/deltaT + K;
        DoubleCsrSparseMatrix kNewSparse = null;
        if (!IsSparseData)
        {
          kNewDense = kGlobalDense + a0 * mGlobalDense;
        }
        else
        {
          kNewSparse = kGlobalSparse + a0 * mGlobalSparse;
        }

        //DoubleMatrix kNew = (kGlobal + a0 * mGlobal);//+ a1 * cGlobal));
        int numberOfSubstepLoadLocal = numberOfSubstepLoad[step - 1];
        for (int substep = 1; substep <= numberOfSubstepLoadLocal; substep++)
        {
          time += deltaTLocal;
          if (IsWriteLogFile)
          {
            IOIGA.WriteLogFileAndConsole(logFile, true, "-----------------------------------------------------------");
            IOIGA.WriteLogFileAndConsole(logFile, true, "Sub step " + substep + " : t=" + time);
          }

          AssemblyTractionVector(out rGlobal, time);
          DoubleVector mRight = a0 * uuu.Col(0) + a2 * uuu.Col(1) + a3 * uuu.Col(2);// DoubleVector(countDOF);
                                                                                    //DoubleVector cRight = a1 * uuu.Column(0) + a4 * uuu.Column(1) + a5 * uuu.Column(2); //DoubleVector(countDOF);
          DoubleMatrix kNewCloneDense = null;
          DoubleCsrSparseMatrix kNewCloneSparse = null;
          if (!IsSparseData)
          {
            rGlobal += MatrixFunctions.Product(mGlobalDense, mRight);// + cGlobal * cRight;
            kNewCloneDense = (DoubleMatrix)kNewDense.Clone();
          }
          else
          {
            rGlobal += MatrixFunctions.Product(mGlobalSparse, mRight);// + cGlobal * cRight;
            kNewCloneSparse = (DoubleCsrSparseMatrix)kNewSparse.Clone();
          }

          ApplyChangeBoundaryConditionAndUpdateIMatrix(time);
          uGlobal = SolveStatic(kNewCloneDense, kNewCloneSparse, ref rGlobal);
          SetUGlobal(uGlobal);
          // Compute velocity, acceleration
          for (int i = 0; i < countDOF; i++)
          {
            double du = a1 * (uGlobal[i] - uuu[i, 0]) - a4 * uuu[i, 1] - a5 * uuu[i, 2];
            double ddu = a0 * (uGlobal[i] - uuu[i, 0]) - a2 * uuu[i, 1] - a3 * uuu[i, 2];
            uuu[i, 0] = uGlobal[i];
            uuu[i, 1] = du;
            uuu[i, 2] = ddu;
          }

          /////////////////////////////////////////////////////////////////
          //ANSYS The generalized HHT-alpha method
          //double gamma = 0.05;
          //double delta = 1.0 / 2.0 + gamma;
          //double alpha = 1.0 / 4.0 * (1.0 + gamma) * (1.0 + gamma);
          //double alphaf = (1.0 - gamma) / 2.0;
          //double alpham = (1.0 - 3.0 * gamma) / 2.0;
          //double a0 = (1.0 - alpham) / (alpha * deltaT * deltaT);
          //double a1 = (1.0 - alphaf) * delta / (alpha * deltaT);
          //double a2 = a0 * deltaT;
          //double a3 = (1.0 - alpham) / (2.0 * alpha) - 1.0;
          //double a4 = (1.0 - alphaf) * delta / alpha - 1.0;
          //double a5 = (1.0 - alphaf) * (delta / (2.0 * alpha) - 1.0) * deltaT;

          //double aa0 = 1.0 / (alpha * deltaT * deltaT);
          //double aa1 = delta / (alpha * deltaT);
          //double aa2 = 1.0 / (alpha * deltaT);
          //double aa3 = 1.0 / (2.0 * alpha) - 1.0;
          //double aa4 = delta / alpha - 1.0;
          //double aa5 = (deltaT / 2.0) * (delta / alpha - 2.0);
          //double aa6 = deltaT * (1 - delta);
          //double aa7 = delta * deltaT;

          //AssemblyStiffnessMatrix(ref kGlobal);
          //AssemblyMassMatrix(ref mGlobal);

          //// kNew = K + b1*M + b4* C
          //DoubleMatrix kNew = (1.0 - alphaf) * kGlobal + a0 * mGlobal;//+ a1 * cGlobal));
          //DoubleMatrix uuu = ApplyInitialBoundaryCondition();// Store u, du, ddu
          //displacementTime = DoubleMatrix(countDOF, numberOfStepSave);
          //DoubleVector time = DoubleVector(numberOfStepSave);
          //int stepSave = 0;
          //DoubleVector rGlobaln = null;
          //AssemblyTractionVector(ref rGlobaln, 0.0);
          //for (int step = 1; step <= numberOfStep; step++)
          //{
          //    DoubleVector rGlobalnn = null;
          //    AssemblyTractionVector(ref rGlobalnn, step * deltaT);
          //    DoubleVector uuuCol0 = uuu.Column(0);
          //    DoubleVector uuuCol1 = uuu.Column(1);
          //    DoubleVector uuuCol2 = uuu.Column(2);
          //    //DoubleVector cRight = a1 * uuu.Column(0) + a4 * uuu.Column(1) + a5 * uuu.Column(2);
          //    rGlobal = (1.0 - alphaf) * rGlobalnn + alphaf * rGlobaln - alphaf * kGlobal * uuuCol0 + mGlobal * (a0 * uuuCol0 + a2 * uuuCol1 + a3 * uuuCol2);
          //    rGlobaln = rGlobalnn;
          //    DoubleMatrix kNewClone = 1.0 * kNew;
          //    ApplyBoundaryCondition(ref kNewClone, ref rGlobal, step * deltaT);
          //    solver = new Solver(kNewClone, rGlobal);
          //    if (IsSparseData)
          //        uGlobal = solver.Solve(SolverType.TFQMR, TOL, maxIter);
          //    else
          //        uGlobal = solver.Solve(SolverType.QR);
          //    // Compute velocity, acceleration
          //    for (int i = 0; i < countDOF; i++)
          //    {
          //        double du = aa1 * (uGlobal[i] - uuu[i, 0]) - aa4 * uuu[i, 1] - aa5 * uuu[i, 2];
          //        double ddu = aa0 * (uGlobal[i] - uuu[i, 0]) - aa2 * uuu[i, 1] - aa3 * uuu[i, 2];
          //        uuu[i, 0] = uGlobal[i];
          //        uuu[i, 1] = du;
          //        uuu[i, 2] = ddu;
          //    }

          int currentStepSave = (stepSave + 1) * (int)Math.Ceiling(Math.Round(((double)numberOfSubstepLoadLocal) / ((double)numberOfStepSaveLocal), 1));

          if (substep == currentStepSave)
          {
            //DisplacementTime.Add(uGlobal);

            if (countData > 0)
            {
              MonitorDataHistorism monitorDataHistorism = GetMonitorDataHistorism();
              GetMonitorDataHistorism().MonitorDataStore.Add(new double[countData]);
              for (int j = 0; j < countData; j++)
              {
                MonitorDataNode nodeData = monitorDataHistorism.GetMonitorDataNode(j);
                double[] xi = nodeData.Xi;
                Result result = nodeData.ResultType;
                AbstractPatch patch = nodeData.Patch;
                switch (result)
                {
                  case Result.UX:
                  case Result.UY:
                  case Result.UZ:
                    double u = 0;
                    if (patch is PatchStructure2D)
                    {
                      u = ((PatchStructure2D)patch).GetApproximateAt(result, xi[0], xi[1]);
                    }
                    else if (patch is PatchStructure3D)
                    {
                      u = ((PatchStructure3D)patch).GetApproximateAt(result, xi[0], xi[1], xi[2]);
                    }
                    GetMonitorDataHistorism().MonitorDataStore[GetMonitorDataHistorism().MonitorDataStore.Count][j] = u;
                    break;
                  case Result.USUM:
                    double uu = 0;
                    if (patch is PatchStructure2D)
                    {
                      double[] temp = new double[] { ((PatchStructure2D)patch).GetApproximateAt(Result.UX, xi[0], xi[1]),
                                    ((PatchStructure2D)patch).GetApproximateAt(Result.UY, xi[0], xi[1])};
                      uu = temp[0] * temp[0] + temp[1] * temp[1];
                    }
                    else if (patch is PatchStructure3D)
                    {
                      double[] temp = new double[] { ((PatchStructure3D)patch).GetApproximateAt(Result.UX, xi[0], xi[1], xi[2]),
                                    ((PatchStructure3D)patch).GetApproximateAt(Result.UY, xi[0], xi[1], xi[2]),
                                    ((PatchStructure3D)patch).GetApproximateAt(Result.UZ, xi[0], xi[1], xi[2])};
                      uu = temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2];
                    }
                    GetMonitorDataHistorism().MonitorDataStore[GetMonitorDataHistorism().MonitorDataStore.Count][j] = Math.Sqrt(uu);
                    break;
                  case Result.SIGMAXX:
                  case Result.SIGMAYY:
                  case Result.SIGMAXY:
                  case Result.SIGMAZZ:
                  case Result.SIGMAYZ:
                  case Result.SIGMAXZ:
                  case Result.SIGMAEQV:
                    if (result != Result.SIGMAEQV)
                    {
                      Extrapolation(result);
                    }
                    else
                    {
                      Extrapolation(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAZZ, Result.SIGMAXY, Result.SIGMAYZ, Result.SIGMAXZ);
                    }
                    if ((patch is PatchStructure3D) || ((patch is PatchStructure2D) && (result != Result.SIGMAYZ && result != Result.SIGMAXZ)))
                    {
                      GetMonitorDataHistorism().MonitorDataStore[GetMonitorDataHistorism().MonitorDataStore.Count][j] = patch.GetApproximateAt(result, xi);
                    }
                    //else if (patch is PatchStructure3D)
                    //{
                    //   MonitorDataStore[stepSave, j] = ((PatchStructure3D)patch).GetStressAt(xi[0], xi[1], xi[2], result);
                    //}
                    break;
                    //case Result.EPSILONXX:
                    //case Result.EPSILONYY:
                    //case Result.EPSILONXY:
                    //case Result.EPSILONZZ:
                    //case Result.EPSILONYZ:
                    //case Result.EPSILONXZ:
                    //case Result.EPSILONEQV:
                    //	model.ExtrapolationStrain();
                    //	if (patch is Structure2DPatch)
                    //	{
                    //		if (result != Result.EPSILONYZ && result != Result.EPSILONXZ)
                    //			data[totalStep - 1, j] = ((Structure2DPatch)patch).GetStrainAt(xi[0], xi[1], result);
                    //	}
                    //	else if (patch is Structure3DPatch)
                    //	{
                    //		data[totalStep - 1, j] = ((Structure3DPatch)patch).GetStrainAt(xi[0], xi[1], xi[2], result);
                    //	}
                    //	break;
                }
              }
            }
            Time.Add(time);
            if (IsSaveStepByStep)
            {
              //IOIGA.WriteListDoubleVectorMessagePackFormatter(FileNameDataTime, DisplacementTime);
              string fileNameDataTime = PathProject + "data_" + step + ".msgpack";
              IOIGA.WriteMessagePackDoubleVector(fileNameDataTime, uGlobal);
              IOIGA.WriteMessagePackData(FileNameTime, Time);
              if (GetMonitorDataHistorism().MonitorDataStore != null)
              {
                IOIGA.WriteExcelResult(PathProject + "result_" + NameProject + ".xls", GetMonitorDataHistorism().MonitorDataStore);
              }
            }
            stepSave++;
          }
        }
      }
      if (!IsSaveStepByStep)
      {
        IOIGA.WriteMessagePackDoubleVector(FileNameDataTime, uGlobal);
        IOIGA.WriteMessagePackData(FileNameTime, Time);
        if (GetMonitorDataHistorism().MonitorDataStore != null)
        {
          IOIGA.WriteExcelResult(PathProject + "result_" + NameProject + ".xls", GetMonitorDataHistorism().MonitorDataStore);
        }
      }
    }
  }
}
