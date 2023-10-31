using System;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Abstract thermal model
  /// </summary>
  public class ModelThermalTransient : AbstractModelThermal
  {
    //private double deltaT;// = 4e-6; step size
    //private int numberOfStep;// = 6000;
    //private int numberOfSubstepLoad;// = 6000;
    //private int numberOfStepSave;// = 100;

    /// <summary>
    /// Constructor class
    /// </summary>
    public ModelThermalTransient(Dimension structureDimension, string pathProject = "temp", string nameProject = "project")
        : base(TypeAnalysisModel.Transient, structureDimension, pathProject, nameProject)
    {
      listInitialConstraint = new System.Collections.Generic.List<AbstractConstraintValue>();
    }

    /// <summary>
    /// Solve problem
    /// </summary>
    public override void Solve()
    {
      if (listPatch.Count == 0)
        throw new NullReferenceException("Model must be assigned Patch");
      WriteInformationProblem("Transient thermal module");

      DoubleMatrix kGlobalDense = null; // stiffness
      DoubleMatrix mGlobalDense = null; // "mass"
      DoubleVector uGlobal = null; // temperature
      DoubleVector rGlobal = null; // right hand side

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



      // Backward Euler for time integration
      double a1 = 1.0 / deltaT[0];
      DoubleMatrix kNewDense = null; // M/deltaT + K;
      DoubleCsrSparseMatrix kNewSparse = null;
      if (!IsSparseData)
      {
        kNewDense = kGlobalDense + a1 * mGlobalDense;
      }
      else
      {
        kNewSparse = kGlobalSparse + a1 * mGlobalSparse;
      }
      DoubleMatrix uuu = ApplyInitialBoundaryCondition();// Store u, du
      //DisplacementTime = new List<DoubleVector>();// (countDOF, numberOfStepSave[0]);
      double[] time = new double[numberOfStepSave[0]];
      int stepSave = 0;

      // Time step
      for (int step = 1; step <= numberOfSubstepLoad[0]; step++)
      {
        if (IsWriteLogFile)
        {
          IOIGA.WriteLogFileAndConsole(logFile, true, "-----------------------------------------------------------");
          IOIGA.WriteLogFileAndConsole(logFile, true, "Load step " + step + " : ");
        }
        AssemblyTractionVector(out rGlobal, step * deltaT[0]);
        AssemplyHeatGenerationSource(ref rGlobal);
        AssemblyHeatTransferConvectionVector(ref rGlobal);

        DoubleVector mRight = a1 * uuu.Col(0); // DoubleVector(countDOF);
                                               //DoubleVector cRight = a1 * uuu.Column(0) + a4 * uuu.Column(1) + a5 * uuu.Column(2); //DoubleVector(countDOF);
        DoubleSymmetricMatrix kNewCloneDense = null;
        DoubleCsrSparseMatrix kNewCloneSparse = null;
        if (!IsSparseData)
        {
          rGlobal += MatrixFunctions.Product(mGlobalDense, mRight);// + cGlobal * cRight;
          kNewCloneDense = (DoubleSymmetricMatrix)kNewDense.Clone();
        }
        else
        {
          rGlobal += MatrixFunctions.Product(mGlobalSparse, mRight);// + cGlobal * cRight;
          kNewCloneSparse = (DoubleCsrSparseMatrix)kNewSparse.Clone();
        }
        ApplyChangeBoundaryConditionAndUpdateIMatrix(step * deltaT[0]);

        uGlobal = SolveStatic(kNewCloneDense, kNewCloneSparse, ref rGlobal);

        // Compute velocity, acceleration
        for (int i = 0; i < countDOF; i++)
        {
          //double du = a1 * (uGlobal[i] - uuu[i, 0]);
          uuu[i, 0] = uGlobal[i]; // update
        }

        int currentStepSave = (stepSave + 1) * (int)Math.Ceiling(Math.Round(((double)numberOfSubstepLoad[0]) / ((double)numberOfStepSave[0]), 1));
        if (step == currentStepSave)
        {
          //DisplacementTime.Add(uGlobal);
          time[stepSave] = step * deltaT[0];
          if (IsSaveStepByStep)
          {
            //IOIGA.WriteListDoubleVectorMessagePackFormatter(FileNameDataTime, DisplacementTime);
            IOIGA.WriteMessagePackData(FileNameTime, Time);
            string fileNameDataTime = PathProject + "data_" + step + ".msgpack";
            IOIGA.WriteMessagePackDoubleVector(fileNameDataTime, uGlobal);
          }
          stepSave++;
        }
        if (IsWriteLogFile)
        {
          IOIGA.WriteLogFileAndConsole(logFile, true, "Elapsed time = " + stopWatchWholeModel.Elapsed);
        }
      }

      if (!IsSaveStepByStep)
      {
        //IOIGA.WriteListDoubleVectorMessagePackFormatter(FileNameDataTime, DisplacementTime);
        IOIGA.WriteMessagePackData(FileNameTime, Time);
        IOIGA.WriteMessagePackDoubleVector(FileNameDataTime, uGlobal);
      }
    }


    //public void AssemblyMassMatrix(out SparseMatrixBuilder<double> mGlobal)
    //{
    //  mGlobal = new SparseMatrixBuilder<double>();
    //  for (int i = 0; i < listPatch.Count; i++)
    //  {
    //    switch (StructureDimension)
    //    {
    //      case Dimension.Plane:
    //        ((PatchThermal2D)listPatch[i]).ComputeMassMatrixPatch(ref mGlobal);
    //        break;
    //      case Dimension.Solid:
    //        ((PatchThermal3D)listPatch[i]).ComputeMassMatrixPatch(ref mGlobal);
    //        break;
    //    }
    //  }
    //}
    //public void AssemblyMassMatrix(out DoubleMatrix mGlobal)
    //{
    //  mGlobal = new DoubleMatrix(countDOF, countDOF);
    //  for (int i = 0; i < listPatch.Count; i++)
    //  {
    //    switch (StructureDimension)
    //    {
    //      case Dimension.Plane:
    //        ((PatchThermal2D)listPatch[i]).ComputeMassMatrixPatch(ref mGlobal);
    //        break;
    //      case Dimension.Solid:
    //        ((PatchThermal3D)listPatch[i]).ComputeMassMatrixPatch(ref mGlobal);
    //        break;
    //    }
    //  }
    //}
  }
}
