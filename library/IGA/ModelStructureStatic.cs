using System;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Structure 2D problem (plane stress, plane strain). Default is plane stress
  /// </summary>
  public class ModelStructureStatic : AbstractModelStructure
  {
    private int numberOfStep;// = 6000;
    private int numberOfSubstepLoad;// = 6000;

    /// <summary>
    /// Constructor class
    /// </summary>
    /// <param name="problem">Define state of stress, plane stress or plane strain</param>
    public ModelStructureStatic(Dimension structureDimension, string pathProject = "temp", string nameProject = "project")
                 : base(TypeAnalysisModel.Static, structureDimension, pathProject, nameProject)
    {
    }

    /// <summary>
    /// Solve problem
    /// </summary>
    public override void Solve()
    {
      WriteInformationProblem("Static structural module");

      //DisplacementTime = new List<DoubleVector>();//(countDOF, 1);
      if (listPatch.Count == 0)
        throw new NullReferenceException("Model must be assigned Patch");

      switch (listPatch[0].Material.TypeMaterialStructure)
      {
        case TypeMaterialStructure.Elasticity:
          DoubleMatrix kGlobalDense = null;
          DoubleVector rGlobal = null;
          DoubleVector uGlobal = null;
          if (!IsSparseData)
          {
            AssemblyStiffnessMatrix(out kGlobalDense);
          }
          else
          {
            AssemblyStiffnessMatrix(ref kGlobalSparse);
          }
          AssemblyTractionVector(out rGlobal);
          uGlobal = SolveStatic(kGlobalDense, kGlobalSparse, ref rGlobal);

          #region store
          //bool a = IsSparseData;
          //if (a)
          //{
          //  if (!IsSparseData)
          //  {
          //    ApplyTAndg(ref kGlobalDense, ref rGlobal, null);
          //    uGlobal = MatrixFunctions.Solve(kGlobalDense, rGlobal);
          //    uGlobal = MatrixFunctions.Product(TDense, uGlobal) + g;
          //  }
          //  else
          //  {
          //    ApplyTAndg(ref kGlobalSparse, ref rGlobal, null);
          //    uGlobal = MatrixFunctions.Solve(kGlobalSparse, rGlobal);
          //    uGlobal = MatrixFunctions.Product(TSparse, uGlobal) + g;
          //  }
          //}
          //else
          //{
          //  int[] tArrayIndependent = this.tArrayIndependent.ToArray();
          //  int[] tArrayDependent = this.tArrayDependent.ToArray();
          //  int[] tArrayUnCoupleAndUnconstraint = this.tArrayUncouple.ToArray();
          //  int[] tArrayIndependentInNew = new int[tArrayIndependent.Length];
          //  int c = 0;
          //  for (int i = tArrayUnCoupleAndUnconstraint.Length; i < tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length; i++)
          //  {
          //    tArrayIndependentInNew[c++] = i;
          //  }

          //  int[] tArrayUncoupleUnconstraintInNew = new int[tArrayUnCoupleAndUnconstraint.Length];
          //  for (int i = 0; i < tArrayUnCoupleAndUnconstraint.Length; i++)
          //  {
          //    tArrayUncoupleUnconstraintInNew[i] = i;
          //  }
          //  DoubleVector gS = null;
          //  if (tArrayDependent.Length != 0)
          //  {
          //    gS = MatrixTool.GetSubVector(g, tArrayDependent);
          //  }
          //  DoubleVector uLocalS = gS;
          //  DoubleVector uLocalUM = null;
          //  DoubleVector uLocalU = null;
          //  DoubleVector uLocalM = null;
          //  DoubleVector Ru = null;
          //  DoubleVector Rm = null;
          //  SparseMatrixBuilder<double> lhsSparse = null;// new SparseMatrixBuilder<double>();
          //  //DoubleCsrSparseMatrix lhsSparse = null;
          //  DoubleMatrix lhsDense = null;
          //  DoubleVector rhs = new DoubleVector(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
          //  DoubleMatrix TsmDense = null;
          //  DoubleCsrSparseMatrix TsmSparse = null;
          //  if (!IsSparseData)
          //  {
          //    DoubleMatrix Kuu = null;
          //    DoubleMatrix Kus = null;
          //    DoubleMatrix Kum = null;
          //    DoubleMatrix Kmm = null;
          //    if (tArrayUnCoupleAndUnconstraint.Length != 0)
          //    {
          //      Kuu = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnCoupleAndUnconstraint, tArrayUnCoupleAndUnconstraint);
          //      Ru = MatrixTool.GetSubVector(rGlobal, tArrayUnCoupleAndUnconstraint);
          //      if (tArrayDependent.Length != 0)
          //      {
          //        Kus = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnCoupleAndUnconstraint, tArrayDependent);
          //        Ru += -(NMathFunctions.Product(Kus, gS));
          //      }
          //    }

          //    if (tArrayDependent.Length != 0 && tArrayIndependent.Length != 0)
          //    {
          //      TsmDense = MatrixTool.GetSubMatrix(TDense, tArrayDependent, tArrayIndependent);
          //      DoubleMatrix TsmTranspose = TsmDense.Transpose();
          //      DoubleMatrix Kms = MatrixTool.GetSubMatrix(kGlobalDense, tArrayIndependent, tArrayDependent);
          //      DoubleMatrix Kss = MatrixTool.GetSubMatrix(kGlobalDense, tArrayDependent, tArrayDependent);
          //      DoubleVector Rs = MatrixTool.GetSubVector(rGlobal, tArrayDependent);
          //      Kmm = MatrixTool.GetSubMatrix(kGlobalDense, tArrayIndependent, tArrayIndependent) + NMathFunctions.Product(TsmTranspose, Kms.Transpose() + NMathFunctions.Product(Kss, TsmDense)) + NMathFunctions.Product(Kms, TsmDense);
          //      Rm = MatrixTool.GetSubVector(rGlobal, tArrayIndependent) + NMathFunctions.Product(TsmTranspose, Rs);
          //      if (tArrayUnCoupleAndUnconstraint.Length != 0)
          //        Kum = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnCoupleAndUnconstraint, tArrayIndependent) + NMathFunctions.Product(Kus, TsmDense);
          //      Rm += -(NMathFunctions.Product(Kms, gS) + NMathFunctions.Product(TsmTranspose, NMathFunctions.Product(Kss, gS)));
          //    }
          //    lhsDense = new DoubleMatrix(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length, tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
          //    //rhs = new DoubleVector(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
          //    if (!AbstractModel.IsParallelProcesing)
          //    {
          //      for (int i = 0; i < tArrayUnCoupleAndUnconstraint.Length; i++)
          //      {
          //        for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
          //        {
          //          lhsDense[i, j] = Kuu[i, j];
          //        }
          //        if (Kum != null)
          //        {
          //          for (int j = 0; j < tArrayIndependent.Length; j++)
          //          {
          //            lhsDense[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
          //          }
          //        }
          //        rhs[i] = Ru[i];
          //      }
          //      for (int i = 0; i < tArrayIndependent.Length; i++)
          //      {
          //        for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
          //        {
          //          lhsDense[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
          //        }

          //        for (int j = 0; j < tArrayIndependent.Length; j++)
          //        {
          //          lhsDense[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
          //        }
          //        rhs[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
          //      }
          //    }
          //    else
          //    {
          //      var degreeOfParallelism = NumberOfCPUs;
          //      var tasks = new Task[degreeOfParallelism];
          //      for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
          //      {
          //        // capturing taskNumber in lambda wouldn't work correctly
          //        int taskNumberCopy = taskNumber;

          //        tasks[taskNumber] = Task.Factory.StartNew(
          //             () =>
          //             {
          //               var max = tArrayUnCoupleAndUnconstraint.Length * (taskNumberCopy + 1) / degreeOfParallelism;
          //               for (int i = tArrayUnCoupleAndUnconstraint.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
          //               {
          //                 for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
          //                 {
          //                   lhsDense[i, j] = Kuu[i, j];
          //                 }
          //                 if (Kum != null)
          //                 {
          //                   for (int j = 0; j < tArrayIndependent.Length; j++)
          //                   {
          //                     lhsDense[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
          //                   }
          //                 }
          //                 rhs[i] = Ru[i];
          //               }
          //             });
          //      }
          //      Task.WaitAll(tasks);

          //      for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
          //      {
          //        // capturing taskNumber in lambda wouldn't work correctly
          //        int taskNumberCopy = taskNumber;

          //        tasks[taskNumber] = Task.Factory.StartNew(
          //             () =>
          //             {
          //               var max = tArrayIndependent.Length * (taskNumberCopy + 1) / degreeOfParallelism;
          //               for (int i = tArrayIndependent.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
          //               {
          //                 for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
          //                 {
          //                   lhsDense[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
          //                 }

          //                 for (int j = 0; j < tArrayIndependent.Length; j++)
          //                 {
          //                   lhsDense[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
          //                 }
          //                 rhs[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
          //               }
          //             });
          //      }
          //      Task.WaitAll(tasks);
          //    }
          //  }
          //  else
          //  {
          //    DoubleCsrSparseMatrix Kuu = null;
          //    DoubleCsrSparseMatrix Kus = null;
          //    DoubleCsrSparseMatrix Kum = null;
          //    DoubleCsrSparseMatrix Kmm = null;
          //    if (tArrayUnCoupleAndUnconstraint.Length != 0)
          //    {
          //      Kuu = MatrixTool.GetSubMatrix(kGlobalSparse, tArrayUnCoupleAndUnconstraint, tArrayUnCoupleAndUnconstraint);
          //      Ru = MatrixTool.GetSubVector(rGlobal, tArrayUnCoupleAndUnconstraint);
          //      if (tArrayDependent.Length != 0)
          //      {
          //        Kus = MatrixTool.GetSubMatrix(kGlobalSparse, tArrayUnCoupleAndUnconstraint, tArrayDependent);
          //      }
          //    }

          //    if (tArrayDependent.Length != 0 && tArrayIndependent.Length != 0)
          //    {
          //      TsmSparse = MatrixTool.GetSubMatrix(TSparse, tArrayDependent, tArrayIndependent);
          //      DoubleCsrSparseMatrix TsmTranspose = MatrixTool.Transpose(TsmSparse);
          //      DoubleCsrSparseMatrix Kms = MatrixTool.GetSubMatrix(kGlobalSparse, tArrayIndependent, tArrayDependent);
          //      DoubleCsrSparseMatrix Kss = MatrixTool.GetSubMatrix(kGlobalSparse, tArrayDependent, tArrayDependent);
          //      DoubleVector Rs = MatrixTool.GetSubVector(rGlobal, tArrayDependent);
          //      DoubleCsrSparseMatrix KsmKss = MatrixTool.Transpose(Kms) + MatrixFunctions.Product(Kss, TsmSparse);
          //      Kmm = MatrixTool.GetSubMatrix(kGlobalSparse, tArrayIndependent, tArrayIndependent) + MatrixFunctions.Product(TsmTranspose, KsmKss) + MatrixFunctions.Product(Kms, TsmSparse);
          //      Rm = MatrixTool.GetSubVector(rGlobal, tArrayIndependent) + MatrixFunctions.Product(TsmTranspose, Rs);
          //      if (tArrayUnCoupleAndUnconstraint.Length != 0)
          //      {
          //        var productKusTsm = MatrixFunctions.Product(Kus, TsmSparse);
          //        DoubleCsrSparseMatrix subSparse = MatrixTool.GetSubMatrix(kGlobalSparse, tArrayUnCoupleAndUnconstraint, tArrayIndependent);
          //        Kum = subSparse + productKusTsm;
          //      }
          //      Rm += -(MatrixFunctions.Product(Kms, gS) + MatrixFunctions.Product(TsmTranspose, MatrixFunctions.Product(Kss, gS)));

          //      if (tArrayUnCoupleAndUnconstraint.Length != 0)
          //      {
          //        Ru += -(MatrixFunctions.Product(Kus, gS));
          //      }
          //    }
          //    lhsSparse = new SparseMatrixBuilder<double>();

          //    if (!AbstractModel.IsParallelProcesing)
          //    {
          //      for (int i = 0; i < tArrayUnCoupleAndUnconstraint.Length; i++)
          //      {
          //        for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
          //        {
          //          lhsSparse[i, j] = Kuu[i, j];
          //        }
          //        if (Kum != null)
          //        {
          //          for (int j = 0; j < tArrayIndependent.Length; j++)
          //          {
          //            lhsSparse[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
          //          }
          //        }
          //        rhs[i] = Ru[i];
          //      }
          //      for (int i = 0; i < tArrayIndependent.Length; i++)
          //      {
          //        for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
          //        {
          //          lhsSparse[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
          //        }

          //        for (int j = 0; j < tArrayIndependent.Length; j++)
          //        {
          //          lhsSparse[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
          //        }
          //        rhs[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
          //      }
          //    }
          //    else
          //    {
          //      var degreeOfParallelism = NumberOfCPUs;
          //      var tasks = new Task[degreeOfParallelism];
          //      for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
          //      {
          //        // capturing taskNumber in lambda wouldn't work correctly
          //        int taskNumberCopy = taskNumber;

          //        tasks[taskNumber] = Task.Factory.StartNew(
          //             () =>
          //             {
          //               var max = tArrayUnCoupleAndUnconstraint.Length * (taskNumberCopy + 1) / degreeOfParallelism;
          //               for (int i = tArrayUnCoupleAndUnconstraint.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
          //               {
          //                 for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
          //                 {
          //                   lhsSparse[i, j] = Kuu[i, j];
          //                 }
          //                 if (Kum != null)
          //                 {
          //                   for (int j = 0; j < tArrayIndependent.Length; j++)
          //                   {
          //                     lhsSparse[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
          //                   }
          //                 }
          //                 rhs[i] = Ru[i];
          //               }
          //             });
          //      }
          //      Task.WaitAll(tasks);

          //      for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
          //      {
          //        // capturing taskNumber in lambda wouldn't work correctly
          //        int taskNumberCopy = taskNumber;

          //        tasks[taskNumber] = Task.Factory.StartNew(
          //             () =>
          //             {
          //               var max = tArrayIndependent.Length * (taskNumberCopy + 1) / degreeOfParallelism;
          //               for (int i = tArrayIndependent.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
          //               {
          //                 for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
          //                 {
          //                   lhsSparse[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
          //                 }

          //                 for (int j = 0; j < tArrayIndependent.Length; j++)
          //                 {
          //                   lhsSparse[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
          //                 }
          //                 rhs[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
          //               }
          //             });
          //      }
          //      Task.WaitAll(tasks);
          //    }
          //  }

          //  if (rGlobal.Length != 0)
          //  {
          //    if (!IsSparseData)
          //      uLocalUM = MatrixFunctions.Solve(lhsDense, rhs);
          //    else
          //    {

          //      uLocalUM = MatrixFunctions.Solve(new DoubleCsrSparseMatrix(lhsSparse, countDOF), rhs);
          //    }
          //    uLocalU = MatrixTool.GetSubVector(uLocalUM, tArrayUncoupleUnconstraintInNew);
          //    uLocalM = MatrixTool.GetSubVector(uLocalUM, tArrayIndependentInNew);
          //    if (uLocalM != null)
          //    {
          //      if (!IsSparseData)
          //        uLocalS += MatrixFunctions.Product(TsmDense, uLocalM);
          //      else
          //        uLocalS += MatrixFunctions.Product(TsmSparse, uLocalM);
          //    }
          //  }

          //  uGlobal = new DoubleVector(countDOF);
          //  if (!AbstractModel.IsParallelProcesing)
          //  {
          //    for (int i = 0; i < tArrayUncoupleUnconstraintInNew.Length; i++)
          //    {
          //      uGlobal[tArrayUnCoupleAndUnconstraint[i]] = uLocalU[i];
          //    }
          //    for (int i = 0; i < tArrayIndependentInNew.Length; i++)
          //    {
          //      uGlobal[tArrayIndependent[i]] = uLocalM[i];
          //    }
          //    for (int i = 0; i < tArrayDependent.Length; i++)
          //    {
          //      uGlobal[tArrayDependent[i]] = uLocalS[i];
          //    }
          //  }
          //  else
          //  {
          //    var degreeOfParallelism = NumberOfCPUs;
          //    var tasks = new Task[degreeOfParallelism];
          //    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
          //    {
          //      // capturing taskNumber in lambda wouldn't work correctly
          //      int taskNumberCopy = taskNumber;

          //      tasks[taskNumber] = Task.Factory.StartNew(
          //           () =>
          //           {
          //             var max = tArrayUncoupleUnconstraintInNew.Length * (taskNumberCopy + 1) / degreeOfParallelism;
          //             for (int i = tArrayUncoupleUnconstraintInNew.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
          //             {
          //               uGlobal[tArrayUnCoupleAndUnconstraint[i]] = uLocalU[i];
          //             }
          //           });
          //    }
          //    Task.WaitAll(tasks);
          //    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
          //    {
          //      // capturing taskNumber in lambda wouldn't work correctly
          //      int taskNumberCopy = taskNumber;

          //      tasks[taskNumber] = Task.Factory.StartNew(
          //           () =>
          //           {
          //             var max = tArrayIndependentInNew.Length * (taskNumberCopy + 1) / degreeOfParallelism;
          //             for (int i = tArrayIndependentInNew.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
          //             {
          //               uGlobal[tArrayIndependent[i]] = uLocalM[i];
          //             }
          //           });
          //    }
          //    Task.WaitAll(tasks);

          //    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
          //    {
          //      // capturing taskNumber in lambda wouldn't work correctly
          //      int taskNumberCopy = taskNumber;

          //      tasks[taskNumber] = Task.Factory.StartNew(
          //           () =>
          //           {
          //             var max = tArrayDependent.Length * (taskNumberCopy + 1) / degreeOfParallelism;
          //             for (int i = tArrayDependent.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
          //             {
          //               uGlobal[tArrayDependent[i]] = uLocalS[i];
          //             }
          //           });
          //    }
          //    Task.WaitAll(tasks);
          //  }
          //}
          #endregion

          //DisplacementTime.Add(uGlobal);
          //IOIGA.WriteListDoubleVectorMessagePackFormatter(FileNameDataTime, DisplacementTime);
          IOIGA.WriteMessagePackDoubleVector(FileNameDataTime, uGlobal);
          break;
        case TypeMaterialStructure.Plasticity:
        case TypeMaterialStructure.Hyperelasticity:
          switch (type)
          {
            case TypeNonlinearSolver.LoadControlledNewtonRaphson:
            case TypeNonlinearSolver.LoadControlledModifiedNewtonRaphson:
              //new SolverLoadControlled(this, type);
              break;
            case TypeNonlinearSolver.DisplacementControlledNewtonRaphson:
            case TypeNonlinearSolver.DisplacementControlledModifiedNewtonRaphson:
              new SolverDisplacementControlled(this, type);
              break;
          }
          break;
      }
    }





    //public void AssemblyInternalForce(ref DoubleVector ResidualGlobal)
    //{
    //  ResidualGlobal = new DoubleVector(countDOF);
    //  for (int i = 0; i < listPatch.Count; i++)
    //  {
    //    AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
    //    switch (StructureDimension)
    //    {
    //      case Dimension.Plane:
    //        ((PatchStructure2D)listPatch[i]).ComputeInternalForcePatch(ref ResidualGlobal);
    //        break;
    //      case Dimension.Solid:
    //        ((PatchStructure3D)listPatch[i]).ComputeInternalForcePatch(ref ResidualGlobal);
    //        break;
    //    }
    //  }
    //}
  }
}
