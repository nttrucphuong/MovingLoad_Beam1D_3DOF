using System;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{

  /// <summary>
  /// Structure 2D problem (plane stress, plane strain). Default is plane stress
  /// </summary>
  public class ModelPiezoelectricStatic : AbstractModelPiezoelectric
  {
    /// <summary>
    /// Constructor class
    /// </summary>
    public ModelPiezoelectricStatic(Dimension dimension, string pathProject = "temp", string nameProject = "project")
         : base(TypeAnalysisModel.Static, dimension, pathProject, nameProject)
    {
    }

    public override void Solve()
    {
      ////////////////////////////////////////////////////////
      ////////// direct couple ////////////////////////
      ////////////////////////////////////////////////////////
      //if (listPatch.Count == 0)
      //    throw new NullReferenceException("Model must be assigned Patch");
      //DoubleMatrix kGlobal = null;
      //DoubleVector uGlobal = null;
      //DoubleVector rGlobal = null;
      //Solver solver = null;
      //string nameFile = null;

      //if (observe != null)
      //    observe.WriteLine("Begin solving piezoelectric analysis");

      //if (observe != null)
      //    observe.WriteLine("Assembly Stiffness matrix");
      //AssemblyStiffnessMatrix(ref kGlobal);
      //if (observe != null)
      //    observe.WriteLine("Assembly Traction vector");
      //AssemblyTractionVector(ref rGlobal);
      //ApplyBoundaryCondition(ref kGlobal, ref rGlobal);
      //if (observe != null)
      //    observe.WriteLine("Solving the equation");
      //solver = new Solver(kGlobal, rGlobal);
      //if (IsUseMatlabSolver)
      //    uGlobal = solver.Solve(matlab);
      //else
      //{
      //    if (IsSparseData)
      //        uGlobal = (DoubleVector)solver.Solve(SolverType.BiCgStab, 1e-6, 100);
      //    else
      //        uGlobal = (DoubleVector)solver.Solve(SolverType.QR);
      //}
      ////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////
      ////////// serial couple ////////////////////////
      ////////////////////////////////////////////////////////
      WriteInformationProblem("Static piezoelectric module");
      //DisplacementTime = new List<DoubleVector>();
      if (listPatch.Count == 0)
        throw new NullReferenceException("Model must be assigned Patch");

      DoubleMatrix kGlobalDense = null;
      DoubleVector u = null;
      DoubleVector uGlobal = null;
      DoubleVector rGlobal = null;

      if (!IsSparseData)
      {
        AssemblyStiffnessMatrix(out kGlobalDense);
      }
      else
      {
        AssemblyStiffnessMatrix(ref kGlobalSparse);
      }

      AssemblyTractionVector(out rGlobal);
      //ApplyBoundaryCondition(ref kGlobal, ref rGlobal);

      if (listActuatorPatch == null && listSensorPatch == null)
      {

        //// |K11  K12  K13  K14| | u  |    |F|
        //// |K21  K22  K23  K24| | u0 | =  | |
        //// |K31  K32  K33  K34| |phi |    |Q|
        //// |K41  K42  K43  K44| |phi0|    | |
        //// (K11-K13*K33^-1*K31) * u =F -K12 * u0 - K14 * phi0 + K13*K33^-1*K32 * u0 + K13*K33^-1*K34 * phi0 - K13*K33^-1*Q
        //DoubleMatrix K11 = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnconstraint[0], tArrayUnconstraint[0]);
        //DoubleMatrix K33 = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnconstraint[1], tArrayUnconstraint[1]);
        //DoubleMatrix K31 = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnconstraint[1], tArrayUnconstraint[0]);
        //DoubleMatrix K13 = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnconstraint[0], tArrayUnconstraint[1]);
        //DoubleVector F = MatrixTool.GetSubVector(rGlobal, tArrayUnconstraint[0]);
        //DoubleVector Q = MatrixTool.GetSubVector(rGlobal, tArrayUnconstraint[1]);
        //DoubleVector u0 = GetConstraint(tArrayConstraint[0]);// GetUConstraint();
        //DoubleVector phi0 = GetConstraint(tArrayConstraint[1]);//GetPhiConstraint();
        //DoubleMatrix K12 = null;
        //DoubleMatrix K14 = null;
        //DoubleMatrix K32 = null;
        //DoubleMatrix K34 = null;
        //if (u0 != null)
        //{
        //  K12 = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnconstraint[0], tArrayConstraint[0]);
        //  K32 = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnconstraint[1], tArrayConstraint[0]);
        //}
        //if (phi0 != null)
        //{
        //  K14 = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnconstraint[0], tArrayConstraint[1]);
        //  K34 = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnconstraint[1], tArrayConstraint[1]);
        //}
        //DoubleMatrix K33Inv = null;
        //DoubleMatrix Knew = null;
        //DoubleVector rNew = null;
        //if (K33 != null)
        //{
        //  K33Inv = MatrixFunctions.PseudoInverse(K33);
        //  //DoubleMatrix K33InvSparse = K33Inv;
        //  Knew = K11 - MatrixFunctions.Product(K13, MatrixFunctions.Product(K33Inv, K31));
        //  rNew = F - MatrixFunctions.Product(MatrixFunctions.Product(K13, K33Inv), Q);
        //  if (u0 != null)
        //  {
        //    rNew += -MatrixFunctions.Product(K12, u0) + MatrixFunctions.Product(MatrixFunctions.Product(K13, MatrixFunctions.Product(K33Inv, K32)), u0);
        //  }
        //  if (phi0 != null)
        //  {
        //    rNew += -MatrixFunctions.Product(K14, phi0) + MatrixFunctions.Product(MatrixFunctions.Product(K13, MatrixFunctions.Product(K33Inv, K34)), phi0);
        //  }
        //}
        //else
        //{
        //  Knew = K11;
        //  rNew = F;
        //}
        //u = MatrixFunctions.Solve(Knew, rNew);
        //DoubleVector phi = null;
        //if (K33 != null)
        //{
        //  Knew = K33;
        //  rNew = Q - MatrixFunctions.Product(K31, u);
        //  if (u0 != null)
        //  {
        //    rNew += -MatrixFunctions.Product(K32, u0);
        //  }
        //  if (phi0 != null)
        //  {
        //    rNew += -MatrixFunctions.Product(K34, phi0);
        //  }
        //  phi = MatrixFunctions.Solve(Knew, rNew);
        //}

        //uGlobal = new DoubleVector(countDOF);
        //for (int i = 0; i < u.Length; i++)
        //{
        //  int tarray = tArrayUnconstraint[0][i];
        //  uGlobal[tarray] = u[i];
        //}
        //if (u0 != null)
        //  for (int i = 0; i < u0.Length; i++)
        //  {
        //    int tarray = tArrayConstraint[0][i];
        //    uGlobal[tarray] = u0[i];
        //  }
        //if (phi != null)
        //{
        //  for (int i = 0; i < phi.Length; i++)
        //  {
        //    int tarray = tArrayUnconstraint[1][i];
        //    uGlobal[tarray] = phi[i];
        //  }
        //}
        //if (phi0 != null)
        //{
        //  for (int i = 0; i < phi0.Length; i++)
        //  {
        //    int tarray = tArrayConstraint[1][i];
        //    uGlobal[tarray] = phi0[i];
        //  }
        //}

        uGlobal = SolveStatic(kGlobalDense, kGlobalSparse, ref rGlobal);
      }
      else
      {
        // Muu * u.. + (Cs + Cu) * u. + Kuu* * u = Fm
        // Cu = Gv * [Kuphi]a * ([Kphiphi]s)^(-1) * [Kphiu]s
        // Kuu* = Kuu + Gd * [Kuphi]a * ([Kphiphi]s)^(-1) * [Kphiu]s

      }
      //DisplacementTime.Add(uGlobal);
      IOIGA.WriteMessagePackDoubleVector(FileNameDataTime, uGlobal);
      //IOIGA.WriteListDoubleVectorMessagePackFormatter(FileNameDataTime, DisplacementTime);
    }
  }
}
