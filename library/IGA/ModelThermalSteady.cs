using System;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Abstract thermal model
  /// </summary>
  public class ModelThermalSteady : AbstractModelThermal
  {
    /// <summary>
    /// Constructor class
    /// </summary>
    public ModelThermalSteady(Dimension structureDimension, string pathProject = "temp", string nameProject = "project")
        : base(TypeAnalysisModel.Static, structureDimension, pathProject, nameProject)
    {
    }

    /// <summary>
    /// Solve problem
    /// </summary>
    public override void Solve()
    {
      if (listPatch.Count == 0)
        throw new NullReferenceException("Model must be assigned Patch");

      WriteInformationProblem("Steady thermal module");
      //DisplacementTime = new List<DoubleVector>();// (countDOF, 1);
      DoubleMatrix kGlobalDense = null; // stiffness
      DoubleVector uGlobal = null; // temperature
      DoubleVector rGlobal = null; // right hand side

      if (!IsSparseData)
      {
        AssemblyStiffnessMatrix(out kGlobalDense);
      }
      else
      {
        AssemblyStiffnessMatrix(ref kGlobalSparse);
      }
      AssemblyTractionVector(out rGlobal);
      AssemplyHeatGenerationSource(ref rGlobal);
      AssemblyHeatTransferConvectionVector(ref rGlobal);

      uGlobal = SolveStatic(kGlobalDense, kGlobalSparse, ref rGlobal);

      //DisplacementTime.Add(uGlobal);
      IOIGA.WriteMessagePackDoubleVector(FileNameDataTime, uGlobal);
    }
  }
}
