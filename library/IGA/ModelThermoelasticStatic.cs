using System;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public class ModelThermoelasticStatic : AbstractModelThermoelastic
  {
    public ModelThermoelasticStatic(Dimension dimension, string pathProject = "temp", string nameProject = "project")
         : base(TypeAnalysisModel.Static, dimension, pathProject, nameProject)
    {
    }

    public override void Solve()
    {
      if (listPatch.Count == 0)
        throw new NullReferenceException("Model must be assigned Patch");
      DoubleMatrix kGlobal = null;
      DoubleVector rGlobal = null;

      AssemblyStiffnessMatrix(out kGlobal);
      AssemblyTractionVector(out rGlobal);
      AssemblyHeatTransferConvectionVector(ref rGlobal);
      ApplyTAndg(ref kGlobal, ref rGlobal, null);
      //DisplacementTime = new List<DoubleVector>();
      DoubleVector u=MatrixFunctions.Solve(kGlobal, rGlobal);
      IOIGA.WriteMessagePackDoubleVector(FileNameDataTime, u);
    }

    public void AssemblyHeatTransferConvectionVector(ref DoubleVector rGlobal)
    {
      for (int i = 0; i < listPatch.Count; i++)
      {
        switch (StructureDimension)
        {
          case Dimension.Plane:
            ((PatchThermoelastic2D)listPatch[i]).ComputeHeatTransferConvectionPatch(ref rGlobal);
            break;
          case Dimension.Solid:
            ((PatchThermoelastic3D)listPatch[i]).ComputeHeatTransferConvectionPatch(ref rGlobal);
            break;
        }
      }
    }
  }
}