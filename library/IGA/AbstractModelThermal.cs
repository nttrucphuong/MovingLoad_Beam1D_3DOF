using System.Collections.Generic;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Abstract thermal model
  /// </summary>
  public abstract class AbstractModelThermal : AbstractModel
  {

    protected List<HeatGenerationSourceBody> heatSource;
    /// <summary>
    /// Constructor class
    /// </summary>
    public AbstractModelThermal(TypeAnalysisModel type, Dimension structureDimension, string pathProject = "temp", string nameProject = "project")
        : base(TypeModelProblem.Thermal, type, structureDimension, pathProject, nameProject)
    {
      heatSource = new List<HeatGenerationSourceBody>();
      listComputeResult.Add(Result.TEMP);
      typeOfFieldsInMultifield.Add(TypeFields.Thermal);
    }

    public void AssemblyHeatTransferConvectionVector(ref DoubleVector rGlobal)
    {
      for (int i = 0; i < listPatch.Count; i++)
      {
        switch (StructureDimension)
        {
          case Dimension.Plane:
            ((PatchThermal2D)listPatch[i]).ComputeHeatTransferConvectionPatch(ref rGlobal);
            break;
          case Dimension.Solid:
            ((PatchThermal3D)listPatch[i]).ComputeHeatTransferConvectionPatch(ref rGlobal);
            break;
        }
      }
    }

    /// <summary>
    /// Add heat source
    /// </summary>
    /// <param name="source">heat generation source</param>
    public void AddHeatSource(HeatGenerationSourceBody source)
    {
      heatSource.Add(source);
    }

    /// <summary>
    /// Get heat source
    /// </summary>
    /// <param name="index">index</param>
    /// <returns></returns>
    public HeatGenerationSourceBody GetHeatSource(int index)
    {
      return heatSource[index];
    }

    /// <summary>
    /// Count number of heat generation source
    /// </summary>
    /// <returns></returns>
    public int CountHeatSource()
    {
      return heatSource.Count;
    }
    public void AssemplyHeatGenerationSource(ref DoubleVector rGlobal)
    {
      for (int i = 0; i < CountHeatSource(); i++)
      {
        HeatGenerationSourceBody load = heatSource[i];
        load.ComputeLocalLoadVector(ref rGlobal);
      }
    }
  }
}
