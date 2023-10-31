using DEMSoft.Common;

namespace DEMSoft.IGA
{

  /// <summary>
  /// Structure 2D problem (plane stress, plane strain). Default is plane stress
  /// </summary>
  public abstract class AbstractModelThermoelastic : AbstractModel, IModelStructure
  {
    public Structure2DState StressState
    { get; set; }
    /// <summary>
    /// Constructor class
    /// </summary>
    /// <param name="problem">Define state of stress, plane stress or plane strain</param>
    public AbstractModelThermoelastic(TypeAnalysisModel type, Dimension dimension, string pathProject = "temp", string nameProject = "project")
                                            : base(TypeModelProblem.StructuralThermal, type, dimension, pathProject, nameProject)
    {
      if (StructureDimension == Dimension.Plane)
      {
        listComputeResult.Add(Result.UX);
        listComputeResult.Add(Result.UY);
        listComputeResult.Add(Result.TEMP);
      }
      else if (StructureDimension == Dimension.Solid)
      {
        listComputeResult.Add(Result.UX);
        listComputeResult.Add(Result.UY);
        listComputeResult.Add(Result.UZ);
        listComputeResult.Add(Result.TEMP);
      }
    }
  }
}