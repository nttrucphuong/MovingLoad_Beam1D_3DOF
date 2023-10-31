using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using DEMSoft.Common;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Abstract patch (as element in IGA)
  /// </summary>
  public abstract class AbstractElement
  {
    /// <summary>
    /// Patch of problem
    /// </summary>
    protected AbstractPatch patch;
    /// <summary>
    /// ID of patch
    /// </summary>
    protected int id;

    /// <summary>
    /// Index of element in patch
    /// </summary>
    protected int[] index;

    protected int numberOfGaussPointOnEachDirection = -1;

    public Material Material
    { get; set; }

    protected DoubleMatrix Ke;
    protected DoubleMatrix KGe;
    protected DoubleMatrix Me;
    protected DoubleVector fi;
    public DoubleMatrix GetStiffnessMatrix()
    { return Ke; }
    public DoubleMatrix GetGeometricStiffnessMatrix()
    { return KGe; }
    public void SetGeometricStiffnessMatrix(DoubleMatrix KGe)
    { this.KGe = KGe; }
    public void SetStiffnessMatrix(DoubleMatrix Ke)
    { this.Ke = Ke; }

    public DoubleMatrix GetMassMatrix()
    { return Me; }

    public void SetMassMatrix(DoubleMatrix Me)
    { this.Me = Me; }

    public DoubleVector GetInternalForce()
    { return fi; }

    public void SetInternalForce(DoubleVector fi)
    { this.fi = fi; }
    /// <summary>
    /// Contructor class
    /// </summary>
    /// <param name="patch">mesh of problem</param>
    /// <param name="idElement">ID of patch</param>
    public AbstractElement(AbstractPatch patch, int idElement)
    {
      this.patch = patch;
      this.id = idElement;
      //this.index = index;
    }

    ///// <summary>
    ///// Compute local stiffness matrix
    ///// </summary>
    ///// <returns></returns>
    //public abstract DoubleMatrix ComputeLocalStiffnessMatrix();

    /// <summary>
    /// Get mesh of problem
    /// </summary>
    /// <returns></returns>
    public AbstractPatch GetPatch()
    {
      return patch;
    }

    /// <summary>
    /// Get index of patch in global coordinate
    /// </summary>
    /// <param name="direction"></param>
    /// <returns></returns>
    public int GetIndexGlobalCoordinate(int direction)
    {
      return index[direction];
    }

    /// <summary>
    /// Get ID of patch
    /// </summary>
    /// <returns></returns>
    public int GetID()
    {
      return id;
    }

    /// <summary>
    /// Get two parametric on end points on patch coressponding index of coordinate
    /// </summary>
    /// <param name="indexCoordinate">index of coordinate : 0 - xi coordinate, 1 - eta coordinate</param>
    /// <returns>two parametric on end points</returns>
    public double[] GetParameterTwoEndElement(int indexCoordinate)
    {
      var basis = patch.GetGeometry(0).Basis;
      var kvNoMulticiply = basis.GetKnotVector(indexCoordinate).GetKnotVectorNoMultiplicity();
      int idx = patch.GetIPN(id, indexCoordinate);
      return new double[] { kvNoMulticiply[idx], kvNoMulticiply[idx + 1] };
    }

    /// <summary>
    /// Check least part of patch is inside region
    /// </summary>
    /// <param name="loc">Location (Region class)</param>
    /// <returns></returns>
    public abstract bool IsInsideRegion(IRegion loc);

    /// <summary>
    /// Check full of patch is inside region
    /// </summary>
    /// <param name="loc">Location (Region class)</param>
    /// <returns></returns>
    public abstract bool IsCompleteInsideRegion(IRegion loc);

    /// <summary>
    /// Get number of dimension
    /// </summary>
    /// <returns>number of dimension</returns>
    public int GetCountDimension()
    {
      return patch.GetCountDimension();
    }

    public int GetNumberOfGaussPointOnEachDirection()
    { return numberOfGaussPointOnEachDirection; }

    public abstract void InitializeGaussPoint();
    public abstract void ComputeMaterialPropertyValueAtGaussPoint(MaterialPropertyName name);
    public abstract void ComputeDrawValueAtGaussPoint(DataInGausspoint name);
    public abstract void ComputeValueAtGaussPoint(DataInGausspoint name);
    public abstract void ComputeStiffnessMatrixElement(out DoubleMatrix Ke);
    public abstract void ComputeGeometricStiffnessMatrixElement(double[,] stressTensor, out DoubleMatrix Ke);
    public abstract void ComputeMassMatrixElement(out DoubleMatrix Me);
    public abstract void ComputeInternalForceElement(out DoubleVector fi);
    public virtual double ComputeSharpCrackSurface()
    { return 0; }
    public abstract int[] GetTArrayGlobal();
    /// <summary>
    /// Get u vector on one element
    /// </summary>
    /// <returns></returns>
    public abstract DoubleVector GetDisplacementLocal();
    public abstract ControlPoint[] GetControlPointsLocal();
    public abstract DoubleVector GetEachVariableLocal(Result variable);
    public abstract DoubleMatrix GradBasisFunction(params double[] xi);
    public abstract DoubleVector ValueBasisFunction(params double[] xi);
    public abstract DoubleMatrix JacobianAt(DoubleMatrix dNdxi);
    internal virtual void UpdateGausspointValue(){ }
  }
}
