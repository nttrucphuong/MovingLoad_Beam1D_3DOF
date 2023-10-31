using DEMSoft.NURBS;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public class PatchThermal2D : AbstractPatch2D
  {
    public PatchThermal2D(Abstract2DParametricGeometry surface)
       : base(surface, TypeStructure.Plane)
    {
    }
    public void ComputeHeatTransferConvectionPatch(ref DoubleVector rGlobal)
    {
      int nel = CountElements();//(n - p) * (m - q);//number of elements
      for (int i = 0; i < nel; i++)//Loop through elements
      {
        ElementThermal2D elem = (ElementThermal2D)listElement[i];
        if (elem.GetHeatTransferConvection() != null)
        {
          elem.ComputeHeatTransferConvectionLoadVectorElement(ref rGlobal);
        }
      }
    }
  }
}
