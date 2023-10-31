using DEMSoft.NURBS;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
    public class PatchThermal3D : AbstractPatch3D
	{
		public PatchThermal3D(NURBSVolume volume)
				  : base(volume, 1)
		{

		}


		public void ComputeHeatTransferConvectionPatch(ref DoubleVector rGlobal)
		{
			int nel = CountElements();//(n - p) * (m - q);//number of elements
												 //if (!IsParallelProcesing)
												 //{
			for (int i = 0; i < nel; i++)//Loop through elements
			{
				ElementThermal3D elem = (ElementThermal3D)listElement[i];
				if (elem.GetHeatTransferConvection() != null)
				{
					elem.ComputeHeatTransferConvectionLoadVectorElement(ref rGlobal);
				}
			}
		}
	}
}


