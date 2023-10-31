using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
    interface IElementStructure
	{
		DoubleVector StressAt(params double[] xi);
		DoubleVector StrainAt(params double[] xi);
		DoubleVector StrainElasticAt(params double[] xi);
		DoubleVector StrainThermoAt(params double[] xi);
	}
}
