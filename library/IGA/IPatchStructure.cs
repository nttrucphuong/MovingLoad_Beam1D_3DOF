using DEMSoft.Common;

namespace DEMSoft.IGA
{
  interface IPatchStructure
    {
        double StressAt(Result re, params double[] xi);
        double StrainAt(Result re, params double[] xi);
        double StressOnGlobalAt(Result re, TypeCoordination coord, params double[] x);
        double StrainOnGlobalAt(Result re, TypeCoordination coord, params double[] x);
        Structure2DState StateStress
        { get; set; }
    }
}
