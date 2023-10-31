using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
    interface IElementThermal
    {
        void ComputeHeatGenerationSourceLoadVectorElement(ref DoubleVector rGlobal, double G);
    }
}
