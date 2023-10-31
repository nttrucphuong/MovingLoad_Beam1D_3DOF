using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
    public class HeatGenerationSourceBody
    {
        private double G;
        private AbstractPatch patch;

        public HeatGenerationSourceBody(AbstractPatch patch, double G)
        {
            this.G = G;
            this.patch = patch;
        }

        public AbstractPatch GetPatch()
        { return patch; }

        public double GetValueHeatGenerationSource()
        { return G; }

        public virtual void ComputeLocalLoadVector(ref DoubleVector rGlobal)
        {
            int nel = patch.CountElements();//(n - p) * (m - q);//number of elements
            //if (!IsParallelProcesing)
            //{
            for (int i = 0; i < nel; i++)//Loop through elements
            {
                IElementThermal elem = (IElementThermal)patch.GetElement(i);
                elem.ComputeHeatGenerationSourceLoadVectorElement(ref rGlobal, G);
            }
            //}
        }
    }
}
