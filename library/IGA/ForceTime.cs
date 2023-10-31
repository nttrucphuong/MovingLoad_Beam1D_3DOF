using DEMSoft.NURBS;
using CenterSpace.NMath.Core;
using DEMSoft.Function;

namespace DEMSoft.IGA
{
    public class ForceTime : AbstractLoad
    {
        private double[] f;
        private FunctionRToR piecewiseLoad;

        public ForceTime(ControlPoint cp, FunctionRToR piecewiseLoad, params double[] f)
                : base(cp)
        {
            this.f = f;
            this.piecewiseLoad = piecewiseLoad;
        }

        public override DoubleVector ComputeLocalLoadVector(double time)
        {
            int countDimensions = GetMeshPart().CountDimension();
            DoubleVector re = new DoubleVector(countDimensions);
            double loadTime = piecewiseLoad.ValueAt(time);
            for (int i = 0; i < countDimensions; i++)
            {
                re[i] = loadTime * f[i];
            }
            return re;
        }

        public ControlPoint GetControlPoint()
        {
            return (ControlPoint)GetMeshPart();
        }

    }
}
