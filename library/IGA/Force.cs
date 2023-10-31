using DEMSoft.NURBS;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
    public class Force : AbstractLoad
    {
        private double[] f;

        public Force(ControlPoint cp, params double[] f)
            : base(cp)
        {
            this.f = f;
        }

        public override DoubleVector ComputeLocalLoadVector(double time = 1)
        {
            //int countDimensions = GetMeshPart().CountDimension();
            //DoubleVector re = new DoubleVector(countDimensions);
            //for (int i = 0; i < countDimensions; i++)
            //{
            //    re[i] = f[i];
            //}
            int countFields = GetMeshPart().GetNumberOfFields();
            DoubleVector re = new DoubleVector(countFields);
            for (int i = 0; i < countFields; i++)
            {
                re[i] = f[i];
            }

            return re;
        }

        public ControlPoint GetControlPoint()
        {
            return (ControlPoint)GetMeshPart();
        }

    }
}
