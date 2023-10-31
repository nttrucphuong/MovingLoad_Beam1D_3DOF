using CenterSpace.NMath.Core;
using DEMSoft.NURBS;

namespace DEMSoft.IGA
{
    public abstract class AbstractLoad
    {
        private AbstractMeshPart mp;

        public AbstractLoad(AbstractMeshPart mp)
        {
            this.mp = mp;
        }

        public abstract DoubleVector ComputeLocalLoadVector(double t);

        public AbstractMeshPart GetMeshPart()
        {
            return mp;
        }
    }
}
