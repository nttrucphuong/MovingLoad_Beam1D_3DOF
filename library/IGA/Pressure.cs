using CenterSpace.NMath.Core;
using DEMSoft.NURBS;

namespace DEMSoft.IGA
{
    /// <summary>
    /// Pressure load on edge of 2D patch
    /// </summary>
    public abstract class Pressure : AbstractLoad
    {
        protected bool isInGlobal;
        protected bool isNaturalCoordinate;

        /// <summary>
        /// Constructor class
        /// </summary>
        /// <param name="mp">Mesh part which be applied load</param>
        /// <param name="isInGlobal">Direction coresponding on global coordinate (x-y) or local coordinate (n-t)</param>
        /// <param name="isNaturalCoordinate">Value coresponding on natural coordinate or physical coordinate</param>
        public Pressure(AbstractMeshPart mp, bool isInGlobal, bool isNaturalCoordinate)
                : base(mp)
        {
            this.isInGlobal = isInGlobal;
            this.isNaturalCoordinate = isNaturalCoordinate;
        }

        /// <summary>
        /// Compute local load vector
        /// </summary>
        /// <returns></returns>
        public abstract override DoubleVector ComputeLocalLoadVector(double time);
    }
}
