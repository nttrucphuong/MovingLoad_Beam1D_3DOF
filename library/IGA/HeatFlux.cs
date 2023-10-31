using CenterSpace.NMath.Core;
using DEMSoft.NURBS;

namespace DEMSoft.IGA
{
    /// <summary>
    /// Pressure load on edge of 2D patch
    /// </summary>
    public abstract class HeatFlux : AbstractLoad
    {
        protected bool isNaturalCoordinate;

        /// <summary>
        /// Constructor class
        /// </summary>
        /// <param name="mp">Mesh part which be applied load</param>
        /// <param name="isInGlobal">Coresponding on global coordinate (x-y) or local coordinate (n-t)</param>
        /// <param name="isNaturalCoordinate">Coresponding on natural coordinate or physical coordinate</param>
        public HeatFlux(AbstractMeshPart mp, bool isNaturalCoordinate)
            : base(mp)
        {
            this.isNaturalCoordinate = isNaturalCoordinate;
        }

        /// <summary>
        /// Compute local load vector
        /// </summary>
        /// <returns></returns>
        public abstract override DoubleVector ComputeLocalLoadVector(double time);
    }
}
