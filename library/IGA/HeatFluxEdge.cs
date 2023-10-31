using System;
using DEMSoft.NURBS;
using DEMSoft.Function;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
    /// <summary>
    /// Pressure load on edge of 2D patch
    /// </summary>
    public class HeatFluxEdge : HeatFlux
    {
        private FunctionRToR press;

        /// <summary>
        /// Constructor class
        /// </summary>
        /// <param name="e">edge of element which be applied load</param>
        /// <param name="isInGlobal">Coresponding on global coordinate (x-y) or local coordinate (n-t)</param>
        /// <param name="p">global coordinate (x-y) : p[0] = px, p[1] = py,  local coordinate (n-t) : p[0] = pn, p[1] = pt</param>
        public HeatFluxEdge(Edge e, FunctionRToR p)
            : base(e, true)
        { press = p; }

        /// <summary>
        /// Constructor class
        /// </summary>
        /// <param name="e">edge of element which be applied load</param>
        /// <param name="isNaturalCoordinate">Coresponding on global coordinate (x-y) or local coordinate (n-t)</param>
        /// <param name="p">global coordinate (x-y) : p[0] = px, p[1] = py,  local coordinate (n-t) : p[0] = pn, p[1] = pt</param>
        public HeatFluxEdge(Edge e, bool isNaturalCoordinate, FunctionRToR p)
            : base(e, isNaturalCoordinate)
        { press = p; }

        /// <summary>
        /// Compute local load vector
        /// </summary>
        /// <returns></returns>
        public override DoubleVector ComputeLocalLoadVector(double time = 1)
        {
            Edge e = (Edge)GetMeshPart();
            int d = 1;
            int p = e.GetDegree(0);
            int n = d * (p + 1);
            DoubleVector re = new DoubleVector(n);

            var patch = e.GetFace().GetElement();
            var mesh = patch.GetPatch();
            var surface = (NURBSSurface)(mesh.GetGeometry(0));
            var basis = surface.Basis;
            var paraEndPatchU = e.GetFace().GetElement().GetParameterTwoEndElement(0);
            var paraEndPatchV = e.GetFace().GetElement().GetParameterTwoEndElement(1);
            double xi = 0;
            double eta = 0;
            double[] paraEndEdge = e.GetParametricEndEdge(0);
            int numGaussPoint = patch.GetNumberOfGaussPointOnEachDirection();
            for (int k = 0; k < numGaussPoint; k++)
            {
                double psi = GaussPoints.GetPoint(numGaussPoint, k);
                double w = GaussPoints.GetWeight(numGaussPoint, k);
                switch (e.GetIndexCoordinate())
                {
                    case 0:
                        xi = 0.5 * ((paraEndPatchU[1] - paraEndPatchU[0]) * psi + paraEndPatchU[1] + paraEndPatchU[0]);
                        eta = paraEndPatchV[e.GetIndexFrontBack()];
                        break;
                    case 1:
                        xi = paraEndPatchU[e.GetIndexFrontBack()];
                        eta = 0.5 * ((paraEndPatchV[1] - paraEndPatchV[0]) * psi + paraEndPatchV[1] + paraEndPatchV[0]);
                        break;
                }
                double paraLoad;
                if (e.GetIndexCoordinate() == 0)
                    paraLoad = xi;
                else
                    paraLoad = eta;
                DoubleVector Ni = new DoubleVector(e.GetBivariateBasisFunctionOnEdge(0, xi, eta));
                double J2 = (e.GetIndexCoordinate() == 0) ? (0.5 * (paraEndPatchU[1] - paraEndPatchU[0]))
                    : (0.5 * (paraEndPatchV[1] - paraEndPatchV[0]));

                DoubleVector dNi = new DoubleVector(e.GetDerivativeBivariateBasisFunctionOnEdge(0, xi, eta));
                double[] J = JacobianAt(dNi).ToArray();
                double J1 = Math.Sqrt(J[0] * J[0] + J[1] * J[1]);

                re += -w * J1 * J2 * Ni * press.ValueAt(paraLoad);//x y
            }
            return re;
        }

        private DoubleVector JacobianAt(DoubleVector dNi)
        {
            int d = GetMeshPart().CountDimension();
            Edge e = (Edge)GetMeshPart();
            int p = e.GetDegree(0);
            int n = p + 1;
            DoubleVector grad = new DoubleVector(d);
            var cpsOnEdge = e.GetControlPointsOnEdge(0);
            for (int j = 0; j < d; j++)
            {
                for (int k = 0; k < n; k++)
                {
                    grad[j] += dNi[k] * cpsOnEdge[k].GetCoordinate(j);
                }
            }
            return grad;
        }
    }
}
