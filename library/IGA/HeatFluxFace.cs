using System;
using DEMSoft.NURBS;
using DEMSoft.Function;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
    public class HeatFuxFace : HeatFlux
    {
        private FunctionR2ToR press;

        /// <summary>
        /// Constructor class
        /// </summary>
        /// <param name="face">face of element which be applied load</param>
        /// <param name="isInGlobal">Coresponding on global coordinate (x-y) or local coordinate (n-t)</param>
        /// <param name="p">global coordinate (x-y) : p[0] = px, p[1] = py,  local coordinate (n-t) : p[0] = pn, p[1] = pt</param>
        public HeatFuxFace(Face face, bool isInGlobal, FunctionR2ToR p)
            : base(face, true)
        { press = p; }

        /// <summary>
        /// Constructor class
        /// </summary>
        /// <param name="face">face of element which be applied load</param>
        /// <param name="isInGlobal">Coresponding on global coordinate (x-y) or local coordinate (n-t)</param>
        /// <param name="p">global coordinate (x-y) : p[0] = px, p[1] = py,  local coordinate (n-t) : p[0] = pn, p[1] = pt</param>
        public HeatFuxFace(Face face, bool isInGlobal, bool isNaturalCoordinate, FunctionR2ToR p)
            : base(face, isNaturalCoordinate)
        { press = p; }

        /// <summary>
        /// Compute local load vector
        /// </summary>
        /// <returns></returns>
        public override DoubleVector ComputeLocalLoadVector(double time = 1)
        {
            int d = 1;
            Face face = (Face)GetMeshPart();
            int p = face.GetDegree(0);
            int q = face.GetDegree(1);
            int n = d * (p + 1) * (q + 1);
            DoubleVector re = new DoubleVector(n);

            var patch = face.GetVolumeBeAttached().GetElement();
            var mesh = patch.GetPatch();
            var vol = (NURBSVolume)(mesh.GetGeometry(0));
            var paraEndPatchU = face.GetVolumeBeAttached().GetElement().GetParameterTwoEndElement(0);
            var paraEndPatchV = face.GetVolumeBeAttached().GetElement().GetParameterTwoEndElement(1);
            var paraEndPatchW = face.GetVolumeBeAttached().GetElement().GetParameterTwoEndElement(2);
            double xi = 0;
            double eta = 0;
            double zeta = 0;
            //double[,] paraEndFace = face.GetParametricEndEdge();
            int numGaussPoint = Math.Max(p + 1, q + 1);
            for (int k = 0; k < numGaussPoint; k++)
            {
                double psi1 = GaussPoints.GetPoint(numGaussPoint, k);
                double w1 = GaussPoints.GetWeight(numGaussPoint, k);
                for (int kk = 0; kk < numGaussPoint; kk++)
                {
                    double psi2 = GaussPoints.GetPoint(numGaussPoint, kk);
                    double w2 = GaussPoints.GetWeight(numGaussPoint, kk);
                    switch (face.GetIndexCoordinate())
                    {
                        case 0:
                            xi = 0.5 * ((paraEndPatchU[1] - paraEndPatchU[0]) * psi1 + paraEndPatchU[1] + paraEndPatchU[0]);
                            eta = 0.5 * ((paraEndPatchV[1] - paraEndPatchV[0]) * psi2 + paraEndPatchV[1] + paraEndPatchV[0]);
                            zeta = paraEndPatchW[face.GetIndexFrontBack()];
                            break;
                        case 1:
                            xi = 0.5 * ((paraEndPatchU[1] - paraEndPatchU[0]) * psi1 + paraEndPatchU[1] + paraEndPatchU[0]);
                            eta = paraEndPatchV[face.GetIndexFrontBack()];
                            zeta = 0.5 * ((paraEndPatchW[1] - paraEndPatchW[0]) * psi2 + paraEndPatchW[1] + paraEndPatchW[0]);
                            break;
                        case 2:
                            xi = paraEndPatchU[face.GetIndexFrontBack()];
                            eta = 0.5 * ((paraEndPatchV[1] - paraEndPatchV[0]) * psi1 + paraEndPatchV[1] + paraEndPatchV[0]);
                            zeta = 0.5 * ((paraEndPatchW[1] - paraEndPatchW[0]) * psi2 + paraEndPatchW[1] + paraEndPatchW[0]);
                            break;
                    }

                    double[] paraLoad = { 0, 0 };

                    if (isNaturalCoordinate)
                    {
                        switch (face.GetIndexCoordinate())
                        {
                            case 0:
                                paraLoad[0] = xi;
                                paraLoad[1] = eta;
                                break;
                            case 1:
                                paraLoad[0] = xi;
                                paraLoad[1] = zeta;
                                break;
                            case 2:
                                paraLoad[0] = eta;
                                paraLoad[1] = zeta;
                                break;
                        }
                    }
                    //else
                    //{
                    //    double[] pp = vol.PointAt(xi, eta, zeta);
                    //    switch (face.GetIndexCoordinate())
                    //    {
                    //        case 0:
                    //            paraLoad[0] = pp[0];
                    //            paraLoad[1] = pp[1];
                    //            break;
                    //        case 1:
                    //            paraLoad[0] = pp[0];
                    //            paraLoad[1] = pp[2];
                    //            break;
                    //        case 2:
                    //            paraLoad[0] = pp[1];
                    //            paraLoad[1] = pp[2];
                    //            break;
                    //    }
                    //}

                    double[,] N = face.GetTrivariateBasisFunctionOnFace(xi, eta, zeta);
                    int num = (N.GetLength(0)) * (N.GetLength(1));
                    DoubleVector Nij = new DoubleVector(num);
                    int dem = 0;
                    for (int j = 0; j < N.GetLength(1); j++)
                        for (int i = 0; i < N.GetLength(0); i++)
                        {
                            Nij[dem] = N[i, j];
                            dem += 1;
                        }

                    double J2 = 0;

                    switch (face.GetIndexCoordinate())
                    {
                        case 0:
                            J2 = 0.25 * (paraEndPatchU[1] - paraEndPatchU[0]) * (paraEndPatchV[1] - paraEndPatchV[0]);
                            break;
                        case 1:
                            J2 = 0.25 * (paraEndPatchU[1] - paraEndPatchU[0]) * (paraEndPatchW[1] - paraEndPatchW[0]);
                            break;
                        case 2:
                            J2 = 0.25 * (paraEndPatchV[1] - paraEndPatchV[0]) * (paraEndPatchW[1] - paraEndPatchW[0]);
                            break;
                    }

                    var elem = face.GetVolumeBeAttached().GetElement();
                    DoubleMatrix gradNE = elem.GradBasisFunction(xi, eta, zeta);
                    DoubleMatrix J = elem.JacobianAt(gradNE);

                    //DoubleMatrix dNij = new DoubleMatrix(face.GetDerivativeTrivariateBasisFunctionOnFace(xi, eta, zeta));
                    //DoubleMatrix J = JacobianAt(dNij);


                    //Jacobian of face mapping
                    double normJ = Math.Abs(normVectorJacobian(J));

                    re += w1 * w2 * normJ * J2 * Nij * press.ValueAt(paraLoad[0], paraLoad[1]);//x y
                }
            }
            return re;
        }

        public double normVectorJacobian(DoubleMatrix J)
        {
            Face face = (Face)GetMeshPart();

            double ee = 0;
            double gg = 0;
            double ff = 0;

            switch (face.GetIndexCoordinate())
            {
                case 0:
                    ee = J[0, 0] * J[0, 0] + J[1, 0] * J[1, 0] + J[2, 0] * J[2, 0];
                    gg = J[0, 1] * J[0, 1] + J[1, 1] * J[1, 1] + J[2, 1] * J[2, 1];
                    ff = J[0, 0] * J[0, 1] + J[1, 0] * J[1, 1] + J[2, 0] * J[2, 1];
                    break;
                case 1:
                    ee = J[0, 0] * J[0, 0] + J[1, 0] * J[1, 0] + J[2, 0] * J[2, 0];
                    gg = J[0, 2] * J[0, 2] + J[1, 2] * J[1, 2] + J[2, 2] * J[2, 2];
                    ff = J[0, 0] * J[0, 2] + J[1, 0] * J[1, 2] + J[2, 0] * J[2, 2];
                    break;
                case 2:
                    ee = J[0, 1] * J[0, 1] + J[1, 1] * J[1, 1] + J[2, 1] * J[2, 1];
                    gg = J[0, 2] * J[0, 2] + J[1, 2] * J[1, 2] + J[2, 2] * J[2, 2];
                    ff = J[0, 1] * J[0, 2] + J[1, 1] * J[1, 2] + J[2, 1] * J[2, 2];
                    break;
            }
            return Math.Sqrt(ee * gg - ff * ff);
        }

        //private DoubleMatrix JacobianAt(DoubleMatrix dNij)
        //{
        //    int d = GetMeshPart().GetNumberOfFields() - 1;
        //    Face face = (Face)GetMeshPart();
        //    int p1 = face.GetDegree(0);
        //    int p2 = face.GetDegree(1);
        //    int nen = (p1 + 1) * (p2 + 1);
        //    DoubleMatrix grad = new DoubleMatrix(d + 1, d);
        //    var cpsOnFace = face.GetControlPointsOnFace();
        //    var inc = face.CreateIndexCoordinate();
        //    for (int j = 0; j < d; j++)
        //        for (int k = 0; k < nen; k++)
        //        {
        //            var controlPointCoord = cpsOnFace[inc[k, 0], inc[k, 1]].GetCoordinate();
        //            for (int i = 0; i <= d; i++)
        //            {
        //                grad[i, j] += dNij[k, j] * controlPointCoord[i];
        //            }
        //        }
        //    return grad;
        //}
    }
}