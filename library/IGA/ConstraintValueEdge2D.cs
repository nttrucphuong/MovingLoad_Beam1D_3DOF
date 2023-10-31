using DEMSoft.Function;
using System;
using DEMSoft.NURBS;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public class ConstraintValueEdge2D : AbstractConstraintValue
    {
        private AbstractPatch2D patch;
        private int indexEdge;
        public ConstraintValueEdge2D(AbstractPatch2D patch, int indexEdge, int fieldID, FunctionRToR piecewiseLoad, double valueConstraint)
           : base(fieldID, piecewiseLoad, valueConstraint)
        {
            this.patch = patch;
            this.indexEdge = indexEdge;
        }

        public ConstraintValueEdge2D(AbstractPatch2D patch, int indexEdge, int fieldID, FunctionRToR piecewiseLoad, FunctionRToR valueConstraint)
           : base(fieldID, piecewiseLoad, valueConstraint)
        {
            this.patch = patch;
            this.indexEdge = indexEdge;
        }
        public override ControlPoint[] GetControlPointsConstrainted()
        {
            return patch.SelectEndPatchControlPoints(indexEdge);
        }
        private int GetNumberOfControlPointsConstrainted()
        {
            return patch.GetSurface().Basis.GetCountBasisFunction(GetIndexCoordinate());
        }
        public override double[] ComputeValueConstraintOnControlPoints(double time)
        {
            double[] value = null;
            if ((GetValueConstraint() is FunctionRToR) && (!(GetValueConstraint() is ConstantFunctionRToR)))
            {
                DoubleMatrix LHS = ComputeLHS();
                DoubleVector RHS = ComputeRHS(time);
                value = MatrixFunctions.Solve(LHS, RHS).ToArray();
            }
            return value;
        }
        private DoubleMatrix ComputeLHS()
        {
            int countCPS = GetNumberOfControlPointsConstrainted();
            ControlPoint[] cps = GetControlPointsConstrainted();
            DoubleMatrix LHS = new DoubleMatrix(countCPS, countCPS);
            NURBSSurface sur = (NURBSSurface)patch.GetGeometry();
            var basis = (BivariateNURBSBasisFunction)sur.Basis;
            int indexCoordinate = GetIndexCoordinate();
            int mod = indexEdge % 2;
            double[] kv1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
            double[] kv2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
            int p = basis.GetDegree(0);
            int q = basis.GetDegree(1);
            int countLocalCPS = basis.GetDegree(indexCoordinate) + 1;
            for (int i = 0; i < kv1.Length - 1; i++)
            {
                DoubleMatrix LHSsmall = new DoubleMatrix(countLocalCPS, countLocalCPS);
                double[] paraEndPatchU = null;// kv1[i];
                double[] paraEndPatchV = null;// kv1[i + 1];
                if (indexCoordinate == 0)
                {
                    paraEndPatchU = new double[] { kv1[i], kv1[i + 1] };
                    if (mod == 0)
                    {
                        paraEndPatchV = new double[] { kv2[0], kv2[1] };
                    }
                    else
                    {
                        paraEndPatchV = new double[] { kv2[kv2.Length - 2], kv2[kv2.Length - 1] };
                    }
                }
                else
                {
                    paraEndPatchV = new double[] { kv2[i], kv2[i + 1] };
                    if (mod == 0)
                    {
                        paraEndPatchU = new double[] { kv1[0], kv1[1] };
                    }
                    else
                    {
                        paraEndPatchU = new double[] { kv1[kv1.Length - 2], kv1[kv1.Length - 1] };
                    }
                }
                double xi = 0;
                double eta = 0;
                int numGaussPoint = basis.GetDegree(indexCoordinate) + 1;
                for (int k = 0; k < numGaussPoint; k++)
                {
                    double psi = GaussPoints.GetPoint(numGaussPoint, k);
                    double w = GaussPoints.GetWeight(numGaussPoint, k);
                    switch (indexCoordinate)
                    {
                        case 0:
                            xi = 0.5 * ((paraEndPatchU[1] - paraEndPatchU[0]) * psi + paraEndPatchU[1] + paraEndPatchU[0]);
                            eta = paraEndPatchV[mod];
                            break;
                        case 1:
                            xi = paraEndPatchU[mod];
                            eta = 0.5 * ((paraEndPatchV[1] - paraEndPatchV[0]) * psi + paraEndPatchV[1] + paraEndPatchV[0]);
                            break;
                    }

                    var Nij = basis.GetValueBivariateBasisFunctions(xi, eta);
                    var dNij = basis.GetDerivativeBivariateBasisFunctions(xi, eta, 1);
                    DoubleVector Ni = new DoubleVector(countLocalCPS);
                    DoubleVector dNi = new DoubleVector(countLocalCPS);
                    for (int j = 0; j < countLocalCPS; j++)
                    {
                        switch (indexCoordinate)
                        {
                            case 0:
                                if (mod == 0)
                                {
                                    Ni[j] = Nij[j, 0];
                                    dNi[j] = dNij[j, 0][1, 0];
                                }
                                else
                                {
                                    Ni[j] = Nij[j, q];
                                    dNi[j] = dNij[j, q][1, 0];
                                }
                                break;
                            case 1:
                                if (mod == 0)
                                {
                                    Ni[j] = Nij[0, j];
                                    dNi[j] = dNij[0, j][0, 1];
                                }
                                else
                                {
                                    Ni[j] = Nij[p, j];
                                    dNi[j] = dNij[p, j][0, 1];
                                }
                                break;
                        }
                    }

                    double J2 = (indexCoordinate == 0) ? (0.5 * (paraEndPatchU[1] - paraEndPatchU[0]))
                          : (0.5 * (paraEndPatchV[1] - paraEndPatchV[0]));

                    int d = patch.GetCountDimension();
                    double[] grad = new double[d];
                    for (int j = 0; j < d; j++)
                    {
                        for (int kk = 0; kk < countLocalCPS; kk++)
                        {
                            grad[j] += dNi[kk] * cps[kk].GetCoordinate(j);
                        }
                    }

                    double J1 = Math.Sqrt(grad[0] * grad[0] + grad[1] * grad[1]);

                    LHSsmall += w * J1 * J2 * MatrixFunctions.OuterProduct(Ni, Ni);
                }
                ///// gan cai nho vo cai bu
                for (int ii = 0; ii < countLocalCPS; ii++)
                    for (int jj = 0; jj < countLocalCPS; jj++)
                    {
                        LHS[i + ii, i + jj] += LHSsmall[ii, jj];
                    }
            }
            return LHS;
        }
        private int GetIndexCoordinate()
        {
            return indexEdge / 2;
        }
        private DoubleVector ComputeRHS(double time)
        {
            int countCPS = GetNumberOfControlPointsConstrainted();
            ControlPoint[] cps = GetControlPointsConstrainted();
            DoubleVector RHS = new DoubleVector(countCPS);
            NURBSSurface sur = (NURBSSurface)patch.GetGeometry();
            var basis = (BivariateNURBSBasisFunction)sur.Basis;
            int indexCoordinate = GetIndexCoordinate();
            int mod = indexEdge % 2;
            double[] kv1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
            double[] kv2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
            int p = basis.GetDegree(0);
            int q = basis.GetDegree(1);
            int countLocalCPS = basis.GetDegree(indexCoordinate) + 1;
            for (int i = 0; i < kv1.Length - 1; i++)
            {
                DoubleVector RHSSmall = new DoubleVector(countLocalCPS);
                double[] paraEndPatchU = null;// kv1[i];
                double[] paraEndPatchV = null;// kv1[i + 1];
                if (indexCoordinate == 0)
                {
                    paraEndPatchU = new double[] { kv1[i], kv1[i + 1] };
                    if (mod == 0)
                    {
                        paraEndPatchV = new double[] { kv2[0], kv2[1] };
                    }
                    else
                    {
                        paraEndPatchV = new double[] { kv2[kv2.Length - 2], kv2[kv2.Length - 1] };
                    }
                }
                else
                {
                    paraEndPatchV = new double[] { kv2[i], kv2[i + 1] };
                    if (mod == 0)
                    {
                        paraEndPatchU = new double[] { kv1[0], kv1[1] };
                    }
                    else
                    {
                        paraEndPatchU = new double[] { kv1[kv1.Length - 2], kv1[kv1.Length - 1] };
                    }
                }
                double xi = 0;
                double eta = 0;
                int numGaussPoint = basis.GetDegree(indexCoordinate) + 1;
                for (int k = 0; k < numGaussPoint; k++)
                {
                    double psi = GaussPoints.GetPoint(numGaussPoint, k);
                    double w = GaussPoints.GetWeight(numGaussPoint, k);
                    double paraLoad = double.NaN;
                    switch (indexCoordinate)
                    {
                        case 0:
                            xi = 0.5 * ((paraEndPatchU[1] - paraEndPatchU[0]) * psi + paraEndPatchU[1] + paraEndPatchU[0]);
                            eta = paraEndPatchV[mod];
                            paraLoad = xi;
                            break;
                        case 1:
                            xi = paraEndPatchU[mod];
                            eta = 0.5 * ((paraEndPatchV[1] - paraEndPatchV[0]) * psi + paraEndPatchV[1] + paraEndPatchV[0]);
                            paraLoad = eta;
                            break;
                    }

                    var Nij = basis.GetValueBivariateBasisFunctions(xi, eta);
                    var dNij = basis.GetDerivativeBivariateBasisFunctions(xi, eta, 1);
                    DoubleVector Ni = new DoubleVector(countLocalCPS);
                    DoubleVector dNi = new DoubleVector(countLocalCPS);
                    for (int j = 0; j < countLocalCPS; j++)
                    {
                        switch (indexCoordinate)
                        {
                            case 0:
                                if (mod == 0)
                                {
                                    Ni[j] = Nij[j, 0];
                                    dNi[j] = dNij[j, 0][1, 0];
                                }
                                else
                                {
                                    Ni[j] = Nij[j, q];
                                    dNi[j] = dNij[j, q][1, 0];
                                }
                                break;
                            case 1:
                                if (mod == 0)
                                {
                                    Ni[j] = Nij[0, j];
                                    dNi[j] = dNij[0, j][0, 1];
                                }
                                else
                                {
                                    Ni[j] = Nij[p, j];
                                    dNi[j] = dNij[p, j][0, 1];
                                }
                                break;
                        }
                    }

                    double J2 = (indexCoordinate == 0) ? (0.5 * (paraEndPatchU[1] - paraEndPatchU[0]))
                          : (0.5 * (paraEndPatchV[1] - paraEndPatchV[0]));

                    int d = patch.GetCountDimension();
                    double[] grad = new double[d];
                    for (int j = 0; j < d; j++)
                    {
                        for (int kk = 0; kk < countLocalCPS; kk++)
                        {
                            grad[j] += dNi[kk] * cps[kk].GetCoordinate(j);
                        }
                    }
                    double J1 = Math.Sqrt(grad[0] * grad[0] + grad[1] * grad[1]);
                    RHSSmall += w * J1 * J2 * Ni * ((FunctionRToR)GetValueConstraint()).ValueAt(paraLoad) * GetPiecewiseLoadTime().ValueAt(time);
                }
                ///// gan cai nho vo cai bu
                for (int ii = 0; ii < countLocalCPS; ii++)
                {
                    RHS[i + ii] += RHSSmall[ii];
                }
            }
            return RHS;
        }
        public AbstractPatch2D GetPatch()
        { return patch; }
        public int GetIndexEdge()
        { return indexEdge; }
    }
}
