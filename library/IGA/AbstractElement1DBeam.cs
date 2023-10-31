using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
    /// <summary>
    /// Abstract 2D element class. It uses only for 2D problem.
    /// </summary>
    public abstract class AbstractElement1DBeam : AbstractElement
    {
        protected Edge edge;
        protected GaussPoints[] gps;

        public AbstractElement1DBeam(AbstractPatch1D patch, int id)
          : base(patch, id)
        {
            edge = new Edge(this);
        }

        /// <summary>
        /// 
        /// </summary>
        public override void InitializeGaussPoint()
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            NURBSBasisFunction basis = (NURBSBasisFunction)((NURBSCurve)(patch.GetGeometry())).Basis;
            var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
            int idx = patch.GetIPN(id, 0);
            var knotVectorNoMultiplicity1 = basis.GetKnotVector().GetKnotVectorNoMultiplicity();
            int nel = knotVectorNoMultiplicity1.Length - 1;
            int nen = patch.GetCountLocalBasisFunctions();// (p + 1) ; // number of local basis functions
            int p = basis.GetDegree(0);

            numberOfGaussPointOnEachDirection = p + 1;//(int)Math.Ceiling((2.0 * Math.Max(p, q) + 1.0) / 2.0);

            gps = new GaussPoints[numberOfGaussPointOnEachDirection];

            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
            {
                double xi = GaussPoints.GetPoint(numberOfGaussPointOnEachDirection, i);
                double w = GaussPoints.GetWeight(numberOfGaussPointOnEachDirection, i);
                xi = 0.5 * ((kvNoMulticiply1[idx + 1] - kvNoMulticiply1[idx]) * xi + kvNoMulticiply1[idx + 1] + kvNoMulticiply1[idx]);
                int span1 = basis.FindSpan(xi, 0);

                GaussPoints gpsij = gps[i] = new GaussPoints(new double[] { xi }, w);
                gpsij.listDataGausspoint = new System.Collections.Generic.List<DataGausspoint>();
                DoubleVector Ni = new DoubleVector(nen);
                DoubleMatrix dNdxi = new DoubleMatrix(nen, 1);
                DoubleMatrix ddNdxi = new DoubleMatrix(nen, 1);

                double[,] gradBasis = basis.GetDerivativeBasisFunctions(xi, 2);
                for (int k = 0; k < nen; k++)
                {
                    int ien = patch.GetIEN(id, k);
                    Ni[k] = gradBasis[0, patch.GetINC(ien, 0) - span1 + p];
                    dNdxi[k, 0] = gradBasis[1, patch.GetINC(ien, 0) - span1 + p];
                    ddNdxi[k, 0] = gradBasis[2, patch.GetINC(ien, 0) - span1 + p];
                }

                DoubleMatrix J = JacobianAt(dNdxi);
                DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
                DoubleMatrix dNdx = MatrixFunctions.Product(dNdxi, invertJ);//dNdxi * invertJ;
                gpsij.SetValue(DataInGausspoint.Ni, Ni);
                gpsij.SetValue(DataInGausspoint.dNdxi, dNdxi);
                double detJ = MatrixFunctions.Determinant(J);
                gpsij.SetValue(DataInGausspoint.detJ, MatrixFunctions.Determinant(J));
                gpsij.SetValue(DataInGausspoint.dNdX, dNdx);
                gpsij.SetValue(DataInGausspoint.ddNdxi, ddNdxi);
                DoubleMatrix jac2 = new DoubleMatrix(1, 1);
                jac2[0, 0] = J[0, 0] * J[0, 0];
                DoubleMatrix jac1 = Jacobian2At(ddNdxi);
                DoubleMatrix term1 = MatrixFunctions.Product(dNdx, jac1.Transpose());
                DoubleMatrix term2 = ddNdxi - term1;
                DoubleMatrix ddNdX = MatrixFunctions.Product(MatrixFunctions.Inverse(jac2), term2.Transpose());
                gpsij.SetValue(DataInGausspoint.ddNdX, ddNdX.Transpose());
                //gps[i, j].ddNdX = ddNdX.Transpose();//bo sung


                if (Material is PlasticityMaterial)
                {
                    gpsij.SetValue(DataInGausspoint.lastAlpha, 0.0);
                    gpsij.SetValue(DataInGausspoint.lastBackStress, new DoubleVector(2));
                    gpsij.SetValue(DataInGausspoint.lastPlasticStrain, new DoubleVector(2));
                    gpsij.SetValue(DataInGausspoint.lastStress, new DoubleVector(2));
                }
                else if (Material is FGMStructureOneVariableMaterial)
                {
                    MaterialProperty matE = Material.GetProperty(MaterialPropertyName.YoungModulus);
                    MaterialProperty matNu = Material.GetProperty(MaterialPropertyName.PoissonRatio);
                    MaterialProperty matK = Material.GetProperty(MaterialPropertyName.IsotropicThermalConductivity);
                    gpsij.SetValue(DataInGausspoint.EModulus, (matE != null) ? ComputeParameterProperty(MaterialPropertyName.YoungModulus, xi) : 0.0);
                    gpsij.SetValue(DataInGausspoint.nu, (matNu != null) ? ComputeParameterProperty(MaterialPropertyName.PoissonRatio, xi) : 0.0);
                    gpsij.SetValue(DataInGausspoint.ThermalConductivity, (matK != null) ? ComputeParameterProperty(MaterialPropertyName.IsotropicThermalConductivity, xi) : 0.0);

                    MaterialProperty matEpsilon = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion);
                    gpsij.SetValue(DataInGausspoint.CoefficientsThermalExpansion, (matEpsilon != null) ? ComputeParameterProperty(MaterialPropertyName.CoefficientThermalExpansion, xi) : 0.0);
                }

            }
        }

        public double ComputeParameterProperty(MaterialPropertyName name, double xi)
        {
            double v;
            if (Material is FGMUserDefinedGradedMaterial)
            {
                double[] point = null;
                if (((FGMUserDefinedGradedMaterial)Material).GetIsGlobal())
                {
                    point = PointAt(xi);
                }
                else
                {
                    point = new double[] { xi };
                }
                if (Material is FGMOneGradedDirectionMaterial)
                {
                    int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
                    v = Material.GetProperty(name).GetValueProperty(point[direction]);
                }
                else
                {
                    v = Material.GetProperty(name).GetValueProperty(point);
                }
            }
            else
            {
                v = Material.GetProperty(name).GetValueProperty();
            }
            return v;
        }
        public override void ComputeDrawValueAtGaussPoint(DataInGausspoint name)
        {
            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
            {
                gps[i].SetValue(DataInGausspoint.DrawValue, gps[i].GetValue(name));
            }
        }
        public override int[] GetTArrayGlobal()
        {
            return edge.GetTArrayGlobal();
        }

        public double[] PointAt(double xi)
        {
            return ((NURBSCurve)(patch.GetGeometry(0))).PointAt(xi);
        }

        public override bool IsInsideRegion(IRegion loc)
        {
            var paramU = GetParameterTwoEndElement(0);
            //bool isInside = false;
            for (int i = 0; i < 2; i++)
            {
                var p = PointAt(paramU[i]);
                if (loc.IsContain(p[0], p[1], p[2]))
                    return true;
            }

            return false;
        }

        public override bool IsCompleteInsideRegion(IRegion loc)
        {
            var paramU = GetParameterTwoEndElement(0);
            //bool isInside = false;
            for (int i = 0; i < 2; i++)
            {
                var p = PointAt(paramU[i]);
                if (!loc.IsContain(p[0], p[1], p[2]))
                    return false;
            }

            return true;
        }

        public Edge GetEdge()
        {
            return edge;
        }

        public void SetEdge(Edge edge)
        {
            this.edge = edge;
        }

        /// <summary>
        /// Compute derivatives wrt local coords
        /// dNdxi
        /// </summary>
        /// <param name="xi1"></param>
        /// <param name="xi2"></param>
        /// <returns></returns>
        public override DoubleMatrix GradBasisFunction(params double[] xi)
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            int d = patch.GetCountDimension();
            NURBSCurve curve = patch.GetCurve();
            var basis = (NURBSBasisFunction)curve.Basis;
            int p = basis.GetDegree(0);
            int nen = patch.GetCountLocalBasisFunctions();// (p + 1) // number of local basis functions
            int span1 = basis.FindSpan(xi[0], 0);
            DoubleMatrix dNdxi = new DoubleMatrix(nen, 1);
            var gradBasis = basis.GetDerivativeBasisFunctions(xi[0], 1);
            for (int i = 0; i < nen; i++)
            {
                int ien = patch.GetIEN(id, i);
                dNdxi[i, 0] = gradBasis[1, patch.GetINC(ien, 0) - span1 + p];
            }
            return dNdxi;
        }

        public override DoubleVector ValueBasisFunction(params double[] xi)
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            NURBSCurve curve = patch.GetCurve();
            var basis = (NURBSBasisFunction)curve.Basis;
            int p = basis.GetDegree(0);
            int nen = patch.GetCountLocalBasisFunctions();// (p + 1) * (q + 1); // number of local basis functions
            int span1 = basis.FindSpan(xi[0], 0);
            DoubleVector Ni = new DoubleVector(nen);
            var valueBasis = basis.GetValueBasisFunctions(xi[0]);
            for (int i = 0; i < nen; i++)
            {
                int ien = patch.GetIEN(id, i);
                Ni[i] = valueBasis[patch.GetINC(ien, 0) - span1 + p];
            }
            return Ni;
        }
        public DoubleMatrix Grad2BasisFunction(params double[] xi)
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            int d = patch.GetCountDimension();
            NURBSCurve curve = patch.GetCurve();
            var basis = (NURBSBasisFunction)curve.Basis;
            int p = basis.GetDegree(0);
            int nen = patch.GetCountLocalBasisFunctions();// (p + 1)  // number of local basis functions
            int span1 = basis.FindSpan(xi[0], 0);
            DoubleMatrix dNdxi = new DoubleMatrix(nen, 2);
            DoubleMatrix ddNdxi = new DoubleMatrix(nen, 3);
            DoubleMatrix reddNdx = new DoubleMatrix(nen, 3);
            var gradBasis = basis.GetDerivativeBasisFunctions(xi[0], 2);
            for (int k = 0; k < nen; k++)
            {
                DoubleMatrix dM = new DoubleMatrix(2, 2);
                int ien = patch.GetIEN(id, k);
                dNdxi[k, 0] = gradBasis[1, patch.GetINC(ien, 0) - span1 + p];
                dNdxi[k, 1] = gradBasis[1, patch.GetINC(ien, 0) - span1 + p];

                ddNdxi[k, 0] = gradBasis[1, patch.GetINC(ien, 0) - span1 + p];
                ddNdxi[k, 1] = gradBasis[1, patch.GetINC(ien, 0) - span1 + p];

                DoubleMatrix J = JacobianAt(dNdxi);
                DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
                DoubleMatrix jac1 = Jacobian2At(ddNdxi);
                DoubleMatrix dNdx = dNdxi * invertJ;
                DoubleMatrix ddNdX = new DoubleMatrix(nen, 2);

                DoubleMatrix jac2 = new DoubleMatrix(2, 2);
                jac2[0, 0] = J[0, 0] * J[0, 0];
                jac2[0, 1] = J[1, 0] * J[1, 0];
                jac2[1, 0] = J[0, 1] * J[0, 1];
                jac2[1, 1] = J[1, 1] * J[1, 1];
                DoubleMatrix term1 = dNdx * jac1.Transpose();
                DoubleMatrix term2 = ddNdxi - term1;
                ddNdX = MatrixFunctions.Inverse(jac2) * term2.Transpose();
                reddNdx = ddNdX.Transpose();
            }
            return reddNdx;
        }
        /// <summary>
        /// Compute the jacobian matrix
        /// dxdxi
        /// </summary>
        /// <param name="dNdxi">Derivatives wrt local coords dNdxi</param>
        /// <returns></returns>
        public override DoubleMatrix JacobianAt(DoubleMatrix dNdxi)
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            int nen = patch.GetCountLocalBasisFunctions();
            int d = patch.GetCountDimension();
            DoubleMatrix dxdxi = new DoubleMatrix(d, d);
            var cps = patch.GetCurve().ControlPoints;
            for (int j = 0; j < d; j++)
            {
                for (int k = 0; k < nen; k++)
                {
                    int ien = patch.GetIEN(id, k);
                    var controlPointCoord = cps[patch.GetINC(ien, 0)].GetCoordinate();
                    for (int i = 0; i < d; i++)
                    {
                        dxdxi[i, j] += dNdxi[k, j] * controlPointCoord[i];
                    }
                }
            }
            return dxdxi;
        }
        protected DoubleMatrix Jacobian2At(DoubleMatrix ddNdxi)
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            int nen = patch.GetCountLocalBasisFunctions();
            int d = patch.GetCountDimension();
            DoubleMatrix ddxdxi = new DoubleMatrix(d, d); //3D: DoubleMatrix(d, 6);//6 là xx,yy,zz,xy,yz,xz //d==3 // 2D: DoubleMatrix(3, d);//d==2
            var cps = patch.GetCurve().ControlPoints;
            for (int j = 0; j < d; j++)
            {
                for (int k = 0; k < nen; k++)
                {
                    int ien = patch.GetIEN(id, k);
                    var controlPointCoord = cps[patch.GetINC(ien, 0)].GetCoordinate();
                    for (int i = 0; i < d; i++)
                    {
                        ddxdxi[i, j] += ddNdxi[k, j] * controlPointCoord[i];
                    }
                }
            }
            return ddxdxi;
        }

        public GaussPoints GetGaussPoint(int i)
        { return gps[i]; }

        public double ComputeLenghtOfElement()
        {
            double lenght = 0;
            var d = GetPatch().GetCountDimension();
            int id = GetID();
            AbstractPatch1D patch = (AbstractPatch1D)GetPatch();
            NURBSCurve curve = patch.GetCurve();
            NURBSBasisFunction basis = (NURBSBasisFunction)((NURBSCurve)(patch.GetGeometry())).Basis;
            int p = basis.GetDegree(0);
            int nen = (p + 1); // number of local basis functions
            var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
            int idx1 = patch.GetIPN(id, 0);
            int numberOfGaussPointOnEachDirection = GetNumberOfGaussPointOnEachDirection();
            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
            {
                double detJbar = 1.0 / 4.0
                       * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1]);
                lenght += GetGaussPoint(i).weight * detJbar * (double)gps[i].GetValue(DataInGausspoint.detJ);//J.Determinant();
            }
            return Math.Abs(lenght);

            /*return face.GetArea();*/

        }

        public override void ComputeMaterialPropertyValueAtGaussPoint(MaterialPropertyName name)
        {
            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
            {
                if (Material is FGMOneGradedDirectionMaterial)
                {
                    double xi1 = gps[i].location[0];
                    double xi2 = gps[i].location[1];
                    double[] point = null;
                    int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
                    if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
                    {
                        point = PointAt(xi1);
                    }
                    else
                    {
                        point = new double[] { xi1, xi2 };
                    }
                    gps[i].SetValue(DataInGausspoint.MaterialPropertyValue, Material.GetProperty(name).GetValueProperty(point[direction]));
                    //gps[i, j].materialPropertyValue = Material.GetProperty(name).GetValueProperty(point[direction]);
                }
                else
                {
                    gps[i].SetValue(DataInGausspoint.MaterialPropertyValue, Material.GetProperty(name).GetValueProperty());
                    //gps[i, j].materialPropertyValue = Material.GetProperty(name).GetValueProperty();
                }
            }
        }
        public override DoubleVector GetDisplacementLocal()
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            var cps = ((NURBSCurve)(patch.GetGeometry())).ControlPoints;
            int d = patch.GetCountField();
            int nen = patch.GetCountLocalBasisFunctions();
            DoubleVector U = new DoubleVector(d * nen);
            for (int i = 0; i < nen; i++)
            {
                var ien = patch.GetIEN(id, i);
                U[d * i] = cps[patch.GetINC(ien, 0)].GetResult(Result.UX);
                U[d * i + 1] = cps[patch.GetINC(ien, 0)].GetResult(Result.UY);
                U[d * i + 2] = cps[patch.GetINC(ien, 0)].GetResult(Result.THETAZ);
                //U[d * i ] = cps[patch.GetINC(ien, 0)].GetResult(Result.UY);
            }
            return U;
        }
        public override ControlPoint[] GetControlPointsLocal()
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            var cps = ((NURBSCurve)(patch.GetGeometry())).ControlPoints;
            int nen = patch.GetCountLocalBasisFunctions();
            ControlPoint[] cpLocal = new ControlPoint[nen];
            for (int i = 0; i < nen; i++)
            {
                var ien = patch.GetIEN(id, i);
                cpLocal[i] = cps[patch.GetINC(ien, 0)];
            }
            return cpLocal;
        }
        public override DoubleVector GetEachVariableLocal(Result variable)
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            var cps = ((NURBSCurve)(patch.GetGeometry())).ControlPoints;
            int nen = patch.GetCountLocalBasisFunctions();
            DoubleVector U = new DoubleVector(nen);
            for (int i = 0; i < nen; i++)
            {
                var ien = patch.GetIEN(id, i);
                U[i] = cps[patch.GetINC(ien, 0)].GetResult(variable);
            }
            return U;
        }
    }
}
