using CenterSpace.NMath.Core;
using DEMSoft.Common;
using DEMSoft.EngineeringData;
using DEMSoft.NURBS;
using System;

namespace DEMSoft.IGA
{
    /// <summary>
    /// Abstract 2D element class. It uses only for 2D problem.
    /// </summary>
    public abstract class AbstractElement1D : AbstractElement
    {
        protected Edge edge;
        protected GaussPoints[] gps;

        public AbstractElement1D(AbstractPatch1D patch, int id)
          : base(patch, id)
        {
            edge = new Edge(this);
        }

        public override void InitializeGaussPoint()
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            NURBSBasisFunction basis = (NURBSBasisFunction)((NURBSCurve)(patch.GetGeometry())).Basis;
            var kvNoMulticiply = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
            int idx = patch.GetIPN(id, 0);
            int nen = patch.GetCountLocalBasisFunctions();// (p + 1); // number of local basis functions
            int p = basis.GetDegree(0);

            numberOfGaussPointOnEachDirection = p + 1;//(int)Math.Ceiling((2.0 * Math.Max(p, q) + 1.0) / 2.0);

            gps = new GaussPoints[numberOfGaussPointOnEachDirection];

            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
            {
                double xi = GaussPoints.GetPoint(numberOfGaussPointOnEachDirection, i);
                double w = GaussPoints.GetWeight(numberOfGaussPointOnEachDirection, i);
                xi = 0.5 * ((kvNoMulticiply[idx + 1] - kvNoMulticiply[idx]) * xi + kvNoMulticiply[idx + 1] + kvNoMulticiply[idx]);
                int span1 = basis.FindSpan(xi, 0);

                GaussPoints gpsij = gps[i] = new GaussPoints(new double[] { xi }, w);
                gpsij.listDataGausspoint = new System.Collections.Generic.List<DataGausspoint>();
                DoubleVector dNdxi = new DoubleVector(nen);
                DoubleVector ddNdxi = new DoubleVector(nen);
                DoubleVector Ni = new DoubleVector(nen);
                if (patch.ExtractionOperator == null)
                {
                    double[,] gradBasis = basis.GetDerivativeBasisFunctions(xi, 2);
                    for (int k = 0; k < nen; k++)
                    {
                        int ien = patch.GetIEN(id, k);
                        Ni[k] = gradBasis[0, patch.GetINC(ien, 0) - span1 + p];
                        dNdxi[k] = gradBasis[1, patch.GetINC(ien, 0) - span1 + p];
                        ddNdxi[k] = gradBasis[2, patch.GetINC(ien, 0) - span1 + p];
                    }
                }
                else
                {
                    //int indexGausspointOnLine = j * numberOfGaussPointOnEachDirection + i;//j + numberOfGaussPointOnEachDirection * i;
                    //DoubleVector Bb = (DoubleVector)patch.GaussPointOnBezierElement[indexGausspointOnLine].GetValue(DataInGausspoint.NBernsteini);//patch.GaussPointOnBezierElement[indexGausspointOnLine].NBernsteini;
                    //                                                                                                                              //int ex = GetIndexGlobalCoordinate(0);
                    //                                                                                                                              //int ey = GetIndexGlobalCoordinate(1);
                    //DoubleVector we = new DoubleVector((p + 1) * (q + 1));
                    //int c = 0;
                    //for (int jj = 0; jj < q + 1; jj++)
                    //    for (int ii = 0; ii < p + 1; ii++)
                    //    {
                    //        we[c] = ((NURBSBasisFunction)patch.GetSurface().Basis).Weights[idx1 + ii, idx2 + jj];
                    //        c++;
                    //    }
                    //DoubleMatrix mCe = new DoubleMatrix(patch.ExtractionOperator[idx1, idx2]);
                    //DoubleVector wb = MatrixFunctions.Product(mCe.Transpose(), we);
                    //double Wb = MatrixFunctions.Dot(Bb, wb);
                    //DoubleMatrix weDiag = new DoubleMatrix(we.Length, we.Length);
                    //for (int iii = 0; iii < we.Length; iii++)
                    //{
                    //    weDiag[iii, iii] = we[iii];
                    //}
                    //DoubleVector Rb = MatrixFunctions.Product(weDiag, MatrixFunctions.Product(mCe, Bb)) / Wb;

                    //////////////////////////;;;;
                    //DoubleVector[] dB = new DoubleVector[2];
                    //DoubleMatrix matdNB = (DoubleMatrix)patch.GaussPointOnBezierElement[indexGausspointOnLine].GetValue(DataInGausspoint.dNBernsteindxi);
                    //dB[0] = matdNB.Col(0);
                    //dB[1] = matdNB.Col(1);
                    //double dWbxi = MatrixFunctions.Dot(dB[0], wb);
                    //double dWbeta = MatrixFunctions.Dot(dB[1], wb);
                    //DoubleVector dRxi = MatrixFunctions.Product(weDiag, MatrixFunctions.Product(mCe, (dB[0] / Wb - dWbxi * Bb / (Wb * Wb)))) / (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1]);
                    //DoubleVector dReta = MatrixFunctions.Product(weDiag, MatrixFunctions.Product(mCe, (dB[1] / Wb - dWbeta * Bb / (Wb * Wb)))) / (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);

                    //Ni = Rb;
                    //for (int k = 0; k < nen; k++)
                    //{
                    //    dNdxi[k, 0] = dRxi[k];
                    //    dNdxi[k, 1] = dReta[k];
                    //}
                }

                DoubleVector J = dNdxi;//[1,1]
                DoubleVector dNdx = dNdxi;//[nen,1]*[1,1]-->[nen,1]
                gpsij.SetValue(DataInGausspoint.Ni, Ni);
                gpsij.SetValue(DataInGausspoint.dNdxi, dNdxi);
                gpsij.SetValue(DataInGausspoint.detJ, J);
                gpsij.SetValue(DataInGausspoint.dNdX, dNdx);

                //if (this is ElementStructurePhaseField1D)
                //{
                //    //DoubleMatrix ddNdX = new DoubleMatrix(nen, 3);
                //    gpsij.SetValue(DataInGausspoint.ddNdxi, ddNdxi);
                //    //gps[i, j].ddNdxi = ddNdxi;//bo sung

                //    DoubleMatrix jac2 = new DoubleMatrix(3, 3);
                //    jac2[0, 0] = J[0, 0] * J[0, 0];
                //    jac2[0, 1] = J[1, 0] * J[1, 0];
                //    jac2[0, 2] = 2 * J[0, 0] * J[1, 0];
                //    jac2[1, 0] = J[0, 1] * J[0, 1];
                //    jac2[1, 1] = J[1, 1] * J[1, 1];
                //    jac2[1, 2] = 2 * J[0, 1] * J[1, 1];
                //    jac2[2, 0] = J[0, 0] * J[0, 1];
                //    jac2[2, 1] = J[1, 0] * J[1, 1];
                //    jac2[2, 2] = J[1, 0] * J[0, 1] + J[0, 0] * J[1, 1];
                //    DoubleMatrix jac1 = Jacobian2At(ddNdxi/*[nen,3]*/); //[3,2]
                //    DoubleMatrix term1 = MatrixFunctions.Product(dNdx, jac1.Transpose());//[nen,2]*[2,3]-->[nen,3]
                //    DoubleMatrix term2 = ddNdxi - term1;//[nen,3]-[nen,3]-->[nen,3]
                //    DoubleMatrix ddNdX = MatrixFunctions.Product(MatrixFunctions.Inverse(jac2), term2.Transpose());//[3,3]*[3,nen]-->[3,nen]
                //    gpsij.SetValue(DataInGausspoint.ddNdX, ddNdX.Transpose());
                //    //gps[i, j].ddNdX = ddNdX.Transpose();//bo sung
                //}
                //if (this is ElementStructurePhaseField2D)
                MaterialProperty Gc = Material.GetProperty(MaterialPropertyName.CriticalEnergyReleaseRate);
                if (Gc != null)
                {
                    gpsij.SetValue(DataInGausspoint.currentPhase, 0.0);
                    gpsij.SetValue(DataInGausspoint.lastPhase, 0.0);
                    gpsij.SetValue(DataInGausspoint.previousPhase, 0.0);
                    gpsij.SetValue(DataInGausspoint.xiEpsilon, 0.0);
                    gpsij.SetValue(DataInGausspoint.currentStress, new DoubleVector(2));
                    gpsij.SetValue(DataInGausspoint.currentStrain, new DoubleVector(2));
                    gpsij.SetValue(DataInGausspoint.lastStress, new DoubleVector(2));
                    gpsij.SetValue(DataInGausspoint.lastStrain, new DoubleVector(2));
                }

                if (Material is PlasticityMaterial)
                {
                    gpsij.SetValue(DataInGausspoint.lastAlpha, 0.0);
                    gpsij.SetValue(DataInGausspoint.lastBackStress, new DoubleVector(2));
                    gpsij.SetValue(DataInGausspoint.lastPlasticStrain, new DoubleVector(2));
                    gpsij.SetValue(DataInGausspoint.lastStress, new DoubleVector(2));

                    //gps[i, j].lastAlpha = 0.0;
                    //gps[i, j].lastBackStress = new DoubleVector(4);
                    //gps[i, j].lastPlasticStrain = new DoubleVector(4);
                    //gps[i, j].lastStress = new DoubleVector(4);
                }
                else if (Material is FGMUserDefinedGradedMaterial)
                {
                    //double[] point = null;
                    //int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
                    //if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
                    //{
                    //  point = PointAt(xi1, xi2);
                    //}
                    //else
                    //{
                    //  point = new double[] { xi1, xi2 };
                    //}
                    //MaterialProperty matE = Material.GetProperty(MaterialPropertyName.YoungModulus);
                    //MaterialProperty matNu = Material.GetProperty(MaterialPropertyName.PoissonRatio);
                    //MaterialProperty matK = Material.GetProperty(MaterialPropertyName.IsotropicThermalConductivity);
                    //gpsij.SetValue(DataInGausspoint.EModulus, (matE != null) ? matE.GetValueProperty(point[direction]) : 0.0);
                    //gpsij.SetValue(DataInGausspoint.nu, (matNu != null) ? matNu.GetValueProperty(point[direction]) : 0.0);
                    //gpsij.SetValue(DataInGausspoint.ThermalConductivity, (matK != null) ? matK.GetValueProperty(point[direction]) : 0.0);


                    MaterialProperty matE = Material.GetProperty(MaterialPropertyName.YoungModulus);
                    MaterialProperty matNu = Material.GetProperty(MaterialPropertyName.PoissonRatio);
                    MaterialProperty matK = Material.GetProperty(MaterialPropertyName.IsotropicThermalConductivity);
                    gpsij.SetValue(DataInGausspoint.EModulus, (matE != null) ? ComputeParameterProperty(MaterialPropertyName.YoungModulus, xi) : 0.0);
                    gpsij.SetValue(DataInGausspoint.nu, (matNu != null) ? ComputeParameterProperty(MaterialPropertyName.PoissonRatio, xi) : 0.0);
                    gpsij.SetValue(DataInGausspoint.ThermalConductivity, (matK != null) ? ComputeParameterProperty(MaterialPropertyName.IsotropicThermalConductivity, xi) : 0.0);

                    MaterialProperty matEpsilon = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion);
                    //gpsij.SetValue(DataInGausspoint.CoefficientsThermalExpansion, (matEpsilon != null) ? matEpsilon.GetValueProperty(point[direction]) : 0.0);
                    gpsij.SetValue(DataInGausspoint.CoefficientsThermalExpansion, (matEpsilon != null) ? ComputeParameterProperty(MaterialPropertyName.CoefficientThermalExpansion, xi) : 0.0);
                }
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
            int nen = patch.GetCountLocalBasisFunctions();// (p + 1)  // number of local basis functions
            int span1 = basis.FindSpan(xi[0], 0);
            DoubleMatrix dNdxi = new DoubleMatrix(nen, 1); //(nen,d) but in 1D = 1
            var gradBasis = basis.GetDerivativeBasisFunctions(xi[0], 1);
            for (int i = 0; i < nen; i++)
            {
                int ien = patch.GetIEN(id, i);
                dNdxi[i, 0] = gradBasis[1, patch.GetINC(ien, 0) - span1 + p];
                //double[,] grad = gradBasis[patch.GetINC(ien, 0) - span1 + p, patch.GetINC(ien, 1) - span2 + q];
            }
            return dNdxi;
        }

        public override DoubleVector ValueBasisFunction(params double[] xi)
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            NURBSCurve curve = patch.GetCurve();
            var basis = (NURBSBasisFunction)curve.Basis;
            int p = basis.GetDegree(0);
            int nen = patch.GetCountLocalBasisFunctions();// (p + 1)  // number of local basis functions
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
        //public DoubleMatrix Grad2BasisFunction(params double[] xi)
        //{
        //  AbstractPatch1D patch = (AbstractPatch1D)this.patch;
        //  int d = patch.GetCountDimension();
        //  NURBSCurve curve = patch.GetSurface();
        //  var basis = (NURBSBasisFunction)curve.Basis;
        //  int p = basis.GetDegree(0);
        //  int q = basis.GetDegree(1);
        //  int nen = patch.GetCountLocalBasisFunctions();// (p + 1) * (q + 1); // number of local basis functions
        //  int span1 = basis.FindSpan(xi[0], 0);
        //  int span2 = basis.FindSpan(xi[1], 1);
        //  DoubleMatrix dNdxi = new DoubleMatrix(nen, 2);
        //  DoubleMatrix ddNdxi = new DoubleMatrix(nen, 3);
        //  DoubleMatrix reddNdx = new DoubleMatrix(nen, 3);
        //  var gradBasis = basis.GetDerivativeBivariateBasisFunctions(xi[0], xi[1], 2);
        //  for (int k = 0; k < nen; k++)
        //  {
        //    DoubleMatrix dM = new DoubleMatrix(2, 2);
        //    int ien = patch.GetIEN(id, k);
        //    dNdxi[k, 0] = gradBasis[patch.GetINC(ien, 0) - span1 + p, patch.GetINC(ien, 1) - span2 + q][1, 0];
        //    dNdxi[k, 1] = gradBasis[patch.GetINC(ien, 0) - span1 + p, patch.GetINC(ien, 1) - span2 + q][0, 1];

        //    ddNdxi[k, 0] = gradBasis[patch.GetINC(ien, 0) - span1 + p, patch.GetINC(ien, 1) - span2 + q][2, 0];
        //    ddNdxi[k, 1] = gradBasis[patch.GetINC(ien, 0) - span1 + p, patch.GetINC(ien, 1) - span2 + q][0, 2];
        //    ddNdxi[k, 2] = gradBasis[patch.GetINC(ien, 0) - span1 + p, patch.GetINC(ien, 1) - span2 + q][1, 1];

        //    DoubleMatrix J = JacobianAt(dNdxi);
        //    DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
        //    DoubleMatrix jac1 = Jacobian2At(ddNdxi);
        //    DoubleMatrix dNdx = dNdxi * invertJ;
        //    DoubleMatrix ddNdX = new DoubleMatrix(nen, 3);

        //    DoubleMatrix jac2 = new DoubleMatrix(3, 3);
        //    jac2[0, 0] = J[0, 0] * J[0, 0];
        //    jac2[0, 1] = J[1, 0] * J[1, 0];
        //    jac2[0, 2] = 2 * J[0, 0] * J[1, 0];
        //    jac2[1, 0] = J[0, 1] * J[0, 1];
        //    jac2[1, 1] = J[1, 1] * J[1, 1];
        //    jac2[1, 2] = 2 * J[0, 1] * J[1, 1];
        //    jac2[2, 0] = J[0, 0] * J[0, 1];
        //    jac2[2, 1] = J[1, 0] * J[1, 1];
        //    jac2[2, 2] = J[1, 0] * J[0, 1] + J[0, 0] * J[1, 1];
        //    DoubleMatrix term1 = dNdx * jac1.Transpose();
        //    DoubleMatrix term2 = ddNdxi - term1;
        //    ddNdX = MatrixFunctions.Inverse(jac2) * term2.Transpose();
        //    reddNdx = ddNdX.Transpose();
        //  }
        //  return reddNdx;
        //}
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
            DoubleMatrix dxdxi = new DoubleMatrix(d, 1);
            var cps = patch.GetCurve().ControlPoints;

            for (int k = 0; k < nen; k++)
            {
                int ien = patch.GetIEN(id, k);
                var controlPointCoord = cps[patch.GetINC(ien, 0)].GetCoordinate();
                for (int i = 0; i < d; i++)
                {
                    dxdxi[i, 0] += dNdxi[k, 0] * controlPointCoord[i];
                }
            }

            return dxdxi;
        }

        protected DoubleVector Jacobian2At(DoubleVector ddNdxi)
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            int nen = patch.GetCountLocalBasisFunctions();
            int d = patch.GetCountDimension();
            DoubleVector ddxdxi = new DoubleVector(3, d);
            var cps = patch.GetCurve().ControlPoints;
            for (int k = 0; k < nen; k++)
            {
                int ien = patch.GetIEN(id, k);
                var controlPointCoord = cps[patch.GetINC(ien, 0)].GetCoordinate();
                ddxdxi[0] += ddNdxi[k] * controlPointCoord[0];
            }
            return ddxdxi;
        }
        public double ComputeAreaOfElement()
        {
            double area = 0;
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
                //DoubleMatrix dNdxi = GetGaussPoint(i, j).dNdxi;
                //DoubleMatrix J = JacobianAt(dNdxi);
                area += GetGaussPoint(i).weight * detJbar * (double)gps[i].GetValue(DataInGausspoint.detJ);//J.Determinant();
            }
            return Math.Abs(area);

            /*return face.GetArea();*/

        }

        public override void ComputeMaterialPropertyValueAtGaussPoint(MaterialPropertyName name)
        {
            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
            {
                if (Material.GetProperty(name) != null)
                {
                    double xi1 = gps[i].location[0];
                    double xi2 = gps[i].location[1];
                    //if (Material is FGMOneGradedDirectionMaterial)
                    //{
                    //  double xi1 = gps[i, j].location[0];
                    //  double xi2 = gps[i, j].location[1];
                    //  double[] point = null;
                    //  int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
                    //  if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
                    //  {
                    //    point = PointAt(xi1, xi2);
                    //  }
                    //  else
                    //  {
                    //    point = new double[] { xi1, xi2 };
                    //  }
                    //  gps[i, j].SetValue(DataInGausspoint.MaterialPropertyValue, Material.GetProperty(name).GetValueProperty(point[direction]));
                    //  //gps[i, j].materialPropertyValue = Material.GetProperty(name).GetValueProperty(point[direction]);
                    //}
                    //else
                    //{
                    //  gps[i, j].SetValue(DataInGausspoint.MaterialPropertyValue, Material.GetProperty(name).GetValueProperty());
                    //  //gps[i, j].materialPropertyValue = Material.GetProperty(name).GetValueProperty();
                    //}
                    gps[i].SetValue(DataInGausspoint.MaterialPropertyValue, ComputeParameterProperty(name, xi1));
                }
            }
        }

        public override void ComputeDrawValueAtGaussPoint(DataInGausspoint name)
        {
            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
            {
                gps[i].SetValue(DataInGausspoint.DrawValue, gps[i].GetValue(name));
                //gps[i, j].listDataGausspoint.Find(x => x.typeData == DataInGausspoint.DrawValue).value = (double)gps[i, j].listDataGausspoint.Find(x => x.typeData == name).value;
                //switch (name)
                //{
                //  case DataInGausspoint.xiEpsilon:
                //    gps[i, j].DrawValue = gps[i, j].xiEpsilon;
                //    break;
                //  case DataInGausspoint.currentPhase:
                //    gps[i, j].DrawValue = gps[i, j].currentPhase;
                //    break;
                //}
            }
        }
        public override DoubleVector GetDisplacementLocal()
        {
            AbstractPatch1D patch = (AbstractPatch1D)this.patch;
            var cps = ((NURBSCurve)(patch.GetGeometry())).ControlPoints;
            int d = patch.GetCountDimension();//2
            int nen = patch.GetCountLocalBasisFunctions();
            DoubleVector U = new DoubleVector(d * nen);
            for (int i = 0; i < nen; i++)
            {
                var ien = patch.GetIEN(id, i);
                ControlPoint cp = cps[patch.GetINC(ien, 0)];
                U[d * i] = cp.GetResult(Result.UX);
                U[d * i + 1] = cp.GetResult(Result.UY);
                U[d * i + 2] = cp.GetResult(Result.THETAZ);
                //U[d * i] = cp.GetResult(Result.UY);
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
        public GaussPoints GetGaussPoint(int i)
        { return gps[i]; }
    }
}
