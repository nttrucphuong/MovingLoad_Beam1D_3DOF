using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
    public abstract class AbstractElementStructure1D : AbstractElement1D, IElementStructure
    {
        public AbstractElementStructure1D(AbstractPatch1D mesh, int id)
              : base(mesh, id)
        {
        }

        public override void ComputeMassMatrixElement(out DoubleMatrix Me)
        {
            PatchStructure1D patch = (PatchStructure1D)this.patch;
            int d = patch.GetCountField();
            NURBSBasisFunction basis = (NURBSBasisFunction)((NURBSCurve)(patch.GetGeometry())).Basis;
            int p = basis.GetDegree(0);
            int nen = (p + 1); // number of local basis functions
            var kvNoMulticiply1 = basis.GetKnotVector().GetKnotVectorNoMultiplicity();
            int idx1 = patch.GetIPN(id, 0);
            Me = new DoubleMatrix(d * nen, d * nen);
            double rho = 0;
            if (!(Material is FGMUserDefinedGradedMaterial))
                rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty();
            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
            {
                GaussPoints gpsi = gps[i];
                double detJbar = 1.0 / 2.0 * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1]);
                DoubleVector Nij = (DoubleVector)gpsi.GetValue(DataInGausspoint.Ni);//gps[i,j].Ni;
                DoubleMatrix Ni = new DoubleMatrix(d * nen, d);
                int count = 0;
                for (int kk = 0; kk < nen; kk++)
                {
                    for (int ii = 0; ii < d; ii++)
                    {
                        Ni[d * count + ii, ii] = Nij[kk];
                    }
                    count++;
                }
                if (Material is FGMUserDefinedGradedMaterial)
                {
                    //double[] point = null;
                    //int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
                    //if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
                    //{
                    //  point = PointAt(gpsij.location[0], gpsij.location[1]);
                    //}
                    //else
                    //{
                    //  point = new double[] { gpsij.location[0], gpsij.location[1] };
                    //}

                    //rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty(point[direction]);
                    rho = ComputeParameterProperty(MaterialPropertyName.Density, gpsi.location[0]);
                    DoubleMatrix NeTNe = rho * MatrixFunctions.Product(Ni, Ni.Transpose());

                    Me += gpsi.weight * detJbar * Math.Abs((double)gpsi.GetValue(DataInGausspoint.detJ)) * NeTNe;
                }
            }
        }

        public override void ComputeGeometricStiffnessMatrixElement(double[,] stressTensor, out DoubleMatrix KGe)
        {
            PatchStructure1D patch = (PatchStructure1D)this.patch;
            int d = patch.GetCountField();
            NURBSBasisFunction basis = (NURBSBasisFunction)((NURBSCurve)(patch.GetGeometry())).Basis;
            int p = basis.GetDegree(0);
            int nen = (p + 1); // number of local basis functions
            var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
            int idx1 = patch.GetIPN(id, 0);
            KGe = new DoubleMatrix(d * nen, d * nen);
            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
            {
                double detJbar = 1.0 / 2.0 * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1]);
                GaussPoints gpsi = gps[i];
                DoubleMatrix dNdX = (DoubleMatrix)gpsi.GetValue(DataInGausspoint.dNdX); //gps[i, j, k].dNdX;
                DoubleMatrix GA = new DoubleMatrix(6, d * nen); ///////////////
                int count = 0;
                for (int ii = 0; ii < nen; ii++)
                {
                    GA[0, 3 * count] = dNdX[ii, 0];
                    GA[1, 3 * count] = dNdX[ii, 1];
                    GA[2, 3 * count] = dNdX[ii, 2];
                    GA[3, 3 * count + 1] = dNdX[ii, 0];
                    GA[4, 3 * count + 1] = dNdX[ii, 1];
                    GA[5, 3 * count + 1] = dNdX[ii, 2];
                    count++;
                }
                DoubleMatrix N0 = new DoubleMatrix(6, 6);
                //double h = 0.025;
                //N0[2, 2] = -1 / h;
                //N0[5, 5] = -1 / h;
                //N0[8, 8] = -1 / h;
                DoubleVector stress = (DoubleVector)gpsi.GetValue(DataInGausspoint.currentStress);
                if (stressTensor == null)
                {
                    N0[0, 0] = stress[0];
                    N0[1, 1] = stress[1];
                    N0[2, 2] = stress[2];
                    N0[0, 1] = N0[1, 0] = stress[3];
                    N0[1, 2] = N0[2, 1] = stress[4];
                    N0[0, 2] = N0[2, 0] = stress[5];

                    N0[3, 3] = stress[0];
                    N0[4, 4] = stress[1];
                    N0[5, 5] = stress[2];
                    N0[3, 4] = N0[4, 3] = stress[3];
                    N0[4, 5] = N0[5, 4] = stress[4];
                    N0[3, 5] = N0[5, 3] = stress[5];
                }
                else
                {
                    int rows = stressTensor.GetLength(0);
                    for (int ii = 0; ii < rows; ii++)
                    {
                        int cols = stressTensor.GetLength(1);
                        for (int jj = 0; jj < cols; jj++)
                        {
                            for (int kk = 0; kk < 1; kk++)
                            {
                                N0[ii + kk * rows, jj + kk * cols] = stressTensor[ii, jj];
                            }
                        }
                    }
                }
                DoubleMatrix GaTN0GA = NMathFunctions.TransposeProduct(GA, MatrixFunctions.Product(-N0, GA));
                KGe += gpsi.weight * detJbar * Math.Abs((double)gpsi.GetValue(DataInGausspoint.detJ)) * GaTN0GA;

                //KGe += a * gps[i, j].weight * detJbar * Math.Abs(gps[i, j].detJ) * NeTNe;
            }
        }
        protected virtual DoubleMatrix CreateMaterialMatrix()
        {
            AnisotropicElasticity prob = (AnisotropicElasticity)Material.GetProperty(MaterialPropertyName.AnisotropicElasticity);
            if (prob != null)
            {
                double[,] ce = prob.GetAnisotropicElasticityMatrix();
                DoubleMatrix C = new DoubleMatrix(4, 4);
                bool isInverse = prob.GetIsInverseMatrix();
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                    {
                        C[i, j] = ce[i, j];
                    }
                if (((PatchStructure1D)patch).StateStress != Structure1DState.Beam)
                {
                    C[2, 2] = 1;
                }
                if (!isInverse)
                    return C;
                else
                    return MatrixFunctions.Inverse(C);
            }
            else
            {
                OrthotropicElasticity prob1 = (OrthotropicElasticity)Material.GetProperty(MaterialPropertyName.OrthotropicElasticity);
                DoubleMatrix C = null;
                if (prob1 == null)
                {
                    var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
                    var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
                    C = new DoubleMatrix(4, 4);
                    switch (((PatchStructure1D)patch).StateStress)
                    {
                        //case Structure1DState.Beam:
                        //    double a = eModulus / (1.0 - nu * nu);
                        //    C[0, 0] = a;
                        //    C[0, 1] = a * nu;
                        //    C[1, 0] = a * nu;
                        //    C[1, 1] = a;
                        //    C[3, 3] = a * (1 - nu) / 2;
                        //    break;
                        //case Structure1DState.Truss:
                        case Structure1DState.Beam:
                            double a = eModulus / ((1.0 + nu) * (1.0 - 2.0 * nu)); //5.32
                            C[0, 0] = a * (1.0 - nu);
                            C[0, 1] = a * nu;
                            C[0, 2] = a * nu;
                            C[1, 0] = a * nu;
                            C[1, 1] = a * (1.0 - nu);
                            C[1, 2] = a * nu;
                            C[2, 0] = a * nu;
                            C[2, 1] = a * nu;
                            C[2, 2] = a * (1.0 - nu);
                            C[3, 3] = a * (1.0 - 2.0 * nu) / 2.0;
                            break;
                    }
                }
                else
                {
                    C = new DoubleMatrix(prob1.GetOrthotropicElasticityMatrix());
                }
                return C;
            }
        }

        /// <summary>
        /// Compute material matrix at 
        /// s point
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        protected DoubleMatrix CreateMaterialMatrix(int i)
        {
            double eModulus = (double)gps[i].GetValue(DataInGausspoint.EModulus);
            double nu = (double)gps[i].GetValue(DataInGausspoint.nu);
            DoubleMatrix C = new DoubleMatrix(4, 4);
            switch (((PatchStructure1D)patch).StateStress)
            {
                case Structure1DState.Beam:
                    double a = eModulus / (1.0 - nu * nu);
                    C[0, 0] = a;
                    C[0, 1] = a * nu;
                    C[1, 0] = a * nu;
                    C[1, 1] = a;
                    C[3, 3] = a * (1 - nu) / 2;
                    //C[3, 3] = a;
                    break;
                    //case Structure2DState.PlaneStrain:
                    //case Structure2DState.Axisymetric:
                    //    a = eModulus / ((1.0 + nu) * (1.0 - 2.0 * nu));
                    //    C[0, 0] = a * (1.0 - nu);
                    //    C[0, 1] = a * nu;
                    //    C[0, 2] = a * nu;
                    //    C[1, 0] = a * nu;
                    //    C[1, 1] = a * (1.0 - nu);
                    //    C[1, 2] = a * nu;
                    //    C[2, 0] = a * nu;
                    //    C[2, 1] = a * nu;
                    //    C[2, 2] = a * (1.0 - nu);
                    //    C[3, 3] = a * (1.0 - 2.0 * nu) / 2.0;
                    //    break;
            }
            return C;
        }

        protected DoubleMatrix CreateMaterialMatrix(double xi)
        {
            //var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty(point[direction]);
            //var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty(point[direction]);
            var eModulus = ComputeParameterProperty(MaterialPropertyName.YoungModulus, xi);
            var nu = ComputeParameterProperty(MaterialPropertyName.PoissonRatio, xi);
            DoubleMatrix C = new DoubleMatrix(4, 4);
            switch (((PatchStructure1D)patch).StateStress)
            {
                case Structure1DState.Beam:
                    double a = eModulus / (1.0 - nu * nu);
                    C[0, 0] = a;
                    C[0, 1] = a * nu;
                    C[1, 0] = a * nu;
                    C[1, 1] = a;
                    C[3, 3] = a * (1 - nu) / 2;
                    //C[3, 3] = a;
                    break;
                    //case Structure2DState.PlaneStrain:
                    //case Structure2DState.Axisymetric:
                    //    a = eModulus / ((1.0 + nu) * (1.0 - 2.0 * nu));
                    //    C[0, 0] = a * (1.0 - nu);
                    //    C[0, 1] = a * nu;
                    //    C[0, 2] = a * nu;
                    //    C[1, 0] = a * nu;
                    //    C[1, 1] = a * (1.0 - nu);
                    //    C[1, 2] = a * nu;
                    //    C[2, 0] = a * nu;
                    //    C[2, 1] = a * nu;
                    //    C[2, 2] = a * (1.0 - nu);
                    //    C[3, 3] = a * (1.0 - 2.0 * nu) / 2.0;
                    //    break;
            }
            return C;
        }
        public abstract DoubleVector StressAt(params double[] xi);
        public abstract DoubleVector StrainAt(params double[] xi);
        public abstract DoubleVector StrainElasticAt(params double[] xi);
        public abstract DoubleVector StrainThermoAt(params double[] xi);
        public override void ComputeValueAtGaussPoint(DataInGausspoint name)
        {
            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
            {
                for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
                {
                    switch (name)
                    {
                        case DataInGausspoint.currentStress:
                            gps[i].SetValue(DataInGausspoint.currentStress, StressAt(gps[i].location));
                            //gps[i, j].currentStress = StressAt(gps[i, j].location);
                            break;
                        case DataInGausspoint.currentStrain:
                            gps[i].SetValue(DataInGausspoint.currentStrain, StrainAt(gps[i].location));
                            //gps[i, j].currentStrain = StrainAt(gps[i, j].location);
                            break;
                    }
                }
            }
        }
    }
}
