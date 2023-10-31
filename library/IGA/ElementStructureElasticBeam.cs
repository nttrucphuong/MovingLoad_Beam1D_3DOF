using CenterSpace.NMath.Core;
using DEMSoft.EngineeringData;
using DEMSoft.Function;
using DEMSoft.NURBS;
using System;

namespace DEMSoft.IGA
{
    public enum TypeBeamFunction
    { Power, MoriTanaka }
    /// <summary>
    /// 
    /// </summary>
    public class ElementStructureElasticBeam : AbstractElementStructureBeam
    {
        public double GetThickness()
        { return ((PatchStructureBeam)patch).Thickness; }
        public TypeBeam GetTypeBeam()
        { return ((PatchStructureBeam)patch).TypeBeam; }
        public FunctionRToR GetKinematicsFunction()
        { return ((PatchStructureBeam)patch).KinematicsFunction; }
        public double getzz()
        { return ((PatchStructureBeam)patch).Getzz; }
        public TypeVFunction GetTypeVFunction()
        { return ((PatchStructureBeam)patch).GetTypeVF; }
        public ElementStructureElasticBeam(AbstractPatch1D patch, int id)
                : base(patch, id)
        {
        }

        private DoubleMatrix CreateMaterialMatrixHSDT()
        {
            double thickness = GetThickness(); //get from h in area
            DoubleMatrix D = new DoubleMatrix(4, 4);
            int numGauss = 7; //tu cho, FGM bac cao se co so diem Gauss nhieu hon
            for (int i = 0; i < numGauss; i++)
            {
                double xi = GaussPoints.GetPoint(numGauss, i);
                double wi = GaussPoints.GetWeight(numGauss, i);
                double zz = 0.5 * thickness * (xi + 1.0) - thickness / 2.0;
                double fz = GetKinematicsFunction().ValueAt(zz);
                double dfz = GetKinematicsFunction().DerivativeAt(zz);
                double detJ = thickness / 2.0;
                var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
                var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();

                double aa = eModulus / (1.0 - nu * nu);
                double bb = eModulus / (2.0 * (1.0 + nu));

                D[0, 0] += detJ * wi * aa;//A
                D[0, 1] += detJ * wi * zz * aa;//B
                D[0, 2] += detJ * wi * fz * aa;//E

                D[1, 1] += detJ * wi * zz * zz * aa;//D
                D[1, 2] += detJ * wi * zz * fz * aa;//F

                D[2, 2] += detJ * wi * fz * fz * aa;//H

                D[3, 3] += detJ * wi * dfz * dfz * bb;//Ds
            }

            for (int ii = 0; ii < 4; ii++)
            {
                for (int jj = ii; jj < 4; jj++)
                {
                    D[jj, ii] = D[ii, jj];
                }
            }
            return D;
        }

        private DoubleMatrix CreateMaterialMatrixFSDTFGM()
        {
            DoubleMatrix D = new DoubleMatrix(3, 3);
            int numGauss = 7; //tu cho, FGM bac cao se co so diem Gauss nhieu hon
            for (int i = 0; i < numGauss; i++)
            {
                double thickness = GetThickness();
                double xi = GaussPoints.GetPoint(numGauss, i);
                double wi = GaussPoints.GetWeight(numGauss, i);
                double zz = 0.5 * (xi + 1.0);
                double fz = GetKinematicsFunction().ValueAt(zz);
                //double dfz = GetKinematicsFunction().DerivativeAt(zz);
                double detJ = thickness / 2.0;
                var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
                var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();

                double aa = eModulus / (1.0 - nu * nu);
                double bb = eModulus / (2.0 * (1.0 + nu));

                D[0, 0] += detJ * wi * aa; //A
                D[0, 1] += detJ * wi * fz * aa;//E

                D[1, 1] += detJ * wi * fz * fz * aa;//H

                D[2, 2] += detJ * wi * (5.0 / 6.0) * bb;  //Ds 
            }

            for (int ii = 0; ii < 3; ii++)
            {
                for (int jj = ii; jj < 3; jj++)
                {
                    D[jj, ii] = D[ii, jj];
                }
            }
            return D;
        }

        private DoubleMatrix CreateMaterialMatrixFSDT()
        {
            DoubleMatrix D = new DoubleMatrix(3, 3);
            int numGauss = 5; //tu cho, FGM bac cao se co so diem Gauss nhieu hon
            double thickness = GetThickness();
            for (int i = 0; i < numGauss; i++)
            {
                double xi = GaussPoints.GetPoint(numGauss, i);
                double wi = GaussPoints.GetWeight(numGauss, i);
                //double zz = 0.5 * (xi + 1.0);
                double zz = 0.5 * thickness * (xi + 1.0) - thickness / 2.0;
                double fz = GetKinematicsFunction().ValueAt(zz);
                double dfz = GetKinematicsFunction().DerivativeAt(zz);
                double detJ = thickness / 2.0;
                var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
                var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();

                double aa = eModulus / (1.0 - nu * nu);
                double bb = eModulus / (2.0 * (1.0 + nu));

                D[0, 0] += detJ * wi * aa;//A
                D[0, 1] += detJ * wi * zz * aa;//B
                D[0, 2] += detJ * wi * zz * aa;//E//fz

                D[1, 1] += detJ * wi * zz * zz * aa;//D
                D[1, 2] += detJ * wi * zz * zz * aa;//F//fz

                D[2, 2] += detJ * wi * zz * zz * aa;//H //fz * fz

                //D[3, 3] += detJ * wi * dfz * dfz * bb;//Ds
            }

            for (int ii = 0; ii < 3; ii++)
            {
                for (int jj = ii; jj < 3; jj++)
                {
                    D[jj, ii] = D[ii, jj];
                }
            }
            return D;
        }

        private double CreateMatrixDs()
        {
            double Ds = 0;
            int numGauss = 5; //tu cho, FGM bac cao se co so diem Gauss nhieu hon
            double thickness = GetThickness();
            for (int i = 0; i < numGauss; i++)
            {
                double xi = GaussPoints.GetPoint(numGauss, i);
                double wi = GaussPoints.GetWeight(numGauss, i);
                //double zz = 0.5 * (xi + 1.0);
                double zz = 0.5 * thickness * (xi + 1.0) - thickness / 2.0;
                double dfz = GetKinematicsFunction().DerivativeAt(zz);
                double detJ = thickness / 2.0;
                var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
                var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
                double bb = eModulus / (2.0 * (1.0 + nu));
                double correctK = Math.Sqrt(5.0 / 6.0);//dfz;
                                                       //if (dfz == 0)
                                                       //  correctK = Math.Sqrt(5.0 / 6.0);
                Ds += detJ * wi * correctK * correctK * bb;
            }
            return Ds;
        }

        /// <summary>
        /// Compute material matrix at gauss point
        /// </summary>
        /// <param name="i"></param>
        /// <param name="j"></param>
        /// <returns></returns>

        public override void ComputeStiffnessMatrixElement(out DoubleMatrix Ke) //start
        {
            PatchStructureBeam patch = (PatchStructureBeam)this.patch;
            int d = patch.GetCountField(); //dimension
            NURBSBasisFunction basis = (NURBSBasisFunction)((NURBSCurve)(patch.GetGeometry())).Basis;
            int p = basis.GetDegree(0); //tra bac tu do theo huong của basic
            int nel = patch.CalculateNumberOfElements();
            int nen = (p + 1); // number of local basis functions
            var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();

            int idx1 = patch.GetIPN(id, 0);

            Ke = new DoubleMatrix(d * nen, d * nen);

            DoubleMatrix B; //6.21
            DoubleMatrix Bs = new DoubleMatrix(1, d * nen);
            double Ds = CreateMatrixDs();

            if (GetTypeBeam() == TypeBeam.FSDT)
            {
                B = new DoubleMatrix(3, d * nen);
            }
            else
            {
                B = new DoubleMatrix(4, d * nen);
            }

            DoubleMatrix D = null;

            if (Material is FGMStructureOneVariableMaterial)
            {
                if (GetTypeBeam() == TypeBeam.HSDT)
                {
                    D = CreateMaterialMatrixHSDT();
                }
                else
                {
                    D = CreateMaterialMatrixFSDTFGM();
                }
            }
            else
            {
                if (GetTypeBeam() == TypeBeam.HSDT)
                {
                    D = CreateMaterialMatrixHSDT();
                }
                else
                {
                    D = CreateMaterialMatrixFSDT();
                }
            }

            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
            {
                double detJbar = 1.0 / 2.0 * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1]);
                DoubleMatrix dNdxi = (DoubleMatrix)gps[i].GetValue(DataInGausspoint.dNdxi);//gps[i].dNdxi;
                DoubleMatrix J = JacobianAt(dNdxi);
                DoubleMatrix dNdx = (DoubleMatrix)gps[i].GetValue(DataInGausspoint.dNdX);
                DoubleMatrix ddNdx = (DoubleMatrix)gps[i].GetValue(DataInGausspoint.ddNdX);//gps[i].ddNdX;
                DoubleVector Ni = (DoubleVector)gps[i].GetValue(DataInGausspoint.Ni);

                for (int k = 0; k < nen; k++)
                {
                    if (GetTypeBeam() == TypeBeam.FSDT)
                    {
                        //B[0, 3 * k] = dNdx[k, 0];
                        //B[1, 3 * k + 2] = dNdx[k, 0];
                        //B[2, 3 * k + 2] = Ni[k];
                        B[0, 3 * k] = dNdx[k, 0];
                        B[1, 3 * k + 1] = -ddNdx[k, 0];
                        B[2, 3 * k + 2] = dNdx[k, 0];
                        Bs[0, 3 * k + 2] = Ni[k];
                    }
                    else
                    {
                        B[0, 3 * k] = dNdx[k, 0];
                        B[1, 3 * k + 1] = -ddNdx[k, 0];
                        B[2, 3 * k + 2] = dNdx[k, 0];
                        B[3, 3 * k + 2] = Ni[k];
                    }
                }
                DoubleMatrix BTDB = NMathFunctions.TransposeProduct(B, NMathFunctions.Product(D, B));
                DoubleMatrix BsTDsBs = NMathFunctions.TransposeProduct(Bs, Ds * Bs);
                Ke += gps[i].weight * detJbar * Math.Abs((double)gps[i].GetValue(DataInGausspoint.detJ)) * (BTDB + BsTDsBs);
            }
            //DoubleMatrix KeSub = new DoubleMatrix(nen, nen);
            //for (int i = 0; i < nen; i++)
            //    for (int j = 0; j < nen; j++)
            //    {
            //        KeSub[i, j] = Ke[3 * i + 1, 3 * j + 1];
            //    }
            //Ke = KeSub;
        }

        public override DoubleVector StrainAt(params double[] xi)
        {
            PatchStructureBeam patch = (PatchStructureBeam)this.patch;
            int d = patch.GetCountField();
            int nen = patch.GetCountLocalBasisFunctions();
            var cps = ((NURBSSurface)(patch.GetGeometry())).ControlPoints;
            DoubleMatrix B;
            if (GetTypeBeam() == TypeBeam.HSDT)
            {
                B = new DoubleMatrix(3, d * nen);
            }
            else
            {
                B = new DoubleMatrix(4, d * nen);
            }

            DoubleMatrix gradNE = GradBasisFunction(xi[0]);
            DoubleMatrix J = JacobianAt(gradNE);
            DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
            DoubleMatrix gradBasis = gradNE * invertJ;
            DoubleVector ValueBasis = ValueBasisFunction(xi[0]);
            DoubleMatrix grad2Basis = Grad2BasisFunction(xi[0]);
            for (int k = 0; k < nen; k++)
            {
                if (GetTypeBeam() == TypeBeam.FSDT)
                {
                    B[0, 3 * k] = gradBasis[k, 0];
                    B[1, 3 * k + 2] = gradBasis[k, 0];
                    B[2, 3 * k + 2] = ValueBasis[k];
                }
                else
                {
                    B[0, 3 * k] = gradBasis[k, 0];
                    B[1, 3 * k + 1] = grad2Basis[k, 0];
                    B[2, 3 * k + 2] = gradBasis[k, 0];
                    B[3, 3 * k + 2] = ValueBasis[k];
                }
            }
            DoubleVector U = GetDisplacementLocal();
            DoubleVector Smatrix = MatrixFunctions.Product(B, U);
            DoubleVector Rmatrix = new DoubleVector(5);
            double zz = getzz();
            return Rmatrix;
        }

        public override DoubleVector StressAt(params double[] xi)
        {
            double zz = getzz();
            DoubleMatrix D = ReMaterial(zz);
            DoubleVector strain = StrainAt(xi);
            return MatrixFunctions.Product(D, strain);
        }

        private DoubleMatrix ReMaterial(double zz)
        {
            DoubleMatrix Q = new DoubleMatrix(3, 3);
            double eModulus = 0;
            double nu = 0;
            if (Material is FGMStructureOneVariableMaterial)
            {
                eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty(zz);
                nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty(zz);
            }
            else
            {
                eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
                nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
            }
            double aa = eModulus / (1 - nu * nu);
            double bb = eModulus / (2.0 * (1 + nu));
            Q[0, 0] = Q[1, 1] = aa;
            Q[0, 1] = Q[1, 0] = aa * nu;
            Q[2, 2] = bb;
            return Q; //5.39
        }


        //public override void UpdateStressGaussPoint()
        //{
        //	for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        //		for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
        //		{
        //			GaussPoints gauss = gps[i, j];
        //			double[] locGauss = gauss.location;
        //			gauss.lastStress = StressAt(locGauss[0], locGauss[1]);
        //		}
        //}

        public override void ComputeInternalForceElement(out DoubleVector fi)
        {
            throw new NotImplementedException();
        }

        private double Ee;
        private double xPhys;

        public void SetCurrentModulus(double Ee)
        {
            this.Ee = Ee;
        }

        public double GetCurrentModulus()
        {
            return Ee;
        }

        public void SetDensityFilter(double xPhys)
        {
            this.xPhys = xPhys;
        }

        public double GetDensityFilter()
        {
            return xPhys;
        }

        public override void ComputeGeometricStiffnessMatrixElement(double[,] stressTensor, out DoubleMatrix Ke)
        {
            throw new NotImplementedException();
        }
    }
}
