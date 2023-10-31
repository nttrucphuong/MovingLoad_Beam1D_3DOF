using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using DEMSoft.Function;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
    public abstract class AbstractElementStructureBeam : AbstractElement1DBeam, IElementStructure
    {
        public AbstractElementStructureBeam(AbstractPatch1D mesh, int id)
              : base(mesh, id)
        {
        }

        public override void ComputeMassMatrixElement(out DoubleMatrix Me)
        {
            PatchStructureBeam patch = (PatchStructureBeam)this.patch;
            int d = patch.GetCountField();
            NURBSBasisFunction basis = (NURBSBasisFunction)((NURBSCurve)(patch.GetGeometry())).Basis;
            int p = basis.GetDegree(0);
            int nen = (p + 1); // number of local basis functions
            var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
            int idx1 = patch.GetIPN(id, 0);
            Me = new DoubleMatrix(d * nen, d * nen);
            DoubleMatrix MM = null;
            if (Material is FGMStructureOneVariableMaterial)
            {
                if (patch.TypeBeam == TypeBeam.FSDT)
                {
                    MM = CreateMatrixMFSDTFGM(patch.Thickness, patch.KinematicsFunction);
                }
                else
                {
                    MM = CreateMatrixMHSDTFGM(patch.Thickness, patch.KinematicsFunction);
                }
            }
            else
            {
                if (patch.TypeBeam == TypeBeam.FSDT)
                {
                    MM = CreateMatrixMFSDT(patch.Thickness, patch.KinematicsFunction);
                }
                else
                {
                    MM = CreateMatrixMHSDT(patch.Thickness, patch.KinematicsFunction);
                }
            }
            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
            {
                for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
                {
                    double detJbar = 1.0 / 4.0 * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1]);
                    //DoubleMatrix dNdxi = gps[i, j].dNdxi;
                    //DoubleMatrix J = JacobianAt(dNdxi);
                    GaussPoints gpsij = gps[i];
                    DoubleMatrix dNdx = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.dNdX);//gps[i, j].dNdX;
                    DoubleVector Ni = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni); //gps[i, j].Ni;
                    DoubleMatrix Ri = null;
                    if (patch.TypeBeam == TypeBeam.FSDT)
                    {
                        Ri = new DoubleMatrix(4, d * nen);
                    }
                    else
                    {
                        Ri = new DoubleMatrix(6, d * nen);
                    }

                    for (int kk = 0; kk < nen; kk++) //7.37
                    {
                        if (patch.TypeBeam == TypeBeam.FSDT)
                        {
                            Ri[0, 3 * kk] = Ni[kk];
                            Ri[0, 3 * kk + 1] = Ni[kk];

                            Ri[1, 3 * kk + 2] = Ni[kk];
                        }
                        else
                        {
                            Ri[0, 3 * kk] = Ni[kk];
                            Ri[0, 3 * kk + 1] = Ni[kk];

                            Ri[1, 3 * kk + 1] = -dNdx[kk, 0];

                            Ri[2, 3 * kk + 2] = Ni[kk];
                        }
                    }
                    DoubleMatrix NeTNe = NMathFunctions.TransposeProduct(Ri, MatrixFunctions.Product(MM, Ri));

                    Me += gpsij.weight * detJbar * Math.Abs((double)gpsij.GetValue(DataInGausspoint.detJ)) * NeTNe;
                }
            }
        }

        private DoubleMatrix CreateMatrixMFSDT(double thickness, FunctionRToR KinematicsFunction)
        {
            DoubleMatrix M = new DoubleMatrix(4, 4);
            int numGauss = 7;
            for (int i = 0; i < numGauss; i++)
            {
                double xi = GaussPoints.GetPoint(numGauss, i);
                double wi = GaussPoints.GetWeight(numGauss, i);
                double zz = 0.5 * thickness * (xi + 1.0) - thickness / 2.0;
                double detJ = thickness / 2.0;
                double fz = KinematicsFunction.ValueAt(zz);
                //var rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty();

                M[0, 0] += detJ * wi;
                M[0, 1] += detJ * wi * zz;
                M[1, 1] += detJ * wi * fz;

                M[2, 2] += detJ * wi;
                M[2, 3] += detJ * wi * zz;
                M[3, 3] += detJ * wi * fz;
            }
            M[1, 0] = M[0, 1];
            M[3, 2] = M[2, 3];

            return M;
        }
        private DoubleMatrix CreateMatrixMFSDTFGM(double thickness, FunctionRToR KinematicsFunction)
        {
            DoubleMatrix M = new DoubleMatrix(4, 4);
            int numGauss = 20;
            for (int i = 0; i < numGauss; i++)
            {
                double xi = GaussPoints.GetPoint(numGauss, i);
                double wi = GaussPoints.GetWeight(numGauss, i);
                double zz = 0.5 * thickness * (xi + 1.0) - thickness / 2.0;
                double detJ = thickness / 2.0;
                double fz = KinematicsFunction.ValueAt(zz);
                //var rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty();

                M[0, 0] += detJ * wi;
                M[0, 1] += detJ * wi * zz;
                M[1, 1] += detJ * wi * fz;

                M[2, 2] += detJ * wi;
                M[2, 3] += detJ * wi * zz;
                M[3, 3] += detJ * wi * fz;
            }
            M[1, 0] = M[0, 1];
            M[3, 2] = M[2, 3];
            return M;
        }
        private DoubleMatrix CreateMatrixMHSDT(double thickness, FunctionRToR KinematicsFunction)
        {
            DoubleMatrix M = new DoubleMatrix(6, 6);
            int numGauss = 7;
            for (int i = 0; i < numGauss; i++)
            {
                double xi = GaussPoints.GetPoint(numGauss, i);
                double wi = GaussPoints.GetWeight(numGauss, i);
                double zz = 0.5 * thickness * (xi + 1.0) - thickness / 2.0;
                double fz = KinematicsFunction.ValueAt(zz);
                double detJ = thickness / 2.0;
                //var rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty();

                M[0, 0] += detJ * wi;
                M[0, 1] += detJ * wi * zz;
                M[0, 2] += detJ * wi * fz;
                M[1, 1] += detJ * wi * zz * zz;
                M[1, 2] += detJ * wi * zz * fz;
                M[2, 2] += detJ * wi * fz * fz;

                M[3, 3] += detJ * wi;
                M[3, 4] += detJ * wi * zz;
                M[3, 5] += detJ * wi * fz;
                M[4, 4] += detJ * wi * zz * zz;
                M[4, 5] += detJ * wi * zz * fz;
                M[5, 5] += detJ * wi * fz * fz;
            }

            M[1, 0] = M[0, 1];
            M[2, 0] = M[0, 2];
            M[2, 1] = M[1, 2];

            M[4, 3] = M[3, 4];
            M[5, 3] = M[3, 5];
            M[5, 4] = M[4, 5];

            return M;
        }
        private DoubleMatrix CreateMatrixMHSDTFGM(double thickness, FunctionRToR KinematicsFunction)
        {
            DoubleMatrix M = new DoubleMatrix(6, 6);
            //int numGauss = 20;
            //for (int i = 0; i < numGauss; i++)
            //{
            //    double xi = GaussPoints.GetPoint(numGauss, i);
            //    double wi = GaussPoints.GetWeight(numGauss, i);
            //    double zz = 0.5 * thickness * (xi + 1.0) - thickness / 2.0;
            //    double fz = KinematicsFunction.ValueAt(zz);
            //    double detJ = thickness / 2.0;
            //    //var rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty();

            //    M[0, 0] += detJ * wi;
            //    M[0, 1] += detJ * wi * zz;
            //    M[0, 2] += detJ * wi * fz;
            //    M[1, 1] += detJ * wi * zz * zz;
            //    M[1, 2] += detJ * wi * zz * fz;
            //    M[2, 2] += detJ * wi * fz * fz;

            //    M[3, 3] += detJ * wi;
            //    M[3, 4] += detJ * wi * zz;
            //    M[3, 5] += detJ * wi * fz;
            //    M[4, 4] += detJ * wi * zz * zz;
            //    M[4, 5] += detJ * wi * zz * fz;
            //    M[5, 5] += detJ * wi * fz * fz;
            //}

            //M[1, 0] = M[0, 1];
            //M[2, 0] = M[0, 2];
            //M[2, 1] = M[1, 2];

            //M[4, 3] = M[3, 4];
            //M[5, 3] = M[3, 5];
            //M[5, 4] = M[4, 5];
            return M;
        }
        public abstract DoubleVector StressAt(params double[] xi);
        public abstract DoubleVector StrainAt(params double[] xi);

        public DoubleVector StrainElasticAt(params double[] xi)
        {
            throw new NotImplementedException();
        }

        public DoubleVector StrainThermoAt(params double[] xi)
        {
            throw new NotImplementedException();
        }
        public override void ComputeValueAtGaussPoint(DataInGausspoint name)
        {
            for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
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
