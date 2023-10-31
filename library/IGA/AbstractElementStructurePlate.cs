using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using DEMSoft.Function;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public abstract class AbstractElementStructurePlate : AbstractElement2DPlate, IElementStructure
  {
    public AbstractElementStructurePlate(AbstractPatch2D mesh, int id)
          : base(mesh, id)
    {
    }

    public override void ComputeMassMatrixElement(out DoubleMatrix Me)
    {
      PatchStructurePlate patch = (PatchStructurePlate)this.patch;
      int d = patch.GetCountField();
      BivariateNURBSBasisFunction basis = (BivariateNURBSBasisFunction)((NURBSSurface)(patch.GetGeometry())).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int nen = (p + 1) * (q + 1); // number of local basis functions
      var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      int idx1 = patch.GetIPN(id, 0);
      int idx2 = patch.GetIPN(id, 1);
      Me = new DoubleMatrix(d * nen, d * nen);
      DoubleMatrix MM = null;
      if (Material is FGMStructureOneVariableMaterial)
      {
        if (patch.TypePlate == TypePlate.FSDT)
        {
          MM = CreateMatrixMFSDTFGM(patch.Thickness);
        }
        else
        {
          MM = CreateMatrixMHSDTFGM(patch.Thickness, patch.KinematicsFunction);
        }
      }
      else
      {
        if (patch.TypePlate == TypePlate.FSDT)
        {
          MM = CreateMatrixMFSDT(patch.Thickness);
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
          double detJbar = 1.0 / 4.0
              * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
              * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          //DoubleMatrix dNdxi = gps[i, j].dNdxi;
          //DoubleMatrix J = JacobianAt(dNdxi);
          GaussPoints gpsij = gps[i, j];
          DoubleMatrix dNdx = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.dNdX);//gps[i, j].dNdX;
          DoubleVector Ni = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni); //gps[i, j].Ni;
          DoubleMatrix Ri = null;
          if (patch.TypePlate == TypePlate.FSDT)
          {
            Ri = new DoubleMatrix(6, d * nen);
          }
          else
          {
            Ri = new DoubleMatrix(9, d * nen);
          }

          for (int kk = 0; kk < nen; kk++)
          {
            if (patch.TypePlate == TypePlate.FSDT)
            {
              Ri[0, 5 * kk] = Ni[kk];
              Ri[1, 5 * kk + 3] = Ni[kk];

              Ri[2, 5 * kk + 1] = Ni[kk];
              Ri[3, 5 * kk + 4] = Ni[kk];

              Ri[4, 5 * kk + 2] = Ni[kk];
            }
            else
            {
              Ri[0, 5 * kk] = Ni[kk];
              Ri[1, 5 * kk + 2] = -dNdx[kk, 0];
              Ri[2, 5 * kk + 3] = Ni[kk];

              Ri[3, 5 * kk + 1] = Ni[kk];
              Ri[4, 5 * kk + 2] = -dNdx[kk, 1];
              Ri[5, 5 * kk + 4] = Ni[kk];

              Ri[6, 5 * kk + 2] = Ni[kk];
            }
          }
          DoubleMatrix NeTNe = NMathFunctions.TransposeProduct(Ri, MatrixFunctions.Product(MM, Ri));

          Me += gpsij.weight * detJbar * Math.Abs((double)gpsij.GetValue(DataInGausspoint.detJ)) * NeTNe;
        }
      }
    }


    private DoubleMatrix CreateMatrixMFSDT(double thickness)
    {
      DoubleMatrix M = new DoubleMatrix(6, 6);
      int numGauss = 7;
      for (int i = 0; i < numGauss; i++)
      {
        double xi = GaussPoints.GetPoint(numGauss, i);
        double wi = GaussPoints.GetWeight(numGauss, i);
        double zz = 0.5 * thickness * (xi + 1.0) - thickness / 2.0;
        double detJ = thickness / 2.0;
        var rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty();

        M[0, 0] += detJ * wi * rho;
        M[0, 1] += detJ * wi * rho * zz;
        M[1, 1] += detJ * wi * rho * zz * zz;

        M[2, 2] += detJ * wi * rho;
        M[2, 3] += detJ * wi * rho * zz;
        M[3, 3] += detJ * wi * rho * zz * zz;

        M[4, 4] += detJ * wi * rho;
        M[4, 5] += detJ * wi * rho * zz;
        M[5, 5] += detJ * wi * rho * zz * zz;
      }
      M[1, 0] = M[0, 1];
      M[3, 2] = M[2, 3];
      M[5, 4] = M[4, 5];
      return M;
    }
    private DoubleMatrix CreateMatrixMFSDTFGM(double thickness)
    {
      DoubleMatrix M = new DoubleMatrix(6, 6);
      int numGauss = 20;
      for (int i = 0; i < numGauss; i++)
      {
        double xi = GaussPoints.GetPoint(numGauss, i);
        double wi = GaussPoints.GetWeight(numGauss, i);
        double zz = 0.5 * thickness * (xi + 1.0) - thickness / 2.0;
        double detJ = thickness / 2.0;
        var rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty(zz);

        M[0, 0] += detJ * wi * rho;
        M[0, 1] += detJ * wi * rho * zz;
        M[1, 1] += detJ * wi * rho * zz * zz;

        M[2, 2] += detJ * wi * rho;
        M[2, 3] += detJ * wi * rho * zz;
        M[3, 3] += detJ * wi * rho * zz * zz;

        M[4, 4] += detJ * wi * rho;
        M[4, 5] += detJ * wi * rho * zz;
        M[5, 5] += detJ * wi * rho * zz * zz;
      }
      M[1, 0] = M[0, 1];
      M[3, 2] = M[2, 3];
      M[5, 4] = M[4, 5];
      return M;
    }
    private DoubleMatrix CreateMatrixMHSDT(double thickness, FunctionRToR KinematicsFunction)
    {
      DoubleMatrix M = new DoubleMatrix(9, 9);
      int numGauss = 7;
      for (int i = 0; i < numGauss; i++)
      {
        double xi = GaussPoints.GetPoint(numGauss, i);
        double wi = GaussPoints.GetWeight(numGauss, i);
        double zz = 0.5 * thickness * (xi + 1.0) - thickness / 2.0;
        double fz = KinematicsFunction.ValueAt(zz);
        double detJ = thickness / 2.0;
        var rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty();

        M[0, 0] += detJ * wi * rho;
        M[0, 1] += detJ * wi * rho * zz;
        M[0, 2] += detJ * wi * rho * fz;
        M[1, 1] += detJ * wi * rho * zz * zz;
        M[1, 2] += detJ * wi * rho * zz * fz;
        M[2, 2] += detJ * wi * rho * fz * fz;

        M[3, 3] += detJ * wi * rho;
        M[3, 4] += detJ * wi * rho * zz;
        M[3, 5] += detJ * wi * rho * fz;
        M[4, 4] += detJ * wi * rho * zz * zz;
        M[4, 5] += detJ * wi * rho * zz * fz;
        M[5, 5] += detJ * wi * rho * fz * fz;

        M[6, 6] += detJ * wi * rho;
        M[6, 7] += detJ * wi * rho * zz;
        M[6, 8] += detJ * wi * rho * fz;
        M[7, 7] += detJ * wi * rho * zz * zz;
        M[7, 8] += detJ * wi * rho * zz * fz;
        M[8, 8] += detJ * wi * rho * fz * fz;
      }
      M[1, 0] = M[0, 1];
      M[2, 0] = M[0, 2];
      M[2, 1] = M[1, 2];

      M[4, 3] = M[3, 4];
      M[5, 3] = M[3, 5];
      M[5, 4] = M[4, 5];

      M[7, 6] = M[6, 7];
      M[8, 6] = M[6, 8];
      M[8, 7] = M[7, 8];
      return M;
    }
    private DoubleMatrix CreateMatrixMHSDTFGM(double thickness, FunctionRToR KinematicsFunction)
    {
      DoubleMatrix M = new DoubleMatrix(9, 9);
      int numGauss = 20;
      for (int i = 0; i < numGauss; i++)
      {
        double xi = GaussPoints.GetPoint(numGauss, i);
        double wi = GaussPoints.GetWeight(numGauss, i);
        double zz = 0.5 * thickness * (xi + 1.0) - thickness / 2.0;
        double fz = KinematicsFunction.ValueAt(zz);
        double detJ = thickness / 2.0;
        var rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty(zz);

        M[0, 0] += detJ * wi * rho;
        M[0, 1] += detJ * wi * rho * zz;
        M[0, 2] += detJ * wi * rho * fz;
        M[1, 1] += detJ * wi * rho * zz * zz;
        M[1, 2] += detJ * wi * rho * zz * fz;
        M[2, 2] += detJ * wi * rho * fz * fz;

        M[3, 3] += detJ * wi * rho;
        M[3, 4] += detJ * wi * rho * zz;
        M[3, 5] += detJ * wi * rho * fz;
        M[4, 4] += detJ * wi * rho * zz * zz;
        M[4, 5] += detJ * wi * rho * zz * fz;
        M[5, 5] += detJ * wi * rho * fz * fz;

        M[6, 6] += detJ * wi * rho;
        M[6, 7] += detJ * wi * rho * zz;
        M[6, 8] += detJ * wi * rho * fz;
        M[7, 7] += detJ * wi * rho * zz * zz;
        M[7, 8] += detJ * wi * rho * zz * fz;
        M[8, 8] += detJ * wi * rho * fz * fz;
      }
      M[1, 0] = M[0, 1];
      M[2, 0] = M[0, 2];
      M[2, 1] = M[1, 2];

      M[4, 3] = M[3, 4];
      M[5, 3] = M[3, 5];
      M[5, 4] = M[4, 5];

      M[7, 6] = M[6, 7];
      M[8, 6] = M[6, 8];
      M[8, 7] = M[7, 8];
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
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          switch (name)
          {
            case DataInGausspoint.currentStress:
              gps[i, j].SetValue(DataInGausspoint.currentStress, StressAt(gps[i, j].location));
              //gps[i, j].currentStress = StressAt(gps[i, j].location);
              break;
            case DataInGausspoint.currentStrain:
              gps[i, j].SetValue(DataInGausspoint.currentStrain, StrainAt(gps[i, j].location));
              //gps[i, j].currentStrain = StrainAt(gps[i, j].location);
              break;
          }
        }
      }
    }
  }
}
