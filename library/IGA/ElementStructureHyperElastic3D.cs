using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public class ElementStructureHyperElastic3D : AbstractElementStructure3D
  {
    public ElementStructureHyperElastic3D(AbstractPatch3D mesh, int id)
        : base(mesh, id)
    {
    }

    private double[,,,] CreateMaterialMatrix(DoubleMatrix b, double detF)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      double mu1 = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
      double K1 = Material.GetProperty(MaterialPropertyName.BulkModulus).GetValueProperty();
      int d = patch.GetCountField();
      double[,,,] C = new double[d, d, d, d];
      DoubleMatrix deltaIJ = new DoubleMatrix(3, 3);
      deltaIJ[0, 0] = deltaIJ[1, 1] = deltaIJ[2, 2] = 1.0;

      double Bqq = b[0, 0] + b[1, 1] + 1;
      if (Material.GetProperty(MaterialPropertyName.NeoHookean) != null)
        for (int i = 0; i < d; i++)
          for (int j = 0; j < d; j++)
            for (int k = 0; k < d; k++)
              for (int l = 0; l < d; l++)
              {
                C[i, j, k, l] = mu1 / Math.Pow(detF, 2.0 / 3.0) *
                    (deltaIJ[i, k] * b[j, l] + b[i, l] * deltaIJ[j, k]
                    - (2.0 / 3.0) * (b[i, j] * deltaIJ[k, l] + deltaIJ[i, j] * b[k, l])
                    + (2.0 / 3.0) * Bqq * deltaIJ[i, j] * deltaIJ[k, l] / 3.0)
                    + K1 * (2.0 * detF - 1) * detF * deltaIJ[i, j] * deltaIJ[k, l];
              }
      if (Material.GetProperty(MaterialPropertyName.MooneyRivlin2Parameter) != null)
        for (int i = 0; i < d; i++)
          for (int j = 0; j < d; j++)
            for (int k = 0; k < d; k++)
              for (int l = 0; l < d; l++)
              {////////////////////////////////////////////////////////
                C[i, j, k, l] = mu1 / Math.Pow(detF, 2.0 / 3.0) *
                    (deltaIJ[i, k] * b[j, l] + b[i, l] * deltaIJ[j, k]
                    - (2.0 / 3.0) * (b[i, j] * deltaIJ[k, l] + deltaIJ[i, j] * b[k, l])
                    + (2.0 / 3.0) * Bqq * deltaIJ[i, j] * deltaIJ[k, l] / 3.0)
                    + K1 * (2.0 * detF - 1) * detF * deltaIJ[i, j] * deltaIJ[k, l];
              }
      return C;
    }

    private DoubleMatrix DeformationGradient(DoubleMatrix dNdx)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      DoubleMatrix F = new DoubleMatrix(3, 3);
      F[0, 0] = 1.0;
      F[1, 1] = 1.0;
      F[2, 2] = 1.0;
      int nen = patch.GetCountLocalBasisFunctions();
      NURBSVolume volume = patch.GetVolume();
      var cps = volume.ControlPoints;
      for (int i = 0; i < patch.GetCountField(); i++)
        for (int j = 0; j < patch.GetCountField(); j++)
          for (int k = 0; k < nen; k++)
          {
            int ien = patch.GetIEN(id, k);
            var currentCp = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(ien, 2)];
            F[i, j] += currentCp.GetULocal(i) * dNdx[k, j];
          }
      return F;
    }

    private double[,] DerivativeCoord(DoubleMatrix dNdx, DoubleMatrix F)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      DoubleMatrix Finv = MatrixFunctions.Inverse(F);
      int nen = patch.GetCountLocalBasisFunctions();
      int d = patch.GetCountField();
      double[,] dNdxs = new double[nen, d];
      for (int k = 0; k < nen; k++)
        for (int i = 0; i < d; i++)
          for (int j = 0; j < d; j++)
          {
            dNdxs[k, i] += dNdx[k, j] * Finv[j, i];
          }
      return dNdxs;
    }

    private DoubleMatrix KirchhoffStress(DoubleMatrix b, double detF)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      int d = patch.GetCountField();
      DoubleMatrix stress = new DoubleMatrix(d, d);
      DoubleMatrix deltaIJ = new DoubleMatrix(3, 3);
      deltaIJ[0, 0] = deltaIJ[1, 1] = deltaIJ[2, 2] = 1.0;
      double Bkk = b[0, 0] + b[1, 1] + b[2, 2];
      if (Material.GetProperty(MaterialPropertyName.NeoHookean) != null)
      {
        double mu1 = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
        double K1 = Material.GetProperty(MaterialPropertyName.BulkModulus).GetValueProperty();
        for (int i = 0; i < d; i++)
          for (int j = 0; j < d; j++)
            stress[i, j] = mu1 * (b[i, j] - Bkk * deltaIJ[i, j] / 3.0) / Math.Pow(detF, 2.0 / 3.0) + K1 * detF * (detF - 1) * deltaIJ[i, j];
      }
      if (Material.GetProperty(MaterialPropertyName.MooneyRivlin2Parameter) != null)
      {
        double mu1 = 2 * Material.GetProperty(MaterialPropertyName.C10).GetValueProperty();
        double mu2 = 2 * Material.GetProperty(MaterialPropertyName.C01).GetValueProperty();
        double K1 = 2.0 / Material.GetProperty(MaterialPropertyName.D1).GetValueProperty();
        for (int i = 0; i < d; i++)
          for (int j = 0; j < d; j++)
            stress[i, j] = mu1 / Math.Pow(detF, 2.0 / 3.0) * (b[i, j] - Bkk * deltaIJ[i, j] / 3.0)
                + mu2 / Math.Pow(detF, 4.0 / 3.0) * (Bkk * b[i, j] - Bkk * Bkk * deltaIJ[i, j] / 3.0
                - (b[i, 0] * b[0, j] + b[i, 1] * b[1, j] + b[i, 2] * b[2, j])
                + (b[0, 0] * b[0, 0] + b[1, 0] * b[0, 1] + b[2, 0] * b[0, 2]
                + b[0, 1] * b[1, 0] + b[1, 1] * b[1, 1] + b[2, 1] * b[1, 2]
                + b[0, 2] * b[2, 0] + b[1, 2] * b[2, 1] + b[2, 2] * b[2, 2]))
                + K1 * detF * (detF - 1) * deltaIJ[i, j];
      }
      return stress;
    }

    public override DoubleVector StressAt(double[] xi)
    {
      DoubleMatrix dNdxi = GradBasisFunction(xi[0], xi[1], xi[2]);
      DoubleMatrix J = JacobianAt(dNdxi);
      DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
      DoubleMatrix dNdx = dNdxi * invertJ;

      DoubleMatrix F = DeformationGradient(dNdx);
      double detF = MatrixFunctions.Determinant(F);
      DoubleMatrix b = F * F.Transpose();//Bbar

      int d = ((PatchStructure3D)patch).GetCountField();
      double[,] stress = new double[d, d];
      double[,] deltaIJ = new double[3, 3];
      deltaIJ[0, 0] = 1.0; deltaIJ[1, 1] = 1.0; deltaIJ[2, 2] = 1.0;
      double Bkk = b[0, 0] + b[1, 1] + b[2, 2];
      if (Material.GetProperty(MaterialPropertyName.NeoHookean) != null)
      {
        double mu1 = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
        double K1 = Material.GetProperty(MaterialPropertyName.BulkModulus).GetValueProperty();
        for (int i = 0; i < d; i++)
          for (int j = 0; j < d; j++)
            stress[i, j] = mu1 * (b[i, j] - Bkk * deltaIJ[i, j] / 3.0) / Math.Pow(detF, 2.0 / 3.0) + K1 * detF * (detF - 1) * deltaIJ[i, j];
      }
      if (Material.GetProperty(MaterialPropertyName.MooneyRivlin2Parameter) != null)
      {
        double mu1 = 2 * Material.GetProperty(MaterialPropertyName.C10).GetValueProperty();
        double mu2 = 2 * Material.GetProperty(MaterialPropertyName.C01).GetValueProperty();
        double K1 = 2.0 / Material.GetProperty(MaterialPropertyName.D1).GetValueProperty();
        for (int i = 0; i < d; i++)
          for (int j = 0; j < d; j++)
            stress[i, j] = mu1 / Math.Pow(detF, 2.0 / 3.0) * (b[i, j] - Bkk * deltaIJ[i, j] / 3.0)
                + mu2 / Math.Pow(detF, 4.0 / 3.0) * (Bkk * b[i, j] - Bkk * Bkk * deltaIJ[i, j] / 3.0
                - (b[i, 0] * b[0, j] + b[i, 1] * b[1, j] + b[i, 2] * b[2, j])
                + (b[0, 0] * b[0, 0] + b[1, 0] * b[0, 1] + b[2, 0] * b[0, 2]
                + b[0, 1] * b[1, 0] + b[1, 1] * b[1, 1] + b[2, 1] * b[1, 2]
                + b[0, 2] * b[2, 0] + b[1, 2] * b[2, 1] + b[2, 2] * b[2, 2]))
                + K1 * detF * (detF - 1) * deltaIJ[i, j];
      }
      return new DoubleVector(new double[] { stress[0, 0] / detF, stress[1, 1] / detF, stress[2, 2] / detF, stress[0, 1] / detF, stress[1, 2] / detF, stress[0, 2] / detF });
    }

    public override void ComputeStiffnessMatrixElement(out DoubleMatrix Ke)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      int d = patch.GetCountField();
      TrivariateNURBSBasisFunction basis = (TrivariateNURBSBasisFunction)((NURBSVolume)(patch.GetGeometry())).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int nen = (p + 1) * (q + 1) * (r + 1); // number of local basis functions
      var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply3 = basis.GetKnotVector(2).GetKnotVectorNoMultiplicity();
      int idx1 = patch.GetIPN(id, 0);
      int idx2 = patch.GetIPN(id, 1);
      int idx3 = patch.GetIPN(id, 2);
      Ke = new DoubleMatrix(d * nen, d * nen);
      //DoubleMatrix B
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {
            double detJbar = 1.0 / 8.0
                * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2])
                * (kvNoMulticiply3[idx3 + 1] - kvNoMulticiply3[idx3]);
            DoubleMatrix dNdxi = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdxi);
            DoubleMatrix J = JacobianAt(dNdxi);
            double detJ = Math.Abs(MatrixFunctions.Determinant(J));
            DoubleMatrix dNdx = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdX);//dNdxi * invertJ;

            DoubleMatrix F = DeformationGradient(dNdx);
            double detF = MatrixFunctions.Determinant(F);
            DoubleMatrix B = F * F.Transpose();//Bbar

            double[,] dNdxs = DerivativeCoord(dNdx, F);
            DoubleMatrix P = KirchhoffStress(B, detF);

            double[,,,] Cijkl = CreateMaterialMatrix(B, detF);

            double w = gps[i, j, k].weight * detJbar * detJ;
            for (int a = 0; a < nen; a++)
              for (int ii = 0; ii < d; ii++)
              {
                int row = d * a + ii;
                for (int b = 0; b < nen; b++)
                  for (int kk = 0; kk < d; kk++)
                  {
                    int col = d * b + kk;
                    for (int jj = 0; jj < d; jj++)
                    {
                      for (int l = 0; l < d; l++)
                        Ke[row, col] += Cijkl[ii, jj, kk, l] * dNdxs[b, l] * dNdxs[a, jj] * w;
                      Ke[row, col] -= P[ii, jj] * dNdxs[a, kk] * dNdxs[b, jj] * w;
                    }
                  }
              }
          }
        }
      }
      //for (int j = 0; j < d * nen; j++)
      //    for (int i = j + 1; i < d * nen; i++)
      //    {
      //        Ke[i, j] = Ke[j, i];
      //    }
    }

    public override void ComputeInternalForceElement(out DoubleVector fi)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      int d = patch.GetCountField();
      TrivariateNURBSBasisFunction basis = (TrivariateNURBSBasisFunction)((NURBSVolume)(patch.GetGeometry())).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int nen = (p + 1) * (q + 1) * (r + 1); // number of local basis functions
      var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply3 = basis.GetKnotVector(2).GetKnotVectorNoMultiplicity();
      int idx1 = patch.GetIPN(id, 0);
      int idx2 = patch.GetIPN(id, 1);
      int idx3 = patch.GetIPN(id, 2);
      fi = new DoubleVector(d * nen);
      //DoubleMatrix B
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {
            double detJbar = 1.0 / 8.0
                * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2])
                * (kvNoMulticiply3[idx3 + 1] - kvNoMulticiply3[idx3]);
            DoubleMatrix dNdxi = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdxi);
            DoubleMatrix J = JacobianAt(dNdxi);
            double detJ = Math.Abs(MatrixFunctions.Determinant(J));
            DoubleMatrix dNdx = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdX);

            DoubleMatrix F = DeformationGradient(dNdx);
            double detF = MatrixFunctions.Determinant(F);
            DoubleMatrix B = F * F.Transpose();//Bbar

            double[,] dNdxs = DerivativeCoord(dNdx, F);
            DoubleMatrix P = KirchhoffStress(B, detF);

            double w = gps[i, j, k].weight * detJbar * detJ;
            for (int a = 0; a < nen; a++)
              for (int ii = 0; ii < d; ii++)
              {
                int row = d * a + ii;
                for (int jj = 0; jj < d; jj++)
                  fi[row] += P[ii, jj] * dNdxs[a, jj] * w;
              }
          }
        }
      }
    }

    public override DoubleVector StrainAt(double[] xi)
    {
      throw new NotImplementedException();
    }

    public override DoubleVector StrainElasticAt(params double[] xi)
    {
      throw new NotImplementedException();
    }

    public override DoubleVector StrainThermoAt(params double[] xi)
    {
      throw new NotImplementedException();
    }
  }
}
