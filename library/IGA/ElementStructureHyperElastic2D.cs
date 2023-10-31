using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public class ElementStructureHyperElastic2D : AbstractElementStructure2D
  {
    public ElementStructureHyperElastic2D(AbstractPatch2D mesh, int id)
        : base(mesh, id)
    {
    }

    private double[,,,] CreateMaterialMatrix(DoubleMatrix b, double detF)
    {
      PatchStructure2D patch = (PatchStructure2D)this.patch;
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
      return C;
    }

    private DoubleMatrix DeformationGradient(DoubleMatrix dNdx)
    {
      PatchStructure2D patch = (PatchStructure2D)this.patch;
      DoubleMatrix F = new DoubleMatrix(3, 3);
      F[0, 0] = F[1, 1] = F[2, 2] = 1.0;
      int nen = patch.GetCountLocalBasisFunctions();
      NURBSSurface surface = patch.GetSurface();
      var cps = surface.ControlPoints;
      for (int i = 0; i < patch.GetCountField(); i++)
        for (int j = 0; j < patch.GetCountField(); j++)
          for (int k = 0; k < nen; k++)
          {
            int ien = patch.GetIEN(id, k);
            var currentCp = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1)];
            F[i, j] += currentCp.GetULocal(i) * dNdx[k, j];
          }
      return F;
    }

    private double[,] DerivativeCoord(DoubleMatrix dNdx, DoubleMatrix F)
    {
      PatchStructure2D patch = (PatchStructure2D)this.patch;
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
      PatchStructure2D patch = (PatchStructure2D)this.patch;
      int d = patch.GetCountField();
      DoubleMatrix stress = new DoubleMatrix(d, d);
      DoubleMatrix deltaIJ = new DoubleMatrix(3, 3);
      deltaIJ[0, 0] = deltaIJ[1, 1] = deltaIJ[2, 2] = 1.0;

      if (Material.GetProperty(MaterialPropertyName.NeoHookean) != null)
      {
        double mu1 = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
        double K1 = Material.GetProperty(MaterialPropertyName.BulkModulus).GetValueProperty();
        double Bkk = b[0, 0] + b[1, 1] + b[2, 2];
        for (int i = 0; i < d; i++)
          for (int j = 0; j < d; j++)
            stress[i, j] = mu1 * (b[i, j] - Bkk * deltaIJ[i, j] / 3.0) / Math.Pow(detF, 2.0 / 3.0) + K1 * detF * (detF - 1) * deltaIJ[i, j];
      }
      if (Material.GetProperty(MaterialPropertyName.MooneyRivlin2Parameter) != null)
      {
        double mu1 = 2.0 * Material.GetProperty(MaterialPropertyName.C10).GetValueProperty();
        double mu2 = 2.0 * Material.GetProperty(MaterialPropertyName.C01).GetValueProperty();
        double K1 = 2.0 / Material.GetProperty(MaterialPropertyName.D1).GetValueProperty();
        //for (int i = 0; i < d; i++)
        //    for (int j = 0; j < d; j++)
        //        stress[i, j] = mu1 / Math.Pow(detF, 2.0 / 3.0) * (b[i, j] - Bkk * deltaIJ[i, j] / 3.0) + K1 * detF * (detF - 1) * deltaIJ[i, j];
      }
      return stress;
    }

    public override DoubleVector StressAt(params double[] xi)
    {
      DoubleMatrix dNdxi = GradBasisFunction(xi[0], xi[1]);
      DoubleMatrix J = JacobianAt(dNdxi);
      DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
      DoubleMatrix dNdx = dNdxi * invertJ;

      DoubleMatrix F = DeformationGradient(dNdx);
      double detF = MatrixFunctions.Determinant(F);
      DoubleMatrix B = F * F.Transpose();//Bbar

      int d = ((PatchStructure2D)patch).GetCountField();
      double[,] stress = new double[d, d];
      double[,] deltaIJ = new double[3, 3];
      deltaIJ[0, 0] = 1.0; deltaIJ[1, 1] = 1.0; deltaIJ[2, 2] = 1.0;
      double mu1 = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
      double K1 = Material.GetProperty(MaterialPropertyName.BulkModulus).GetValueProperty();
      double Bkk = B[0, 0] + B[1, 1] + B[2, 2];
      if (Material.GetProperty(MaterialPropertyName.NeoHookean) != null)
        for (int i = 0; i < d; i++)
          for (int j = 0; j < d; j++)
            stress[i, j] = mu1 * (B[i, j] - Bkk * deltaIJ[i, j] / 3.0) / Math.Pow(detF, 2.0 / 3.0) + K1 * detF * (detF - 1) * deltaIJ[i, j];

      return new DoubleVector(new double[] { stress[0, 0] / detF, stress[1, 1] / detF, stress[0, 1] / detF });
    }

    public override void ComputeStiffnessMatrixElement(out DoubleMatrix Ke)
    {
      PatchStructure2D patch = (PatchStructure2D)this.patch;
      int d = patch.GetCountField();
      BivariateNURBSBasisFunction basis = (BivariateNURBSBasisFunction)((NURBSSurface)(patch.GetGeometry())).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int nen = (p + 1) * (q + 1); // number of local basis functions
      var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      int idx1 = patch.GetIPN(id, 0);
      int idx2 = patch.GetIPN(id, 1);
      Ke = new DoubleMatrix(d * nen, d * nen);
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          double detJbar = 1.0 / 4.0
              * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
              * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          //DoubleMatrix dNdxi = gps[i, j].dNdxi;
          //DoubleMatrix J = JacobianAt(dNdxi);
          //DoubleMatrix invertJ = (DoubleMatrix)J.Inverse();
          double detJ = (double)gps[i, j].GetValue(DataInGausspoint.detJ);
          DoubleMatrix dNdx = (DoubleMatrix)gps[i, j].GetValue(DataInGausspoint.dNdX);//dNdxi * invertJ;

          DoubleMatrix F = DeformationGradient(dNdx);
          double detF = MatrixFunctions.Determinant(F);
          DoubleMatrix B = F * F.Transpose();//Bbar

          double[,] dNdxs = DerivativeCoord(dNdx, F);
          DoubleMatrix P = KirchhoffStress(B, detF);

          double[,,,] Cijkl = CreateMaterialMatrix(B, detF);

          double w = gps[i, j].weight * detJbar * Math.Abs(detJ);
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
                    for (int ll = 0; ll < d; ll++)
                      Ke[row, col] += Cijkl[ii, jj, kk, ll] * dNdxs[b, ll] * dNdxs[a, jj] * w;
                    Ke[row, col] -= P[ii, jj] * dNdxs[a, kk] * dNdxs[b, jj] * w;
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
      PatchStructure2D patch = (PatchStructure2D)this.patch;
      int d = patch.GetCountField();
      BivariateNURBSBasisFunction basis = (BivariateNURBSBasisFunction)((NURBSSurface)(patch.GetGeometry())).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int nen = (p + 1) * (q + 1); // number of local basis functions
      var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      int idx1 = patch.GetIPN(id, 0);
      int idx2 = patch.GetIPN(id, 1);
      fi = new DoubleVector(d * nen);
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          double detJbar = 1.0 / 4.0
              * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
              * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          //DoubleMatrix dNdxi = gps[i, j].dNdxi;
          //DoubleMatrix J = JacobianAt(dNdxi);
          //DoubleMatrix invertJ = (DoubleMatrix)J.Inverse();
          double detJ = Math.Abs((double)gps[i, j].GetValue(DataInGausspoint.detJ));
          DoubleMatrix dNdx = (DoubleMatrix)gps[i, j].GetValue(DataInGausspoint.dNdX);//dNdxi * invertJ;

          DoubleMatrix F = DeformationGradient(dNdx);
          double detF = MatrixFunctions.Determinant(F);
          DoubleMatrix B = F * (DoubleMatrix)F.Transpose();//Bbar

          double[,] dNdxs = DerivativeCoord(dNdx, F);
          DoubleMatrix P = KirchhoffStress(B, detF);

          double w = gps[i, j].weight * detJbar * detJ;
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

    public override DoubleVector StrainAt(params double[] xi)
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
    //public override void UpdateStressGaussPoint()
    //{
    //	throw new NotImplementedException();
    //}
  }
}
