using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public class ElementPiezoelectric3D : ElementStructureElastic3D
  {
    public ElementPiezoelectric3D(AbstractPatch3D patch, int id)
          : base(patch, id)
    {
    }

    public override void ComputeStiffnessMatrixElement(out DoubleMatrix Ke)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      int d = patch.GetCountField();
      int dimension = patch.GetCountDimension();
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
      DoubleMatrix Keuu = new DoubleMatrix(dimension * nen, dimension * nen);
      DoubleMatrix Kett = new DoubleMatrix(nen, nen);
      DoubleMatrix Keut = new DoubleMatrix(dimension * nen, nen);
      DoubleMatrix Bu = new DoubleMatrix(6, dimension * nen);
      DoubleMatrix Bt = new DoubleMatrix(3, nen);
      DoubleMatrix D = null;
      DoubleMatrix e = null;
      DoubleMatrix epsS = null;
      if (!(Material is FGMStructureOneVariableMaterial))
      {
        D = CreateMaterialMatrix();
        e = CreatePiezoElectricStressMatrix(D);
        epsS = CreateDielectricMatrix(D);
      }
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
            DoubleMatrix dNdx = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdX);
            double detJ = Math.Abs((double)gps[i, j, k].GetValue(DataInGausspoint.detJ));
            for (int kk = 0; kk < nen; kk++)
            {
              ////xx yy zz yz xz xy
              //Bu[0, 3 * kk] = dNdx[kk, 0];
              //Bu[1, 3 * kk + 1] = dNdx[kk, 1];
              //Bu[2, 3 * kk + 2] = dNdx[kk, 2];
              //Bu[3, 3 * kk + 1] = dNdx[kk, 2];
              //Bu[3, 3 * kk + 2] = dNdx[kk, 1];
              //Bu[4, 3 * kk] = dNdx[kk, 2];
              //Bu[4, 3 * kk + 2] = dNdx[kk, 0];
              //Bu[5, 3 * kk] = dNdx[kk, 1];
              //Bu[5, 3 * kk + 1] = dNdx[kk, 0];

              //xx yy zz xy yz xz
              Bu[0, 3 * kk] = dNdx[kk, 0];
              Bu[1, 3 * kk + 1] = dNdx[kk, 1];
              Bu[2, 3 * kk + 2] = dNdx[kk, 2];
              Bu[3, 3 * kk] = dNdx[kk, 1];
              Bu[3, 3 * kk + 1] = dNdx[kk, 0];
              Bu[4, 3 * kk + 1] = dNdx[kk, 2];
              Bu[4, 3 * kk + 2] = dNdx[kk, 1];
              Bu[5, 3 * kk] = dNdx[kk, 2];
              Bu[5, 3 * kk + 2] = dNdx[kk, 0];

              Bt[0, kk] = dNdx[kk, 0];
              Bt[1, kk] = dNdx[kk, 1];
              Bt[2, kk] = dNdx[kk, 2];
            }
            if (Material is FGMStructureOneVariableMaterial)
            {
              D = CreateMaterialMatrix(i, j, k);
              e = CreatePiezoElectricStressMatrix(D);
              epsS = CreateDielectricMatrix(D);
            }
            DoubleMatrix BuTDBu = NMathFunctions.TransposeProduct(Bu, MatrixFunctions.Product(D, Bu));
            DoubleMatrix BuTeBt = MatrixFunctions.TransposeProduct(Bu, MatrixFunctions.Product(e, Bt));
            DoubleMatrix BtTepsSBt = NMathFunctions.TransposeProduct(Bt, MatrixFunctions.Product(epsS, Bt));
            Keuu += gps[i, j, k].weight * detJbar * detJ * BuTDBu;
            Keut += gps[i, j, k].weight * detJbar * detJ * BuTeBt;
            Kett -= gps[i, j, k].weight * detJbar * detJ * BtTepsSBt;
          }
        }
      }
      for (int j = 0; j < nen; j++)
        for (int i = 0; i < nen; i++)
        {
          for (int jj = 0; jj < 3; jj++)
          {
            for (int ii = 0; ii < 3; ii++)
            {
              Ke[i * 4 + ii, j * 4 + jj] = Keuu[i * 3 + ii, j * 3 + jj];
            }
            Ke[i * 4 + 3, j * 4 + jj] = Ke[j * 4 + jj, i * 4 + 3] = Keut[j * 3 + jj, i];
          }
          Ke[i * 4 + 3, j * 4 + 3] = Kett[i, j];
        }
    }

    private DoubleMatrix CreateDielectricMatrix(DoubleMatrix C)
    {
      AnisotropicDielectricMatrix prob = (AnisotropicDielectricMatrix)Material.GetProperty(MaterialPropertyName.AnisotropicDielectricMatrix);
      DoubleMatrix epsS = new DoubleMatrix(3, 3);
      if (prob != null)
      {
        double[,] epsSTemp = prob.GetMatrix();
        bool isStrain = prob.GetIsStrain();
        epsS[0, 0] = epsSTemp[0, 0];
        epsS[1, 1] = epsSTemp[1, 1];
        epsS[2, 2] = epsSTemp[2, 2];
        if (!isStrain)
        {
          DoubleMatrix e = CreatePiezoElectricStressMatrix(C);
          DoubleMatrix d = CreatePiezoElectricStrainMatrix(C);
          epsS = epsS - NMathFunctions.TransposeProduct(e, d);
        }
      }
      return epsS;
    }

    /// <summary>
    /// e
    /// </summary>
    /// <returns></returns>
    private DoubleMatrix CreatePiezoElectricStressMatrix(DoubleMatrix C)
    {
      PiezoelectricMatrix prob = (PiezoelectricMatrix)Material.GetProperty(MaterialPropertyName.PiezoelectricMatrix);
      DoubleMatrix e = new DoubleMatrix(6, 3);
      if (prob != null)
      {
        double[,] eTemp = prob.GetMatrix();
        bool isStrain = prob.GetIsStrain();
        for (int i = 0; i < 6; i++)
          for (int j = 0; j < 3; j++)
          {
            e[i, j] = eTemp[i, j];
          }
        if (isStrain)
        {
          //DoubleMatrix C = CreateMaterialMatrix();
          e = NMathFunctions.Product(C, e);
        }
      }
      return e;
    }

    /// <summary>
    /// d
    /// </summary>
    /// <returns></returns>
    private DoubleMatrix CreatePiezoElectricStrainMatrix(DoubleMatrix C)
    {
      PiezoelectricMatrix prob = (PiezoelectricMatrix)Material.GetProperty(MaterialPropertyName.PiezoelectricMatrix);
      DoubleMatrix d = new DoubleMatrix(6, 3);
      if (prob != null)
      {
        double[,] dTemp = prob.GetMatrix();
        bool isStrain = prob.GetIsStrain();
        for (int i = 0; i < 6; i++)
          for (int j = 0; j < 3; j++)
          {
            d[i, j] = dTemp[i, j];
          }
        if (!isStrain)
        {
          //DoubleMatrix C = CreateMaterialMatrix();
          d = NMathFunctions.Product(MatrixFunctions.Inverse(C), d);
        }
      }
      return d;
    }

    public override DoubleVector StressAt(params double[] xi)
    {
      DoubleMatrix D = null;
      if (!(Material is FGMStructureOneVariableMaterial))
        D = CreateMaterialMatrix();
      else
        D = CreateMaterialMatrix(xi[0], xi[1], xi[2]);
      DoubleMatrix e = CreatePiezoElectricStressMatrix(D);
      return MatrixFunctions.Product(D, StrainAt(xi)) - MatrixFunctions.Product(e, ElectricFieldPotential(xi));
    }

    private DoubleVector ElectricFieldPotential(params double[] xi)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      int d = 1;
      int nen = patch.GetCountLocalBasisFunctions();
      var cps = ((NURBSVolume)(patch.GetGeometry())).ControlPoints;
      DoubleMatrix B = new DoubleMatrix(3, d * nen);

      DoubleMatrix gradNE = GradBasisFunction(xi);
      DoubleMatrix J = JacobianAt(gradNE);
      DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
      DoubleMatrix gradBasis = MatrixFunctions.Product(gradNE, invertJ);
      for (int k = 0; k < nen; k++)
      {
        B[0, k] = gradBasis[k, 0];
        B[1, k] = gradBasis[k, 1];
        B[2, k] = gradBasis[k, 2];
      }
      DoubleVector U = GetElectricLocal();
      return -MatrixFunctions.Product(B, U);
    }
    private DoubleVector GetElectricLocal()
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      var cps = ((NURBSVolume)(patch.GetGeometry())).ControlPoints;
      int d = 1;
      int nen = patch.GetCountLocalBasisFunctions();
      DoubleVector U = new DoubleVector(d * nen);
      for (int i = 0; i < nen; i++)
      {
        var ien = patch.GetIEN(id, i);
        U[i] = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(ien, 2)].GetResult(Result.VOLT);
      }
      return U;
    }
  }
}
