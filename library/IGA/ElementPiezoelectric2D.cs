using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public class ElementPiezoelectric2D : ElementStructureElastic2D
  {
    public ElementPiezoelectric2D(AbstractPatch2D patch, int id)
          : base(patch, id)
    {
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
      DoubleMatrix Keuu = new DoubleMatrix(2 * nen, 2 * nen);
      DoubleMatrix Kett = new DoubleMatrix(nen, nen);
      DoubleMatrix Keut = new DoubleMatrix(2 * nen, nen);
      DoubleMatrix Bu = new DoubleMatrix(4, 2 * nen);
      DoubleMatrix Bt = new DoubleMatrix(3, nen);
      DoubleMatrix D = CreateMaterialMatrix();
      DoubleMatrix e = CreatePiezoElectricStressMatrix();
      DoubleMatrix epsS = CreateDielectricMatrix();
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          int c = 0;
          double x1 = 1;
          double a = 1;
          GaussPoints gpsij = gps[i, j];
          if (patch.StateStress == Structure2DState.Axisymetric)
          {
            c = 1;
            x1 = PointAt(gpsij.location[0], gpsij.location[1])[0];
            a = 2 * Math.PI * x1;
          }
          else if (patch.StateStress == Structure2DState.PlaneStress)
          {
            a = patch.Thickness;
            if (a == 0)
              throw new InvalidArgumentException("Plane stress condition must had thickness");
          }
          double detJbar = 1.0 / 4.0 * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1]) * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          DoubleMatrix dNdx = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.dNdX);
          DoubleVector N = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);
          double detJ = Math.Abs((double)gpsij.GetValue(DataInGausspoint.detJ));
          for (int k = 0; k < nen; k++)
          {
            //xx yy zz xy
            Bu[0, 2 * k] = dNdx[k, 0];
            Bu[2, 2 * k] = c * N[k] / x1;
            Bu[3, 2 * k] = dNdx[k, 1];
            Bu[1, 2 * k + 1] = dNdx[k, 1];
            Bu[3, 2 * k + 1] = dNdx[k, 0];

            //x y z (z is rotate)
            Bt[0, k] = dNdx[k, 0];
            Bt[1, k] = dNdx[k, 1];
            Bt[2, k] = c * N[k] / x1;
          }
          DoubleMatrix BuTDBu = NMathFunctions.TransposeProduct(Bu, MatrixFunctions.Product(D, Bu));
          DoubleMatrix BuTeBt = MatrixFunctions.TransposeProduct(Bu, MatrixFunctions.Product(e, Bt));
          DoubleMatrix BtTepsSBt = NMathFunctions.TransposeProduct(Bt, MatrixFunctions.Product(epsS, Bt));
          Keuu += a * gpsij.weight * detJbar * detJ * BuTDBu;
          Keut += a * gpsij.weight * detJbar * detJ * BuTeBt;
          Kett -= a * gpsij.weight * detJbar * detJ * BtTepsSBt;
        }
      }
      for (int j = 0; j < nen; j++)
        for (int i = 0; i < nen; i++)
        {
          for (int jj = 0; jj < 2; jj++)
          {
            for (int ii = 0; ii < 2; ii++)
            {
              Ke[i * 3 + ii, j * 3 + jj] = Keuu[i * 2 + ii, j * 2 + jj];
            }
            Ke[i * 3 + 2, j * 3 + jj] = Ke[j * 3 + jj, i * 3 + 2] = Keut[j * 2 + jj, i];
          }
          Ke[i * 3 + 2, j * 3 + 2] = Kett[i, j];
        }
    }

    private DoubleMatrix CreateDielectricMatrix()
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
          //DoubleMatrix C = CreateMaterialMatrix();
          DoubleMatrix e = CreatePiezoElectricStressMatrix();
          DoubleMatrix d = CreatePiezoElectricStrainMatrix();
          epsS = epsS - MatrixFunctions.Product(e.Transpose(), d);
        }
      }
      return epsS;
    }

    /// <summary>
    /// e
    /// </summary>
    /// <returns></returns>
    private DoubleMatrix CreatePiezoElectricStressMatrix()
    {
      PiezoelectricMatrix prob = (PiezoelectricMatrix)Material.GetProperty(MaterialPropertyName.PiezoelectricMatrix);
      DoubleMatrix e = new DoubleMatrix(4, 3);
      if (prob != null)
      {
        double[,] eTemp = prob.GetMatrix();
        bool isStrain = prob.GetIsStrain();
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 3; j++)
          {
            e[i, j] = eTemp[i, j];
          }
        if (isStrain)
        {
          DoubleMatrix C = CreateMaterialMatrix();
          e = MatrixFunctions.Product(C, e);
        }
      }
      return e;
    }

    /// <summary>
    /// d
    /// </summary>
    /// <returns></returns>
    private DoubleMatrix CreatePiezoElectricStrainMatrix()
    {
      PiezoelectricMatrix prob = (PiezoelectricMatrix)Material.GetProperty(MaterialPropertyName.PiezoelectricMatrix);
      DoubleMatrix d = new DoubleMatrix(4, 3);
      if (prob != null)
      {
        double[,] dTemp = prob.GetMatrix();
        bool isStrain = prob.GetIsStrain();
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 3; j++)
          {
            d[i, j] = dTemp[i, j];
          }
        if (!isStrain)
        {
          DoubleMatrix C = CreateMaterialMatrix();
          d = MatrixFunctions.Product(MatrixFunctions.Inverse(C), d);
        }
      }
      return d;
    }

    public override DoubleVector StressAt(params double[] xi)
    {
      DoubleMatrix D = CreateMaterialMatrix();
      DoubleMatrix e = CreatePiezoElectricStressMatrix();
      return MatrixFunctions.Product(D, StrainAt(xi)) - MatrixFunctions.Product(e, ElectricFieldPotential(xi));
    }

    private DoubleVector ElectricFieldPotential(params double[] xi)
    {
      PatchStructure2D patch = (PatchStructure2D)this.patch;
      int d = 1;// patch.GetCountField();
      int nen = patch.GetCountLocalBasisFunctions();
      var cps = ((NURBSSurface)(patch.GetGeometry())).ControlPoints;
      DoubleMatrix B = new DoubleMatrix(3, d * nen);

      DoubleMatrix gradNE = GradBasisFunction(xi[0], xi[1]);
      DoubleMatrix J = JacobianAt(gradNE);
      DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
      DoubleMatrix gradBasis = MatrixFunctions.Product(gradNE, invertJ);
      DoubleVector N = ValueBasisFunction(xi[0], xi[1]);
      int c = 0;
      double x1 = 1;
      if (patch.StateStress == Structure2DState.Axisymetric)
      {
        c = 1;
        x1 = PointAt(xi[0], xi[1])[0];
      }
      for (int k = 0; k < nen; k++)
      {
        B[0, k] = gradBasis[k, 0];
        B[1, k] = gradBasis[k, 1];
        B[2, k] = c * N[k] / x1;
      }
      DoubleVector U = GetElectricLocal();
      return -MatrixFunctions.Product(B, U);
    }
    private DoubleVector GetElectricLocal()
    {
      AbstractPatch2D patch = (AbstractPatch2D)this.patch;
      var cps = ((NURBSSurface)(patch.GetGeometry())).ControlPoints;
      int d = 1;// patch.GetCountField();
      int nen = patch.GetCountLocalBasisFunctions();
      DoubleVector U = new DoubleVector(d * nen);
      for (int i = 0; i < nen; i++)
      {
        var ien = patch.GetIEN(id, i);
        U[d * i] = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1)].GetResult(Result.VOLT);
      }
      return U;
    }
  }
}
