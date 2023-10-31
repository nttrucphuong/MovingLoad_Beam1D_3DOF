using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using DEMSoft.Function;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public enum TypeVFunction
  { Power, MoriTanaka }
  public class ElementStructureElasticPlate : AbstractElementStructurePlate
  {
    public double GetThickness()
    { return ((PatchStructurePlate)patch).Thickness; }
    public TypePlate GetTypePlate()
    { return ((PatchStructurePlate)patch).TypePlate; }
    public FunctionRToR GetKinematicsFunction()
    { return ((PatchStructurePlate)patch).KinematicsFunction; }
    public double getzz()
    { return ((PatchStructurePlate)patch).Getzz; }
    public TypeVFunction GetTypeVFunction()
    { return ((PatchStructurePlate)patch).GetTypeVF; }
    public ElementStructureElasticPlate(AbstractPatch2D patch, int id)
            : base(patch, id)
    {
    }

    private DoubleMatrix CreateMaterialMatrixHSDT()
    {
      double thickness = GetThickness();
      DoubleMatrix D = new DoubleMatrix(11, 11);
      int numGauss = 7;
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

        D[0, 0] += detJ * wi * aa;
        D[0, 1] += detJ * wi * aa * nu;
        D[2, 2] += detJ * wi * aa * (1.0 - nu) / 2.0;

        D[0, 3] += detJ * wi * zz * aa;
        D[0, 4] += detJ * wi * zz * aa * nu;
        D[2, 5] += detJ * wi * zz * aa * (1.0 - nu) / 2.0;

        D[0, 6] += detJ * wi * fz * aa;
        D[0, 7] += detJ * wi * fz * aa * nu;
        D[2, 8] += detJ * wi * fz * aa * (1.0 - nu) / 2.0;

        D[3, 3] += detJ * wi * zz * zz * aa;
        D[3, 4] += detJ * wi * zz * zz * aa * nu;
        D[5, 5] += detJ * wi * zz * zz * aa * (1.0 - nu) / 2.0;

        D[3, 6] += detJ * wi * zz * fz * aa;
        D[3, 7] += detJ * wi * zz * fz * aa * nu;
        D[5, 8] += detJ * wi * zz * fz * aa * (1.0 - nu) / 2.0;

        D[6, 6] += detJ * wi * fz * fz * aa;
        D[6, 7] += detJ * wi * fz * fz * aa * nu;
        D[8, 8] += detJ * wi * fz * fz * aa * (1.0 - nu) / 2.0;

        D[9, 9] += detJ * wi * dfz * dfz * bb;
      }
      D[1, 1] = D[0, 0];
      D[1, 3] = D[0, 4];
      D[1, 4] = D[0, 3];
      D[1, 6] = D[0, 7];
      D[1, 7] = D[0, 6];
      D[4, 4] = D[3, 3];
      D[4, 6] = D[3, 7];
      D[4, 7] = D[3, 6];
      D[7, 7] = D[6, 6];
      D[10, 10] = D[9, 9];

      for (int ii = 0; ii < 11; ii++)
      {
        for (int jj = ii; jj < 11; jj++)
        {
          D[jj, ii] = D[ii, jj];
        }
      }
      return D;
    }

    /// <summary>
    /// Compute material matrix at gauss point
    /// </summary>
    /// <param name="i"></param>
    /// <param name="j"></param>
    /// <returns></returns>
    private DoubleMatrix CreateMaterialMatrixHSDTFGM()
    {
      double thickness = GetThickness();
      DoubleMatrix D = new DoubleMatrix(11, 11);
      int numGauss = 20;
      double eModulus = 0;
      double nu = 0;
      if (Material is FGMStructureOneVariableMaterial)
      {
        eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
        nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
      }
      for (int i = 0; i < numGauss; i++)
      {
        double xi = GaussPoints.GetPoint(numGauss, i);
        double wi = GaussPoints.GetWeight(numGauss, i);
        double zz = 0.5 * thickness * (xi + 1.0) - thickness / 2.0;
        double fz = GetKinematicsFunction().ValueAt(zz);
        double dfz = GetKinematicsFunction().DerivativeAt(zz);
        double detJ = thickness / 2.0;
        if (Material is FGMStructureOneVariableMaterial)
        {
          FunctionRToR E = ((FunctionRToR)Material.GetProperty(MaterialPropertyName.YoungModulus).GetNameClassObjectOfValueProperty());
          FunctionRToR VC = E.GetSubFunction();
          var emModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty(-thickness / 2.0);
          var ecModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty(thickness / 2.0);
          var nuym = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty(-thickness / 2.0);
          var nuyc = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty(thickness / 2.0);
          double kmModulus = emModulus / (3.0 * (1.0 - 2.0 * nuym));
          double kcModulus = ecModulus / (3.0 * (1.0 - 2.0 * nuyc));
          double gmModulus = emModulus / (2.0 * (1.0 + nuym));
          double gcModulus = ecModulus / (2.0 * (1.0 + nuyc));
          double fx = gmModulus * (9.0 * kmModulus + 8.0 * gmModulus) / (6.0 * (kmModulus + 2.0 * gmModulus));
          double kModulus = (kcModulus - kmModulus) * VC.ValueAt(zz) / (1.0 + (1.0 - VC.ValueAt(zz)) * (kcModulus - kmModulus) / (kmModulus + (4.0 / 3.0) * gmModulus)) + kmModulus;
          double gModulus = (gcModulus - gmModulus) * VC.ValueAt(zz) / (1.0 + (1.0 - VC.ValueAt(zz)) * (gcModulus - gmModulus) / (gmModulus + fx)) + gmModulus;
          eModulus = 9.0 * kModulus * gModulus / (3.0 * kModulus + gModulus);
          nu = (3.0 * kModulus - 2.0 * gModulus) / (2.0 * (3.0 * kModulus + gModulus));
        }
        if (GetTypeVFunction() == TypeVFunction.Power)
        {
          eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty(zz);
          nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty(zz);
        }
        //7.17
        double aa = eModulus / (1.0 - nu * nu);
        double bb = eModulus / (2.0 * (1.0 + nu));
        //A
        D[0, 0] += detJ * wi * aa;
        D[0, 1] += detJ * wi * aa * nu;
        D[2, 2] += detJ * wi * aa * (1.0 - nu) / 2.0;
        //B
        D[0, 3] += detJ * wi * zz * aa;
        D[0, 4] += detJ * wi * zz * aa * nu;
        D[2, 5] += detJ * wi * zz * aa * (1.0 - nu) / 2.0;
        //E
        D[0, 6] += detJ * wi * fz * aa;
        D[0, 7] += detJ * wi * fz * aa * nu;
        D[2, 8] += detJ * wi * fz * aa * (1.0 - nu) / 2.0;
        //D
        D[3, 3] += detJ * wi * zz * zz * aa;
        D[3, 4] += detJ * wi * zz * zz * aa * nu;
        D[5, 5] += detJ * wi * zz * zz * aa * (1.0 - nu) / 2.0;
        //F
        D[3, 6] += detJ * wi * zz * fz * aa;
        D[3, 7] += detJ * wi * zz * fz * aa * nu;
        D[5, 8] += detJ * wi * zz * fz * aa * (1.0 - nu) / 2.0;
        //H
        D[6, 6] += detJ * wi * fz * fz * aa;
        D[6, 7] += detJ * wi * fz * fz * aa * nu;
        D[8, 8] += detJ * wi * fz * fz * aa * (1.0 - nu) / 2.0;
        //Ds
        D[9, 9] += detJ * wi * dfz * dfz * bb;
      }
      D[1, 1] = D[0, 0];
      D[1, 3] = D[0, 4];
      D[1, 4] = D[0, 3];
      D[1, 6] = D[0, 7];
      D[1, 7] = D[0, 6];
      D[4, 4] = D[3, 3];
      D[4, 6] = D[3, 7];
      D[4, 7] = D[3, 6];
      D[7, 7] = D[6, 6];
      D[10, 10] = D[9, 9];

      for (int ii = 0; ii < 11; ii++)
      {
        for (int jj = ii; jj < 11; jj++)
        {
          D[jj, ii] = D[ii, jj];
        }
      }
      return D;
    }
    private DoubleMatrix CreateMaterialMatrixFSDT()
    {
      double thickness = GetThickness();
      DoubleMatrix D = new DoubleMatrix(8, 8);
      int numGauss = 5;
      for (int i = 0; i < numGauss; i++)
      {
        double xi = GaussPoints.GetPoint(numGauss, i);
        double wi = GaussPoints.GetWeight(numGauss, i);
        double zz = 0.5 * thickness * (xi + 1) - thickness / 2.0;
        double detJ = thickness / 2.0;
        var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
        var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();

        double aa = eModulus / (1.0 - nu * nu);
        double bb = eModulus / (2.0 * (1.0 + nu));
        //A
        D[0, 0] += detJ * wi * aa;
        D[0, 1] += detJ * wi * aa * nu;
        D[2, 2] += detJ * wi * aa * (1.0 - nu) / 2.0;
        //B
        D[0, 3] += detJ * wi * zz * aa;
        D[0, 4] += detJ * wi * zz * aa * nu;
        D[2, 5] += detJ * wi * zz * aa * (1.0 - nu) / 2.0;
        //D
        D[3, 3] += detJ * wi * zz * zz * aa;
        D[3, 4] += detJ * wi * zz * zz * aa * nu;
        D[5, 5] += detJ * wi * zz * zz * aa * (1.0 - nu) / 2.0;
        //Ds - lay f'z = 5/6
        D[6, 6] += detJ * wi * (5.0 / 6.0) * bb;
      }
      D[1, 1] = D[0, 0];
      D[1, 3] = D[0, 4];
      D[1, 4] = D[0, 3];
      D[4, 4] = D[3, 3];
      D[7, 7] = D[6, 6];

      for (int ii = 0; ii < 8; ii++)
      {
        for (int jj = ii; jj < 8; jj++)
        {
          D[jj, ii] = D[ii, jj];
        }
      }
      return D;
    }
    private DoubleMatrix CreateMaterialMatrixFSDTFGM()
    {
      double thickness = GetThickness();
      DoubleMatrix D = new DoubleMatrix(8, 8);
      int numGauss = 20;
      FunctionRToR E = ((FunctionRToR)Material.GetProperty(MaterialPropertyName.YoungModulus).GetNameClassObjectOfValueProperty());
      FunctionRToR VC = E.GetSubFunction();
      var emModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty(-thickness / 2.0);
      var ecModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty(thickness / 2.0);
      var nuym = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty(-thickness / 2.0);
      var nuyc = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty(thickness / 2.0);
      double kmModulus = emModulus / (3.0 * (1.0 - 2.0 * nuym));
      double kcModulus = ecModulus / (3.0 * (1.0 - 2.0 * nuyc));
      double gmModulus = emModulus / (2.0 * (1.0 + nuym));
      double gcModulus = ecModulus / (2.0 * (1.0 + nuyc));
      double fx = gmModulus * (9.0 * kmModulus + 8.0 * gmModulus) / (6.0 * (kmModulus + 2.0 * gmModulus));
      for (int i = 0; i < numGauss; i++)
      {
        double xi = GaussPoints.GetPoint(numGauss, i);
        double wi = GaussPoints.GetWeight(numGauss, i);
        double zz = 0.5 * thickness * (xi + 1) - thickness / 2.0;
        double detJ = thickness / 2.0;
        double kModulus = (kcModulus - kmModulus) * VC.ValueAt(zz) / (1.0 + (1.0 - VC.ValueAt(zz)) * (kcModulus - kmModulus) / (kmModulus + (4.0 / 3.0) * gmModulus)) + kmModulus;
        double gModulus = (gcModulus - gmModulus) * VC.ValueAt(zz) / (1.0 + (1.0 - VC.ValueAt(zz)) * (gcModulus - gmModulus) / (gmModulus + fx)) + gmModulus;
        double eModulus = 9.0 * kModulus * gModulus / (3.0 * kModulus + gModulus);
        double nu = (3.0 * kModulus - 2.0 * gModulus) / (2.0 * (3.0 * kModulus + gModulus));
        if (GetTypeVFunction() == TypeVFunction.Power)
        {
          eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty(zz);
          nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty(zz);
        }
        double aa = eModulus / (1.0 - nu * nu);
        double bb = eModulus / (2.0 * (1.0 + nu));

        D[0, 0] += detJ * wi * aa;
        D[0, 1] += detJ * wi * aa * nu;
        D[2, 2] += detJ * wi * aa * (1.0 - nu) / 2.0;

        D[0, 3] += detJ * wi * zz * aa;
        D[0, 4] += detJ * wi * zz * aa * nu;
        D[2, 5] += detJ * wi * zz * aa * (1.0 - nu) / 2.0;

        D[3, 3] += detJ * wi * zz * zz * aa;
        D[3, 4] += detJ * wi * zz * zz * aa * nu;
        D[5, 5] += detJ * wi * zz * zz * aa * (1.0 - nu) / 2.0;

        D[6, 6] += detJ * wi * (5.0 / 6.0) * bb;
      }
      D[1, 1] = D[0, 0];
      D[1, 3] = D[0, 4];
      D[1, 4] = D[0, 3];
      D[4, 4] = D[3, 3];
      D[7, 7] = D[6, 6];

      for (int ii = 0; ii < 8; ii++)
      {
        for (int jj = ii; jj < 8; jj++)
        {
          D[jj, ii] = D[ii, jj];
        }
      }
      return D;
    }


    public override void ComputeStiffnessMatrixElement(out DoubleMatrix Ke)
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
      Ke = new DoubleMatrix(d * nen, d * nen);
      DoubleMatrix B;
      if (GetTypePlate() == TypePlate.FSDT)
      {
        B = new DoubleMatrix(8, d * nen);
      }
      else
      {
        B = new DoubleMatrix(11, d * nen);
      }
      DoubleMatrix D = null;

      if (Material is FGMStructureOneVariableMaterial)
      {
        if (GetTypePlate() == TypePlate.FSDT)
        {
          D = CreateMaterialMatrixFSDTFGM();
        }
        else
        {
          D = CreateMaterialMatrixHSDTFGM();
        }
      }
      else
      {
        if (GetTypePlate() == TypePlate.FSDT)
        {
          D = CreateMaterialMatrixFSDT();
        }
        else
        {
          D = CreateMaterialMatrixHSDT();
        }
      }

      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          double detJbar = 1.0 / 4.0
            * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
            * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          DoubleMatrix dNdxi = (DoubleMatrix)gps[i, j].GetValue(DataInGausspoint.dNdxi);//gps[i, j].dNdxi;
          DoubleMatrix J = JacobianAt(dNdxi);
          DoubleMatrix dNdx = (DoubleMatrix)gps[i, j].GetValue(DataInGausspoint.dNdX);
          DoubleMatrix ddNdX = (DoubleMatrix)gps[i, j].GetValue(DataInGausspoint.ddNdX);//gps[i, j].ddNdX;
          DoubleVector Ni = (DoubleVector)gps[i, j].GetValue(DataInGausspoint.Ni);
          for (int k = 0; k < nen; k++)
          {
            if (GetTypePlate() == TypePlate.FSDT)
            {
                //BAm
              B[0, 5 * k] = dNdx[k, 0];
              B[1, 5 * k + 1] = dNdx[k, 1];
              B[2, 5 * k] = dNdx[k, 1];
              B[2, 5 * k + 1] = dNdx[k, 0];
                //BAb2
              B[3, 5 * k + 3] = dNdx[k, 0];
              B[4, 5 * k + 4] = dNdx[k, 1];
              B[5, 5 * k + 3] = dNdx[k, 1];
              B[5, 5 * k + 4] = dNdx[k, 0];
                //BAs+BAg
              B[6, 5 * k + 2] = dNdx[k, 0];
              B[6, 5 * k + 3] = Ni[k];
              B[7, 5 * k + 2] = dNdx[k, 1];
              B[7, 5 * k + 4] = Ni[k];
            }
            else
            {
                            //BAm
              B[0, 5 * k] = dNdx[k, 0];
              B[1, 5 * k + 1] = dNdx[k, 1];
              B[2, 5 * k] = dNdx[k, 1];
              B[2, 5 * k + 1] = dNdx[k, 0];
                            //BAb1
              B[3, 5 * k + 2] = -ddNdX[k, 0];
              B[4, 5 * k + 2] = -ddNdX[k, 1];
              B[5, 5 * k + 2] = -2.0 * ddNdX[k, 2];
                            //BAb2
              B[6, 5 * k + 3] = dNdx[k, 0];
              B[7, 5 * k + 4] = dNdx[k, 1];
              B[8, 5 * k + 3] = dNdx[k, 1];
              B[8, 5 * k + 4] = dNdx[k, 0];
                            //BAs
              B[9, 5 * k + 3] = Ni[k];
              B[10, 5 * k + 4] = Ni[k];
            }
          }

          DoubleMatrix BTDB = NMathFunctions.TransposeProduct(B, MatrixFunctions.Product(D, B));
          Ke += gps[i, j].weight * detJbar * Math.Abs((double)gps[i, j].GetValue(DataInGausspoint.detJ)) * BTDB;
        }
      }
    }

    public override DoubleVector StrainAt(params double[] xi)
    {
      PatchStructurePlate patch = (PatchStructurePlate)this.patch;
      int d = patch.GetCountField();
      int nen = patch.GetCountLocalBasisFunctions();
      var cps = ((NURBSSurface)(patch.GetGeometry())).ControlPoints;
      DoubleMatrix B;
      if (GetTypePlate() == TypePlate.FSDT)
      {
        B = new DoubleMatrix(8, d * nen);
      }
      else
      {
        B = new DoubleMatrix(11, d * nen);
      }

      DoubleMatrix gradNE = GradBasisFunction(xi[0], xi[1]);
      DoubleMatrix J = JacobianAt(gradNE);
      DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
      DoubleMatrix gradBasis = gradNE * invertJ;
      DoubleVector ValueBasis = ValueBasisFunction(xi[0], xi[1]);
      DoubleMatrix grad2Basis = Grad2BasisFunction(xi[0], xi[1]);
      for (int k = 0; k < nen; k++)
      {
        if (GetTypePlate() == TypePlate.FSDT)
        {
          B[0, 5 * k] = gradBasis[k, 0];
          B[1, 5 * k + 1] = gradBasis[k, 1];
          B[2, 5 * k] = gradBasis[k, 1];
          B[2, 5 * k + 1] = gradBasis[k, 0];

          B[3, 5 * k + 3] = gradBasis[k, 0];
          B[4, 5 * k + 4] = gradBasis[k, 1];
          B[5, 5 * k + 3] = gradBasis[k, 1];
          B[5, 5 * k + 4] = gradBasis[k, 0];

          B[6, 5 * k + 2] = gradBasis[k, 0];
          B[6, 5 * k + 3] = ValueBasis[k];
          B[7, 5 * k + 2] = gradBasis[k, 1];
          B[7, 5 * k + 4] = ValueBasis[k];
        }
        else
        {//7.30
          B[0, 5 * k] = gradBasis[k, 0];
          B[1, 5 * k + 1] = gradBasis[k, 1];
          B[2, 5 * k] = gradBasis[k, 1];
          B[2, 5 * k + 1] = gradBasis[k, 0];

          B[3, 5 * k + 2] = -grad2Basis[k, 0];
          B[4, 5 * k + 2] = -grad2Basis[k, 1];
          B[5, 5 * k + 2] = -2 * grad2Basis[k, 2];

          B[6, 5 * k + 3] = gradBasis[k, 0];
          B[7, 5 * k + 4] = gradBasis[k, 1];
          B[8, 5 * k + 3] = gradBasis[k, 1];
          B[8, 5 * k + 4] = gradBasis[k, 0];

          B[9, 5 * k + 3] = ValueBasis[k];
          B[10, 5 * k + 4] = ValueBasis[k];
        }
      }
      DoubleVector U = GetDisplacementLocal();
      DoubleVector Smatrix = MatrixFunctions.Product(B, U);
      DoubleVector Rmatrix = new DoubleVector(5);
      double zz = getzz();
      if (GetTypePlate() == TypePlate.FSDT)
      {
        Rmatrix[0] = Smatrix[0] + zz * Smatrix[3];
        Rmatrix[1] = Smatrix[1] + zz * Smatrix[4];
        Rmatrix[2] = Smatrix[2] + zz * Smatrix[5];
        Rmatrix[3] = Smatrix[6];
        Rmatrix[4] = Smatrix[7];
      }
      else
      {
        FunctionRToR fz = GetKinematicsFunction();
        Rmatrix[0] = Smatrix[0] + zz * Smatrix[3] + fz.ValueAt(zz) * Smatrix[6];
        Rmatrix[1] = Smatrix[1] + zz * Smatrix[4] + fz.ValueAt(zz) * Smatrix[7];
        Rmatrix[2] = Smatrix[2] + zz * Smatrix[5] + fz.ValueAt(zz) * Smatrix[8];
        Rmatrix[3] = Smatrix[9];
        Rmatrix[4] = Smatrix[10];
      }
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
      DoubleMatrix Q = new DoubleMatrix(5, 5);
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
      Q[2, 2] = Q[3, 3] = Q[4, 4] = bb;
      return Q;
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
    //private double Getfz(double zz)
    //{
    //    TypePlate type = GetTypePlate();
    //    double fz = 0;
    //    double h = GetThickness();
    //    switch (type)
    //    {
    //        case TypePlate.HSDT1:
    //            fz = 5.0 * zz / 4.0 - 5.0 * Math.Pow(zz, 3.0) / (3.0 * h * h);
    //            break;
    //        case TypePlate.HSDT2:
    //            fz = zz - 4.0 * Math.Pow(zz, 3.0) / (3.0 * h * h);
    //            break;
    //        case TypePlate.HSDT3:
    //            fz = zz * Math.Exp(-2.0 * zz * zz / (h * h));
    //            break;
    //        case TypePlate.HSDT4:
    //            fz = 7.0 * zz / 8.0 - 2.0 * Math.Pow(zz, 3.0) / (h * h) + 2.0 * Math.Pow(zz, 5.0) / Math.Pow(h, 4.0);
    //            break;
    //        case TypePlate.HSDT5:
    //            fz = h * Math.Atan(2.0 * zz / h) - zz;
    //            break;
    //        case TypePlate.HSDT6:
    //            fz = Math.Atan(Math.Sin(Math.PI * zz / h));
    //            break;
    //        case TypePlate.HSDT7:
    //            fz = zz;
    //            break;
    //    }

    //    return fz;
    //}
    //private double Getdfz(double zz)
    //{
    //    TypePlate type = GetTypePlate();
    //    double dfz = 0;
    //    double h = GetThickness();
    //    switch (type)
    //    {
    //        case TypePlate.HSDT1:
    //            dfz = 5.0 / 4.0 - 5.0 * zz * zz / (h * h);
    //            break;
    //        case TypePlate.HSDT2:
    //            dfz = 1.0 - 4.0 * zz * zz / (h * h);
    //            break;
    //        case TypePlate.HSDT3:
    //            dfz = (1.0 - 4.0 * zz * zz / (h * h)) * Math.Exp(-2.0 * zz * zz / (h * h));
    //            break;
    //        case TypePlate.HSDT4:
    //            dfz = 7.0 / 8.0 - 6.0 * zz * zz / (h * h) + 10.0 * Math.Pow(zz, 4.0) / Math.Pow(h, 4.0);
    //            break;
    //        case TypePlate.HSDT5:
    //            dfz = (h * h - 4.0 * zz * zz) / (h * h + 4.0 * zz * zz);
    //            break;
    //        case TypePlate.HSDT6:
    //            dfz = (Math.PI / h) * Math.Cos(Math.PI * zz / h) / (1.0 + Math.Sin(Math.PI * zz / h) * Math.Sin(Math.PI * zz / h));
    //            break;
    //        case TypePlate.HSDT7:
    //            dfz = Math.Sqrt(5.0 / 6.0);
    //            break;
    //    }

    //    return dfz;
    //}
  }
}
