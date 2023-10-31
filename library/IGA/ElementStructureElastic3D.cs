using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public class ElementStructureElastic3D : AbstractElementStructure3D
  {
    public ElementStructureElastic3D(AbstractPatch3D mesh, int id)
         : base(mesh, id)
    {
    }
    protected DoubleMatrix CreateMaterialMatrixLamda()
    {
      DoubleMatrix D = new DoubleMatrix(6, 6);
      D[0, 0] = coefficientLamda(1, 1, 1, 1);
      D[0, 1] = D[1, 0] = coefficientLamda(1, 1, 2, 2);
      D[0, 2] = D[2, 0] = coefficientLamda(1, 1, 3, 3);
      D[0, 3] = D[3, 0] = coefficientLamda(1, 1, 2, 3);
      D[0, 4] = D[4, 0] = coefficientLamda(1, 1, 1, 3);
      D[0, 5] = D[5, 0] = coefficientLamda(1, 1, 1, 2);
      D[1, 1] = coefficientLamda(2, 2, 2, 2);
      D[1, 2] = D[2, 1] = coefficientLamda(2, 2, 3, 3);
      D[1, 3] = D[3, 1] = coefficientLamda(2, 2, 2, 3);
      D[1, 4] = D[4, 1] = coefficientLamda(2, 2, 1, 3);
      D[1, 5] = D[5, 1] = coefficientLamda(2, 2, 1, 2);
      D[2, 2] = coefficientLamda(3, 3, 3, 3);
      D[2, 3] = D[3, 2] = coefficientLamda(3, 3, 2, 3);
      D[2, 4] = D[4, 2] = coefficientLamda(3, 3, 1, 3);
      D[2, 5] = D[5, 2] = coefficientLamda(3, 3, 1, 2);
      D[3, 3] = coefficientLamda(2, 3, 3, 2);
      D[3, 4] = D[4, 3] = coefficientLamda(2, 3, 1, 3);
      D[3, 5] = D[5, 3] = coefficientLamda(2, 3, 1, 2);
      D[4, 4] = coefficientLamda(1, 3, 3, 1);
      D[4, 5] = D[5, 4] = coefficientLamda(1, 3, 1, 2);
      D[5, 5] = coefficientLamda(1, 2, 2, 1);
      return D;
    }

    protected DoubleMatrix CreateMaterialMatrixMu()
    {
      DoubleMatrix D = new DoubleMatrix(6, 6);
      D[0, 0] = coefficientMu(1, 1, 1, 1);
      D[0, 1] = D[1, 0] = coefficientMu(1, 1, 2, 2);
      D[0, 2] = D[2, 0] = coefficientMu(1, 1, 3, 3);
      D[0, 3] = D[3, 0] = coefficientMu(1, 1, 2, 3);
      D[0, 4] = D[4, 0] = coefficientMu(1, 1, 1, 3);
      D[0, 5] = D[5, 0] = coefficientMu(1, 1, 1, 2);
      D[1, 1] = coefficientMu(2, 2, 2, 2);
      D[1, 2] = D[2, 1] = coefficientMu(2, 2, 3, 3);
      D[1, 3] = D[3, 1] = coefficientMu(2, 2, 2, 3);
      D[1, 4] = D[4, 1] = coefficientMu(2, 2, 1, 3);
      D[1, 5] = D[5, 1] = coefficientMu(2, 2, 1, 2);
      D[2, 2] = coefficientMu(3, 3, 3, 3);
      D[2, 3] = D[3, 2] = coefficientMu(3, 3, 2, 3);
      D[2, 4] = D[4, 2] = coefficientMu(3, 3, 1, 3);
      D[2, 5] = D[5, 2] = coefficientMu(3, 3, 1, 2);
      D[3, 3] = coefficientMu(2, 3, 3, 2);
      D[3, 4] = D[4, 3] = coefficientMu(2, 3, 1, 3);
      D[3, 5] = D[5, 3] = coefficientMu(2, 3, 1, 2);
      D[4, 4] = coefficientMu(1, 3, 3, 1);
      D[4, 5] = D[5, 4] = coefficientMu(1, 3, 1, 2);
      D[5, 5] = coefficientMu(1, 2, 2, 1);
      return D;
    }
    protected DoubleMatrix CreateMaterialMatrix()
    {
      AnisotropicElasticity prob = (AnisotropicElasticity)Material.GetProperty(MaterialPropertyName.AnisotropicElasticity);
      if (prob != null)
      {
        double[,] ce = prob.GetAnisotropicElasticityMatrix();
        DoubleMatrix C = new DoubleMatrix(6, 6);
        bool isInverse = prob.GetIsInverseMatrix();
        for (int i = 0; i < 6; i++)
          for (int j = 0; j < 6; j++)
          {
            C[i, j] = ce[i, j];
          }
        if (!isInverse)
          return C;
        else
          return MatrixFunctions.Inverse(C);
      }
      else
      {
        DoubleMatrix C = null;
        OrthotropicElasticity prob1 = (OrthotropicElasticity)Material.GetProperty(MaterialPropertyName.OrthotropicElasticity);
        if (prob1 == null)
        {
          var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
          var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
          C = new DoubleMatrix(6, 6);
          //Xem lai thu tu xx yy zz xy yz xz
          //C[0, 0] = coefficient(1, 1, 1, 1);
          //C[0, 1] = C[1, 0] = coefficient(1, 1, 2, 2);
          //C[0, 2] = C[2, 0] = coefficient(1, 1, 3, 3);
          //C[0, 3] = C[3, 0] = coefficient(1, 1, 2, 3);
          //C[0, 4] = C[4, 0] = coefficient(1, 1, 1, 3);
          //C[0, 5] = C[5, 0] = coefficient(1, 1, 1, 2);
          //C[1, 1] = coefficient(2, 2, 2, 2);
          //C[1, 2] = C[2, 1] = coefficient(2, 2, 3, 3);
          //C[1, 3] = C[3, 1] = coefficient(2, 2, 2, 3);
          //C[1, 4] = C[4, 1] = coefficient(2, 2, 1, 3);
          //C[1, 5] = C[5, 1] = coefficient(2, 2, 1, 2);
          //C[2, 2] = coefficient(3, 3, 3, 3);
          //C[2, 3] = C[3, 2] = coefficient(3, 3, 2, 3);
          //C[2, 4] = C[4, 2] = coefficient(3, 3, 1, 3);
          //C[2, 5] = C[5, 2] = coefficient(3, 3, 1, 2);
          //C[3, 3] = coefficient(2, 3, 3, 2);
          //C[3, 4] = C[4, 3] = coefficient(2, 3, 1, 3);
          //C[3, 5] = C[5, 3] = coefficient(2, 3, 1, 2);
          //C[4, 4] = coefficient(1, 3, 3, 1);
          //C[4, 5] = C[5, 4] = coefficient(1, 3, 1, 2);
          //C[5, 5] = coefficient(1, 2, 2, 1);
          double a = eModulus / ((1 + nu) * (1 - 2 * nu));
          C[0, 0] = a * (1 - nu);
          C[0, 1] = C[1, 0] = a * nu;
          C[0, 2] = C[2, 0] = a * nu;
          C[1, 2] = C[2, 1] = a * nu;
          C[1, 1] = a * (1 - nu);
          C[2, 2] = a * (1 - nu);
          C[3, 3] = C[4, 4] = C[5, 5] = a * (1 - 2 * nu) / 2.0;
        }
        else
        {
          C = new DoubleMatrix(prob1.GetOrthotropicElasticityMatrix());
        }
        return C;
      }
    }

    private double coefficient(int i, int j, int k, int l)
    {
      var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
      var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
      double G = eModulus / (2 * (1 + nu));
      double lamda = eModulus * nu / ((1 + nu) * (1 - 2 * nu));
      return lamda * KroneckerDelta(i, j) * KroneckerDelta(k, l) + G * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
    }

    private double coefficientLamda(int i, int j, int k, int l)
    {
      var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
      var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
      double lamda = eModulus * nu / ((1 + nu) * (1 - 2 * nu));
      return lamda * KroneckerDelta(i, j) * KroneckerDelta(k, l);
    }
    private double coefficientMu(int i, int j, int k, int l)
    {
      var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
      var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
      double G = eModulus / (2 * (1 + nu));
      return G * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
    }

    protected double KroneckerDelta(int i, int j)
    {
      return (i == j) ? 1.0 : 0.0;
    }

    /// <summary>
    /// Material matrix at gauss point
    /// </summary>
    /// <param name="i">Index of gp at xi direction</param>
    /// <param name="j">Index of gp at eta direction</param>
    /// <param name="k">Index of gp at zeta direction</param>
    /// <returns></returns>
    protected DoubleMatrix CreateMaterialMatrix(int i, int j, int k)
    {
      // xx yy zz yz xz xy
      DoubleMatrix D = new DoubleMatrix(6, 6);
      D[0, 0] = coefficient(1, 1, 1, 1, i, j, k);
      D[0, 1] = D[1, 0] = coefficient(1, 1, 2, 2, i, j, k);
      D[0, 2] = D[2, 0] = coefficient(1, 1, 3, 3, i, j, k);
      D[0, 3] = D[3, 0] = coefficient(1, 1, 2, 3, i, j, k);
      D[0, 4] = D[4, 0] = coefficient(1, 1, 1, 3, i, j, k);
      D[0, 5] = D[5, 0] = coefficient(1, 1, 1, 2, i, j, k);
      D[1, 1] = coefficient(2, 2, 2, 2, i, j, k);
      D[1, 2] = D[2, 1] = coefficient(2, 2, 3, 3, i, j, k);
      D[1, 3] = D[3, 1] = coefficient(2, 2, 2, 3, i, j, k);
      D[1, 4] = D[4, 1] = coefficient(2, 2, 1, 3, i, j, k);
      D[1, 5] = D[5, 1] = coefficient(2, 2, 1, 2, i, j, k);
      D[2, 2] = coefficient(3, 3, 3, 3, i, j, k);
      D[2, 3] = D[3, 2] = coefficient(3, 3, 2, 3, i, j, k);
      D[2, 4] = D[4, 2] = coefficient(3, 3, 1, 3, i, j, k);
      D[2, 5] = D[5, 2] = coefficient(3, 3, 1, 2, i, j, k);
      D[3, 3] = coefficient(2, 3, 3, 2, i, j, k);
      D[3, 4] = D[4, 3] = coefficient(2, 3, 1, 3, i, j, k);
      D[3, 5] = D[5, 3] = coefficient(2, 3, 1, 2, i, j, k);
      D[4, 4] = coefficient(1, 3, 3, 1, i, j, k);
      D[4, 5] = D[5, 4] = coefficient(1, 3, 1, 2, i, j, k);
      D[5, 5] = coefficient(1, 2, 2, 1, i, j, k);
      return D;
    }

    protected DoubleMatrix CreateMaterialMatrixLamda(int i, int j, int k)
    {
      // xx yy zz yz xz xy
      DoubleMatrix D = new DoubleMatrix(6, 6);
      D[0, 0] = coefficientLamda(1, 1, 1, 1, i, j, k);
      D[0, 1] = D[1, 0] = coefficientLamda(1, 1, 2, 2, i, j, k);
      D[0, 2] = D[2, 0] = coefficientLamda(1, 1, 3, 3, i, j, k);
      D[0, 3] = D[3, 0] = coefficientLamda(1, 1, 2, 3, i, j, k);
      D[0, 4] = D[4, 0] = coefficientLamda(1, 1, 1, 3, i, j, k);
      D[0, 5] = D[5, 0] = coefficientLamda(1, 1, 1, 2, i, j, k);
      D[1, 1] = coefficientLamda(2, 2, 2, 2, i, j, k);
      D[1, 2] = D[2, 1] = coefficientLamda(2, 2, 3, 3, i, j, k);
      D[1, 3] = D[3, 1] = coefficientLamda(2, 2, 2, 3, i, j, k);
      D[1, 4] = D[4, 1] = coefficientLamda(2, 2, 1, 3, i, j, k);
      D[1, 5] = D[5, 1] = coefficientLamda(2, 2, 1, 2, i, j, k);
      D[2, 2] = coefficientLamda(3, 3, 3, 3, i, j, k);
      D[2, 3] = D[3, 2] = coefficientLamda(3, 3, 2, 3, i, j, k);
      D[2, 4] = D[4, 2] = coefficientLamda(3, 3, 1, 3, i, j, k);
      D[2, 5] = D[5, 2] = coefficientLamda(3, 3, 1, 2, i, j, k);
      D[3, 3] = coefficientLamda(2, 3, 3, 2, i, j, k);
      D[3, 4] = D[4, 3] = coefficientLamda(2, 3, 1, 3, i, j, k);
      D[3, 5] = D[5, 3] = coefficientLamda(2, 3, 1, 2, i, j, k);
      D[4, 4] = coefficientLamda(1, 3, 3, 1, i, j, k);
      D[4, 5] = D[5, 4] = coefficientLamda(1, 3, 1, 2, i, j, k);
      D[5, 5] = coefficientLamda(1, 2, 2, 1, i, j, k);
      return D;
    }

    protected DoubleMatrix CreateMaterialMatrixMu(int i, int j, int k)
    {
      // xx yy zz yz xz xy
      DoubleMatrix D = new DoubleMatrix(6, 6);
      D[0, 0] = coefficientMu(1, 1, 1, 1, i, j, k);
      D[0, 1] = D[1, 0] = coefficientMu(1, 1, 2, 2, i, j, k);
      D[0, 2] = D[2, 0] = coefficientMu(1, 1, 3, 3, i, j, k);
      D[0, 3] = D[3, 0] = coefficientMu(1, 1, 2, 3, i, j, k);
      D[0, 4] = D[4, 0] = coefficientMu(1, 1, 1, 3, i, j, k);
      D[0, 5] = D[5, 0] = coefficientMu(1, 1, 1, 2, i, j, k);
      D[1, 1] = coefficientMu(2, 2, 2, 2, i, j, k);
      D[1, 2] = D[2, 1] = coefficientMu(2, 2, 3, 3, i, j, k);
      D[1, 3] = D[3, 1] = coefficientMu(2, 2, 2, 3, i, j, k);
      D[1, 4] = D[4, 1] = coefficientMu(2, 2, 1, 3, i, j, k);
      D[1, 5] = D[5, 1] = coefficientMu(2, 2, 1, 2, i, j, k);
      D[2, 2] = coefficientMu(3, 3, 3, 3, i, j, k);
      D[2, 3] = D[3, 2] = coefficientMu(3, 3, 2, 3, i, j, k);
      D[2, 4] = D[4, 2] = coefficientMu(3, 3, 1, 3, i, j, k);
      D[2, 5] = D[5, 2] = coefficientMu(3, 3, 1, 2, i, j, k);
      D[3, 3] = coefficientMu(2, 3, 3, 2, i, j, k);
      D[3, 4] = D[4, 3] = coefficientMu(2, 3, 1, 3, i, j, k);
      D[3, 5] = D[5, 3] = coefficientMu(2, 3, 1, 2, i, j, k);
      D[4, 4] = coefficientMu(1, 3, 3, 1, i, j, k);
      D[4, 5] = D[5, 4] = coefficientMu(1, 3, 1, 2, i, j, k);
      D[5, 5] = coefficientMu(1, 2, 2, 1, i, j, k);
      return D;
    }

    private double coefficient(int i, int j, int k, int l, int ii, int jj, int kk)
    {
      double eModulus = (double)gps[ii, jj, kk].GetValue(DataInGausspoint.EModulus);
      var nu = (double)gps[ii, jj, kk].GetValue(DataInGausspoint.nu);
      double G = eModulus / (2 * (1 + nu));
      double lamda = eModulus * nu / ((1 + nu) * (1 - 2 * nu));
      //MaterialProperty matGeneralNonlocal1 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter1);
      //MaterialProperty matGeneralNonlocal2 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter2);
      //if (matGeneralNonlocal1 != null && matGeneralNonlocal2 != null)
      //{
      //  double nuy1 = (double)gps[ii, jj, kk].GetValue(DataInGausspoint.GeneralNonLocalParameter1);
      //  double nuy2 = (double)gps[ii, jj, kk].GetValue(DataInGausspoint.GeneralNonLocalParameter2);

      //  //nuy1 = ComputeParameterProperty(MaterialPropertyName.GeneralNonLocalParameter1, gps[i, j, k].location[0], gps[i, j, k].location[1], gps[i, j, k].location[2]);
      //  //nuy2 = ComputeParameterProperty(MaterialPropertyName.GeneralNonLocalParameter2, gps[i, j, k].location[0], gps[i, j, k].location[1], gps[i, j, k].location[2]);

      //  return lamda * nuy1 * KroneckerDelta(i, j) * KroneckerDelta(k, l) + G * nuy2 * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
      //}
      //else
      //{
      return lamda * KroneckerDelta(i, j) * KroneckerDelta(k, l) + G * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
      //}
    }

    private double coefficientLamda(int i, int j, int k, int l, int ii, int jj, int kk)
    {
      double eModulus = (double)gps[ii, jj, kk].GetValue(DataInGausspoint.EModulus);
      var nu = (double)gps[ii, jj, kk].GetValue(DataInGausspoint.nu);
      double lamda = eModulus * nu / ((1 + nu) * (1 - 2 * nu));
      return lamda * KroneckerDelta(i, j) * KroneckerDelta(k, l);
    }

    private double coefficientMu(int i, int j, int k, int l, int ii, int jj, int kk)
    {
      double eModulus = (double)gps[ii, jj, kk].GetValue(DataInGausspoint.EModulus);
      var nu = (double)gps[ii, jj, kk].GetValue(DataInGausspoint.nu);
      double G = eModulus / (2 * (1 + nu));
      return G * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
    }

    protected DoubleMatrix CreateMaterialMatrix(double xi, double eta, double zeta)
    {
      // xx yy zz yz xz xy
      DoubleMatrix D = new DoubleMatrix(6, 6);
      D[0, 0] = coefficient(1, 1, 1, 1, xi, eta, zeta);
      D[0, 1] = D[1, 0] = coefficient(1, 1, 2, 2, xi, eta, zeta);
      D[0, 2] = D[2, 0] = coefficient(1, 1, 3, 3, xi, eta, zeta);
      D[0, 3] = D[3, 0] = coefficient(1, 1, 2, 3, xi, eta, zeta);
      D[0, 4] = D[4, 0] = coefficient(1, 1, 1, 3, xi, eta, zeta);
      D[0, 5] = D[5, 0] = coefficient(1, 1, 1, 2, xi, eta, zeta);
      D[1, 1] = coefficient(2, 2, 2, 2, xi, eta, zeta);
      D[1, 2] = D[2, 1] = coefficient(2, 2, 3, 3, xi, eta, zeta);
      D[1, 3] = D[3, 1] = coefficient(2, 2, 2, 3, xi, eta, zeta);
      D[1, 4] = D[4, 1] = coefficient(2, 2, 1, 3, xi, eta, zeta);
      D[1, 5] = D[5, 1] = coefficient(2, 2, 1, 2, xi, eta, zeta);
      D[2, 2] = coefficient(3, 3, 3, 3, xi, eta, zeta);
      D[2, 3] = D[3, 2] = coefficient(3, 3, 2, 3, xi, eta, zeta);
      D[2, 4] = D[4, 2] = coefficient(3, 3, 1, 3, xi, eta, zeta);
      D[2, 5] = D[5, 2] = coefficient(3, 3, 1, 2, xi, eta, zeta);
      D[3, 3] = coefficient(2, 3, 3, 2, xi, eta, zeta);
      D[3, 4] = D[4, 3] = coefficient(2, 3, 1, 3, xi, eta, zeta);
      D[3, 5] = D[5, 3] = coefficient(2, 3, 1, 2, xi, eta, zeta);
      D[4, 4] = coefficient(1, 3, 3, 1, xi, eta, zeta);
      D[4, 5] = D[5, 4] = coefficient(1, 3, 1, 2, xi, eta, zeta);
      D[5, 5] = coefficient(1, 2, 2, 1, xi, eta, zeta);
      return D;
    }

    private double coefficient(int i, int j, int k, int l, double xi, double eta, double zeta)
    {
      //double[] point = null;
      //int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
      //if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
      //{
      //  point = PointAt(xi, eta, zeta);
      //}
      //else
      //{
      //  point = new double[] { xi, eta, zeta };
      //}
      //var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty(point[direction]);
      //var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty(point[direction]);
      var eModulus = ComputeParameterProperty(MaterialPropertyName.YoungModulus, xi, eta, zeta);
      var nu = ComputeParameterProperty(MaterialPropertyName.PoissonRatio, xi, eta, zeta);
      double G = eModulus / (2 * (1 + nu));
      double lamda = eModulus * nu / ((1 + nu) * (1 - 2 * nu));

      //MaterialProperty matGeneralNonlocal1 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter1);
      //MaterialProperty matGeneralNonlocal2 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter2);
      //if (matGeneralNonlocal1 != null && matGeneralNonlocal2 != null)
      //{
      //  double nuy1 = ComputeParameterProperty(MaterialPropertyName.GeneralNonLocalParameter1, xi, eta, zeta);
      //  double nuy2 = ComputeParameterProperty(MaterialPropertyName.GeneralNonLocalParameter2, xi, eta, zeta);
      //  return lamda * nuy1 * KroneckerDelta(i, j) * KroneckerDelta(k, l) + G * nuy2 * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
      //}
      //else
      //{
      return lamda * KroneckerDelta(i, j) * KroneckerDelta(k, l) + G * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
      //}
    }

    public override void ComputeStiffnessMatrixElement(out DoubleMatrix Ke)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      int d = GetCountDimension();
      TrivariateNURBSBasisFunction basis = (TrivariateNURBSBasisFunction)patch.GetVolume().Basis;
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
      DoubleMatrix B = new DoubleMatrix(6, d * nen);
      DoubleMatrix D = null;
      DoubleMatrix DLamda = null;
      DoubleMatrix DMu = null;
      MaterialProperty matGeneralNonlocal1 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter1);
      MaterialProperty matGeneralNonlocal2 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter2);
      if (!(Material is FGMUserDefinedGradedMaterial))
      {
        D = CreateMaterialMatrix();
        if (matGeneralNonlocal1 != null && matGeneralNonlocal2 != null)
        {
          DLamda = CreateMaterialMatrixLamda();
          DMu = CreateMaterialMatrixMu();
        }
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
            double detJ = Math.Abs((double)gps[i, j, k].GetValue(DataInGausspoint.detJ));
            DoubleMatrix dNdx = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdX);
            for (int kk = 0; kk < nen; kk++)
            {
              //                     //xx yy zz yz xz xy
              //B[0, 3 * kk] = dNdx[kk, 0];
              //B[1, 3 * kk + 1] = dNdx[kk, 1];
              //B[2, 3 * kk + 2] = dNdx[kk, 2];
              //B[3, 3 * kk + 1] = dNdx[kk, 2];
              //B[3, 3 * kk + 2] = dNdx[kk, 1];
              //B[4, 3 * kk] = dNdx[kk, 2];
              //B[4, 3 * kk + 2] = dNdx[kk, 0];
              //B[5, 3 * kk] = dNdx[kk, 1];
              //B[5, 3 * kk + 1] = dNdx[kk, 0];

              //xx yy zz xy yz xz
              B[0, 3 * kk] = dNdx[kk, 0];
              B[1, 3 * kk + 1] = dNdx[kk, 1];
              B[2, 3 * kk + 2] = dNdx[kk, 2];
              B[3, 3 * kk] = dNdx[kk, 1];
              B[3, 3 * kk + 1] = dNdx[kk, 0];
              B[4, 3 * kk + 1] = dNdx[kk, 2];
              B[4, 3 * kk + 2] = dNdx[kk, 1];
              B[5, 3 * kk] = dNdx[kk, 2];
              B[5, 3 * kk + 2] = dNdx[kk, 0];
            }
            if (Material is FGMUserDefinedGradedMaterial)
              D = CreateMaterialMatrix(i, j, k);
            DoubleMatrix BTDB = null;
            double nuy1 = 0;
            double nuy2 = 0;
            if (matGeneralNonlocal1 != null && matGeneralNonlocal2 != null)
            {
              //nuy1 = (double)gps[i, j, k].GetValue(DataInGausspoint.GeneralNonLocalParameter1);
              //nuy2 = (double)gps[i, j, k].GetValue(DataInGausspoint.GeneralNonLocalParameter2);
              nuy1 = ComputeParameterProperty(MaterialPropertyName.GeneralNonLocalParameter1, gps[i, j, k].location[0], gps[i, j, k].location[1], gps[i, j, k].location[2]);
              nuy2 = ComputeParameterProperty(MaterialPropertyName.GeneralNonLocalParameter2, gps[i, j, k].location[0], gps[i, j, k].location[1], gps[i, j, k].location[2]);
              if (Material is FGMUserDefinedGradedMaterial)
              {
                DLamda = CreateMaterialMatrixLamda(i, j, k);
                DMu = CreateMaterialMatrixMu(i, j, k);
              }
              DoubleMatrix d3Ndx = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.d3NdX);
              DoubleMatrix BGradient2 = new DoubleMatrix(6, d * nen);
              for (int kk = 0; kk < nen; kk++)
              {
                //xx yy zz xy yz xz
                BGradient2[0, 3 * kk] = d3Ndx[kk, 0] + d3Ndx[kk, 4] + d3Ndx[kk, 8];
                BGradient2[1, 3 * kk + 1] = d3Ndx[kk, 3] + d3Ndx[kk, 1] + d3Ndx[kk, 6];
                BGradient2[2, 3 * kk + 2] = d3Ndx[kk, 7] + d3Ndx[kk, 5] + d3Ndx[kk, 2];
                BGradient2[3, 3 * kk] = d3Ndx[kk, 3] + d3Ndx[kk, 1] + d3Ndx[kk, 6];
                BGradient2[3, 3 * kk + 1] = d3Ndx[kk, 0] + d3Ndx[kk, 4] + d3Ndx[kk, 8];
                BGradient2[4, 3 * kk + 1] = d3Ndx[kk, 7] + d3Ndx[kk, 5] + d3Ndx[kk, 2];
                BGradient2[4, 3 * kk + 2] = d3Ndx[kk, 3] + d3Ndx[kk, 1] + d3Ndx[kk, 6];
                BGradient2[5, 3 * kk] = d3Ndx[kk, 7] + d3Ndx[kk, 5] + d3Ndx[kk, 2];
                BGradient2[5, 3 * kk + 2] = d3Ndx[kk, 0] + d3Ndx[kk, 4] + d3Ndx[kk, 8];
              }
              //BTDB = NMathFunctions.Product(B.Transpose(), MatrixFunctions.Product(D, B) - nuy1 * MatrixFunctions.Product(DLamda, BGradient2) - nuy2 * MatrixFunctions.Product(DMu, BGradient2));
              //BTDB = NMathFunctions.Product(MatrixFunctions.Product(B.Transpose(), D) - nuy1 * MatrixFunctions.Product(BGradient2.Transpose(), DLamda) - nuy2 * MatrixFunctions.Product(BGradient2.Transpose(), DMu), B);
              
              
              BTDB = NMathFunctions.Product(MatrixFunctions.Product(B.Transpose(), D) - MatrixFunctions.Product(BGradient2.Transpose(), nuy1 * DLamda + nuy2 * DMu), B);
            }
            else
            {
              BTDB = NMathFunctions.Product(B.Transpose(), MatrixFunctions.Product(D, B));
            }
            Ke += gps[i, j, k].weight * detJbar * detJ * BTDB;
          }
        }
      }
    }
    public override DoubleVector StrainAt(params double[] xi)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      int d = patch.GetCountDimension();
      int nen = patch.GetCountLocalBasisFunctions();
      var cps = patch.GetVolume().ControlPoints;
      DoubleMatrix B = new DoubleMatrix(6, d * nen);
      DoubleMatrix gradNE = GradBasisFunction(xi);
      DoubleMatrix J = JacobianAt(gradNE);
      DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
      DoubleMatrix gradBasis = MatrixFunctions.Product(gradNE, invertJ);
      for (int kk = 0; kk < nen; kk++)
      {
        //            //xx yy zz yz xz xy
        //            B[0, 3 * kk] = gradBasis[kk, 0];
        //B[1, 3 * kk + 1] = gradBasis[kk, 1];
        //B[2, 3 * kk + 2] = gradBasis[kk, 2];
        //B[3, 3 * kk + 1] = gradBasis[kk, 2];
        //B[3, 3 * kk + 2] = gradBasis[kk, 1];
        //B[4, 3 * kk] = gradBasis[kk, 2];
        //B[4, 3 * kk + 2] = gradBasis[kk, 0];
        //B[5, 3 * kk] = gradBasis[kk, 1];
        //B[5, 3 * kk + 1] = gradBasis[kk, 0];

        //xx yy zz xy yz xz
        B[0, 3 * kk] = gradBasis[kk, 0];
        B[1, 3 * kk + 1] = gradBasis[kk, 1];
        B[2, 3 * kk + 2] = gradBasis[kk, 2];
        B[3, 3 * kk] = gradBasis[kk, 1];
        B[3, 3 * kk + 1] = gradBasis[kk, 0];
        B[4, 3 * kk + 1] = gradBasis[kk, 2];
        B[4, 3 * kk + 2] = gradBasis[kk, 1];
        B[5, 3 * kk] = gradBasis[kk, 2];
        B[5, 3 * kk + 2] = gradBasis[kk, 0];
      }
      DoubleVector U = GetDisplacementLocal();
      return MatrixFunctions.Product(B, U);
    }

    public override DoubleVector StressAt(params double[] xi)
    {
      DoubleMatrix D = null;
      if (!(Material is FGMStructureOneVariableMaterial))
        D = CreateMaterialMatrix();
      else
        D = CreateMaterialMatrix(xi[0], xi[1], xi[2]);
      DoubleVector strain = StrainAt(xi);
      return MatrixFunctions.Product(D, strain);
    }

    public override void ComputeInternalForceElement(out DoubleVector fi)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      int d = GetCountDimension();
      TrivariateNURBSBasisFunction basis = (TrivariateNURBSBasisFunction)patch.GetVolume().Basis;
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
      DoubleMatrix Bu = new DoubleMatrix(6, 3 * nen);

      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {
            GaussPoints gpsijk = gps[i, j,k];
            double detJbar = 1.0 / 8.0
                 * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                 * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2])
                 * (kvNoMulticiply3[idx3 + 1] - kvNoMulticiply3[idx3]);
            DoubleMatrix dNdx = (DoubleMatrix)gpsijk.GetValue(DataInGausspoint.dNdX);//gps[i, j].dNdX;
            DoubleVector N = (DoubleVector)gpsijk.GetValue(DataInGausspoint.Ni);//gps[i, j].Ni;
            double detJ = Math.Abs((double)gpsijk.GetValue(DataInGausspoint.detJ));//Math.Abs(gps[i, j].detJ);
            for (int kk = 0; kk < nen; kk++)
            {
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
            }

            DoubleVector stress = StressAt(gpsijk.location);
            fi += gpsijk.weight * detJbar * detJ * NMathFunctions.Product(Bu.Transpose(), stress);
          }
        }
      }
    }

    public override DoubleVector StrainElasticAt(params double[] xi)
    {
      return StrainAt(xi);
    }

    public override DoubleVector StrainThermoAt(params double[] xi)
    {
      throw new NotImplementedException();
    }


  }
}
