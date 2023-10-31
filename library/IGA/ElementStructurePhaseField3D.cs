using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;
using DEMSoft.Function;
using System.Linq;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public class ElementStructurePhaseField3D : AbstractElementStructure3D
  {
    public ElementStructurePhaseField3D(AbstractPatch3D mesh, int id)
         : base(mesh, id)
    {
    }

    private DoubleMatrix JTensor()
    {
      DoubleMatrix J = new DoubleMatrix(6, 6);
      J[0, 0] = ItemJTensor(1, 1, 1, 1);
      J[0, 1] = J[1, 0] = ItemJTensor(1, 1, 2, 2);
      J[0, 2] = J[2, 0] = ItemJTensor(1, 1, 3, 3);
      J[0, 3] = J[3, 0] = ItemJTensor(1, 1, 2, 3);
      J[0, 4] = J[4, 0] = ItemJTensor(1, 1, 1, 3);
      J[0, 5] = J[5, 0] = ItemJTensor(1, 1, 1, 2);
      J[1, 1] = ItemJTensor(2, 2, 2, 2);
      J[1, 2] = J[2, 1] = ItemJTensor(2, 2, 3, 3);
      J[1, 3] = J[3, 1] = ItemJTensor(2, 2, 2, 3);
      J[1, 4] = J[4, 1] = ItemJTensor(2, 2, 1, 3);
      J[1, 5] = J[5, 1] = ItemJTensor(2, 2, 1, 2);
      J[2, 2] = ItemJTensor(3, 3, 3, 3);
      J[2, 3] = J[3, 2] = ItemJTensor(3, 3, 2, 3);
      J[2, 4] = J[4, 2] = ItemJTensor(3, 3, 1, 3);
      J[2, 5] = J[5, 2] = ItemJTensor(3, 3, 1, 2);
      J[3, 3] = ItemJTensor(2, 3, 3, 2);
      J[3, 4] = J[4, 3] = ItemJTensor(2, 3, 1, 3);
      J[3, 5] = J[5, 3] = ItemJTensor(2, 3, 1, 2);
      J[4, 4] = ItemJTensor(1, 3, 3, 1);
      J[4, 5] = J[5, 4] = ItemJTensor(1, 3, 1, 2);
      J[5, 5] = ItemJTensor(1, 2, 2, 1);
      return J;
    }
    private double ItemJTensor(int i, int j, int k, int l)
    {
      //return 1.0 / 2.0 * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
      //return KroneckerDelta(i, k) * KroneckerDelta(j, l);//Liu 2016
      return KroneckerDelta(i, j) * KroneckerDelta(k, l);//Formular
    }
    //protected double KroneckerDelta(int i, int j)
    //{
    //  return (i == j) ? 1.0 : 0.0;
    //}

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
      int d = patch.GetCountField();
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

      DoubleMatrix Keuu = null;
      DoubleMatrix Kepp = null;
      DoubleMatrix Keup = null;
      DoubleMatrix Kepu = null;
      DoubleMatrix Bu = new DoubleMatrix(6, 3 * nen);
      DoubleMatrix Bp = new DoubleMatrix(3, nen);
      switch (AbstractModel.TypeSchemeNonlinearSolver)
      {
        case TypeSchemeNonlinearSolver.Monolithic:
          Keuu = new DoubleMatrix(3 * nen, 3 * nen);
          Kepp = new DoubleMatrix(nen, nen);
          Keup = new DoubleMatrix(3 * nen, nen);
          Kepu = new DoubleMatrix(nen, 3 * nen);
          break;
        case TypeSchemeNonlinearSolver.SingleStaggered:
          Keuu = new DoubleMatrix(3 * nen, 3 * nen);
          Kepp = new DoubleMatrix(nen, nen);
          break;
        case TypeSchemeNonlinearSolver.Staggered:
          if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
          {
            Keuu = new DoubleMatrix(3 * nen, 3 * nen);
          }
          else
          {
            Kepp = new DoubleMatrix(nen, nen);
          }
          break;
      }
      double l0 = ModelStructurePhaseFieldStatic.l0;
      //double lch = ModelStructurePhaseFieldStatic.lch;
      double k1 = ModelStructurePhaseFieldStatic.k;
      DoubleMatrix BuT = null;
      DoubleMatrix BpT = null;
      double previousPhase = 0;
      double currentPhase = 0;
      DoubleVector ue = null;////////////
      DoubleVector strain = null;
      DoubleMatrix DsigmaDEpsilon = null;

      DoubleVector principleStrain;
      DoubleVector[] eigenVectors;
      DoubleMatrix BuTDBu = null;
      DoubleMatrix kepp = null;
      double xiEpsilon = 0;
      double xiEpsilonPos = 0;
      double xiEpsilonNeg = 0;
      double lamda = 0, mu = 0, gc = 0, ft = 0, lch = 0;
      if (!(Material is FGMUserDefinedGradedMaterial))
      {
        lamda = Material.GetProperty(MaterialPropertyName.LameParameter).GetValueProperty();
        mu = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
        gc = Material.GetProperty(MaterialPropertyName.CriticalEnergyReleaseRate).GetValueProperty();
        ft = Material.GetProperty(MaterialPropertyName.TensileYieldStrength).GetValueProperty();
      }

      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {
            GaussPoints gpsijk = gps[i, j, k];
            double detJbar = 1.0 / 8.0
                 * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                 * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2])
                 * (kvNoMulticiply3[idx3 + 1] - kvNoMulticiply3[idx3]);
            double detJ = Math.Abs((double)gpsijk.GetValue(DataInGausspoint.detJ));
            DoubleVector N = (DoubleVector)gpsijk.GetValue(DataInGausspoint.Ni);//gps[i, j].Ni;
            DoubleMatrix dNdx = (DoubleMatrix)gpsijk.GetValue(DataInGausspoint.dNdX);
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

              //x y z
              Bp[0, kk] = dNdx[kk, 0];
              Bp[1, kk] = dNdx[kk, 1];
              Bp[2, kk] = dNdx[kk, 2];
            }
            if (Material is FGMUserDefinedGradedMaterial)
            {
              lamda = ComputeParameterProperty(MaterialPropertyName.LameParameter, gpsijk.location[0], gpsijk.location[1], gpsijk.location[2]);
              mu = ComputeParameterProperty(MaterialPropertyName.ShearModulus, gpsijk.location[0], gpsijk.location[1], gpsijk.location[2]);
              gc = ComputeParameterProperty(MaterialPropertyName.CriticalEnergyReleaseRate, gpsijk.location[0], gpsijk.location[1], gpsijk.location[2]);
              ft = ComputeParameterProperty(MaterialPropertyName.TensileYieldStrength, gpsijk.location[0], gpsijk.location[1], gpsijk.location[2]);
            }
            double E = mu * (3 * lamda + 2 * mu) / (lamda + mu);
            lch = E * gc / (ft * ft);

            DoubleMatrix NTN = NMathFunctions.OuterProduct(N, N);
            double valDDOmegaD, valCAlpha, valDDAlphaD;
            switch (AbstractModel.TypeSchemeNonlinearSolver)
            {
              case TypeSchemeNonlinearSolver.Monolithic:
                #region monolithic
                BuT = Bu.Transpose();
                BpT = Bp.Transpose();
                currentPhase = (double)gpsijk.GetValue(DataInGausspoint.currentPhase);

                strain = (DoubleVector)gpsijk.GetValue(DataInGausspoint.currentStrain);
                ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
                DsigmaDEpsilon = CreateMaterialMatrixIsotropicAnisotropic(lamda, mu, k1, l0, lch, currentPhase, principleStrain, eigenVectors);
                BuTDBu = NMathFunctions.Product(BuT, NMathFunctions.Product(DsigmaDEpsilon, Bu));
                Keuu += gpsijk.weight * detJbar * detJ * BuTDBu;

                xiEpsilon = (double)gpsijk.GetValue(DataInGausspoint.xiEpsilon);
                valDDOmegaD = ddOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                valDDAlphaD = ddAlphaD(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                if (ModelStructurePhaseFieldStatic.typeOrderPhaseField == TypeOrderPhaseField.SecondOrder)
                {
                  kepp = 2.0 * gc * l0 / valCAlpha * NMathFunctions.Product(BpT, Bp) + (valDDAlphaD * gc / (valCAlpha * l0) + valDDOmegaD * xiEpsilon) * NTN;
                }
                else
                {
                  DoubleMatrix ddNdX = (DoubleMatrix)gpsijk.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
                  DoubleVector BB = new DoubleVector(nen);
                  for (int kk = 0; kk < nen; kk++)
                  {
                    BB[kk] = ddNdX[kk, 0] + ddNdX[kk, 1] + ddNdX[kk, 2];
                  }
                  DoubleMatrix F4 = gc * Math.Pow(l0, 3) / (8.0 * valCAlpha) * NMathFunctions.OuterProduct(BB, BB);
                  kepp = gc * l0 / valCAlpha * NMathFunctions.Product(BpT, Bp) + (valDDAlphaD * gc / (l0 * valCAlpha) + valDDOmegaD * xiEpsilon) * NTN + F4;
                }
                Kepp += gpsijk.weight * detJbar * detJ * kepp;

                double valDOmegaD = dOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                DoubleVector dSdphi = valDOmegaD * ComputeDStressDPhi(lamda, mu, principleStrain, eigenVectors);
                Keup += gpsijk.weight * detJbar * detJ * NMathFunctions.Product(BuT, NMathFunctions.OuterProduct(dSdphi, N));
                DoubleVector dHDEpsilon = valDOmegaD * ComputeDHDEpsilon(lamda, mu, principleStrain, eigenVectors);
                Kepu += gpsijk.weight * detJbar * detJ * NMathFunctions.Product(NMathFunctions.OuterProduct(N, dHDEpsilon), Bu);
                #endregion
                break;
              case TypeSchemeNonlinearSolver.SingleStaggered:
                #region singlestaggered
                BuT = Bu.Transpose();
                BpT = Bp.Transpose();
                currentPhase = (double)gpsijk.GetValue(DataInGausspoint.currentPhase);

                strain = (DoubleVector)gpsijk.GetValue(DataInGausspoint.currentStrain);
                ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
                DsigmaDEpsilon = CreateMaterialMatrixIsotropicAnisotropic(lamda, mu, k1, l0, lch, currentPhase, principleStrain, eigenVectors);
                BuTDBu = NMathFunctions.Product(BuT, NMathFunctions.Product(DsigmaDEpsilon, Bu));
                Keuu += gpsijk.weight * detJbar * detJ * BuTDBu;

                xiEpsilon = (double)gpsijk.GetValue(DataInGausspoint.xiEpsilon);
                valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);//pi
                valDDAlphaD = ddAlphaD(ModelStructurePhaseFieldStatic.typePhaseFieldModel);//-2
                valDDOmegaD = ddOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                if (ModelStructurePhaseFieldStatic.typeOrderPhaseField == TypeOrderPhaseField.SecondOrder)
                {
                  kepp = 2.0 * gc * l0 / valCAlpha * NMathFunctions.Product(BpT, Bp) + (valDDAlphaD * gc / (l0 * valCAlpha) + valDDOmegaD * xiEpsilon) * NTN;
                }
                else
                {
                  DoubleMatrix ddNdX = (DoubleMatrix)gpsijk.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
                  DoubleVector BB = new DoubleVector(nen);
                  for (int kk = 0; kk < nen; kk++)
                  {
                    BB[kk] = ddNdX[kk, 0] + ddNdX[kk, 1] + ddNdX[kk, 2];
                  }
                  DoubleMatrix F4 = gc * Math.Pow(l0, 3) / (8.0 * valCAlpha) * NMathFunctions.OuterProduct(BB, BB);////
                  kepp = gc * l0 / valCAlpha * NMathFunctions.Product(BpT, Bp) + (valDDAlphaD * gc / (l0 * valCAlpha) + valDDOmegaD * xiEpsilon) * NTN + F4;
                }
                Kepp += gpsijk.weight * detJbar * detJ * kepp;
                #endregion
                break;
              case TypeSchemeNonlinearSolver.Staggered:
                #region staggered
                if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
                {
                  //displacement
                  BuT = Bu.Transpose();
                  currentPhase = (double)gpsijk.GetValue(DataInGausspoint.currentPhase);
                  strain = (DoubleVector)gpsijk.GetValue(DataInGausspoint.currentStrain);
                  ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
                  DsigmaDEpsilon = CreateMaterialMatrixIsotropicAnisotropic(lamda, mu, k1, l0, lch, currentPhase, principleStrain, eigenVectors);
                  BuTDBu = NMathFunctions.Product(BuT, NMathFunctions.Product(DsigmaDEpsilon, Bu));
                  Keuu += gpsijk.weight * detJbar * detJ * BuTDBu;
                }
                else
                {
                  BpT = Bp.Transpose();
                  xiEpsilon = (double)gpsijk.GetValue(DataInGausspoint.xiEpsilon);
                  currentPhase = (double)gpsijk.GetValue(DataInGausspoint.currentPhase);
                  valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);//pi
                  valDDAlphaD = ddAlphaD(ModelStructurePhaseFieldStatic.typePhaseFieldModel);//-2
                  valDDOmegaD = ddOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                  if (ModelStructurePhaseFieldStatic.typeOrderPhaseField == TypeOrderPhaseField.SecondOrder)
                  {
                    // l0 biểu diễn độ rộng vết nứt
                    //Singh2018
                    kepp = 2.0 * gc * l0 / valCAlpha * NMathFunctions.Product(BpT, Bp) + (valDDAlphaD * gc / (l0 * valCAlpha) + valDDOmegaD * xiEpsilon) * NTN;
                  }
                  else
                  {
                    DoubleMatrix ddNdX = (DoubleMatrix)gpsijk.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
                    DoubleVector BB = new DoubleVector(nen);
                    for (int kk = 0; kk < nen; kk++)
                    {
                      BB[kk] = ddNdX[kk, 0] + ddNdX[kk, 1] + ddNdX[kk, 2];
                    }
                    // l0 biểu diễn độ rộng vết nứt
                    DoubleMatrix F4 = gc * Math.Pow(l0, 3) / (8.0 * valCAlpha) * NMathFunctions.OuterProduct(BB, BB);////
                    kepp = gc * l0 / valCAlpha * NMathFunctions.Product(BpT, Bp) + (valDDAlphaD * gc / (l0 * valCAlpha) + valDDOmegaD * xiEpsilon) * NTN + F4;
                  }
                  Kepp += gpsijk.weight * detJbar * detJ * kepp;
                }
                #endregion
                break;
            }
          }
        }
      }

      for (int j = 0; j < nen; j++)
        for (int i = 0; i < nen; i++)
        {
          switch (AbstractModel.TypeSchemeNonlinearSolver)
          {
            case TypeSchemeNonlinearSolver.Monolithic:
              for (int jj = 0; jj < 3; jj++)
              {
                for (int ii = 0; ii < 3; ii++)
                {
                  Ke[i * 4 + ii, j * 4 + jj] = Keuu[i * 3 + ii, j * 3 + jj];
                  Ke[i * 4 + ii, j * 4 + 3] = Keup[i * 3 + ii, j];
                }
                Ke[i * 4 + 3, j * 4 + jj] = Kepu[i, j * 3 + jj];
              }
              Ke[i * 4 + 3, j * 4 + 3] = Kepp[i, j];
              break;
            case TypeSchemeNonlinearSolver.SingleStaggered:
              for (int jj = 0; jj < 3; jj++)
              {
                for (int ii = 0; ii < 3; ii++)
                {
                  Ke[i * 4 + ii, j * 4 + jj] = Keuu[i * 3 + ii, j * 3 + jj];
                }
              }
              Ke[i * 4 + 3, j * 4 + 3] = Kepp[i, j];
              break;
            case TypeSchemeNonlinearSolver.Staggered:
              if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
              {
                for (int jj = 0; jj < 3; jj++)
                {
                  for (int ii = 0; ii < 3; ii++)
                  {
                    Ke[i * 4 + ii, j * 4 + jj] = Keuu[i * 3 + ii, j * 3 + jj];
                  }
                }
              }
              else
              {
                Ke[i * 4 + 3, j * 4 + 3] = Kepp[i, j];
              }
              break;
          }
        }
    }

    private DoubleMatrix CreateMaterialMatrixIsotropicAnisotropic(double lamda, double mu, double k1, double l0, double lch, double phiCurrent, DoubleVector principleStrain, DoubleVector[] eigenVectors)
    {
      DoubleMatrix DsigmaDEpsilon = null;
      switch (ModelStructurePhaseFieldStatic.typeDegradationModel)
      {
        case TypeDegradationEnergy.Isotropic:
        case TypeDegradationEnergy.Hybrid:
          DsigmaDEpsilon = CreateMaterialMatrixIsotropic(phiCurrent, k1, l0, lch, lamda, mu);
          break;
        case TypeDegradationEnergy.Anisotropic:
          DsigmaDEpsilon = CreateMaterialMatrixAnisotropic(principleStrain, eigenVectors, phiCurrent, k1, l0, lch, lamda, mu);
          break;
      }
      return DsigmaDEpsilon;
    }

    private DoubleMatrix CreateMaterialMatrixIsotropic(double phi, double k, double l0, double lch, double lamda, double mu)
    {
      AnisotropicElasticity prob = (AnisotropicElasticity)Material.GetProperty(MaterialPropertyName.AnisotropicElasticity);
      DoubleMatrix C = null;
      if (prob != null)
      {
        double[,] ce = prob.GetAnisotropicElasticityMatrix();
        C = new DoubleMatrix(6, 6);
        bool isInverse = prob.GetIsInverseMatrix();
        for (int i = 0; i < 6; i++)
          for (int j = 0; j < 6; j++)
          {
            C[i, j] = ce[i, j];
          }
        if (isInverse)
          C = NMathFunctions.Inverse(C);
      }
      else
      {
        OrthotropicElasticity prob1 = (OrthotropicElasticity)Material.GetProperty(MaterialPropertyName.OrthotropicElasticity);
        if (prob1 == null)
        {
          C = new DoubleMatrix(6, 6);
          //Xem lai thu tu xx yy zz yz xz xy
          C[0, 0] = coefficient(1, 1, 1, 1, lamda, mu);
          C[0, 1] = C[1, 0] = coefficient(1, 1, 2, 2, lamda, mu);
          C[0, 2] = C[2, 0] = coefficient(1, 1, 3, 3, lamda, mu);
          C[0, 3] = C[3, 0] = coefficient(1, 1, 2, 3, lamda, mu);
          C[0, 4] = C[4, 0] = coefficient(1, 1, 1, 3, lamda, mu);
          C[0, 5] = C[5, 0] = coefficient(1, 1, 1, 2, lamda, mu);
          C[1, 1] = coefficient(2, 2, 2, 2, lamda, mu);
          C[1, 2] = C[2, 1] = coefficient(2, 2, 3, 3, lamda, mu);
          C[1, 3] = C[3, 1] = coefficient(2, 2, 2, 3, lamda, mu);
          C[1, 4] = C[4, 1] = coefficient(2, 2, 1, 3, lamda, mu);
          C[1, 5] = C[5, 1] = coefficient(2, 2, 1, 2, lamda, mu);
          C[2, 2] = coefficient(3, 3, 3, 3, lamda, mu);
          C[2, 3] = C[3, 2] = coefficient(3, 3, 2, 3, lamda, mu);
          C[2, 4] = C[4, 2] = coefficient(3, 3, 1, 3, lamda, mu);
          C[2, 5] = C[5, 2] = coefficient(3, 3, 1, 2, lamda, mu);
          C[3, 3] = coefficient(2, 3, 3, 2, lamda, mu);
          C[3, 4] = C[4, 3] = coefficient(2, 3, 1, 3, lamda, mu);
          C[3, 5] = C[5, 3] = coefficient(2, 3, 1, 2, lamda, mu);
          C[4, 4] = coefficient(1, 3, 3, 1, lamda, mu);
          C[4, 5] = C[5, 4] = coefficient(1, 3, 1, 2, lamda, mu);
          C[5, 5] = coefficient(1, 2, 2, 1, lamda, mu);
        }
        else
        {
          C = new DoubleMatrix(prob1.GetOrthotropicElasticityMatrix());
        }
      }
      return (OmegaD(phi, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel) + k) * C;
    }

    private double coefficient(int i, int j, int k, int l, double lamda, double mu)
    {
      //https://en.wikipedia.org/wiki/Linear_elasticity
      //https://en.wikipedia.org/wiki/Hooke%27s_law#Isotropic_materials
      return lamda * KroneckerDelta(i, j) * KroneckerDelta(k, l) + mu * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
      //return K * KroneckerDelta(i, j) * KroneckerDelta(k, l) + mu * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k) - 2.0 / 3.0 * KroneckerDelta(i, j) * KroneckerDelta(k, l));
    }

    private DoubleMatrix CreateMaterialMatrixAnisotropic(DoubleVector principleStrain, DoubleVector[] eigenVectors, double phi, double k, double l0, double lch, double lamda, double mu)
    {
      //if (principleStrain[0] != 0 || principleStrain[1] != 0 || principleStrain[2] != 0)
      //{
      double trStrain = principleStrain[0] + principleStrain[1] + principleStrain[2];
      DoubleMatrix PPositive = new DoubleMatrix(6, 6);
      DoubleMatrix PNegative = new DoubleMatrix(6, 6);
      ComputeAnisotropicProjectionTensor(principleStrain, eigenVectors, out PPositive, out PNegative);
      DoubleMatrix J = JTensor();
      DoubleMatrix C = (OmegaD(phi, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel) + k) * (lamda * HeavisideFunction(trStrain) * J + 2 * mu * PPositive)
        + (lamda * HeavisideFunction(-trStrain) * J + 2 * mu * PNegative);
      return C;
      //}
      //else
      //{
      //  return (Math.Pow(1 - phi, 2) + k) * CreateMaterialMatrix();
      //}
    }

    private void ComputeAnisotropicProjectionTensor(DoubleVector principleStrain, DoubleVector[] eigenVectors, out DoubleMatrix PPositive, out DoubleMatrix PNegative)
    {
      DoubleMatrix[] Ma = ComputeAnisotropicMa(eigenVectors);
      PPositive = new DoubleMatrix(6, 6);
      PNegative = new DoubleMatrix(6, 6);
      for (int a = 0; a < 3; a++)
      {
        DoubleMatrix secondTermPos = new DoubleMatrix(6, 6);
        DoubleMatrix secondTermNeg = new DoubleMatrix(6, 6);
        for (int b = 0; b < 3; b++)
        {
          if (b != a)
          {
            DoubleMatrix temp = 1.0 / 4.0 * (ComputeAnisotropicGabMatrix(Ma[a], Ma[b]) + ComputeAnisotropicGabMatrix(Ma[b], Ma[a]));//current
            double dab = principleStrain[a] - principleStrain[b];
            secondTermPos += HeavisideFunction(dab) * temp;
            secondTermNeg += HeavisideFunction(-dab) * temp;
            /////////////////////////////////////////////////////////////////
            //DoubleMatrix temp = 1.0 / 2.0 * (ComputeAnisotropicGabMatrix(Ma[a], Ma[b]) + ComputeAnisotropicGabMatrix(Ma[b], Ma[a]));//current
            //double dab = principleStrain[a] - principleStrain[b];
            //secondTermPos += HeavisideFunction(dab) * temp;
            //secondTermNeg += HeavisideFunction(-dab) * temp;
            ////////////////////////////////////////////////////////////////////////
            //DoubleMatrix temp = 1.0 / 2.0 * (ComputeAnisotropicGabMatrix(Ma[a], Ma[b]) + ComputeAnisotropicGabMatrix(Ma[b], Ma[a]));//current
            //double dab = principleStrain[a] - principleStrain[b];
            //if (Math.Abs(dab) <= 1e-4)
            //  dab = Math.Abs(dab);
            //secondTermPos += HeavisideFunction(dab) * temp;
            //secondTermNeg += HeavisideFunction(-dab) * temp;
            /////////////////////////////////////////////////////////////////
            //DoubleMatrix temp = 1.0 / 2.0 * (ComputeAnisotropicGabMatrix(Ma[a], Ma[b]) + ComputeAnisotropicGabMatrix(Ma[b], Ma[a]));
            //secondTermPos += (PositiveRampFunction(prinTemp[a]) - PositiveRampFunction(prinTemp[b])) / (prinTemp[a] - prinTemp[b]) * temp;
            //secondTermNeg += (NegativeRampFunction(prinTemp[a]) - NegativeRampFunction(prinTemp[b])) / (prinTemp[a] - prinTemp[b]) * temp;
          }
        }

        DoubleMatrix Qa = ComputeAnisotropicQaMatrix(Ma[a]);
        PPositive += HeavisideFunction(principleStrain[a]) * Qa + secondTermPos;
        PNegative += HeavisideFunction(-principleStrain[a]) * Qa + secondTermNeg;
      }
    }
    private DoubleMatrix[] ComputeAnisotropicMa(DoubleVector[] eigenVectors)
    {

      DoubleMatrix[] Ma = new DoubleMatrix[3];
      Ma[0] = MatrixFunctions.OuterProduct(eigenVectors[0], eigenVectors[0]);
      Ma[1] = MatrixFunctions.OuterProduct(eigenVectors[1], eigenVectors[1]);
      Ma[2] = MatrixFunctions.OuterProduct(eigenVectors[2], eigenVectors[2]);
      return Ma;
    }
    private double ComputeAnisotropicQa(DoubleMatrix Ma, int i, int j, int k, int l)
    {
      return Ma[i - 1, j - 1] * Ma[k - 1, l - 1];
    }
    private DoubleMatrix ComputeAnisotropicQaMatrix(DoubleMatrix Ma)
    {
      DoubleMatrix Qa = new DoubleMatrix(6, 6);
      Qa[0, 0] = ComputeAnisotropicQa(Ma, 1, 1, 1, 1);
      Qa[0, 1] = Qa[1, 0] = ComputeAnisotropicQa(Ma, 1, 1, 2, 2);
      Qa[0, 2] = Qa[2, 0] = ComputeAnisotropicQa(Ma, 1, 1, 3, 3);
      Qa[0, 3] = Qa[3, 0] = ComputeAnisotropicQa(Ma, 1, 1, 2, 3);
      Qa[0, 4] = Qa[4, 0] = ComputeAnisotropicQa(Ma, 1, 1, 1, 3);
      Qa[0, 5] = Qa[5, 0] = ComputeAnisotropicQa(Ma, 1, 1, 1, 2);
      Qa[1, 1] = ComputeAnisotropicQa(Ma, 2, 2, 2, 2);
      Qa[1, 2] = Qa[2, 1] = ComputeAnisotropicQa(Ma, 2, 2, 3, 3);
      Qa[1, 3] = Qa[3, 1] = ComputeAnisotropicQa(Ma, 2, 2, 2, 3);
      Qa[1, 4] = Qa[4, 1] = ComputeAnisotropicQa(Ma, 2, 2, 1, 3);
      Qa[1, 5] = Qa[5, 1] = ComputeAnisotropicQa(Ma, 2, 2, 1, 2);
      Qa[2, 2] = ComputeAnisotropicQa(Ma, 3, 3, 3, 3);
      Qa[2, 3] = Qa[3, 2] = ComputeAnisotropicQa(Ma, 3, 3, 2, 3);
      Qa[2, 4] = Qa[4, 2] = ComputeAnisotropicQa(Ma, 3, 3, 1, 3);
      Qa[2, 5] = Qa[5, 2] = ComputeAnisotropicQa(Ma, 3, 3, 1, 2);
      Qa[3, 3] = ComputeAnisotropicQa(Ma, 2, 3, 3, 2);
      Qa[3, 4] = Qa[4, 3] = ComputeAnisotropicQa(Ma, 2, 3, 1, 3);
      Qa[3, 5] = Qa[5, 3] = ComputeAnisotropicQa(Ma, 2, 3, 1, 2);
      Qa[4, 4] = ComputeAnisotropicQa(Ma, 1, 3, 3, 1);
      Qa[4, 5] = Qa[5, 4] = ComputeAnisotropicQa(Ma, 1, 3, 1, 2);
      Qa[5, 5] = ComputeAnisotropicQa(Ma, 1, 2, 2, 1);
      return Qa;
    }
    private double ComputeAnisotropicGab(DoubleMatrix Ma, DoubleMatrix Mb, int i, int j, int k, int l)
    {
      return Ma[i - 1, k - 1] * Mb[j - 1, l - 1] + Ma[i - 1, l - 1] * Mb[j - 1, k - 1];
    }
    private DoubleMatrix ComputeAnisotropicGabMatrix(DoubleMatrix Ma, DoubleMatrix Mb)
    {
      DoubleMatrix Gab = new DoubleMatrix(6, 6);
      Gab[0, 0] = ComputeAnisotropicGab(Ma, Mb, 1, 1, 1, 1);
      Gab[0, 1] = Gab[1, 0] = ComputeAnisotropicGab(Ma, Mb, 1, 1, 2, 2);
      Gab[0, 2] = Gab[2, 0] = ComputeAnisotropicGab(Ma, Mb, 1, 1, 3, 3);
      Gab[0, 3] = Gab[3, 0] = ComputeAnisotropicGab(Ma, Mb, 1, 1, 2, 3);
      Gab[0, 4] = Gab[4, 0] = ComputeAnisotropicGab(Ma, Mb, 1, 1, 1, 3);
      Gab[0, 5] = Gab[5, 0] = ComputeAnisotropicGab(Ma, Mb, 1, 1, 1, 2);
      Gab[1, 1] = ComputeAnisotropicGab(Ma, Mb, 2, 2, 2, 2);
      Gab[1, 2] = Gab[2, 1] = ComputeAnisotropicGab(Ma, Mb, 2, 2, 3, 3);
      Gab[1, 3] = Gab[3, 1] = ComputeAnisotropicGab(Ma, Mb, 2, 2, 2, 3);
      Gab[1, 4] = Gab[4, 1] = ComputeAnisotropicGab(Ma, Mb, 2, 2, 1, 3);
      Gab[1, 5] = Gab[5, 1] = ComputeAnisotropicGab(Ma, Mb, 2, 2, 1, 2);
      Gab[2, 2] = ComputeAnisotropicGab(Ma, Mb, 3, 3, 3, 3);
      Gab[2, 3] = Gab[3, 2] = ComputeAnisotropicGab(Ma, Mb, 3, 3, 2, 3);
      Gab[2, 4] = Gab[4, 2] = ComputeAnisotropicGab(Ma, Mb, 3, 3, 1, 3);
      Gab[2, 5] = Gab[5, 2] = ComputeAnisotropicGab(Ma, Mb, 3, 3, 1, 2);
      Gab[3, 3] = ComputeAnisotropicGab(Ma, Mb, 2, 3, 3, 2);
      Gab[3, 4] = Gab[4, 3] = ComputeAnisotropicGab(Ma, Mb, 2, 3, 1, 3);
      Gab[3, 5] = Gab[5, 3] = ComputeAnisotropicGab(Ma, Mb, 2, 3, 1, 2);
      Gab[4, 4] = ComputeAnisotropicGab(Ma, Mb, 1, 3, 3, 1);
      Gab[4, 5] = Gab[5, 4] = ComputeAnisotropicGab(Ma, Mb, 1, 3, 1, 2);
      Gab[5, 5] = ComputeAnisotropicGab(Ma, Mb, 1, 2, 2, 1);
      return Gab;
    }
    private double HeavisideFunction(double a)
    {
      if (a > 0)
        return 1;
      else if (a == 0)
        return 0.5;
      else
        return 0;
      //if (a > 0)
      //  return 1;
      //else
      //  return 0;
    }
    private void ComputePrincipleStrains(DoubleVector strain, out DoubleVector eigenValues, out DoubleVector[] eigenVectors)
    {
      DoubleMatrix strainTensor = new DoubleMatrix(3, 3);
      strainTensor[0, 0] = strain[0];
      strainTensor[1, 1] = strain[1];
      strainTensor[2, 2] = strain[2];
      strainTensor[0, 1] = strainTensor[1, 0] = strain[3] / 2.0;
      strainTensor[1, 2] = strainTensor[2, 1] = strain[4] / 2.0;
      strainTensor[0, 2] = strainTensor[2, 0] = strain[5] / 2.0;
      var eigDecomp = new DoubleSymEigDecomp(strainTensor);
      eigenValues = eigDecomp.EigenValues;
      eigenVectors = new DoubleVector[3];
      eigenVectors[0] = eigDecomp.EigenVector(0);
      eigenVectors[1] = eigDecomp.EigenVector(1);
      eigenVectors[2] = eigDecomp.EigenVector(2);
    }
    private DoubleVector ComputeDStressDPhi(double lamda, double mu, DoubleVector principleStrain, DoubleVector[] eigenVectors)
    {
      DoubleVector stress = null;
      switch (ModelStructurePhaseFieldStatic.typeDegradationModel)
      {
        case TypeDegradationEnergy.Isotropic:
        case TypeDegradationEnergy.Hybrid:
          DoubleVector stress0 = null;
          ComputeStressIsotropic(principleStrain, eigenVectors, lamda, mu, out stress0);
          stress = stress0;
          break;
        case TypeDegradationEnergy.Anisotropic:
          DoubleVector stress0Pos = null;
          DoubleVector stress0Neg = null;
          ComputeStressAnisotropic(principleStrain, eigenVectors, lamda, mu, out stress0Pos, out stress0Neg);
          stress = stress0Pos;
          break;
      }

      return stress;
    }

    private void ComputeStressIsotropic(DoubleVector principleStrain, DoubleVector[] eigenVectors, double lamda, double mu, out DoubleVector stress0)
    {
      var prinVector1 = eigenVectors[0];
      var prinVector2 = eigenVectors[1];
      var prinVector3 = eigenVectors[2];

      double trEpsilon = principleStrain[0] + principleStrain[1] + principleStrain[2];
      DoubleMatrix stressTensor = (lamda * trEpsilon + 2 * mu * principleStrain[0]) * NMathFunctions.OuterProduct(prinVector1, prinVector1)
        + (lamda * trEpsilon + 2 * mu * principleStrain[1]) * NMathFunctions.OuterProduct(prinVector2, prinVector2)
        + (lamda * trEpsilon + 2 * mu * principleStrain[2]) * NMathFunctions.OuterProduct(prinVector3, prinVector3);
      stress0 = new DoubleVector(6);
      stress0[0] = stressTensor[0, 0];
      stress0[1] = stressTensor[1, 1];
      stress0[2] = stressTensor[2, 2];
      stress0[3] = stressTensor[0, 1];
      stress0[4] = stressTensor[1, 2];
      stress0[5] = stressTensor[0, 2];
    }

    private void ComputeStressAnisotropic(DoubleVector principleStrain, DoubleVector[] eigenVectors, double lamda, double mu, out DoubleVector stress0Pos, out DoubleVector stress0Neg)
    {
      var prinVector1 = eigenVectors[0];
      var prinVector2 = eigenVectors[1];
      var prinVector3 = eigenVectors[2];

      double trEpsilon = principleStrain[0] + principleStrain[1] + principleStrain[2];
      DoubleMatrix unit = new DoubleMatrix(3, 3);
      unit[0, 0] = unit[1, 1] = unit[2, 2] = 1.0;
      //DoubleMatrix stressTensorPos = (lamda * PositiveRampFunction(trEpsilon) + 2 * mu * PositiveRampFunction(principleStrain[0])) * NMathFunctions.OuterProduct(prinVector1, prinVector1)
      //+ (lamda * PositiveRampFunction(trEpsilon) + 2 * mu * PositiveRampFunction(principleStrain[1])) * NMathFunctions.OuterProduct(prinVector2, prinVector2)
      //+ (lamda * PositiveRampFunction(trEpsilon) + 2 * mu * PositiveRampFunction(principleStrain[2])) * NMathFunctions.OuterProduct(prinVector3, prinVector3);
      //DoubleMatrix stressTensorNeg = (lamda * NegativeRampFunction(trEpsilon) + 2 * mu * NegativeRampFunction(principleStrain[0])) * NMathFunctions.OuterProduct(prinVector1, prinVector1)
      //+ (lamda * NegativeRampFunction(trEpsilon) + 2 * mu * NegativeRampFunction(principleStrain[1])) * NMathFunctions.OuterProduct(prinVector2, prinVector2)
      //+ (lamda * NegativeRampFunction(trEpsilon) + 2 * mu * NegativeRampFunction(principleStrain[2])) * NMathFunctions.OuterProduct(prinVector3, prinVector3);

      DoubleMatrix stressTensorPos = lamda * PositiveRampFunction(trEpsilon) * unit +
      2 * mu * (PositiveRampFunction(principleStrain[0]) * NMathFunctions.OuterProduct(prinVector1, prinVector1)
      + PositiveRampFunction(principleStrain[1]) * NMathFunctions.OuterProduct(prinVector2, prinVector2)
      + PositiveRampFunction(principleStrain[2]) * NMathFunctions.OuterProduct(prinVector3, prinVector3));
      DoubleMatrix stressTensorNeg = lamda * NegativeRampFunction(trEpsilon) * unit + 2 * mu * (NegativeRampFunction(principleStrain[0]) * NMathFunctions.OuterProduct(prinVector1, prinVector1)
        + NegativeRampFunction(principleStrain[1]) * NMathFunctions.OuterProduct(prinVector2, prinVector2)
        + NegativeRampFunction(principleStrain[2]) * NMathFunctions.OuterProduct(prinVector3, prinVector3));

      stress0Pos = new DoubleVector(6);
      stress0Neg = new DoubleVector(6);
      stress0Pos[0] = stressTensorPos[0, 0];
      stress0Pos[1] = stressTensorPos[1, 1];
      stress0Pos[2] = stressTensorPos[2, 2];
      stress0Pos[3] = stressTensorPos[0, 1];
      stress0Pos[4] = stressTensorPos[1, 2];
      stress0Pos[5] = stressTensorPos[0, 2];

      stress0Neg[0] = stressTensorNeg[0, 0];
      stress0Neg[1] = stressTensorNeg[1, 1];
      stress0Neg[2] = stressTensorNeg[2, 2];
      stress0Neg[3] = stressTensorNeg[0, 1];
      stress0Neg[4] = stressTensorNeg[1, 2];
      stress0Neg[5] = stressTensorNeg[0, 2];
    }
    /// <summary>
    /// Miehe 2010 A phase field model for rate-independent crack propagation: Robust algorithmic implementation based on operator splits
    /// </summary>
    /// <param name="x"></param>
    /// <returns></returns>
    private double PositiveRampFunction(double x)
    {
      //return (Math.Abs(x) + x) / 2.0;
      if (x > 0)
        return x;
      else
        return 0;
    }
    private double NegativeRampFunction(double x)
    {
      //return (Math.Abs(x) - x) / 2.0;
      if (x < 0)
        return x;
      else
        return 0;
    }
    private DoubleVector ComputeDHDEpsilon(double lamda, double mu, DoubleVector principleStrain, DoubleVector[] eigenVectors)
    {
      DoubleVector stress = null;
      switch (ModelStructurePhaseFieldStatic.typeDegradationModel)
      {
        case TypeDegradationEnergy.Isotropic:
          DoubleVector stress0 = null;
          ComputeStressIsotropic(principleStrain, eigenVectors, lamda, mu, out stress0);
          stress = stress0;
          break;
        case TypeDegradationEnergy.Hybrid:
        case TypeDegradationEnergy.Anisotropic:
          DoubleVector stress0Pos = null;
          DoubleVector stress0Neg = null;
          ComputeStressAnisotropic(principleStrain, eigenVectors, lamda, mu, out stress0Pos, out stress0Neg);
          stress = stress0Pos;
          break;
      }

      return stress;
    }

    /// <summary>
    /// The energetic degradation function
    /// </summary>
    private double OmegaD(double d, double l0, double lch, TypeModelPhasefield type)
    {
      double val = 0;
      switch (type)
      {
        case TypeModelPhasefield.AT1:
        case TypeModelPhasefield.AT2:
          val = Math.Pow(1.0 - d, 2);
          break;
        case TypeModelPhasefield.Cubic:
          double s = 1e-4;
          val = s * (Math.Pow(1.0 - d, 3) - Math.Pow(1.0 - d, 2)) + 3 * Math.Pow(1.0 - d, 2) - 2 * Math.Pow(1.0 - d, 3);
          break;
        case TypeModelPhasefield.PFCZM_Linear:
                case TypeModelPhasefield.PFCZM_BLinear:
                case TypeModelPhasefield.PFCZM_Exponential:
        case TypeModelPhasefield.PFCZM_Hyperbolic:
        case TypeModelPhasefield.PFCZM_Cornelissen:
          double a1 = 4.0 * lch / (Math.PI * l0);
          DetermineCoefficientCZM(type, out double p, out double a2, out double a3);
          double fac1 = Math.Pow(1 - d, p);
          double fac2 = fac1 + a1 * d + a1 * a2 * d * d + a1 * a2 * a3 * d * d * d;

          val = fac1 / fac2;
          break;
      }
      return val;
    }
    private void DetermineCoefficientCZM(TypeModelPhasefield type, out double p, out double a2, out double a3)
    {
      switch (type)
      {
        case TypeModelPhasefield.PFCZM_Linear:
          p = 2;
          a2 = -0.5;
          a3 = 0;
          break;
        case TypeModelPhasefield.PFCZM_BLinear:
          p = 2;
          a2 = 0.03687;
          a3 = 20.8343;
          break;
        case TypeModelPhasefield.PFCZM_Exponential:
          p = 2.5;
          a2 = Math.Pow(2.0, 5.0 / 3.0) - 3.0;
          a3 = 0;
          break;
        case TypeModelPhasefield.PFCZM_Hyperbolic:
          p = 4.0;
          a2 = Math.Pow(2.0, 7.0 / 3.0) - 4.5;
          a3 = 0;
          break;
        case TypeModelPhasefield.PFCZM_Cornelissen:
          p = 2.0;
          a2 = 1.3868;
          a3 = 0.6567;//0.9106;
          break;
        default:
          p = 0;
          a2 = 0;
          a3 = 0;
          break;
      }
    }
    private double dOmegaD(double d, double l0, double lch, TypeModelPhasefield type)
    {
      double val = 0;
      switch (type)
      {
        case TypeModelPhasefield.AT1:
        case TypeModelPhasefield.AT2:
          val = -2.0 * (1.0 - d);
          break;
        case TypeModelPhasefield.Cubic:
          double s = 1e-4;
          val = s * (-3 * Math.Pow(1.0 - d, 2) + 2 * (1.0 - d)) - 6 * (1.0 - d) + 6 * Math.Pow(1.0 - d, 2);
          break;
        case TypeModelPhasefield.PFCZM_Linear:
                case TypeModelPhasefield.PFCZM_BLinear:
                case TypeModelPhasefield.PFCZM_Exponential:
        case TypeModelPhasefield.PFCZM_Hyperbolic:
        case TypeModelPhasefield.PFCZM_Cornelissen:
          double a1 = 4.0 * lch / (Math.PI * l0);
          DetermineCoefficientCZM(type, out double p, out double a2, out double a3);
          double fac1 = Math.Pow(1 - d, p);
          double dfac1 = -p * Math.Pow(1 - d, p - 1);

          double fac2 = fac1 + a1 * d + a1 * a2 * d * d + a1 * a2 * a3 * d * d * d;
          double dfac2 = dfac1 + a1 + 2.0 * a1 * a2 * d + 3.0 * a1 * a2 * a3 * d * d;

          val = (dfac1 * fac2 - fac1 * dfac2) / (fac2 * fac2);
          break;
      }
      return val;
    }
    private double ddOmegaD(double d, double l0, double lch, TypeModelPhasefield type)
    {
      double val = 0;
      switch (type)
      {
        case TypeModelPhasefield.AT1:
        case TypeModelPhasefield.AT2:
          val = 2.0;
          break;
        case TypeModelPhasefield.Cubic:
          double s = 1e-4;
          val = s * (6 * (1.0 - d) - 2.0) + 6.0 - 12 * (1.0 - d);
          break;
        case TypeModelPhasefield.PFCZM_Linear:
                case TypeModelPhasefield.PFCZM_BLinear:
                case TypeModelPhasefield.PFCZM_Exponential:
        case TypeModelPhasefield.PFCZM_Hyperbolic:
        case TypeModelPhasefield.PFCZM_Cornelissen:
          double a1 = 4.0 * lch / (Math.PI * l0);
          DetermineCoefficientCZM(type, out double p, out double a2, out double a3);
          double fac1 = Math.Pow(1 - d, p);
          double dfac1 = -p * Math.Pow(1 - d, p - 1);
          double ddfac1 = p * (p - 1) * Math.Pow(1 - d, p - 2);

          double fac2 = fac1 + a1 * d + a1 * a2 * d * d + a1 * a2 * a3 * d * d * d;
          double dfac2 = dfac1 + a1 + 2.0 * a1 * a2 * d + 3.0 * a1 * a2 * a3 * d * d;
          double ddfac2 = ddfac1 + 2.0 * a1 * a2 + 6.0 * a1 * a2 * a3 * d;

          val = ((ddfac1 * fac2 - fac1 * ddfac2) * fac2 - 2 * (dfac1 * fac2 - fac1 * dfac2) * dfac2) / (fac2 * fac2 * fac2);
          break;
      }
      return val;
    }

    private double alphaD(double d, TypeModelPhasefield type)
    {
      double xi = TypeOfAlphaFunction(type);
      return xi * d + (1.0 - xi) * d * d;
    }
    private double dAlphaD(double d, TypeModelPhasefield type)
    {
      double xi = TypeOfAlphaFunction(type);
      return xi + 2.0 * (1.0 - xi) * d;
    }
    private double ddAlphaD(TypeModelPhasefield type)
    {
      double xi = TypeOfAlphaFunction(type);
      return 2.0 * (1.0 - xi);
    }
    private double cAlpha(TypeModelPhasefield type)
    {
      double xi = TypeOfAlphaFunction(type);
      FunctionRToR alpha = new PolynomialFunctionRToR(0, xi, 1 - xi);

      FunctionRToR SqrtAlpha = new PowerFunctionRToR(4.0, 0.5, alpha);
      return SqrtAlpha.Integration(0, 1, 20);
    }
    private static double TypeOfAlphaFunction(TypeModelPhasefield type)
    {
      double xi = 0;
      switch (type)
      {
        case TypeModelPhasefield.AT1:
          xi = 1.0;
          break;
        case TypeModelPhasefield.AT2:
        case TypeModelPhasefield.Cubic:
          xi = 0;
          break;
        case TypeModelPhasefield.PFCZM_Linear:
                case TypeModelPhasefield.PFCZM_BLinear:
                case TypeModelPhasefield.PFCZM_Exponential:
        case TypeModelPhasefield.PFCZM_Hyperbolic:
        case TypeModelPhasefield.PFCZM_Cornelissen:
          xi = 2.0;
          break;
      }
      return xi;
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
      int d = patch.GetCountField();
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
      DoubleVector fiU = null;
      DoubleVector fiP = null;
      DoubleMatrix Bu = new DoubleMatrix(6, 3 * nen);
      DoubleMatrix Bp = new DoubleMatrix(3, nen);

      switch (AbstractModel.TypeSchemeNonlinearSolver)
      {
        case TypeSchemeNonlinearSolver.Monolithic:
        case TypeSchemeNonlinearSolver.SingleStaggered:
          fiU = new DoubleVector(3 * nen);
          fiP = new DoubleVector(nen);
          break;
        case TypeSchemeNonlinearSolver.Staggered:
          if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
          {
            fiU = new DoubleVector(3 * nen);
          }
          else
          {
            fiP = new DoubleVector(nen);
          }
          break;
      }

      double lamda = 0, mu = 0, gc = 0, ft = 0, lch = 0;
      if (!(Material is FGMUserDefinedGradedMaterial))
      {
        lamda = Material.GetProperty(MaterialPropertyName.LameParameter).GetValueProperty();
        mu = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
        gc = Material.GetProperty(MaterialPropertyName.CriticalEnergyReleaseRate).GetValueProperty();
        ft = Material.GetProperty(MaterialPropertyName.TensileYieldStrength).GetValueProperty();
      }
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {
            GaussPoints gpsijk = gps[i, j, k];
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

              //x y z
              Bp[0, kk] = dNdx[kk, 0];
              Bp[1, kk] = dNdx[kk, 1];
              Bp[2, kk] = dNdx[kk, 2];
            }
            if (Material is FGMUserDefinedGradedMaterial)
            {
              lamda = ComputeParameterProperty(MaterialPropertyName.LameParameter, gpsijk.location[0], gpsijk.location[1], gpsijk.location[2]);
              mu = ComputeParameterProperty(MaterialPropertyName.ShearModulus, gpsijk.location[0], gpsijk.location[1], gpsijk.location[2]);
              gc = ComputeParameterProperty(MaterialPropertyName.CriticalEnergyReleaseRate, gpsijk.location[0], gpsijk.location[1], gpsijk.location[2]);
              ft = ComputeParameterProperty(MaterialPropertyName.TensileYieldStrength, gpsijk.location[0], gpsijk.location[1], gpsijk.location[2]);
            }
            double E = mu * (3 * lamda + 2 * mu) / (lamda + mu);
            lch = E * gc / (ft * ft);
            //double lch = ModelStructurePhaseFieldStatic.lch;
            double l0 = ModelStructurePhaseFieldStatic.l0;
            double k1 = ModelStructurePhaseFieldStatic.k;
            double previousPhase;
            double currentPhase;
            DoubleVector stress = null;
            DoubleVector phie = null;
            double xiEpsilon = 0;
            DoubleMatrix F1 = null;
            DoubleVector F3 = null;
            DoubleVector F2 = null;
            DoubleVector ue = null;
            DoubleVector strain = null;
            DoubleVector principleStrain;
            DoubleVector[] eigenVectors;
            double valDOmegaD;
            double valCAlpha;
            double valDAlphaD;

            switch (AbstractModel.TypeSchemeNonlinearSolver)
            {
              case TypeSchemeNonlinearSolver.Monolithic:
              case TypeSchemeNonlinearSolver.SingleStaggered:////////////thay doi Wu degradation
                #region singlestaggered
                phie = GetEachVariableLocal(Result.PHASEFIELD);
                currentPhase = (double)gpsijk.GetValue(DataInGausspoint.currentPhase);
                stress = (DoubleVector)gpsijk.GetValue(DataInGausspoint.currentStress);
                fiU += gpsijk.weight * detJbar * detJ * NMathFunctions.Product(Bu.Transpose(), stress);
                //double xiEpsilonPos, xiEpsilonNeg;
                xiEpsilon = (double)gpsijk.GetValue(DataInGausspoint.xiEpsilon);
                //xiEpsilon = gps[i, j].xiEpsilon;
                valDOmegaD = dOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                valDAlphaD = dAlphaD(currentPhase, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                if (ModelStructurePhaseFieldStatic.typeOrderPhaseField == TypeOrderPhaseField.SecondOrder)
                {
                  // l0 là kích thước của vết nứt
                  F3 = (valDAlphaD * gc / (l0 * valCAlpha) + valDOmegaD * xiEpsilon) * N;
                  F1 = 2.0 * gc * l0 / valCAlpha * NMathFunctions.Product(Bp.Transpose(), Bp);
                }
                else
                {
                  //Borden(2014)
                  F3 = (valDAlphaD * gc / (l0 * valCAlpha) + valDOmegaD * xiEpsilon) * N;
                  DoubleMatrix ddNdX = (DoubleMatrix)gpsijk.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
                  DoubleVector BB = new DoubleVector(nen);
                  for (int kk = 0; kk < nen; kk++)
                  {
                    BB[k] = ddNdX[kk, 0] + ddNdX[kk, 1] + ddNdX[kk, 2];
                  }
                  //Weinberg (2015) == Borden (2014) divide by 2*l0/gc
                  F1 = gc * l0 / valCAlpha * NMathFunctions.Product(Bp.Transpose(), Bp) + gc * Math.Pow(l0, 3) / (8.0 * valCAlpha) * NMathFunctions.OuterProduct(BB, BB);
                }
                F2 = NMathFunctions.Product(F1, phie);
                fiP += gpsijk.weight * detJbar * detJ * (F2 + F3);
                #endregion
                break;
              case TypeSchemeNonlinearSolver.Staggered:
                #region staggered
                if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
                {
                  //displacement
                  stress = (DoubleVector)gpsijk.GetValue(DataInGausspoint.currentStress);
                  fiU += gpsijk.weight * detJbar * detJ * NMathFunctions.Product(Bu.Transpose(), stress);
                }
                else
                {
                  //phasefield             
                  phie = GetEachVariableLocal(Result.PHASEFIELD);
                  currentPhase = (double)gpsijk.GetValue(DataInGausspoint.currentPhase);
                  xiEpsilon = (double)gpsijk.GetValue(DataInGausspoint.xiEpsilon);

                  valDOmegaD = dOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                  valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                  valDAlphaD = dAlphaD(currentPhase, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                  if (ModelStructurePhaseFieldStatic.typeOrderPhaseField == TypeOrderPhaseField.SecondOrder)
                  {
                    //// l0 là kích thước của vết nứt
                    F3 = (valDAlphaD * gc / (l0 * valCAlpha) + valDOmegaD * xiEpsilon) * N;
                    F1 = 2.0 * gc * l0 / valCAlpha * NMathFunctions.Product(Bp.Transpose(), Bp);
                  }
                  else
                  {
                    //Borden(2014)
                    F3 = (valDAlphaD * gc / (l0 * valCAlpha) + valDOmegaD * xiEpsilon) * N;
                    DoubleMatrix ddNdX = (DoubleMatrix)gpsijk.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
                    DoubleVector BB = new DoubleVector(nen);
                    for (int kk = 0; kk < nen; kk++)
                    {
                      BB[kk] = ddNdX[kk, 0] + ddNdX[kk, 1] + ddNdX[kk, 2];
                    }
                    //Weinberg (2015) == Borden (2014) divide by 2*l0/gc
                    F1 = gc * l0 / valCAlpha * NMathFunctions.Product(Bp.Transpose(), Bp) + gc * Math.Pow(l0, 3) / (8.0 * valCAlpha) * NMathFunctions.OuterProduct(BB, BB);
                  }
                  F2 = NMathFunctions.Product(F1, phie);
                  fiP += gpsijk.weight * detJbar * detJ * (F2 + F3);
                }
                #endregion
                break;
            }
          }
        }
      }

      for (int i = 0; i < nen; i++)
      {
        switch (AbstractModel.TypeSchemeNonlinearSolver)
        {
          case TypeSchemeNonlinearSolver.Monolithic:
          case TypeSchemeNonlinearSolver.SingleStaggered:
            for (int ii = 0; ii < 3; ii++)//x y
            {
              fi[i * 4 + ii] = fiU[i * 3 + ii];
            }
            fi[i * 4 + 3] = fiP[i];//phi
            break;
          case TypeSchemeNonlinearSolver.Staggered:
            if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
            {
              for (int ii = 0; ii < 3; ii++)//x y
              {
                fi[i * 4 + ii] = fiU[i * 3 + ii];
              }
            }
            else
            {
              fi[i * 4 + 3] = fiP[i];//phi
            }
            break;
        }
      }
    }

    internal override void UpdateGausspointValue()
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      int d = patch.GetCountField();
      TrivariateNURBSBasisFunction basis = (TrivariateNURBSBasisFunction)patch.GetVolume().Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int nen = (p + 1) * (q + 1) * (r + 1); // number of local basis functions
      DoubleMatrix Bu = new DoubleMatrix(6, 3 * nen);
      double lamda = 0, mu = 0, gc = 0, ft = 0;
      if (!(Material is FGMUserDefinedGradedMaterial))
      {
        lamda = Material.GetProperty(MaterialPropertyName.LameParameter).GetValueProperty();
        mu = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
        gc = Material.GetProperty(MaterialPropertyName.CriticalEnergyReleaseRate).GetValueProperty();
        ft = Material.GetProperty(MaterialPropertyName.TensileYieldStrength).GetValueProperty();
      }
      DoubleMatrix dNdx = null;
      DoubleVector N = null;
      DoubleVector ue = null;////////////
      DoubleVector strain = null;
      double previousPhase;
      double currentPhase;
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {
            GaussPoints gpsijk = gps[i, j, k];
            if (Material is FGMUserDefinedGradedMaterial)
            {
              lamda = ComputeParameterProperty(MaterialPropertyName.LameParameter, gpsijk.location[0], gpsijk.location[1], gpsijk.location[2]);
              mu = ComputeParameterProperty(MaterialPropertyName.ShearModulus, gpsijk.location[0], gpsijk.location[1], gpsijk.location[2]);
              gc = ComputeParameterProperty(MaterialPropertyName.CriticalEnergyReleaseRate, gpsijk.location[0], gpsijk.location[1], gpsijk.location[2]);
              ft = ComputeParameterProperty(MaterialPropertyName.TensileYieldStrength, gpsijk.location[0], gpsijk.location[1], gpsijk.location[2]);
            }
            double E = mu * (3 * lamda + 2 * mu) / (lamda + mu);
            double lch = E * gc / (ft * ft);
            double xiEpsilon = 0;
            double xiEpsilonPos = 0;
            double xiEpsilonNeg = 0;
            DoubleVector principleStrain = null;
            DoubleVector[] eigenVectors = null;
            DoubleVector phie = null;
            double phiCurrent;
            switch (AbstractModel.TypeSchemeNonlinearSolver)
            {
              case TypeSchemeNonlinearSolver.Monolithic:
              case TypeSchemeNonlinearSolver.SingleStaggered:
              case TypeSchemeNonlinearSolver.Staggered:
                ////update after displacement step
                dNdx = (DoubleMatrix)gpsijk.GetValue(DataInGausspoint.dNdX);//gps[i, j].dNdX;
                N = (DoubleVector)gpsijk.GetValue(DataInGausspoint.Ni);//gps[i, j].Ni;
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

                strain = ComputeStrainGauss(patch, Bu);
                gpsijk.SetValue(DataInGausspoint.currentStrain, strain);
                ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);

                switch (ModelStructurePhaseFieldStatic.typeDegradationModel)
                {
                  case TypeDegradationEnergy.Isotropic:
                    xiEpsilon = ComputeXiEpsilonIsotropic(principleStrain, eigenVectors, lamda, mu);
                    break;
                  case TypeDegradationEnergy.Hybrid:
                  case TypeDegradationEnergy.Anisotropic:
                    ComputeXiEpsilonAnisotropic(principleStrain, eigenVectors, lamda, mu, out xiEpsilonPos, out xiEpsilonNeg, ft);
                    xiEpsilon = xiEpsilonPos;
                    break;
                }
                if (xiEpsilon > (double)gpsijk.GetValue(DataInGausspoint.xiEpsilon))
                {
                  gpsijk.SetValue(DataInGausspoint.xiEpsilon, xiEpsilon);
                }
                phie = GetEachVariableLocal(Result.PHASEFIELD);
                phiCurrent = NMathFunctions.Dot(N, phie);
                if (phiCurrent < 0.0)
                  phiCurrent = 0.0;
                if (phiCurrent > 1.0)
                  phiCurrent = 1.0;

                gpsijk.SetValue(DataInGausspoint.currentPhase, phiCurrent);
                currentPhase = phiCurrent;//(double)gpsij.GetValue(DataInGausspoint.currentPhase);
                DoubleVector stress = ComputeStress(lamda, mu, lch, currentPhase, principleStrain, eigenVectors);
                gpsijk.SetValue(DataInGausspoint.currentStress, stress);
                break;
            }
          }
        }
      }
    }
    private DoubleVector ComputeStress(double lamda, double mu, double lch, double phiCurrent, DoubleVector principleStrain, DoubleVector[] eigenVectors)
    {
      double k1 = ModelStructurePhaseFieldStatic.k;
      double l0 = ModelStructurePhaseFieldStatic.l0;

      DoubleVector stress = null;
      double valOmega = OmegaD(phiCurrent, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel) + k1;
      switch (ModelStructurePhaseFieldStatic.typeDegradationModel)
      {
        case TypeDegradationEnergy.Isotropic:
        case TypeDegradationEnergy.Hybrid:
          DoubleVector stress0 = null;
          ComputeStressIsotropic(principleStrain, eigenVectors, lamda, mu, out stress0);
          stress = valOmega * stress0;
          break;
        case TypeDegradationEnergy.Anisotropic:
          DoubleVector stress0Pos = null;
          DoubleVector stress0Neg = null;
          ComputeStressAnisotropic(principleStrain, eigenVectors, lamda, mu, out stress0Pos, out stress0Neg);
          stress = valOmega * stress0Pos + stress0Neg;
          break;
      }

      return stress;
    }
    private double ComputeXiEpsilonIsotropic(DoubleVector principleStrain, DoubleVector[] eigenVectors, double lamda, double mu)
    {
      return 1.0 / 2.0 * lamda * Math.Pow(principleStrain[0] + principleStrain[1] + principleStrain[2], 2) + mu * (Math.Pow(principleStrain[0], 2) + Math.Pow(principleStrain[1], 2) + Math.Pow(principleStrain[2], 2));
    }

    private void ComputeXiEpsilonAnisotropic(DoubleVector principleStrain, DoubleVector[] eigenVectors, double lamda, double mu, out double xiEpsilonPositive, out double xiEpsilonNegative, double ft)
    {
      xiEpsilonPositive = xiEpsilonNegative = 0;
      switch (ModelStructurePhaseFieldStatic.typePhaseFieldModel)
      {
        case TypeModelPhasefield.AT1:
        case TypeModelPhasefield.AT2:
        case TypeModelPhasefield.Cubic:
          double trEpsilon = principleStrain[0] + principleStrain[1] + principleStrain[2];

          DoubleMatrix epPos = PositiveRampFunction(principleStrain[0]) * NMathFunctions.OuterProduct(eigenVectors[0], eigenVectors[0])
            + PositiveRampFunction(principleStrain[1]) * NMathFunctions.OuterProduct(eigenVectors[1], eigenVectors[1])
            + PositiveRampFunction(principleStrain[2]) * NMathFunctions.OuterProduct(eigenVectors[2], eigenVectors[2]);
          DoubleMatrix epNeg = NegativeRampFunction(principleStrain[0]) * NMathFunctions.OuterProduct(eigenVectors[0], eigenVectors[0])
            + NegativeRampFunction(principleStrain[1]) * NMathFunctions.OuterProduct(eigenVectors[1], eigenVectors[1])
            + NegativeRampFunction(principleStrain[2]) * NMathFunctions.OuterProduct(eigenVectors[2], eigenVectors[2]);

          DoubleMatrix epPosSquare = NMathFunctions.Product(epPos, epPos);
          DoubleMatrix epNegSquare = NMathFunctions.Product(epNeg, epNeg);

          xiEpsilonPositive = lamda / 2.0 * Math.Pow(PositiveRampFunction(trEpsilon), 2) + mu * (epPosSquare[0, 0] + epPosSquare[1, 1] + epPosSquare[2, 2]);
          xiEpsilonNegative = lamda / 2.0 * Math.Pow(NegativeRampFunction(trEpsilon), 2) + mu * (epNegSquare[0, 0] + epNegSquare[1, 1] + epNegSquare[2, 2]);
          break;
        case TypeModelPhasefield.PFCZM_Linear:
                case TypeModelPhasefield.PFCZM_BLinear:
                case TypeModelPhasefield.PFCZM_Exponential:
        case TypeModelPhasefield.PFCZM_Hyperbolic:
        case TypeModelPhasefield.PFCZM_Cornelissen:
          double E = mu * (3.0 * lamda + 2.0 * mu) / (lamda + mu);
          DoubleVector stress = ComputeEfficientStress(lamda, mu, principleStrain, eigenVectors);
          double smax = FailureCriterion(ModelStructurePhaseFieldStatic.typeFailureCriterion, stress);
          xiEpsilonPositive = Math.Max(0.5 * Math.Pow(ft, 2.0) / E, 0.5 * Math.Pow(smax, 2.0) / E);
          break;
      }
    }

    private double FailureCriterion(TypeFailureCriterion type, DoubleVector stress)
    {
      double sigmaEq = 0;
      switch (type)
      {
        case TypeFailureCriterion.Rankine:
          ComputePrincipleStress(stress, out DoubleVector eigValue, out DoubleVector[] eigVector);
          sigmaEq = PositiveRampFunction(eigValue.Max());
          break;
        case TypeFailureCriterion.ModifiedVonMises:
          double rhoC = 10.0;
          double I1 = stress[0] + stress[1] + stress[2];
          double J2 = 1.0 / 6.0 * (Math.Pow(stress[0] - stress[1], 2) + Math.Pow(stress[1] - stress[2], 2) + Math.Pow(stress[2] - stress[0], 2)) + Math.Pow(stress[3], 2);
          sigmaEq = (rhoC - 1.0) / (2.0 * rhoC) * I1 + 1.0 / (2.0 * rhoC) * Math.Sqrt(Math.Pow(rhoC - 1, 2) * I1 * I1 + 12.0 * rhoC * J2);
          break;
        case TypeFailureCriterion.Wu2017:
          double beta = 10.0 - 1.0;//fc/ft=10.0
          J2 = 1.0 / 6.0 * (Math.Pow(stress[0] - stress[1], 2) + Math.Pow(stress[1] - stress[2], 2) + Math.Pow(stress[2] - stress[0], 2)) + Math.Pow(stress[3], 2);
          ComputePrincipleStress(stress, out DoubleVector eigValue1, out DoubleVector[] eigVector1);
          double sigmaMax = PositiveRampFunction(eigValue1.Max());
          sigmaEq = 1.0 / (1.0 + beta) * (beta * sigmaMax + Math.Sqrt(3.0 * J2));
          break;
      }
      return sigmaEq;
    }
    private void ComputePrincipleStress(DoubleVector stress, out DoubleVector eigenValues, out DoubleVector[] eigenVectors)
    {
      DoubleMatrix stressTensor = new DoubleMatrix(3, 3);
      stressTensor[0, 0] = stress[0];
      stressTensor[1, 1] = stress[1];
      stressTensor[2, 2] = stress[2];
      stressTensor[0, 1] = stressTensor[1, 0] = stress[3];
      stressTensor[1, 2] = stressTensor[2, 1] = stress[4];
      stressTensor[0, 2] = stressTensor[2, 0] = stress[5];
      var eigDecomp = new DoubleSymEigDecomp(stressTensor);
      eigenValues = eigDecomp.EigenValues;
      eigenVectors = new DoubleVector[3];
      eigenVectors[0] = eigDecomp.EigenVector(0);
      eigenVectors[1] = eigDecomp.EigenVector(1);
      eigenVectors[2] = eigDecomp.EigenVector(2);
    }
    private DoubleVector ComputeEfficientStress(double lamda, double mu, DoubleVector principleStrain, DoubleVector[] eigenVectors)
    {
      DoubleVector stress = null;
      switch (ModelStructurePhaseFieldStatic.typeDegradationModel)
      {
        case TypeDegradationEnergy.Isotropic:
        case TypeDegradationEnergy.Hybrid:
          DoubleVector stress0 = null;
          ComputeStressIsotropic(principleStrain, eigenVectors, lamda, mu, out stress0);
          stress = stress0;
          break;
        case TypeDegradationEnergy.Anisotropic:
          DoubleVector stress0Pos = null;
          DoubleVector stress0Neg = null;
          ComputeStressAnisotropic(principleStrain, eigenVectors, lamda, mu, out stress0Pos, out stress0Neg);
          stress = stress0Pos;
          break;
      }
      return stress;
    }
    private DoubleVector ComputeStrainGauss(PatchStructure3D patch, DoubleMatrix Bu)
    {
      DoubleVector ue = GetDisplacementLocal();////////////
      DoubleVector strain = NMathFunctions.Product(Bu, ue);
      return strain;
    }

    public override double ComputeSharpCrackSurface()
    {
      double Ad = 0;
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      int d = patch.GetCountField();
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
      DoubleMatrix Bp = new DoubleMatrix(3, nen);
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {
            GaussPoints gpsijk = gps[i, j, k];
            double detJbar = 1.0 / 8.0
                 * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                 * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2])
                 * (kvNoMulticiply3[idx3 + 1] - kvNoMulticiply3[idx3]);
            DoubleMatrix dNdx = (DoubleMatrix)gpsijk.GetValue(DataInGausspoint.dNdX);//gps[i, j].dNdX;
            DoubleVector N = (DoubleVector)gpsijk.GetValue(DataInGausspoint.Ni);//gps[i, j].Ni;
            double detJ = Math.Abs((double)gpsijk.GetValue(DataInGausspoint.detJ));//Math.Abs(gps[i, j].detJ);
            for (int kk = 0; kk < nen; kk++)
            {
              //x y z
              Bp[0, kk] = dNdx[kk, 0];
              Bp[1, kk] = dNdx[kk, 1];
              Bp[2, kk] = dNdx[kk, 2];
            }
            double l0 = ModelStructurePhaseFieldStatic.l0;
            double currentPhase = (double)gpsijk.GetValue(DataInGausspoint.currentPhase);
            DoubleVector phie = GetEachVariableLocal(Result.PHASEFIELD);
            double valA = 0;
            if (ModelStructurePhaseFieldStatic.typeOrderPhaseField == TypeOrderPhaseField.SecondOrder)
            {
              //// l0 là kích thước của vết nứt
              double valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
              double valAlphaD = alphaD(currentPhase, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
              DoubleVector Bphi = NMathFunctions.Product(Bp, phie);
              valA = 1.0 / valCAlpha * (1.0 / l0 * valAlphaD + l0 * Math.Pow(Bphi.TwoNormSquared(), 2));
            }
            else
            {

            }
            Ad += gpsijk.weight * detJbar * detJ * valA;
          }
        }
      }
      return Ad;
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
