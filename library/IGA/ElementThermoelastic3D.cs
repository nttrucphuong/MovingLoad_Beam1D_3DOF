using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public class ElementThermoelastic3D : AbstractElementStructure3D
  {
    public ElementThermoelastic3D(AbstractPatch3D patch, int id)
         : base(patch, id)
    {
    }

    ////////////////// Doi Luu
    private double[] convectionHeatTransfer;//The convective heat transfer boundary on 4 edges
    private double[] temperatureHeatTransferConvection;//temperature on egde has convective heat transfer
    public double[] GetHeatTransferConvection()
    {
      return convectionHeatTransfer;
    }

    public void SetHeatTransferConvection(double hValue, double tValue, int indexFace)
    {
      if (convectionHeatTransfer == null)
      {
        convectionHeatTransfer = new double[6];
        temperatureHeatTransferConvection = new double[6];
      }
      convectionHeatTransfer[indexFace] = hValue;
      temperatureHeatTransferConvection[indexFace] = tValue;
    }

    public void ComputeHeatTransferConvectionStiffnessMatrixElement(ref DoubleMatrix Ke)
    {
      for (int iii = 0; iii < 6; iii++)
      {
        if (convectionHeatTransfer[iii] != 0)
        {
          Face face = GetVolume().GetFace(iii);
          var paraEndPatchU = GetParameterTwoEndElement(0);
          var paraEndPatchV = GetParameterTwoEndElement(1);
          var paraEndPatchW = GetParameterTwoEndElement(2);
          double xi = 0;
          double eta = 0;
          double zeta = 0;
          //int numGaussPoint = NotFiniteNumberException;//patch.GetNumberOfGaussPoint();
          int nen = face.CountControlPointOnFace();
          DoubleMatrix keSub = new DoubleMatrix(nen, nen);
          //this.KeHeatTransfer = new DoubleMatrix(Ke.ColumnCount);

          for (int kk = 0; kk < numberOfGaussPointOnEachDirection; kk++)
          {
            double psi1 = GaussPoints.GetPoint(numberOfGaussPointOnEachDirection, kk);
            double w1 = GaussPoints.GetWeight(numberOfGaussPointOnEachDirection, kk);
            for (int kkk = 0; kkk < numberOfGaussPointOnEachDirection; kkk++)
            {
              double psi2 = GaussPoints.GetPoint(numberOfGaussPointOnEachDirection, kkk);
              double w2 = GaussPoints.GetWeight(numberOfGaussPointOnEachDirection, kkk);
              switch (face.GetIndexCoordinate())
              {
                case 0:
                  xi = 0.5 * ((paraEndPatchU[1] - paraEndPatchU[0]) * psi1 + paraEndPatchU[1] + paraEndPatchU[0]);
                  eta = 0.5 * ((paraEndPatchV[1] - paraEndPatchV[0]) * psi2 + paraEndPatchV[1] + paraEndPatchV[0]);
                  zeta = paraEndPatchW[face.GetIndexFrontBack()];
                  break;
                case 1:
                  xi = 0.5 * ((paraEndPatchU[1] - paraEndPatchU[0]) * psi1 + paraEndPatchU[1] + paraEndPatchU[0]);
                  eta = paraEndPatchV[face.GetIndexFrontBack()];
                  zeta = 0.5 * ((paraEndPatchW[1] - paraEndPatchW[0]) * psi2 + paraEndPatchW[1] + paraEndPatchW[0]);
                  break;
                case 2:
                  xi = paraEndPatchU[face.GetIndexFrontBack()];
                  eta = 0.5 * ((paraEndPatchV[1] - paraEndPatchV[0]) * psi1 + paraEndPatchV[1] + paraEndPatchV[0]);
                  zeta = 0.5 * ((paraEndPatchW[1] - paraEndPatchW[0]) * psi2 + paraEndPatchW[1] + paraEndPatchW[0]);
                  break;
              }

              double[,] N = face.GetTrivariateBasisFunctionOnFace(xi, eta, zeta);
              int num = (N.GetLength(0)) * (N.GetLength(1));
              DoubleVector Nij = new DoubleVector(num);
              int dem = 0;
              for (int j = 0; j < N.GetLength(1); j++)
                for (int i = 0; i < N.GetLength(0); i++)
                {
                  Nij[dem] = N[i, j];
                  dem += 1;
                }

              double J2 = 0;
              switch (face.GetIndexCoordinate())
              {
                case 0:
                  J2 = 0.25 * (paraEndPatchU[1] - paraEndPatchU[0]) * (paraEndPatchV[1] - paraEndPatchV[0]);
                  break;
                case 1:
                  J2 = 0.25 * (paraEndPatchU[1] - paraEndPatchU[0]) * (paraEndPatchW[1] - paraEndPatchW[0]);
                  break;
                case 2:
                  J2 = 0.25 * (paraEndPatchV[1] - paraEndPatchV[0]) * (paraEndPatchW[1] - paraEndPatchW[0]);
                  break;
              }

              var elem = face.GetVolumeBeAttached().GetElement();
              DoubleMatrix gradNE = elem.GradBasisFunction(xi, eta, zeta);
              DoubleMatrix J = elem.JacobianAt(gradNE);
              //Jacobian of face mapping
              double normJ = Math.Abs(normVectorJacobian(J));

              DoubleMatrix NijTNij = MatrixFunctions.OuterProduct(Nij, Nij);
              keSub += w1 * w2 * normJ * J2 * NijTNij * convectionHeatTransfer[iii];
            }
          }
          int[] tArray = face.GetTArray();
          int c2 = 0;

          for (int j = 3; j < face.CountDofOnFace(); j += 4)
          {
            int c1 = 0;
            for (int i = j; i < face.CountDofOnFace(); i += 4)
            {
              if (tArray[i] != -1 && tArray[j] != -1)
              {
                int indexCP1 = FindEnumeratePatch(tArray[i]);
                int indexCP2 = FindEnumeratePatch(tArray[j]);
                int m = FindIEN(indexCP1);
                int n = FindIEN(indexCP2);

                Ke[m * 4 + 3, n * 4 + 3] += keSub[c1, c2];
                Ke[n * 4 + 3, m * 4 + 3] += keSub[c1, c2];
              }
              c1++;
            }
            c2++;
          }
        }
      }
    }





    public void ComputeHeatTransferConvectionLoadVectorElement(ref DoubleVector rGlobal)
    {
      //DoubleVector rLocal = null;
      for (int iii = 0; iii < 6; iii++)
      {
        if (convectionHeatTransfer[iii] != 0)
        {
          Face face = GetVolume().GetFace(iii);
          var paraEndPatchU = GetParameterTwoEndElement(0);
          var paraEndPatchV = GetParameterTwoEndElement(1);
          var paraEndPatchW = GetParameterTwoEndElement(2);
          double xi = 0;
          double eta = 0;
          double zeta = 0;
          int numGaussPoint = numberOfGaussPointOnEachDirection;
          DoubleVector rLocal = new DoubleVector(face.CountControlPointOnFace());
          for (int kk = 0; kk < numGaussPoint; kk++)
          {
            double psi1 = GaussPoints.GetPoint(numGaussPoint, kk);
            double w1 = GaussPoints.GetWeight(numGaussPoint, kk);
            for (int kkk = 0; kkk < numGaussPoint; kkk++)
            {
              double psi2 = GaussPoints.GetPoint(numGaussPoint, kkk);
              double w2 = GaussPoints.GetWeight(numGaussPoint, kkk);
              switch (face.GetIndexCoordinate())
              {
                case 0:
                  xi = 0.5 * ((paraEndPatchU[1] - paraEndPatchU[0]) * psi1 + paraEndPatchU[1] + paraEndPatchU[0]);
                  eta = 0.5 * ((paraEndPatchV[1] - paraEndPatchV[0]) * psi2 + paraEndPatchV[1] + paraEndPatchV[0]);
                  zeta = paraEndPatchW[face.GetIndexFrontBack()];
                  break;
                case 1:
                  xi = 0.5 * ((paraEndPatchU[1] - paraEndPatchU[0]) * psi1 + paraEndPatchU[1] + paraEndPatchU[0]);
                  eta = paraEndPatchV[face.GetIndexFrontBack()];
                  zeta = 0.5 * ((paraEndPatchW[1] - paraEndPatchW[0]) * psi2 + paraEndPatchW[1] + paraEndPatchW[0]);
                  break;
                case 2:
                  xi = paraEndPatchU[face.GetIndexFrontBack()];
                  eta = 0.5 * ((paraEndPatchV[1] - paraEndPatchV[0]) * psi1 + paraEndPatchV[1] + paraEndPatchV[0]);
                  zeta = 0.5 * ((paraEndPatchW[1] - paraEndPatchW[0]) * psi2 + paraEndPatchW[1] + paraEndPatchW[0]);
                  break;
              }
              double[,] N = face.GetTrivariateBasisFunctionOnFace(xi, eta, zeta);
              int num = (N.GetLength(0)) * (N.GetLength(1));
              DoubleVector Nij = new DoubleVector(num);
              int dem = 0;
              for (int j = 0; j < N.GetLength(1); j++)
                for (int i = 0; i < N.GetLength(0); i++)
                {
                  Nij[dem] = N[i, j];
                  dem += 1;
                }
              double J2 = 0;
              switch (face.GetIndexCoordinate())
              {
                case 0:
                  J2 = 0.25 * (paraEndPatchU[1] - paraEndPatchU[0]) * (paraEndPatchV[1] - paraEndPatchV[0]);
                  break;
                case 1:
                  J2 = 0.25 * (paraEndPatchU[1] - paraEndPatchU[0]) * (paraEndPatchW[1] - paraEndPatchW[0]);
                  break;
                case 2:
                  J2 = 0.25 * (paraEndPatchV[1] - paraEndPatchV[0]) * (paraEndPatchW[1] - paraEndPatchW[0]);
                  break;
              }

              var elem = face.GetVolumeBeAttached().GetElement();
              DoubleMatrix gradNE = elem.GradBasisFunction(xi, eta, zeta);
              DoubleMatrix J = elem.JacobianAt(gradNE);
              //Jacobian of face mapping
              double normJ = Math.Abs(normVectorJacobian(J));

              rLocal += w1 * w2 * normJ * J2 * Nij * convectionHeatTransfer[iii] * temperatureHeatTransferConvection[iii];
            }
          }
          int c = 0;
          int[] tArrayGlobal = face.GetTArrayGlobal();
          for (int i = 3; i < face.CountDofOnFace(); i += 4)
          {
            int m = tArrayGlobal[i];
            rGlobal[m] += rLocal[c++];
          }
        }
      }
    }





    private double normVectorJacobian(DoubleMatrix J)
    {
      double ee = 0;
      double gg = 0;
      double ff = 0;
      for (int iii = 0; iii < 6; iii++)
      {
        if (convectionHeatTransfer[iii] != 0)
        {
          Face face = GetVolume().GetFace(iii);
          switch (face.GetIndexCoordinate())
          {
            case 0:
              ee = J[0, 0] * J[0, 0] + J[1, 0] * J[1, 0] + J[2, 0] * J[2, 0];
              gg = J[0, 1] * J[0, 1] + J[1, 1] * J[1, 1] + J[2, 1] * J[2, 1];
              ff = J[0, 0] * J[0, 1] + J[1, 0] * J[1, 1] + J[2, 0] * J[2, 1];
              break;
            case 1:
              ee = J[0, 0] * J[0, 0] + J[1, 0] * J[1, 0] + J[2, 0] * J[2, 0];
              gg = J[0, 2] * J[0, 2] + J[1, 2] * J[1, 2] + J[2, 2] * J[2, 2];
              ff = J[0, 0] * J[0, 2] + J[1, 0] * J[1, 2] + J[2, 0] * J[2, 2];
              break;
            case 2:
              ee = J[0, 1] * J[0, 1] + J[1, 1] * J[1, 1] + J[2, 1] * J[2, 1];
              gg = J[0, 2] * J[0, 2] + J[1, 2] * J[1, 2] + J[2, 2] * J[2, 2];
              ff = J[0, 1] * J[0, 2] + J[1, 1] * J[1, 2] + J[2, 1] * J[2, 2];
              break;
          }
        }
      }

      return Math.Sqrt(ee * gg - ff * ff);
    }


    public int FindIEN(int numberOfGlobalBasis)
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      for (int i = 0; i < patch.GetCountLocalBasisFunctions(); i++)
      {
        if (patch.GetIEN(id, i) == numberOfGlobalBasis)
          return i;
      }
      return -1;
    }

    public int FindEnumeratePatch(int tArray)
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      for (int i = 0; i < patch.GetCountGlobalBasisFunctions(); i++)
      {
        if (patch.GetEnumerateInPatch(0, 3, i) == tArray)
          return i;
      }
      return -1;
    }

    ///////////////////


    private DoubleMatrix CreateMaterialMatrix()
    {
      DoubleMatrix D = new DoubleMatrix(6, 6);
      D[0, 0] = coefficient(1, 1, 1, 1);
      D[0, 1] = D[1, 0] = coefficient(1, 1, 2, 2);
      D[0, 2] = D[2, 0] = coefficient(1, 1, 3, 3);
      D[0, 3] = D[3, 0] = coefficient(1, 1, 2, 3);
      D[0, 4] = D[4, 0] = coefficient(1, 1, 1, 3);
      D[0, 5] = D[5, 0] = coefficient(1, 1, 1, 2);
      D[1, 1] = coefficient(2, 2, 2, 2);
      D[1, 2] = D[2, 1] = coefficient(2, 2, 3, 3);
      D[1, 3] = D[3, 1] = coefficient(2, 2, 2, 3);
      D[1, 4] = D[4, 1] = coefficient(2, 2, 1, 3);
      D[1, 5] = D[5, 1] = coefficient(2, 2, 1, 2);
      D[2, 2] = coefficient(3, 3, 3, 3);
      D[2, 3] = D[3, 2] = coefficient(3, 3, 2, 3);
      D[2, 4] = D[4, 2] = coefficient(3, 3, 1, 3);
      D[2, 5] = D[5, 2] = coefficient(3, 3, 1, 2);
      D[3, 3] = coefficient(2, 3, 3, 2);
      D[3, 4] = D[4, 3] = coefficient(2, 3, 1, 3);
      D[3, 5] = D[5, 3] = coefficient(2, 3, 1, 2);
      D[4, 4] = coefficient(1, 3, 3, 1);
      D[4, 5] = D[5, 4] = coefficient(1, 3, 1, 2);
      D[5, 5] = coefficient(1, 2, 2, 1);
      return D;
    }

    private double coefficient(int i, int j, int k, int l)
    {
      var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
      var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
      double G = eModulus / (2.0 * (1.0 + nu));
      double lamda = eModulus * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
      return lamda * KroneckerDelta(i, j) * KroneckerDelta(k, l) + G * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
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
    private DoubleMatrix CreateMaterialMatrix(int i, int j, int k)
    {
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

    private double coefficient(int i, int j, int k, int l, int ii, int jj, int kk)
    {
      double eModulus = (double)gps[ii, jj, kk].GetValue(DataInGausspoint.EModulus);
      double nu = (double)gps[ii, jj, kk].GetValue(DataInGausspoint.nu);
      double G = eModulus / (2.0 * (1.0 + nu));
      double lamda = eModulus * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
      return lamda * KroneckerDelta(i, j) * KroneckerDelta(k, l) + G * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
    }

    private DoubleMatrix CreateMaterialMatrix(double xi, double eta, double zeta)
    {
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
      double[] point = null;
      int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
      if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
      {
        point = PointAt(xi, eta, zeta);
      }
      else
      {
        point = new double[] { xi, eta, zeta };
      }
      var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty(point[direction]);
      var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty(point[direction]);
      double G = eModulus / (2.0 * (1.0 + nu));
      double lamda = eModulus * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
      return lamda * KroneckerDelta(i, j) * KroneckerDelta(k, l) + G * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
    }



    private DoubleMatrix CreateMaterialThermalMatrix()
    {
      var thermalConductivity = Material.GetProperty(MaterialPropertyName.IsotropicThermalConductivity).GetValueProperty();

      DoubleMatrix C = new DoubleMatrix(3, 3);
      C[0, 0] = C[1, 1] = C[2, 2] = thermalConductivity;
      return C;
    }

    private DoubleMatrix CreateMaterialThermalMatrix(int i, int j, int k)
    {
      var thermalConductivity = (double)gps[i, j, k].GetValue(DataInGausspoint.ThermalConductivity);

      DoubleMatrix C = new DoubleMatrix(3, 3);
      C[0, 0] = C[1, 1] = C[2, 2] = thermalConductivity;
      return C;
    }

    private DoubleMatrix CreateMaterialThermalMatrix(double xi, double eta, double zeta)
    {
      double[] point = null;
      int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
      if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
      {
        point = PointAt(xi, eta, zeta);
      }
      else
      {
        point = new double[] { xi, eta, zeta };
      }
      var thermalConductivity = Material.GetProperty(MaterialPropertyName.FGMThermalConductivity).GetValueProperty(point[direction]);

      // sua lai
      DoubleMatrix C = new DoubleMatrix(3, 3);
      C[0, 0] = C[1, 1] = C[2, 2] = thermalConductivity;
      return C;
    }


    private DoubleVector CreateVectorCoefficientsExpansion()
    {
      var coefficientsthermalexpansion = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion).GetValueProperty();
      DoubleVector A = new DoubleVector(6);
      A[0] = A[1] = A[2] = coefficientsthermalexpansion;
      return A;
    }

    private DoubleVector CreateVectorCoefficientsExpansion(int i, int j, int k)
    {
      double coefficientsthermalexpansion = (double)gps[i, j, k].GetValue(DataInGausspoint.CoefficientsThermalExpansion);
      DoubleVector A = new DoubleVector(6);
      A[0] = A[1] = A[2] = coefficientsthermalexpansion;
      return A;
    }

    private DoubleVector CreateVectorCoefficientsExpansion(double xi, double eta, double zeta)
    {
      double[] point = null;
      int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
      if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
      {
        point = PointAt(xi, eta, zeta);
      }
      else
      {
        point = new double[] { xi, eta, zeta };
      }
      var coefficientsthermalexpansion = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion).GetValueProperty(point[direction]);
      DoubleVector A = new DoubleVector(6);
      A[0] = A[1] = A[2] = coefficientsthermalexpansion;
      return A;
    }
    private DoubleVector CrateVectorThermoelasticCoefficent()
    {
      DoubleMatrix D = CreateMaterialMatrix();
      DoubleVector Anpha = CreateVectorCoefficientsExpansion();
      DoubleVector Beta = MatrixFunctions.Product(D, Anpha);
      return Beta;
    }

    private DoubleVector CrateVectorThermoelasticCoefficent(int i, int j, int k)
    {
      DoubleMatrix D = CreateMaterialMatrix(i, j, k);
      DoubleVector Anpha = CreateVectorCoefficientsExpansion(i, j, k);
      DoubleVector Beta = MatrixFunctions.Product(D, Anpha);
      return Beta;
    }

    private DoubleVector CrateVectorThermoelasticCoefficent(double xi, double eta, double zeta)
    {

      DoubleMatrix D = CreateMaterialMatrix(xi, eta, zeta);
      DoubleVector Anpha = CreateVectorCoefficientsExpansion(xi, eta, zeta);
      DoubleVector Beta = MatrixFunctions.Product(D, Anpha);
      return Beta;
    }


    public override void ComputeStiffnessMatrixElement(out DoubleMatrix Ke)
    {
      PatchThermoelastic3D patch = (PatchThermoelastic3D)this.patch;
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

      DoubleMatrix Keuu = new DoubleMatrix(3 * nen, 3 * nen);
      DoubleMatrix Kett = new DoubleMatrix(nen, nen);
      DoubleMatrix Keut = new DoubleMatrix(3 * nen, nen);
      DoubleMatrix Bu = new DoubleMatrix(6, 3 * nen);
      DoubleMatrix Bt = new DoubleMatrix(3, nen);
      DoubleMatrix Du = null;
      DoubleMatrix Dt = null;
      DoubleVector Beta = null;
      if (!(Material is FGMStructureOneVariableMaterial))
      {
        Du = CreateMaterialMatrix();
        Dt = CreateMaterialThermalMatrix();
        Beta = CrateVectorThermoelasticCoefficent();
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
            DoubleVector Ni = (DoubleVector)gps[i, j, k].GetValue(DataInGausspoint.Ni);//gps[i, j, k].Ni;
            DoubleMatrix dNdx = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdX);
            for (int kk = 0; kk < nen; kk++)
            {
              Bu[0, 3 * kk] = dNdx[kk, 0];
              Bu[1, 3 * kk + 1] = dNdx[kk, 1];
              Bu[2, 3 * kk + 2] = dNdx[kk, 2];
              Bu[3, 3 * kk + 1] = dNdx[kk, 2];
              Bu[3, 3 * kk + 2] = dNdx[kk, 1];
              Bu[4, 3 * kk] = dNdx[kk, 2];
              Bu[4, 3 * kk + 2] = dNdx[kk, 0];
              Bu[5, 3 * kk] = dNdx[kk, 1];
              Bu[5, 3 * kk + 1] = dNdx[kk, 0];

              Bt[0, kk] = dNdx[kk, 0];
              Bt[1, kk] = dNdx[kk, 1];
              Bt[2, kk] = dNdx[kk, 2];
            }
            if (Material is FGMStructureOneVariableMaterial)
            {
              Du = CreateMaterialMatrix(i, j, k);
              Dt = CreateMaterialThermalMatrix(i, j, k);
              Beta = CrateVectorThermoelasticCoefficent(i, j, k);
            }
            DoubleMatrix BuTDuBu = NMathFunctions.TransposeProduct(Bu, MatrixFunctions.Product(Du, Bu));
            DoubleMatrix BtTDtBt = NMathFunctions.TransposeProduct(Bt, MatrixFunctions.Product(Dt, Bt));
            DoubleMatrix BuTBetaNT = Bu.Transpose() * MatrixFunctions.OuterProduct(Beta, Ni);
            double detJ = Math.Abs((double)gps[i, j, k].GetValue(DataInGausspoint.detJ));
            Keuu += gps[i, j, k].weight * detJbar * detJ * BuTDuBu;
            Kett += gps[i, j, k].weight * detJbar * detJ * BtTDtBt;
            Keut -= gps[i, j, k].weight * detJbar * detJ * BuTBetaNT;
          }
        }
      }
      for (int j = 0; j < nen; j++)
        for (int i = j; i < nen; i++)
        {
          for (int jj = 0; jj < 3; jj++)
          {
            for (int ii = 0; ii < 3; ii++)
            {
              Ke[i * 4 + ii, j * 4 + jj] = Keuu[i * 3 + ii, j * 3 + jj];
            }
            Ke[j * 4 + jj, i * 4 + 3] = Ke[i * 4 + 3, j * 4 + jj] = Keut[j * 3 + jj, i];
          }
          Ke[i * 4 + 3, j * 4 + 3] = Kett[i, j];
        }
    }

    public override DoubleVector StrainAt(params double[] xi)
    {
      PatchThermoelastic3D patch = (PatchThermoelastic3D)this.patch;
      int d = GetCountDimension();
      int nen = patch.GetCountLocalBasisFunctions();
      var cps = patch.GetVolume().ControlPoints;
      DoubleMatrix Bu = new DoubleMatrix(6, d * nen);


      DoubleMatrix gradNE = GradBasisFunction(xi[0], xi[1], xi[2]);
      DoubleMatrix J = JacobianAt(gradNE);
      DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
      DoubleMatrix gradBasis = gradNE * invertJ;
      for (int kk = 0; kk < nen; kk++)
      {
        Bu[0, 3 * kk] = gradBasis[kk, 0];
        Bu[1, 3 * kk + 1] = gradBasis[kk, 1];
        Bu[2, 3 * kk + 2] = gradBasis[kk, 2];
        Bu[3, 3 * kk + 1] = gradBasis[kk, 2];
        Bu[3, 3 * kk + 2] = gradBasis[kk, 1];
        Bu[4, 3 * kk] = gradBasis[kk, 2];
        Bu[4, 3 * kk + 2] = gradBasis[kk, 0];
        Bu[5, 3 * kk] = gradBasis[kk, 1];
        Bu[5, 3 * kk + 1] = gradBasis[kk, 0];
      }
      DoubleVector U = GetDisplacementLocal();
      return MatrixFunctions.Product(Bu, U);
    }

    public override DoubleVector StressAt(params double[] xi)
    {
      DoubleMatrix Du = null;
      if (!(Material is FGMStructureOneVariableMaterial))
      {
        Du = CreateMaterialMatrix();
      }
      else
      {
        Du = CreateMaterialMatrix(xi[0], xi[1], xi[2]);
      }
      return MatrixFunctions.Product(Du, (StrainAt(xi) - StrainThermoAt(xi)));
    }

    public DoubleVector GetDisplacementLocal()
    {
      PatchThermoelastic3D patch = (PatchThermoelastic3D)this.patch;
      var cps = patch.GetVolume().ControlPoints;
      int d = GetCountDimension();
      int nen = patch.GetCountLocalBasisFunctions();
      DoubleVector U = new DoubleVector(d * nen);
      for (int i = 0; i < nen; i++)
      {
        var ien = patch.GetIEN(id, i);
        U[d * i] = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(0, ien, 2)].GetResult(Result.UX);
        U[d * i + 1] = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(0, ien, 2)].GetResult(Result.UY);
        U[d * i + 2] = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(0, ien, 2)].GetResult(Result.UZ);
      }
      return U;
    }

    public override void ComputeInternalForceElement(out DoubleVector fi)
    {
      throw new NotImplementedException();
    }

    public override DoubleVector StrainElasticAt(params double[] xi)
    {
      return StrainAt(xi) - StrainThermoAt(xi);
    }

    public override DoubleVector StrainThermoAt(params double[] xi)
    {
      PatchThermoelastic3D patch = (PatchThermoelastic3D)this.patch;
      double Tref = Material.GetProperty(MaterialPropertyName.ReferenceTemperature).GetValueProperty();
      double T = patch.GetApproximateAt(Result.TEMP, xi[0], xi[1], xi[2]);
      DoubleVector Anpha = null;
      if (!(Material is FGMStructureOneVariableMaterial))
        Anpha = CreateVectorCoefficientsExpansion();
      else
        Anpha = CreateVectorCoefficientsExpansion(xi[0], xi[1], xi[2]);
      return Anpha * (T - Tref);
    }
  }
}
