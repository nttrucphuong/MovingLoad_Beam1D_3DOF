using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public class ElementThermal3D : AbstractElement3D
  {
    private double[] convectionHeatTransfer;
    private double[] temperatureHeatTransferConvection;

    public ElementThermal3D(AbstractPatch3D mesh, int id)
      : base(mesh, id)
    {
    }
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
          DoubleSymmetricMatrix keSub = new DoubleSymmetricMatrix(face.CountDofOnFace());
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

              DoubleSymmetricMatrix NijTNij = new DoubleSymmetricMatrix(MatrixFunctions.OuterProduct(Nij, Nij));
              keSub += w1 * w2 * normJ * J2 * NijTNij * convectionHeatTransfer[iii];
            }
          }
          int[] tArray = face.GetTArray();
          for (int j = 0; j < face.CountDofOnFace(); j++)
            for (int i = j; i < face.CountDofOnFace(); i++)
            {
              int m = FindIEN(tArray[i]);
              int n = FindIEN(tArray[j]);
              Ke[m, n] += keSub[i, j];
              Ke[n, m] += keSub[i, j];
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
          DoubleVector rLocal = new DoubleVector(face.CountDofOnFace());
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
          int[] tArrayGlobal = face.GetTArrayGlobal();
          for (int i = 0; i < face.CountDofOnFace(); i++)
          {
            int m = tArrayGlobal[i];
            rGlobal[m] += rLocal[i];
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

    public void ComputeHeatGenerationSourceLoadVectorElement(ref DoubleVector rGlobal, double G)
    {
      PatchThermal3D patch = (PatchThermal3D)this.patch;
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
      DoubleVector rLocal = new DoubleVector(d * nen);
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

            DoubleVector Ni = (DoubleVector)gps[i, j, k].GetValue(DataInGausspoint.Ni);
            rLocal += gps[i, j, k].weight * detJbar * Math.Abs((double)gps[i, j, k].GetValue(DataInGausspoint.detJ)) * G * Ni;
          }
        }
      }
      int[] tArrayGlobal = GetVolume().GetTArrayGlobal();
      for (int i = 0; i < GetVolume().CountDof(); i++)
      {
        int m = tArrayGlobal[i];
        rGlobal[m] += rLocal[i];
      }
    }

    public double[] GetTLocal()
    {
      PatchThermal3D patch = (PatchThermal3D)this.patch;
      var cps = patch.GetVolume().ControlPoints;
      int d = patch.GetCountField();
      int nen = patch.GetCountLocalBasisFunctions();
      double[] U = new double[d * nen];
      for (int i = 0; i < nen; i++)
      {
        var ien = patch.GetIEN(id, i);
        U[d * i] = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(ien, 2)].GetULocal(0);
      }
      return U;
    }


    private DoubleMatrix CreateMaterialMatrix()
    {
      var thermalConductivity = Material.GetProperty(MaterialPropertyName.IsotropicThermalConductivity).GetValueProperty();

      DoubleMatrix C = new DoubleMatrix(3, 3);
      C[0, 0] = C[1, 1] = C[2, 2] = thermalConductivity;
      return C;
    }

    private DoubleMatrix CreateMaterialMatrix(int i, int j, int k)
    {
      double thermalConductivity = (double)gps[i, j, k].GetValue(DataInGausspoint.ThermalConductivity);

      DoubleMatrix C = new DoubleMatrix(3, 3);
      C[0, 0] = C[1, 1] = C[2, 2] = thermalConductivity;
      return C;
    }

    private DoubleMatrix CreateMaterialMatrix(double xi, double eta, double zeta)
    {
      double[] point = null;
      int direction = ((FGMThermalOneGradedDirectionMaterial)Material).GetDirectionGraded();
      if (((FGMThermalOneGradedDirectionMaterial)Material).GetIsGlobal())
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


    public override void ComputeStiffnessMatrixElement(out DoubleMatrix Ke)
    {
      PatchThermal3D patch = (PatchThermal3D)this.patch;
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
      DoubleMatrix B = new DoubleMatrix(GetCountDimension(), d * nen);
      DoubleMatrix D = null;
      if (!(Material is FGMThermalOneGradedDirectionMaterial))
      {
        D = CreateMaterialMatrix();
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
              * (kvNoMulticiply3[idx3 + 1] - kvNoMulticiply2[idx3]);

            DoubleMatrix dNdx = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdX);// gps[i, j, k].dNdX;//dNdxi * invertJ;
            for (int l = 0; l < nen; l++)
            {
              B[0, d * l] = dNdx[l, 0];
              B[1, d * l] = dNdx[l, 1];
              B[2, d * l] = dNdx[l, 2];
            }
            if (Material is FGMThermalOneGradedDirectionMaterial)
              D = CreateMaterialMatrix(i, j, k);
            DoubleMatrix BTDB = NMathFunctions.TransposeProduct(B, MatrixFunctions.Product(D, B));
            Ke += gps[i, j, k].weight * detJbar * Math.Abs((double)gps[i, j, k].GetValue(DataInGausspoint.detJ)) * BTDB;
          }
        }
      }
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

    public override void ComputeMassMatrixElement(out DoubleMatrix Me)
    {
      PatchThermal3D patch = (PatchThermal3D)this.patch;
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
      Me = new DoubleMatrix(nen, nen);
      double rho = 0;
      double cc = 0;
      if (!(Material is FGMThermalOneGradedDirectionMaterial))
      {
        rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty();
        cc = Material.GetProperty(MaterialPropertyName.SpecificHeatCapacity).GetValueProperty();
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
            DoubleVector Ni = (DoubleVector)gps[i, j, k].GetValue(DataInGausspoint.Ni);
            //DoubleMatrix Ni = DoubleMatrix(nen, 1);
            //int count = 0;
            //for (int kk = 0; kk < nen; kk++)
            //{
            //    Ni[count, 0] = Nij[kk];
            //    count++;
            //}
            if (Material is FGMThermalOneGradedDirectionMaterial)
            {
              double[] point = null;
              int direction = ((FGMThermalOneGradedDirectionMaterial)Material).GetDirectionGraded();
              if (((FGMThermalOneGradedDirectionMaterial)Material).GetIsGlobal())
              {
                point = PointAt(gps[i, j, k].location[0], gps[i, j, k].location[1], gps[i, j, k].location[2]);
              }
              else
              {
                point = new double[] { gps[i, j, k].location[0], gps[i, j, k].location[1] };
              }
              rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty(point[direction]);
              cc = Material.GetProperty(MaterialPropertyName.SpecificHeatCapacity).GetValueProperty(point[direction]);
            }
            DoubleMatrix NeTNe = rho * cc * MatrixFunctions.OuterProduct(Ni, Ni);

            Me += gps[i, j, k].weight * detJbar * Math.Abs((double)gps[i, j, k].GetValue(DataInGausspoint.detJ)) * NeTNe;
          }
        }
      }
    }

    public override void ComputeInternalForceElement(out DoubleVector fi)
    {
      throw new NotImplementedException();
    }

    public override void ComputeGeometricStiffnessMatrixElement(double[,] stressTensor, out DoubleMatrix Ke)
    {
      throw new NotImplementedException();
    }

    public override void ComputeValueAtGaussPoint(DataInGausspoint name)
    {
      throw new NotImplementedException();
    }
  }
}
