using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public class ElementThermoelastic2D : AbstractElementStructure2D
  {
    public ElementThermoelastic2D(AbstractPatch2D patch, int id)
          : base(patch, id)
    {
    }

    /////////////////////////// Doi Luu
    private double[] convectionHeatTransfer;//The convective heat transfer boundary on 4 edges
    private double[] temperatureHeatTransferConvection;//temperature on egde has convective heat transfer
    public double[] GetHeatTransferConvection()
    {
      return convectionHeatTransfer;
    }

    public void SetHeatTransferConvection(double hValue, double tValue, int indexEdge)
    {
      if (convectionHeatTransfer == null)
      {
        convectionHeatTransfer = new double[4];
        temperatureHeatTransferConvection = new double[4];
      }
      convectionHeatTransfer[indexEdge] = hValue;
      temperatureHeatTransferConvection[indexEdge] = tValue;
    }

    public void ComputeHeatTransferConvectionStiffnessMatrixElement(ref DoubleMatrix Ke)
    {
      for (int iii = 0; iii < 4; iii++)
      {
        if (convectionHeatTransfer[iii] != 0)
        {
          Edge e = GetFace().GetEdge(iii);
          var paraEndPatchU = GetParameterTwoEndElement(0);
          var paraEndPatchV = GetParameterTwoEndElement(1);
          double xi = 0;
          double eta = 0;
          //int numGaussPoint = NotFiniteNumberException;//patch.GetNumberOfGaussPoint();
          int nen = e.CountControlPointOnEdge();
          DoubleMatrix keSub = new DoubleMatrix(nen, nen);
          //this.KeHeatTransfer = new DoubleMatrix(Ke.ColumnCount);
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {
            double psi = GaussPoints.GetPoint(numberOfGaussPointOnEachDirection, k);
            double w = GaussPoints.GetWeight(numberOfGaussPointOnEachDirection, k);
            switch (e.GetIndexCoordinate())
            {
              case 0:
                xi = 0.5 * ((paraEndPatchU[1] - paraEndPatchU[0]) * psi + paraEndPatchU[1] + paraEndPatchU[0]);
                eta = paraEndPatchV[e.GetIndexFrontBack()];
                break;
              case 1:
                xi = paraEndPatchU[e.GetIndexFrontBack()];
                eta = 0.5 * ((paraEndPatchV[1] - paraEndPatchV[0]) * psi + paraEndPatchV[1] + paraEndPatchV[0]);
                break;
            }



            DoubleVector Ni = new DoubleVector(e.GetBivariateBasisFunctionOnEdge(0, xi, eta));
            double J2 = (e.GetIndexCoordinate() == 0) ? (0.5 * (paraEndPatchU[1] - paraEndPatchU[0]))
                  : (0.5 * (paraEndPatchV[1] - paraEndPatchV[0]));

            DoubleVector dNi = new DoubleVector(e.GetDerivativeBivariateBasisFunctionOnEdge(0, xi, eta));
            double[] J = JacobianOnEdge(dNi, iii).ToArray();
            double J1 = Math.Sqrt(J[0] * J[0] + J[1] * J[1]);

            DoubleMatrix NiTNi = MatrixFunctions.OuterProduct(Ni, Ni);
            keSub += w * J1 * J2 * NiTNi * convectionHeatTransfer[iii];
          }

          int[] tArray = e.GetTArray();
          int c2 = 0;

          for (int j = 2; j < e.CountDofOnEdge(); j += 3)
          {
            int c1 = 0;
            for (int i = j; i < e.CountDofOnEdge(); i += 3)
            {
              if (tArray[i] != -1 && tArray[j] != -1)
              {
                int indexCP1 = FindEnumeratePatch(tArray[i]);
                int indexCP2 = FindEnumeratePatch(tArray[j]);
                int m = FindIEN(indexCP1);
                int n = FindIEN(indexCP2);

                Ke[m * 3 + 2, n * 3 + 2] += keSub[c1, c2];
                Ke[n * 3 + 2, m * 3 + 2] += keSub[c1, c2];
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
      for (int iii = 0; iii < 4; iii++)
      {
        if (convectionHeatTransfer[iii] != 0)
        {
          Edge e = GetFace().GetEdge(iii);
          var paraEndPatchU = GetParameterTwoEndElement(0);
          var paraEndPatchV = GetParameterTwoEndElement(1);
          double xi = 0;
          double eta = 0;
          int numGaussPoint = numberOfGaussPointOnEachDirection;
          DoubleVector rLocal = new DoubleVector(e.CountControlPointOnEdge());
          for (int k = 0; k < numGaussPoint; k++)
          {
            double psi = GaussPoints.GetPoint(numGaussPoint, k);
            double w = GaussPoints.GetWeight(numGaussPoint, k);
            switch (e.GetIndexCoordinate())
            {
              case 0:
                xi = 0.5 * ((paraEndPatchU[1] - paraEndPatchU[0]) * psi + paraEndPatchU[1] + paraEndPatchU[0]);
                eta = paraEndPatchV[e.GetIndexFrontBack()];
                break;
              case 1:
                xi = paraEndPatchU[e.GetIndexFrontBack()];
                eta = 0.5 * ((paraEndPatchV[1] - paraEndPatchV[0]) * psi + paraEndPatchV[1] + paraEndPatchV[0]);
                break;
            }

            DoubleVector Ni = new DoubleVector(e.GetBivariateBasisFunctionOnEdge(0, xi, eta));
            double J2 = (e.GetIndexCoordinate() == 0) ? (0.5 * (paraEndPatchU[1] - paraEndPatchU[0]))
                  : (0.5 * (paraEndPatchV[1] - paraEndPatchV[0]));

            DoubleVector dNi = new DoubleVector(e.GetDerivativeBivariateBasisFunctionOnEdge(0, xi, eta));
            double[] J = JacobianOnEdge(dNi, iii).ToArray();
            double J1 = Math.Sqrt(J[0] * J[0] + J[1] * J[1]);

            rLocal += w * J1 * J2 * Ni * convectionHeatTransfer[iii] * temperatureHeatTransferConvection[iii];
          }
          int c = 0;
          int[] tArrayGlobal = e.GetTArrayGlobal();
          for (int i = 2; i < e.CountDofOnEdge(); i += 3)
          {
            int m = tArrayGlobal[i];
            rGlobal[m] += rLocal[c++];
          }
        }
      }
    }

    private DoubleVector JacobianOnEdge(DoubleVector dNi, int indexEdge)
    {
      Edge e = GetFace().GetEdge(indexEdge);
      int d = e.CountDimension();
      int p = e.GetDegree(0);
      int n = p + 1;
      DoubleVector grad = new DoubleVector(d);
      var cpsOnEdge = e.GetControlPointsOnEdge(0);
      for (int j = 0; j < d; j++)
      {
        for (int k = 0; k < n; k++)
        {
          grad[j] += dNi[k] * cpsOnEdge[k].GetCoordinate(j);
        }
      }
      return grad;
    }

    public int FindIEN(int numberOfGlobalBasis)
    {
      AbstractPatch2D patch = (AbstractPatch2D)this.patch;
      for (int i = 0; i < patch.GetCountLocalBasisFunctions(); i++)
      {
        if (patch.GetIEN(id, i) == numberOfGlobalBasis)
          return i;
      }
      return -1;
    }
    public int FindEnumeratePatch(int tArray)
    {
      AbstractPatch2D patch = (AbstractPatch2D)this.patch;
      for (int i = 0; i < patch.GetCountGlobalBasisFunctions(); i++)
      {
        if (patch.GetEnumerateInPatch(0, 2, i) == tArray)
          return i;
      }
      return -1;
    }
    ///////////////////////////////


    public override void ComputeMassMatrixElement(out DoubleMatrix Me)
    {
      PatchThermoelastic2D patch = (PatchThermoelastic2D)this.patch;
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
      double rho = 0;
      if (!(Material is FGMStructureOneVariableMaterial))
        rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty();
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          double detJbar = 1.0 / 4.0
                * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          //DoubleMatrix dNdxi = gps[i, j].dNdxi;
          //DoubleMatrix J = JacobianAt(dNdxi);
          DoubleVector Nij = (DoubleVector)gps[i, j].GetValue(DataInGausspoint.Ni);
          DoubleMatrix Ni = new DoubleMatrix(d * nen, d);
          int count = 0;
          for (int kk = 0; kk < nen; kk++)
          {
            Ni[d * count, 0] = Nij[kk];
            Ni[d * count + 1, 1] = Nij[kk];
            count++;
          }
          if (Material is FGMStructureOneVariableMaterial)
          {
            double[] point = null;
            int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
            if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
            {
              point = PointAt(gps[i, j].location[0], gps[i, j].location[1]);
            }
            else
            {
              point = new double[] { gps[i, j].location[0], gps[i, j].location[1] };
            }

            rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty(point[direction]);
          }
          DoubleMatrix NeTNe = rho * MatrixFunctions.Product(Ni, Ni.Transpose());
          //gps[i, j].detJ
          Me += gps[i, j].weight * detJbar * Math.Abs((double)gps[i, j].GetValue(DataInGausspoint.detJ)) * NeTNe;
        }
      }

      //for (int j = 0; j < d * nen; j++)
      //  for (int i = j + 1; i < d * nen; i++)
      //  {
      //    Me[i, j] = Me[j, i];
      //  }
    }

    private DoubleMatrix CreateMaterialThermalMatrix()
    {
      var thermalConductivity = Material.GetProperty(MaterialPropertyName.IsotropicThermalConductivity).GetValueProperty();

      DoubleMatrix C = new DoubleMatrix(3, 3);
      C[0, 0] = C[1, 1] = C[2, 2] = thermalConductivity;
      return C;
    }

    private DoubleMatrix CreateMaterialThermalMatrix(int i, int j)
    {
      var thermalConductivity = (double)gps[i, j].GetValue(DataInGausspoint.ThermalConductivity);// gps[i, j].ThermalConductivity;

      DoubleMatrix C = new DoubleMatrix(3, 3);
      C[0, 0] = C[1, 1] = C[2, 2] = thermalConductivity;
      return C;
    }

    private DoubleMatrix CreateMaterialThermalMatrix(double xi, double eta)
    {
      double[] point = null;
      int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
      if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
      {
        point = PointAt(xi, eta);
      }
      else
      {
        point = new double[] { xi, eta };
      }
      var thermalConductivity = Material.GetProperty(MaterialPropertyName.FGMThermalConductivity).GetValueProperty(point[direction]);
      DoubleMatrix C = new DoubleMatrix(3, 3);
      C[0, 0] = C[1, 1] = C[2, 2] = thermalConductivity;
      return C;
    }

    //////////////////////////////// Tien
    private DoubleVector CreateVectorCoefficientsExpansion()
    {
      var coefficientsthermalexpansion = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion).GetValueProperty();
      DoubleVector A = new DoubleVector(4);
      A[0] = A[1] = A[2] = coefficientsthermalexpansion;
      return A;
    }

    private DoubleVector CreateVectorCoefficientsExpansion(int i, int j)
    {
      var coefficientsthermalexpansion = (double)gps[i, j].GetValue(DataInGausspoint.CoefficientsThermalExpansion);//gps[i, j].coefficientsthermalexpansion;
      DoubleVector A = new DoubleVector(4);
      A[0] = A[1] = A[2] = coefficientsthermalexpansion;
      return A;
    }

    private DoubleVector CreateVectorCoefficientsExpansion(double xi, double eta)
    {
      double[] point = null;
      int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
      if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
      {
        point = PointAt(xi, eta);
      }
      else
      {
        point = new double[] { xi, eta };
      }
      var coefficientsthermalexpansion = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion).GetValueProperty(point[direction]);
      DoubleVector A = new DoubleVector(4);
      A[0] = A[1] = A[2] = coefficientsthermalexpansion;
      return A;
    }
    private DoubleVector CreateVectorThermoelasticCoefficent()
    {
      DoubleMatrix D = CreateMaterialMatrix();
      DoubleVector Anpha = CreateVectorCoefficientsExpansion();
      DoubleVector Beta = MatrixFunctions.Product(D, Anpha);
      return Beta;
    }

    private DoubleVector CreateVectorThermoelasticCoefficent(int i, int j)
    {
      DoubleMatrix D = CreateMaterialMatrix(i, j);
      DoubleVector Anpha = CreateVectorCoefficientsExpansion(i, j);
      DoubleVector Beta = MatrixFunctions.Product(D, Anpha);
      return Beta;
    }

    private DoubleVector CrateVectorThermoelasticCoefficent(double xi, double eta)
    {

      DoubleMatrix D = CreateMaterialMatrix(xi, eta);
      DoubleVector Anpha = CreateVectorCoefficientsExpansion(xi, eta);
      DoubleVector Beta = MatrixFunctions.Product(D, Anpha);
      return Beta;
    }


    /////////////////////////////////////////////////
    public override void ComputeStiffnessMatrixElement(out DoubleMatrix Ke)
    {
      PatchThermoelastic2D patch = (PatchThermoelastic2D)this.patch;
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
      DoubleSymmetricMatrix Keuu = new DoubleSymmetricMatrix(2 * nen);
      DoubleSymmetricMatrix Kett = new DoubleSymmetricMatrix(nen);
      DoubleMatrix Keut = new DoubleMatrix(2 * nen, nen);
      DoubleMatrix Bu = new DoubleMatrix(4, 2 * nen);
      DoubleMatrix Bt = new DoubleMatrix(3, nen);
      DoubleMatrix Du = null;
      DoubleMatrix Dt = null;
      DoubleVector Beta = null;


      if (!(Material is FGMStructureOneVariableMaterial))
      {
        Du = CreateMaterialMatrix();
        Dt = CreateMaterialThermalMatrix();
        Beta = CreateVectorThermoelasticCoefficent();
      }
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          int c = 0;
          double x1 = 1;
          double a = 1;
          if (patch.StateStress == Structure2DState.Axisymetric)
          {
            c = 1;
            x1 = PointAt(gps[i, j].location[0], gps[i, j].location[1])[0];
            a = 2 * Math.PI * x1;
          }
          else if (patch.StateStress == Structure2DState.PlaneStress)
          {
            a = patch.Thickness;
            if (a == 0)
              throw new InvalidArgumentException("Plane stress condition must had thickness");
          }
          double detJbar = 1.0 / 4.0
                * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);

          DoubleMatrix dNdx = (DoubleMatrix)gps[i, j].GetValue(DataInGausspoint.dNdX);// gps[i, j].dNdX;//dNdxi * invertJ;
          DoubleVector N = (DoubleVector)gps[i, j].GetValue(DataInGausspoint.Ni);
          for (int k = 0; k < nen; k++)
          {
            Bu[0, 2 * k] = dNdx[k, 0];
            Bu[2, 2 * k] = c * N[k] / x1;
            Bu[3, 2 * k] = dNdx[k, 1];
            Bu[1, 2 * k + 1] = dNdx[k, 1];
            Bu[3, 2 * k + 1] = dNdx[k, 0];
            Bt[0, k] = dNdx[k, 0];
            Bt[1, k] = dNdx[k, 1];
            Bt[2, k] = c * N[k] / x1;
          }
          if (Material is FGMStructureOneVariableMaterial)
          {
            Du = CreateMaterialMatrix(i, j);
            Dt = CreateMaterialThermalMatrix(i, j);
            Beta = CreateVectorThermoelasticCoefficent(i, j);
          }
          DoubleSymmetricMatrix BuTDuBu = new DoubleSymmetricMatrix(NMathFunctions.TransposeProduct(Bu, MatrixFunctions.Product(Du, Bu)));
          DoubleSymmetricMatrix BtTDtBt = new DoubleSymmetricMatrix(NMathFunctions.TransposeProduct(Bt, MatrixFunctions.Product(Dt, Bt)));
          DoubleMatrix BuTBetaNT = Bu.Transpose() * MatrixFunctions.OuterProduct(Beta, N);
          double detJ = Math.Abs((double)gps[i, j].GetValue(DataInGausspoint.detJ));
          Keuu += a * gps[i, j].weight * detJbar * detJ * BuTDuBu;
          Kett += a * gps[i, j].weight * detJbar * detJ * BtTDtBt;
          Keut -= a * gps[i, j].weight * detJbar * detJ * BuTBetaNT;
        }
      }

      for (int j = 0; j < nen; j++)
        for (int i = j; i < nen; i++)
        {
          for (int jj = 0; jj < 2; jj++)
          {
            for (int ii = 0; ii < 2; ii++)
            {
              Ke[i * 3 + ii, j * 3 + jj] = Keuu[i * 2 + ii, j * 2 + jj];
            }
            Ke[j * 3 + jj, i * 3 + 2] = Ke[i * 3 + 2, j * 3 + jj] = Keut[j * 2 + jj, i];
          }
          Ke[i * 3 + 2, j * 3 + 2] = Kett[i, j];
        }
    }




    public override DoubleVector StrainElasticAt(params double[] xi)
    {
      PatchThermoelastic2D patch = (PatchThermoelastic2D)this.patch;
      int d = 2;
      int nen = patch.GetCountLocalBasisFunctions();
      var cps = ((NURBSSurface)(patch.GetGeometry())).ControlPoints;
      DoubleMatrix B = new DoubleMatrix(4, d * nen);

      DoubleMatrix gradNE = GradBasisFunction(xi[0], xi[1]);
      DoubleMatrix J = JacobianAt(gradNE);
      DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
      DoubleMatrix gradBasis = gradNE * invertJ;
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
        B[0, 2 * k] = gradBasis[k, 0];
        B[2, 2 * k] = c * N[k] / x1;
        B[3, 2 * k] = gradBasis[k, 1];
        B[1, 2 * k + 1] = gradBasis[k, 1];
        B[3, 2 * k + 1] = gradBasis[k, 0];
      }
      DoubleVector U = GetDisplacementLocal();
      return MatrixFunctions.Product(B, U);
    }

    public override DoubleVector StrainThermoAt(params double[] xi)
    {
      double Tref = Material.GetProperty(MaterialPropertyName.ReferenceTemperature).GetValueProperty();
      PatchThermoelastic2D patch = (PatchThermoelastic2D)this.patch;
      double T = patch.GetApproximateAt(Result.TEMP, xi[0], xi[1]);
      DoubleVector Anpha = null;
      if (!(Material is FGMStructureOneVariableMaterial))
        Anpha = CreateVectorCoefficientsExpansion();
      else
        Anpha = CreateVectorCoefficientsExpansion(xi[0], xi[1]);
      return Anpha * (T - Tref);

    }

    public override DoubleVector StrainAt(params double[] xi)
    {
      DoubleVector StrainElastic = StrainElasticAt(xi[0], xi[1]);
      DoubleVector StrainThermo = StrainThermoAt(xi[0], xi[1]);
      return StrainElastic + StrainThermo;

    }

    private DoubleVector CreatStressSupportPlaneStrain()
    {
      var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
      var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
      var coefficientsthermalexpansion = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion).GetValueProperty();
      DoubleVector AnphaE = new DoubleVector(4);
      AnphaE[0] = AnphaE[1] = AnphaE[2] = (eModulus / (1 - 2 * nu)) * coefficientsthermalexpansion;
      return AnphaE;
    }

    private DoubleVector CreatStressSupportPlaneStrain(int i, int j)
    {
      var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
      var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
      var coefficientsthermalexpansion = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion).GetValueProperty();
      DoubleVector AnphaE = new DoubleVector(4);
      AnphaE[0] = AnphaE[1] = AnphaE[2] = (eModulus / (1 - 2 * nu)) * coefficientsthermalexpansion;
      return AnphaE;
    }

    private DoubleVector CreatStressSupportPlaneStrain(double xi, double eta)
    {
      double[] point = null;
      int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
      if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
      {
        point = PointAt(xi, eta);
      }
      else
      {
        point = new double[] { xi, eta };
      }
      var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty(point[direction]);
      var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty(point[direction]);
      var coefficientsthermalexpansion = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion).GetValueProperty(point[direction]);
      DoubleVector AnphaE = new DoubleVector(4);
      AnphaE[0] = AnphaE[1] = AnphaE[2] = (eModulus / (1 - 2 * nu)) * coefficientsthermalexpansion;
      return AnphaE;
    }

    private DoubleVector CreatStressSupportPlaneStress()
    {
      var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
      var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
      var coefficientsthermalexpansion = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion).GetValueProperty();
      DoubleVector AnphaE = new DoubleVector(4);
      AnphaE[0] = AnphaE[1] = AnphaE[2] = (eModulus / (1 - nu)) * coefficientsthermalexpansion;
      return AnphaE;
    }

    private DoubleVector CreatStressSupportPlaneStress(int i, int j)
    {
      double eModulus = (double)gps[i, j].GetValue(DataInGausspoint.EModulus);
      var nu = (double)gps[i, j].GetValue(DataInGausspoint.nu);
      var coefficientsthermalexpansion = (double)gps[i, j].GetValue(DataInGausspoint.CoefficientsThermalExpansion);//gps[i, j].coefficientsthermalexpansion;
      DoubleVector AnphaE = new DoubleVector(4);
      AnphaE[0] = AnphaE[1] = AnphaE[2] = (eModulus / (1 - nu)) * coefficientsthermalexpansion;
      return AnphaE;
    }

    private DoubleVector CreatStressSupportPlaneStress(double xi, double eta)
    {
      double[] point = null;
      int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
      if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
      {
        point = PointAt(xi, eta);
      }
      else
      {
        point = new double[] { xi, eta };
      }
      var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty(point[direction]);
      var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty(point[direction]);
      var coefficientsthermalexpansion = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion).GetValueProperty(point[direction]);
      DoubleVector AnphaE = new DoubleVector(4);
      AnphaE[0] = AnphaE[1] = AnphaE[2] = (eModulus / (1 - nu)) * coefficientsthermalexpansion;
      return AnphaE;
    }
    public override DoubleVector StressAt(params double[] xi)
    {
      DoubleMatrix Du = null;
      //DoubleVector Anpha = null;
      DoubleVector AnphaE = null;
      double Tref = Material.GetProperty(MaterialPropertyName.ReferenceTemperature).GetValueProperty();
      DoubleVector Strain = StrainElasticAt(xi[0], xi[1]);
      PatchThermoelastic2D patch = (PatchThermoelastic2D)this.patch;
      double T = patch.GetApproximateAt(Result.TEMP, xi[0], xi[1]);

      switch (((IPatchStructure)patch).StateStress)
      {
        case Structure2DState.PlaneStress:
          {
            if (!(Material is FGMStructureOneVariableMaterial))
            {
              Du = CreateMaterialMatrix();
              //Anpha= CreateVectorCoefficientsExpansion();
              AnphaE = CreatStressSupportPlaneStress();
            }
            else
            {
              Du = CreateMaterialMatrix(xi[0], xi[1]);
              //Anpha = CreateVectorCoefficientsExpansion(xi[0], xi[1]);
              AnphaE = CreatStressSupportPlaneStress(xi[0], xi[1]);
            }
          }
          break;
        case Structure2DState.PlaneStrain:
          {
            if (!(Material is FGMStructureOneVariableMaterial))
            {
              Du = CreateMaterialMatrix();
              //Anpha= CreateVectorCoefficientsExpansion();
              AnphaE = CreatStressSupportPlaneStrain();
            }
            else
            {
              Du = CreateMaterialMatrix(xi[0], xi[1]);
              //Anpha = CreateVectorCoefficientsExpansion(xi[0], xi[1]);
              AnphaE = CreatStressSupportPlaneStrain(xi[0], xi[1]);
            }
          }
          break;
      }
      return MatrixFunctions.Product(Du, Strain) - AnphaE * (T - Tref);
    }

    //public override DoubleVector StressAt(params double[] xi)
    //{
    //    DoubleMatrix Du = null;
    //    DoubleVector Beta = null;
    //    double Tref = 25;
    //    DoubleVector Strain = StrainAt(xi[0], xi[1]);
    //    PatchThermoelastic2D patch = (PatchThermoelastic2D)this.patch;
    //    double T = patch.GetApproximateAt(Result.TEMP, xi[0], xi[1]);
    //    if (!(Material is FGMStructureOneVariableMaterial))
    //    {
    //        Du = CreateMaterialMatrix();
    //        Beta = CreateVectorThermoelasticCoefficent();
    //    }
    //    else
    //    {
    //        Du = CreateMaterialMatrix(xi[0], xi[1]);
    //        Beta = CrateVectorThermoelasticCoefficent(xi[0], xi[1]);

    //    }
    //    return Du * Strain - Beta * (T - Tref);
    //}

    public override void ComputeInternalForceElement(out DoubleVector fi)
    {
      throw new NotImplementedException();
    }
  }
}
