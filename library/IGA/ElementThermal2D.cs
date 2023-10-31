using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public class ElementThermal2D : AbstractElement2D, IElementThermal
  {
    private double[] convectionHeatTransfer;//The convective heat transfer boundary on 4 edges
    private double[] temperatureHeatTransferConvection;//temperature on egde has convective heat transfer


    public ElementThermal2D(AbstractPatch2D mesh, int id)
          : base(mesh, id)
    {
    }

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
          DoubleSymmetricMatrix keSub = new DoubleSymmetricMatrix(e.CountDofOnEdge());
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

            //double[,] Nuv = basis.GetValueBivariateBasisFunctions(xi, eta);
            //DoubleVector Ni = new DoubleVector(((Thermal2DPatch)this.patch).GetCountLocalBasisFunctions());

            //int count = 0;
            //for (int j = 0; j <= q; j++)
            //    for (int i = 0; i <= p; i++)
            //    {
            //        Ni[count++] = Nuv[i, j];
            //    }

            DoubleVector Ni = new DoubleVector(e.GetBivariateBasisFunctionOnEdge(0, xi, eta));
            double J2 = (e.GetIndexCoordinate() == 0) ? (0.5 * (paraEndPatchU[1] - paraEndPatchU[0]))
                  : (0.5 * (paraEndPatchV[1] - paraEndPatchV[0]));

            DoubleVector dNi = new DoubleVector(e.GetDerivativeBivariateBasisFunctionOnEdge(0, xi, eta));
            double[] J = JacobianOnEdge(dNi, iii).ToArray();
            double J1 = Math.Sqrt(J[0] * J[0] + J[1] * J[1]);

            keSub += new DoubleSymmetricMatrix(w * J1 * J2 * MatrixFunctions.OuterProduct(Ni, Ni) * convectionHeatTransfer[iii]);
          }

          int[] tArray = e.GetTArray();
          for (int j = 0; j < e.CountDofOnEdge(); j++)
            for (int i = j; i < e.CountDofOnEdge(); i++)
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
          DoubleVector rLocal = new DoubleVector(e.CountDofOnEdge());
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
          int[] tArrayGlobal = e.GetTArrayGlobal();
          for (int i = 0; i < e.CountDofOnEdge(); i++)
          {
            int m = tArrayGlobal[i];
            rGlobal[m] += rLocal[i];
          }
        }
      }
    }

    public void ComputeHeatGenerationSourceLoadVectorElement(ref DoubleVector rGlobal, double G)
    {
      PatchThermal2D patch = (PatchThermal2D)this.patch;
      int d = patch.GetCountField();
      BivariateNURBSBasisFunction basis = (BivariateNURBSBasisFunction)((NURBSSurface)(patch.GetGeometry())).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int nen = (p + 1) * (q + 1); // number of local basis functions
      var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      int idx1 = patch.GetIPN(id, 0);
      int idx2 = patch.GetIPN(id, 1);
      DoubleVector rLocal = new DoubleVector(d * nen);
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          double detJbar = 1.0 / 4.0
                * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          //DoubleMatrix dNdxi = gps[i, j].dNdxi;
          //DoubleMatrix J = JacobianAt(dNdxi);

          GaussPoints gpsij = gps[i, j];
          DoubleVector Ni = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);//gps[i, j].Ni;
          rLocal += gpsij.weight * detJbar * Math.Abs((double)gpsij.GetValue(DataInGausspoint.detJ)) * G * Ni;
        }
      }
      int[] tArrayGlobal = GetFace().GetTArrayGlobal();
      for (int i = 0; i < GetFace().CountDof(); i++)
      {
        int m = tArrayGlobal[i];
        rGlobal[m] += rLocal[i];
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

    private DoubleMatrix CreateMaterialMatrix()
    {
      var thermalConductivity = Material.GetProperty(MaterialPropertyName.IsotropicThermalConductivity).GetValueProperty();

      DoubleMatrix C = new DoubleMatrix(2, 2);
      C[0, 0] = C[1, 1] = thermalConductivity;
      return C;
    }

    private DoubleMatrix CreateMaterialMatrix(int i, int j)
    {
      double thermalConductivity = (double)gps[i, j].GetValue(DataInGausspoint.ThermalConductivity);

      DoubleMatrix C = new DoubleMatrix(2, 2);
      C[0, 0] = C[1, 1] = thermalConductivity;
      return C;
    }

    private DoubleMatrix CreateMaterialMatrix(double xi, double eta)
    {
      double[] point = null;
      int direction = ((FGMThermalOneGradedDirectionMaterial)Material).GetDirectionGraded();
      if (((FGMThermalOneGradedDirectionMaterial)Material).GetIsGlobal())
      {
        point = PointAt(xi, eta);
      }
      else
      {
        point = new double[] { xi, eta };
      }
      var thermalConductivity = Material.GetProperty(MaterialPropertyName.FGMThermalConductivity).GetValueProperty(point[direction]);


      DoubleMatrix C = new DoubleMatrix(2, 2);
      C[0, 0] = C[1, 1] = thermalConductivity;
      return C;
    }
    public override void ComputeStiffnessMatrixElement(out DoubleMatrix Ke)
    {
      PatchThermal2D patch = (PatchThermal2D)this.patch;
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
          double detJbar = 1.0 / 4.0
                * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          //DoubleMatrix dNdxi = gps[i, j].dNdxi;
          //DoubleMatrix J = JacobianAt(dNdxi);
          //DoubleMatrix invertJ = (DoubleMatrix)J.Inverse();
          GaussPoints gpsij = gps[i, j];
          DoubleMatrix dNdx = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.dNdX);//gps[i, j].dNdX;//dNdxi * invertJ;
          for (int k = 0; k < nen; k++)
          {
            B[0, d * k] = dNdx[k, 0];
            B[1, d * k] = dNdx[k, 1];
          }
          if (Material is FGMThermalOneGradedDirectionMaterial)
            D = CreateMaterialMatrix(i, j);
          DoubleMatrix BTDB = NMathFunctions.TransposeProduct(B, MatrixFunctions.Product(D, B));
          Ke += gpsij.weight * detJbar * Math.Abs((double)gpsij.GetValue(DataInGausspoint.detJ)) * BTDB;
        }
      }
    }
    /////////////////////////////////

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

    public override void ComputeMassMatrixElement(out DoubleMatrix Me)
    {
      PatchThermal2D patch = (PatchThermal2D)this.patch;
      BivariateNURBSBasisFunction basis = (BivariateNURBSBasisFunction)((NURBSSurface)(patch.GetGeometry())).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int nen = (p + 1) * (q + 1); // number of local basis functions
      var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      int idx1 = patch.GetIPN(id, 0);
      int idx2 = patch.GetIPN(id, 1);
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
          double detJbar = 1.0 / 4.0
              * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
              * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          //DoubleMatrix dNdxi = gps[i, j].dNdxi;
          //DoubleMatrix J = JacobianAt(dNdxi);
          GaussPoints gpsij = gps[i, j];
          DoubleVector Ni = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);//gps[i, j].Ni;
          //DoubleMatrix Ni = new DoubleMatrix(nen, 1);
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
              point = PointAt(gpsij.location[0], gpsij.location[1]);
            }
            else
            {
              point = new double[] { gpsij.location[0], gpsij.location[1] };
            }
            rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty(point[direction]);
            cc = Material.GetProperty(MaterialPropertyName.SpecificHeatCapacity).GetValueProperty(point[direction]);
          }
          DoubleMatrix NeTNe = rho * cc * MatrixFunctions.OuterProduct(Ni, Ni);

          Me += gpsij.weight * detJbar * Math.Abs((double)gpsij.GetValue(DataInGausspoint.detJ)) * NeTNe;
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
