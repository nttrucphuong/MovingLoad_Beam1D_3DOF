using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public abstract class AbstractElementStructure3D : AbstractElement3D, IElementStructure
  {
    public AbstractElementStructure3D(AbstractPatch3D mesh, int id)
        : base(mesh, id)
    {
    }

    public override void ComputeMassMatrixElement(out DoubleMatrix Me)
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
      Me = new DoubleMatrix(d * nen, d * nen);
      double rho = 0;
      double nuy = 0;
      if (!(Material is FGMUserDefinedGradedMaterial))
      {
        rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty();
        if (Material.GetProperty(MaterialPropertyName.NonLocalParameter) != null)
          nuy = Material.GetProperty(MaterialPropertyName.NonLocalParameter).GetValueProperty();
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
            //DoubleMatrix dNdxi = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdxi);// gps[i, j, k].dNdxi;
            //DoubleMatrix J = JacobianAt(dNdxi);
            DoubleVector Nijk = (DoubleVector)gps[i, j, k].GetValue(DataInGausspoint.Ni);// basis.GetValueTrivariateBasisFunctions(xi1, xi2, xi3);
            DoubleMatrix Ni = new DoubleMatrix(d * nen, d);
            DoubleMatrix dNi = new DoubleMatrix(d * nen, d);
            DoubleMatrix ddNi = new DoubleMatrix(d * nen, d);
            DoubleMatrix d3Ni = new DoubleMatrix(d * nen, d);
            DoubleMatrix d4Ni = new DoubleMatrix(d * nen, d);
            DoubleMatrix Grad2Ni = null;
            DoubleMatrix Grad0Ni = null;
            DoubleMatrix Grad1Ni = null;
            DoubleMatrix Grad3Ni = null;
            DoubleMatrix Grad4Ni = null;
            DoubleMatrix dNdX = null;
            DoubleMatrix ddNdX = null;
            DoubleMatrix d3NdX = null;
            DoubleMatrix d4NdX = null;
            MaterialProperty matNonlocal = Material.GetProperty(MaterialPropertyName.NonLocalParameter);
            MaterialProperty matGeneralNonlocal1 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter1);
            MaterialProperty matGeneralNonlocal2 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter2);

            double nuy1 = 0;
            double nuy2 = 0;
            if (matNonlocal != null)
            {
              //nuy = ComputeParameterProperty(MaterialPropertyName.NonLocalParameter, gps[i, j, k].location[0], gps[i, j, k].location[1], gps[i, j, k].location[2]);
              if ((Material is FGMUserDefinedGradedMaterial) || (Material is FGMOneGradedDirectionMaterial))
                nuy = (double)gps[i, j, k].GetValue(DataInGausspoint.NonLocalParameter);
              Grad2Ni = new DoubleMatrix(d * nen, d);
              ddNdX = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.ddNdX);
            }
            if (matGeneralNonlocal1 != null && matGeneralNonlocal2 != null)
            {
              nuy1 = ComputeParameterProperty(MaterialPropertyName.GeneralNonLocalParameter1, gps[i, j, k].location[0], gps[i, j, k].location[1], gps[i, j, k].location[2]);
              nuy2 = ComputeParameterProperty(MaterialPropertyName.GeneralNonLocalParameter2, gps[i, j, k].location[0], gps[i, j, k].location[1], gps[i, j, k].location[2]);
              Grad0Ni = new DoubleMatrix(d * nen, d);
              Grad1Ni = new DoubleMatrix(d * nen, d);
              Grad2Ni = new DoubleMatrix(d * nen, d);
              Grad3Ni = new DoubleMatrix(d * nen, d);
              Grad4Ni = new DoubleMatrix(d * nen, d);

              dNdX = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdX);
              ddNdX = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.ddNdX);
              d3NdX = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.d3NdX);
              d4NdX = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.d4NdX);
            }
            int count = 0;
            for (int ii = 0; ii < nen; ii++)
            {
              double N = Nijk[ii];
              Ni[d * count, 0] = N;
              Ni[d * count + 1, 1] = N;
              Ni[d * count + 2, 2] = N;

              if (matNonlocal != null)
              {
                double gama2 = (ddNdX[ii, 0] + ddNdX[ii, 1] + ddNdX[ii, 2]);
                //double gama2 = ddNdX[ii, 0] + ddNdX[ii, 1];
                double Nnuy = N - nuy * gama2;
                Grad2Ni[d * count, 0] = Nnuy;
                Grad2Ni[d * count + 1, 1] = Nnuy;
                Grad2Ni[d * count + 2, 2] = Nnuy;
              }
              if (matGeneralNonlocal1 != null && matGeneralNonlocal2 != null)
              {
                double Gama2 = (ddNdX[ii, 0] + ddNdX[ii, 1] + ddNdX[ii, 2]);
                double Gama4 = (d4NdX[ii, 0] + d4NdX[ii, 1] + d4NdX[ii, 2] + 2 * d4NdX[ii, 4] + 2 * d4NdX[ii, 7] + 2 * d4NdX[ii, 10]);
                double N4nuy = N - (nuy1 + nuy2) * Gama2 + nuy1 * nuy2 * Gama4;
                Grad4Ni[d * count, 0] = N4nuy;
                Grad4Ni[d * count + 1, 1] = N4nuy;
                Grad4Ni[d * count + 2, 2] = N4nuy;

                double Gama0 = N;
                double Gama1 = dNdX[ii, 0] + dNdX[ii, 1] + dNdX[ii, 2];
                double Gama3 = d3NdX[ii, 0] + d3NdX[ii, 1] + d3NdX[ii, 2] + d3NdX[ii, 3] + d3NdX[ii, 4] + d3NdX[ii, 5] + d3NdX[ii, 6] + d3NdX[ii, 7] + d3NdX[ii, 8];
                double N3nuy = -2 * (nuy1 + nuy2) * Gama1 + 4 * nuy1 * nuy2 * Gama3;
                Grad3Ni[d * count, 0] = N3nuy;
                Grad3Ni[d * count + 1, 1] = N3nuy;
                Grad3Ni[d * count + 2, 2] = N3nuy;

                double N2nuy = -(nuy1 + nuy2) * Gama0 + 6 * nuy1 * nuy2 * Gama2;
                Grad2Ni[d * count, 0] = N2nuy;
                Grad2Ni[d * count + 1, 1] = N2nuy;
                Grad2Ni[d * count + 2, 2] = N2nuy;

                double N1nuy = 4 * nuy1 * nuy2 * Gama1;
                Grad1Ni[d * count, 0] = N1nuy;
                Grad1Ni[d * count + 1, 1] = N1nuy;
                Grad1Ni[d * count + 2, 2] = N1nuy;

                double N0nuy = nuy1 * nuy2 * Gama0;
                Grad0Ni[d * count, 0] = N0nuy;
                Grad0Ni[d * count + 1, 1] = N0nuy;
                Grad0Ni[d * count + 2, 2] = N0nuy;

                //dNi[d * count, 0] = dNdX[ii, 0];
                //dNi[d * count + 1, 1] = dNdX[ii, 1];
                //dNi[d * count + 2, 2] = dNdX[ii, 2];

                //ddNi[d * count, 0] = ddNdX[ii, 0];
                //ddNi[d * count + 1, 1] = ddNdX[ii, 1];
                //ddNi[d * count + 2, 2] = ddNdX[ii, 2];

                //d3Ni[d * count, 0] = d3NdX[ii, 0];
                //d3Ni[d * count + 1, 1] = d3NdX[ii, 1];
                //d3Ni[d * count + 2, 2] = d3NdX[ii, 2];

                //d4Ni[d * count, 0] = d4NdX[ii, 0];
                //d4Ni[d * count + 1, 1] = d4NdX[ii, 1];
                //d4Ni[d * count + 2, 2] = d4NdX[ii, 2];
              }
              count++;
            }

            if (Material is FGMUserDefinedGradedMaterial)
            {
              //double[] point = null;
              //int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
              //if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
              //{
              //  point = PointAt(gps[i, j, k].location[0], gps[i, j, k].location[1], gps[i, j, k].location[2]);
              //}
              //else
              //{
              //  point = new double[] { gps[i, j, k].location[0], gps[i, j, k].location[1], gps[i, j, k].location[2] };
              //}

              //rho = Material.GetProperty(MaterialPropertyName.Density).GetValueProperty(point[direction]);
              //rho = ComputeParameterProperty(MaterialPropertyName.Density, gps[i, j, k].location[0], gps[i, j, k].location[1], gps[i, j, k].location[2]);
              rho = (double)gps[i, j, k].GetValue(DataInGausspoint.Density);
            }
            DoubleMatrix NeTNe = null;
            //DoubleMatrix NeTNe = rho * MatrixFunctions.Product(Ni, NiGrad2Ni.Transpose());
            if (matNonlocal == null && matGeneralNonlocal1 == null && matGeneralNonlocal2 == null)
            {
              NeTNe = rho * MatrixFunctions.Product(Ni, Ni.Transpose());
            }
            if (matNonlocal != null)
            {
              NeTNe = rho * MatrixFunctions.Product(Ni, Grad2Ni.Transpose());
              //NeTNe = rho * MatrixFunctions.Product(Grad2Ni, Ni.Transpose());
            }
            if (matGeneralNonlocal1 != null && matGeneralNonlocal2 != null)
            {
              //NeTNe = rho * (MatrixFunctions.Product(Grad4Ni, Ni.Transpose()) + MatrixFunctions.Product(Grad3Ni, dNi.Transpose()) + MatrixFunctions.Product(Grad2Ni, ddNi.Transpose()) + MatrixFunctions.Product(Grad1Ni, d3Ni.Transpose()) + MatrixFunctions.Product(Grad0Ni, d4Ni.Transpose()));
              NeTNe = rho * (MatrixFunctions.Product(Grad4Ni, Grad0Ni.Transpose()) + MatrixFunctions.Product(Grad3Ni, Grad1Ni.Transpose()) + MatrixFunctions.Product(Grad2Ni, Grad2Ni.Transpose()) + MatrixFunctions.Product(Grad1Ni, Grad3Ni.Transpose()) + MatrixFunctions.Product(Grad0Ni, Grad4Ni.Transpose()));
            }
            Me += gps[i, j, k].weight * detJbar * Math.Abs((double)gps[i, j, k].GetValue(DataInGausspoint.detJ)) * NeTNe;
          }
        }
      }
    }

    public override void ComputeGeometricStiffnessMatrixElement(double[,] stressTensor, out DoubleMatrix KGe)
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
      KGe = new DoubleMatrix(d * nen, d * nen);
      double nuy = 0;
      if (!(Material is FGMUserDefinedGradedMaterial))
      {
        if (Material.GetProperty(MaterialPropertyName.NonLocalParameter) != null)
          nuy = Material.GetProperty(MaterialPropertyName.NonLocalParameter).GetValueProperty();
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
            GaussPoints gpsijk = gps[i, j, k];
            DoubleMatrix dNdX = (DoubleMatrix)gpsijk.GetValue(DataInGausspoint.dNdX); //gps[i, j, k].dNdX;
            DoubleMatrix d3NdX = null;
            DoubleMatrix GA = new DoubleMatrix(9, d * nen);
            DoubleMatrix Grad2GA = null;

            MaterialProperty matNonlocal = Material.GetProperty(MaterialPropertyName.NonLocalParameter);
            MaterialProperty matGeneralNonlocal1 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter1);
            MaterialProperty matGeneralNonlocal2 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter2);
            if (matNonlocal != null)
            {
              //nuy = ComputeParameterProperty(MaterialPropertyName.NonLocalParameter, gps[i, j, k].location[0], gps[i, j, k].location[1], gps[i, j, k].location[2]);
              if ((Material is FGMUserDefinedGradedMaterial) || (Material is FGMOneGradedDirectionMaterial))
                nuy = (double)gps[i, j, k].GetValue(DataInGausspoint.NonLocalParameter);
              Grad2GA = new DoubleMatrix(9, d * nen);
              d3NdX = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.d3NdX);
            }
            int count = 0;
            for (int ii = 0; ii < nen; ii++)
            {
              GA[0, 3 * count] = dNdX[ii, 0];
              GA[1, 3 * count] = dNdX[ii, 1];
              GA[2, 3 * count] = dNdX[ii, 2];
              GA[3, 3 * count + 1] = dNdX[ii, 0];
              GA[4, 3 * count + 1] = dNdX[ii, 1];
              GA[5, 3 * count + 1] = dNdX[ii, 2];
              GA[6, 3 * count + 2] = dNdX[ii, 0];
              GA[7, 3 * count + 2] = dNdX[ii, 1];
              GA[8, 3 * count + 2] = dNdX[ii, 2];

              if (matNonlocal != null)
              {
                double grad2GAx = dNdX[ii, 0] - nuy * (d3NdX[ii, 0] + d3NdX[ii, 4] + d3NdX[ii, 8]);
                double grad2GAy = dNdX[ii, 1] - nuy * (d3NdX[ii, 3] + d3NdX[ii, 1] + d3NdX[ii, 6]);
                double grad2GAz = dNdX[ii, 2] - nuy * (d3NdX[ii, 7] + d3NdX[ii, 5] + d3NdX[ii, 2]);
                Grad2GA[0, 3 * count] = grad2GAx;
                Grad2GA[1, 3 * count] = grad2GAy;
                Grad2GA[2, 3 * count] = grad2GAz;
                Grad2GA[3, 3 * count + 1] = grad2GAx;
                Grad2GA[4, 3 * count + 1] = grad2GAy;
                Grad2GA[5, 3 * count + 1] = grad2GAz;
                Grad2GA[6, 3 * count + 2] = grad2GAx;
                Grad2GA[7, 3 * count + 2] = grad2GAy;
                Grad2GA[8, 3 * count + 2] = grad2GAz;
              }

              count++;
            }
            DoubleMatrix N0 = new DoubleMatrix(9, 9);
            //double h = 0.025;
            //N0[2, 2] = -1 / h;
            //N0[5, 5] = -1 / h;
            //N0[8, 8] = -1 / h;
            DoubleVector stress = (DoubleVector)gpsijk.GetValue(DataInGausspoint.currentStress);
            if (stressTensor == null)
            {
              N0[0, 0] = stress[0];
              N0[1, 1] = stress[1];
              N0[2, 2] = stress[2];
              N0[0, 1] = N0[1, 0] = stress[3];
              N0[1, 2] = N0[2, 1] = stress[4];
              N0[0, 2] = N0[2, 0] = stress[5];

              N0[3, 3] = stress[0];
              N0[4, 4] = stress[1];
              N0[5, 5] = stress[2];
              N0[3, 4] = N0[4, 3] = stress[3];
              N0[4, 5] = N0[5, 4] = stress[4];
              N0[3, 5] = N0[5, 3] = stress[5];

              N0[6, 6] = stress[0];
              N0[7, 7] = stress[1];
              N0[8, 8] = stress[2];
              N0[6, 7] = N0[7, 6] = stress[3];
              N0[7, 8] = N0[8, 7] = stress[4];
              N0[6, 8] = N0[8, 6] = stress[5];
            }
            else
            {
              int rows = stressTensor.GetLength(0);
              for (int ii = 0; ii < rows; ii++)
              {
                int cols = stressTensor.GetLength(1);
                for (int jj = 0; jj < cols; jj++)
                {
                  for (int kk = 0; kk < 3; kk++)
                  {
                    N0[ii + kk * rows, jj + kk * cols] = stressTensor[ii, jj];
                  }
                }
              }
            }

            DoubleMatrix GaTN0GA = null;

            if (matNonlocal == null && matGeneralNonlocal1 == null && matGeneralNonlocal2 == null)
            {
              GaTN0GA = NMathFunctions.TransposeProduct(GA, MatrixFunctions.Product(-N0, GA));
            }
            if (matNonlocal != null)
            {
              GaTN0GA = NMathFunctions.TransposeProduct(GA, MatrixFunctions.Product(-N0, Grad2GA));
              //GaTN0GA = NMathFunctions.TransposeProduct(Grad2GA, MatrixFunctions.Product(-N0, GA));
            }
            KGe += gpsijk.weight * detJbar * Math.Abs((double)gpsijk.GetValue(DataInGausspoint.detJ)) * GaTN0GA;
          }
        }
      }
    }
    public abstract DoubleVector StressAt(params double[] xi);
    public abstract DoubleVector StrainAt(params double[] xi);
    public abstract DoubleVector StrainElasticAt(params double[] xi);
    public abstract DoubleVector StrainThermoAt(params double[] xi);
    public override void ComputeValueAtGaussPoint(DataInGausspoint name)
    {
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {
            switch (name)
            {
              case DataInGausspoint.currentStress:
                gps[i, j, k].SetValue(DataInGausspoint.currentStress, StressAt(gps[i, j, k].location));
                //gps[i, j, k].currentStress = StressAt(gps[i, j, k].location);
                break;
              case DataInGausspoint.currentStrain:
                gps[i, j, k].SetValue(DataInGausspoint.currentStrain, StrainAt(gps[i, j, k].location));
                //gps[i, j, k].currentStrain = StrainAt(gps[i, j, k].location);
                break;
            }
          }
        }
      }
    }
  }
}
