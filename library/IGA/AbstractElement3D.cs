using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public abstract class AbstractElement3D : AbstractElement
  {
    protected Volume volume;//2 face on 2 field ux,uy,theta, pressure
    protected GaussPoints[,,] gps;

    public AbstractElement3D(AbstractPatch3D mesh, int id)
        : base(mesh, id)
    {
      volume = new Volume(this);
    }

    public override void InitializeGaussPoint()
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      TrivariateNURBSBasisFunction basis = (TrivariateNURBSBasisFunction)((NURBSVolume)(patch.GetGeometry())).Basis;
      var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply3 = basis.GetKnotVector(2).GetKnotVectorNoMultiplicity();
      int idx1 = patch.GetIPN(id, 0);
      int idx2 = patch.GetIPN(id, 1);
      int idx3 = patch.GetIPN(id, 2);
      int nen = patch.GetCountLocalBasisFunctions();// number of local basis functions
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      numberOfGaussPointOnEachDirection = Math.Max(Math.Max(p, q), r) + 1;//(int)Math.Ceiling((2 * Math.Max(Math.Max(p, q), r) + 1) / 2.0);

      gps = new GaussPoints[numberOfGaussPointOnEachDirection, numberOfGaussPointOnEachDirection, numberOfGaussPointOnEachDirection];
      MaterialProperty matNuy = Material.GetProperty(MaterialPropertyName.NonLocalParameter);
      MaterialProperty matNuy1 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter1);
      MaterialProperty matNuy2 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter2);

      bool isNeed_d2NdX = false;
      bool isNeed_d3NdX = false;
      bool isNeed_d4NdX = false;
      if (matNuy != null)
      {
        isNeed_d2NdX = true;
        if (AbstractModel.TypeAnalysisModel is TypeAnalysisModel.Buckling)
          isNeed_d3NdX = true;
      }
      if (matNuy1 != null && matNuy2 != null)
      {
        isNeed_d2NdX = true;
        isNeed_d3NdX = true;
        isNeed_d4NdX = true;
      }
      if (this is ElementStructurePhaseField3D)
      {
        isNeed_d2NdX = true;
      }

        for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        double xi1 = GaussPoints.GetPoint(numberOfGaussPointOnEachDirection, i);
        double w1 = GaussPoints.GetWeight(numberOfGaussPointOnEachDirection, i);
        xi1 = 0.5 * ((kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1]) * xi1 + kvNoMulticiply1[idx1 + 1] + kvNoMulticiply1[idx1]);
        int span1 = basis.FindSpan(xi1, 0);
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          double xi2 = GaussPoints.GetPoint(numberOfGaussPointOnEachDirection, j);
          double w2 = GaussPoints.GetWeight(numberOfGaussPointOnEachDirection, j);
          xi2 = 0.5 * ((kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]) * xi2 + kvNoMulticiply2[idx2 + 1] + kvNoMulticiply2[idx2]);
          int span2 = basis.FindSpan(xi2, 1);
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {
            double xi3 = GaussPoints.GetPoint(numberOfGaussPointOnEachDirection, k);
            double w3 = GaussPoints.GetWeight(numberOfGaussPointOnEachDirection, k);
            xi3 = 0.5 * ((kvNoMulticiply3[idx3 + 1] - kvNoMulticiply3[idx3]) * xi3 + kvNoMulticiply3[idx3 + 1] + kvNoMulticiply3[idx3]);
            int span3 = basis.FindSpan(xi3, 2);

            GaussPoints gpsijk = gps[i, j, k] = new GaussPoints(new double[] { xi1, xi2, xi3 }, w1 * w2 * w3);
            gpsijk.listDataGausspoint = new System.Collections.Generic.List<DataGausspoint>();
            DoubleMatrix dNdxi = new DoubleMatrix(nen, 3);
            DoubleVector Ni = new DoubleVector(nen);
            DoubleMatrix ddNdxi = null;
            DoubleMatrix d3Ndxi = null;
            DoubleMatrix d4Ndxi = null;
            //if (matNuy != null || (matNuy1 != null && matNuy2 != null))
            if (isNeed_d2NdX)
            {
              ddNdxi = new DoubleMatrix(nen, 6);
            }
            //if (matNuy1 != null && matNuy2 != null)
            //{
            //  d3Ndxi = new DoubleMatrix(nen, 10);
            //  d4Ndxi = new DoubleMatrix(nen, 15);
            //}
            if (isNeed_d3NdX)
              d3Ndxi = new DoubleMatrix(nen, 10);
            if (isNeed_d4NdX)
              d4Ndxi = new DoubleMatrix(nen, 15);

            if (patch.ExtractionOperator == null)
            {
              double[,,][,,] gradBasis = null;
              //if (matNuy != null/*|| this is ElementStructurePhaseField3D*/)
              if (isNeed_d4NdX)
              {
                gradBasis = basis.GetDerivativeTrivariateBasisFunctions(xi1, xi2, xi3, 4);
              }
              else if (isNeed_d3NdX)
              {
                gradBasis = basis.GetDerivativeTrivariateBasisFunctions(xi1, xi2, xi3, 3);
              }
              else if (isNeed_d2NdX)
              {
                gradBasis = basis.GetDerivativeTrivariateBasisFunctions(xi1, xi2, xi3, 2);
              }
              //else if (matNuy1 != null && matNuy2 != null)
              //{
              //  gradBasis = basis.GetDerivativeTrivariateBasisFunctions(xi1, xi2, xi3, 4);
              //}
              else
                gradBasis = basis.GetDerivativeTrivariateBasisFunctions(xi1, xi2, xi3, 1);
              for (int ii = 0; ii < nen; ii++)
              {
                int ien = patch.GetIEN(id, ii);
                int v1 = patch.GetINC(0, ien, 0) - span1 + p;
                int v2 = patch.GetINC(0, ien, 1) - span2 + q;
                int v3 = patch.GetINC(0, ien, 2) - span3 + r;
                Ni[ii] = gradBasis[v1, v2, v3][0, 0, 0];
                dNdxi[ii, 0] = gradBasis[v1, v2, v3][1, 0, 0];
                dNdxi[ii, 1] = gradBasis[v1, v2, v3][0, 1, 0];
                dNdxi[ii, 2] = gradBasis[v1, v2, v3][0, 0, 1];
                //if (matNuy != null || (matNuy1 != null && matNuy2 != null)/*|| this is ElementStructurePhaseField3D*/)
                if (isNeed_d2NdX)
                {
                  ddNdxi[ii, 0] = gradBasis[v1, v2, v3][2, 0, 0];
                  ddNdxi[ii, 1] = gradBasis[v1, v2, v3][0, 2, 0];
                  ddNdxi[ii, 2] = gradBasis[v1, v2, v3][0, 0, 2];
                  ddNdxi[ii, 3] = gradBasis[v1, v2, v3][1, 1, 0];
                  ddNdxi[ii, 4] = gradBasis[v1, v2, v3][0, 1, 1];
                  ddNdxi[ii, 5] = gradBasis[v1, v2, v3][1, 0, 1];
                }
                //if (matNuy1 != null && matNuy2 != null)
                //{
                //  d3Ndxi[ii, 0] = gradBasis[v1, v2, v3][3, 0, 0];
                //  d3Ndxi[ii, 1] = gradBasis[v1, v2, v3][0, 3, 0];
                //  d3Ndxi[ii, 2] = gradBasis[v1, v2, v3][0, 0, 3];
                //  d3Ndxi[ii, 3] = gradBasis[v1, v2, v3][2, 1, 0];
                //  d3Ndxi[ii, 4] = gradBasis[v1, v2, v3][1, 2, 0];
                //  d3Ndxi[ii, 5] = gradBasis[v1, v2, v3][0, 2, 1];
                //  d3Ndxi[ii, 6] = gradBasis[v1, v2, v3][0, 1, 2];
                //  d3Ndxi[ii, 7] = gradBasis[v1, v2, v3][2, 0, 1];
                //  d3Ndxi[ii, 8] = gradBasis[v1, v2, v3][1, 0, 2];
                //  d3Ndxi[ii, 9] = gradBasis[v1, v2, v3][1, 1, 1];

                //  d4Ndxi[ii, 0] = gradBasis[v1, v2, v3][4, 0, 0];
                //  d4Ndxi[ii, 1] = gradBasis[v1, v2, v3][0, 4, 0];
                //  d4Ndxi[ii, 2] = gradBasis[v1, v2, v3][0, 0, 4];
                //  d4Ndxi[ii, 3] = gradBasis[v1, v2, v3][3, 1, 0];
                //  d4Ndxi[ii, 4] = gradBasis[v1, v2, v3][2, 2, 0];
                //  d4Ndxi[ii, 5] = gradBasis[v1, v2, v3][1, 3, 0];
                //  d4Ndxi[ii, 6] = gradBasis[v1, v2, v3][0, 3, 1];
                //  d4Ndxi[ii, 7] = gradBasis[v1, v2, v3][0, 2, 2];
                //  d4Ndxi[ii, 8] = gradBasis[v1, v2, v3][0, 1, 3];
                //  d4Ndxi[ii, 9] = gradBasis[v1, v2, v3][3, 0, 1];
                //  d4Ndxi[ii, 10] = gradBasis[v1, v2, v3][2, 0, 2];
                //  d4Ndxi[ii, 11] = gradBasis[v1, v2, v3][1, 0, 3];
                //  d4Ndxi[ii, 12] = gradBasis[v1, v2, v3][2, 1, 1];
                //  d4Ndxi[ii, 13] = gradBasis[v1, v2, v3][1, 2, 1];
                //  d4Ndxi[ii, 14] = gradBasis[v1, v2, v3][1, 1, 2];
                //}
                if (isNeed_d3NdX)
                {
                  d3Ndxi[ii, 0] = gradBasis[v1, v2, v3][3, 0, 0];
                  d3Ndxi[ii, 1] = gradBasis[v1, v2, v3][0, 3, 0];
                  d3Ndxi[ii, 2] = gradBasis[v1, v2, v3][0, 0, 3];
                  d3Ndxi[ii, 3] = gradBasis[v1, v2, v3][2, 1, 0];
                  d3Ndxi[ii, 4] = gradBasis[v1, v2, v3][1, 2, 0];
                  d3Ndxi[ii, 5] = gradBasis[v1, v2, v3][0, 2, 1];
                  d3Ndxi[ii, 6] = gradBasis[v1, v2, v3][0, 1, 2];
                  d3Ndxi[ii, 7] = gradBasis[v1, v2, v3][2, 0, 1];
                  d3Ndxi[ii, 8] = gradBasis[v1, v2, v3][1, 0, 2];
                  d3Ndxi[ii, 9] = gradBasis[v1, v2, v3][1, 1, 1];
                }
                if (isNeed_d4NdX)
                {
                  d4Ndxi[ii, 0] = gradBasis[v1, v2, v3][4, 0, 0];
                  d4Ndxi[ii, 1] = gradBasis[v1, v2, v3][0, 4, 0];
                  d4Ndxi[ii, 2] = gradBasis[v1, v2, v3][0, 0, 4];
                  d4Ndxi[ii, 3] = gradBasis[v1, v2, v3][3, 1, 0];
                  d4Ndxi[ii, 4] = gradBasis[v1, v2, v3][2, 2, 0];
                  d4Ndxi[ii, 5] = gradBasis[v1, v2, v3][1, 3, 0];
                  d4Ndxi[ii, 6] = gradBasis[v1, v2, v3][0, 3, 1];
                  d4Ndxi[ii, 7] = gradBasis[v1, v2, v3][0, 2, 2];
                  d4Ndxi[ii, 8] = gradBasis[v1, v2, v3][0, 1, 3];
                  d4Ndxi[ii, 9] = gradBasis[v1, v2, v3][3, 0, 1];
                  d4Ndxi[ii, 10] = gradBasis[v1, v2, v3][2, 0, 2];
                  d4Ndxi[ii, 11] = gradBasis[v1, v2, v3][1, 0, 3];
                  d4Ndxi[ii, 12] = gradBasis[v1, v2, v3][2, 1, 1];
                  d4Ndxi[ii, 13] = gradBasis[v1, v2, v3][1, 2, 1];
                  d4Ndxi[ii, 14] = gradBasis[v1, v2, v3][1, 1, 2];
                }
              }
            }
            else
            {
              int indexGausspointOnLine = k * numberOfGaussPointOnEachDirection * numberOfGaussPointOnEachDirection + j * numberOfGaussPointOnEachDirection + i;
              DoubleVector Bb = (DoubleVector)patch.GaussPointOnBezierElement[indexGausspointOnLine].GetValue(DataInGausspoint.NBernsteini);//patch.GaussPointOnBezierElement[indexGausspointOnLine].NBernsteini;
              //int ex = GetIndexGlobalCoordinate(0);
              //int ey = GetIndexGlobalCoordinate(1);
              DoubleVector we = new DoubleVector((p + 1) * (q + 1) * (r + 1));
              int c = 0;
              for (int kk = 0; kk < r + 1; kk++)
                for (int jj = 0; jj < q + 1; jj++)
                  for (int ii = 0; ii < p + 1; ii++)
                  {
                    we[c] = ((TrivariateNURBSBasisFunction)patch.GetVolume().Basis).Weights[idx1 + ii, idx2 + jj, idx3 + kk];
                    c++;
                  }
              DoubleMatrix mCe = new DoubleMatrix(patch.ExtractionOperator[idx1, idx2, idx3]);
              DoubleVector wb = MatrixFunctions.Product(mCe.Transpose(), we);
              double Wb = MatrixFunctions.Dot(Bb, wb);
              DoubleMatrix weMat = new DoubleMatrix(we.Length, we.Length);
              for (int iii = 0; iii < we.Length; iii++)
                weMat[iii, iii] = we[iii];
              DoubleVector Rb = MatrixFunctions.Product(weMat, MatrixFunctions.Product(mCe, Bb)) / Wb;

              ////////////////////////;;;;
              DoubleVector[] dB = new DoubleVector[3];
              DoubleMatrix matdNB = (DoubleMatrix)patch.GaussPointOnBezierElement[indexGausspointOnLine].GetValue(DataInGausspoint.dNBernsteindxi);
              dB[0] = matdNB.Col(0);
              dB[1] = matdNB.Col(1);
              dB[2] = matdNB.Col(2);
              //dB[0] = patch.GaussPointOnBezierElement[indexGausspointOnLine].dNBernsteindxi.Col(0);
              //dB[1] = patch.GaussPointOnBezierElement[indexGausspointOnLine].dNBernsteindxi.Col(1);
              //dB[2] = patch.GaussPointOnBezierElement[indexGausspointOnLine].dNBernsteindxi.Col(2);
              double dWbxi = MatrixFunctions.Dot(dB[0], wb);
              double dWbeta = MatrixFunctions.Dot(dB[1], wb);
              double dWbzeta = MatrixFunctions.Dot(dB[2], wb);
              DoubleVector dRxi = MatrixFunctions.Product(weMat, MatrixFunctions.Product(mCe, (dB[0] / Wb - dWbxi * Bb / (Wb * Wb)))) / (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1]);
              DoubleVector dReta = MatrixFunctions.Product(weMat, MatrixFunctions.Product(mCe, (dB[1] / Wb - dWbeta * Bb / (Wb * Wb)))) / (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
              DoubleVector dRzeta = MatrixFunctions.Product(weMat, MatrixFunctions.Product(mCe, (dB[2] / Wb - dWbzeta * Bb / (Wb * Wb)))) / (kvNoMulticiply3[idx3 + 1] - kvNoMulticiply3[idx3]);
              Ni = Rb;
              for (int t = 0; t < nen; t++)
              {
                dNdxi[t, 0] = dRxi[t];
                dNdxi[t, 1] = dReta[t];
                dNdxi[t, 2] = dRzeta[t];
              }
            }

            DoubleMatrix dxdxi = JacobianAt(dNdxi);//[3,3][x,xi]
            DoubleMatrix invertJ = MatrixFunctions.Inverse(dxdxi);
            DoubleMatrix dNdx = MatrixFunctions.Product(dNdxi, invertJ);
            gpsijk.SetValue(DataInGausspoint.Ni, Ni);
            gpsijk.SetValue(DataInGausspoint.dNdxi, dNdxi);
            gpsijk.SetValue(DataInGausspoint.detJ, MatrixFunctions.Determinant(dxdxi));
            gpsijk.SetValue(DataInGausspoint.dNdX, dNdx);
            DoubleMatrix ddxdxi = null;
            DoubleMatrix d3xdxi = null;
            DoubleMatrix ddNdX = null;
            DoubleMatrix d3NdX = null;
            DoubleMatrix d4NdX = null;

            //if (matNuy != null || (matNuy1 != null && matNuy2 != null) /*|| this is ElementStructurePhaseField3D*/)
            if (isNeed_d2NdX)
            {
              gpsijk.SetValue(DataInGausspoint.ddNdxi, ddNdxi);
              DoubleMatrix jac2 = new DoubleMatrix(6, 6);
              jac2[0, 0] = dxdxi[0, 0] * dxdxi[0, 0];
              jac2[0, 1] = dxdxi[1, 0] * dxdxi[1, 0];
              jac2[0, 2] = dxdxi[2, 0] * dxdxi[2, 0];
              jac2[0, 3] = 2 * dxdxi[0, 0] * dxdxi[1, 0];
              jac2[0, 4] = 2 * dxdxi[1, 0] * dxdxi[2, 0];
              jac2[0, 5] = 2 * dxdxi[0, 0] * dxdxi[2, 0];
              jac2[1, 0] = dxdxi[0, 1] * dxdxi[0, 1];
              jac2[1, 1] = dxdxi[1, 1] * dxdxi[1, 1];
              jac2[1, 2] = dxdxi[2, 1] * dxdxi[2, 1];
              jac2[1, 3] = 2 * dxdxi[0, 1] * dxdxi[1, 1];
              jac2[1, 4] = 2 * dxdxi[1, 1] * dxdxi[2, 1];
              jac2[1, 5] = 2 * dxdxi[0, 1] * dxdxi[2, 1];
              jac2[2, 0] = dxdxi[0, 2] * dxdxi[0, 2];
              jac2[2, 1] = dxdxi[1, 2] * dxdxi[1, 2];
              jac2[2, 2] = dxdxi[2, 2] * dxdxi[2, 2];
              jac2[2, 3] = 2 * dxdxi[0, 2] * dxdxi[1, 2];
              jac2[2, 4] = 2 * dxdxi[1, 2] * dxdxi[2, 2];
              jac2[2, 5] = 2 * dxdxi[0, 2] * dxdxi[2, 2];
              jac2[3, 0] = dxdxi[0, 0] * dxdxi[0, 1];
              jac2[3, 1] = dxdxi[1, 0] * dxdxi[1, 1];
              jac2[3, 2] = dxdxi[2, 0] * dxdxi[2, 1];
              jac2[3, 3] = dxdxi[1, 0] * dxdxi[0, 1] + dxdxi[0, 0] * dxdxi[1, 1];
              jac2[3, 4] = dxdxi[2, 0] * dxdxi[1, 1] + dxdxi[1, 0] * dxdxi[2, 1];
              jac2[3, 5] = dxdxi[2, 0] * dxdxi[0, 1] + dxdxi[0, 0] * dxdxi[2, 1];
              jac2[4, 0] = dxdxi[0, 1] * dxdxi[0, 2];
              jac2[4, 1] = dxdxi[1, 1] * dxdxi[1, 2];
              jac2[4, 2] = dxdxi[2, 1] * dxdxi[2, 2];
              jac2[4, 3] = dxdxi[1, 1] * dxdxi[0, 2] + dxdxi[0, 1] * dxdxi[1, 2];
              jac2[4, 4] = dxdxi[2, 1] * dxdxi[1, 2] + dxdxi[1, 1] * dxdxi[2, 2];
              jac2[4, 5] = dxdxi[2, 1] * dxdxi[0, 2] + dxdxi[0, 1] * dxdxi[2, 2];
              jac2[5, 0] = dxdxi[0, 0] * dxdxi[0, 2];
              jac2[5, 1] = dxdxi[1, 0] * dxdxi[1, 2];
              jac2[5, 2] = dxdxi[2, 0] * dxdxi[2, 2];
              jac2[5, 3] = dxdxi[1, 0] * dxdxi[0, 2] + dxdxi[0, 0] * dxdxi[1, 2];
              jac2[5, 4] = dxdxi[2, 0] * dxdxi[1, 2] + dxdxi[1, 0] * dxdxi[2, 2];
              jac2[5, 5] = dxdxi[2, 0] * dxdxi[0, 2] + dxdxi[0, 0] * dxdxi[2, 2];

              ddxdxi = Jacobian2At(ddNdxi);
              DoubleMatrix term12 = MatrixFunctions.Product(dNdx, ddxdxi);
              DoubleMatrix term2 = ddNdxi - term12;
              ddNdX = MatrixFunctions.Product(MatrixFunctions.Inverse(jac2), term2.Transpose()).Transpose();
              gpsijk.SetValue(DataInGausspoint.ddNdX, ddNdX);
            }

            //if (matNuy1 != null && matNuy2 != null)
            //{
            //  DoubleMatrix d3xdxi = Jacobian13At(d3Ndxi);
            //  DoubleMatrix term13 = MatrixFunctions.Product(dNdx, d3xdxi);
            //  DoubleMatrix jac23 = new DoubleMatrix(10, 6);

            //  jac23[0, 0] = 3 * dxdxi[0, 0] * ddxdxi[0, 0];
            //  jac23[0, 1] = 3 * dxdxi[1, 0] * ddxdxi[1, 0];
            //  jac23[0, 2] = 3 * dxdxi[2, 0] * ddxdxi[2, 0];
            //  jac23[0, 3] = 3 * (dxdxi[0, 0] * ddxdxi[1, 0] + dxdxi[1, 0] * ddxdxi[0, 0]);
            //  jac23[0, 4] = 3 * (dxdxi[1, 0] * ddxdxi[2, 0] + dxdxi[2, 0] * ddxdxi[1, 0]);
            //  jac23[0, 5] = 3 * (dxdxi[0, 0] * ddxdxi[2, 0] + dxdxi[2, 0] * ddxdxi[0, 0]);
            //  jac23[1, 0] = 3 * dxdxi[0, 1] * ddxdxi[0, 1];
            //  jac23[1, 1] = 3 * dxdxi[1, 1] * ddxdxi[1, 1];
            //  jac23[1, 2] = 3 * dxdxi[2, 1] * ddxdxi[2, 1];
            //  jac23[1, 3] = 3 * (dxdxi[0, 1] * ddxdxi[1, 1] + dxdxi[1, 1] * ddxdxi[0, 1]);
            //  jac23[1, 4] = 3 * (dxdxi[1, 1] * ddxdxi[2, 1] + dxdxi[2, 1] * ddxdxi[1, 1]);
            //  jac23[1, 5] = 3 * (dxdxi[0, 1] * ddxdxi[2, 1] + dxdxi[2, 1] * ddxdxi[0, 1]);
            //  jac23[2, 0] = 3 * dxdxi[0, 2] * ddxdxi[0, 2];
            //  jac23[2, 1] = 3 * dxdxi[1, 2] * ddxdxi[1, 2];
            //  jac23[2, 2] = 3 * dxdxi[2, 2] * ddxdxi[2, 2];
            //  jac23[2, 3] = 3 * (dxdxi[0, 2] * ddxdxi[1, 2] + dxdxi[1, 2] * ddxdxi[0, 2]);
            //  jac23[2, 4] = 3 * (dxdxi[1, 2] * ddxdxi[2, 2] + dxdxi[2, 2] * ddxdxi[1, 2]);
            //  jac23[2, 5] = 3 * (dxdxi[0, 2] * ddxdxi[2, 2] + dxdxi[2, 2] * ddxdxi[0, 2]);
            //  jac23[3, 0] = dxdxi[0, 1] * ddxdxi[0, 0] + 2 * dxdxi[0, 0] * ddxdxi[0, 3];
            //  jac23[3, 1] = dxdxi[1, 1] * ddxdxi[1, 0] + 2 * dxdxi[1, 0] * ddxdxi[1, 3];
            //  jac23[3, 2] = dxdxi[2, 1] * ddxdxi[2, 0] + 2 * dxdxi[2, 0] * ddxdxi[2, 3];
            //  jac23[3, 3] = 2 * (dxdxi[0, 0] * ddxdxi[1, 3] + dxdxi[1, 0] * ddxdxi[0, 3]) + dxdxi[0, 1] * ddxdxi[1, 0] + dxdxi[1, 1] * ddxdxi[0, 0];
            //  jac23[3, 4] = 2 * (dxdxi[1, 0] * ddxdxi[2, 3] + dxdxi[2, 0] * ddxdxi[1, 3]) + dxdxi[1, 1] * ddxdxi[2, 0] + dxdxi[2, 1] * ddxdxi[1, 0];
            //  jac23[3, 5] = 2 * (dxdxi[0, 0] * ddxdxi[2, 3] + dxdxi[2, 0] * ddxdxi[0, 3]) + dxdxi[0, 1] * ddxdxi[2, 0] + dxdxi[2, 1] * ddxdxi[0, 0];
            //  jac23[4, 0] = dxdxi[0, 0] * ddxdxi[0, 1] + 2 * dxdxi[0, 1] * ddxdxi[0, 3];
            //  jac23[4, 1] = dxdxi[1, 0] * ddxdxi[1, 1] + 2 * dxdxi[1, 1] * ddxdxi[1, 3];
            //  jac23[4, 2] = dxdxi[2, 0] * ddxdxi[2, 1] + 2 * dxdxi[2, 1] * ddxdxi[2, 3];
            //  jac23[4, 3] = 2 * (dxdxi[0, 1] * ddxdxi[1, 3] + dxdxi[1, 1] * ddxdxi[0, 3]) + dxdxi[0, 0] * ddxdxi[1, 1] + dxdxi[1, 0] * ddxdxi[0, 1];
            //  jac23[4, 4] = 2 * (dxdxi[1, 1] * ddxdxi[2, 3] + dxdxi[2, 1] * ddxdxi[1, 3]) + dxdxi[1, 0] * ddxdxi[2, 1] + dxdxi[2, 0] * ddxdxi[1, 1];
            //  jac23[4, 5] = 2 * (dxdxi[0, 1] * ddxdxi[2, 3] + dxdxi[2, 1] * ddxdxi[0, 3]) + dxdxi[0, 0] * ddxdxi[2, 1] + dxdxi[2, 0] * ddxdxi[0, 1];
            //  jac23[5, 0] = dxdxi[0, 2] * ddxdxi[0, 1] + 2 * dxdxi[0, 1] * ddxdxi[0, 4];
            //  jac23[5, 1] = dxdxi[1, 2] * ddxdxi[1, 1] + 2 * dxdxi[1, 1] * ddxdxi[1, 4];
            //  jac23[5, 2] = dxdxi[2, 2] * ddxdxi[2, 1] + 2 * dxdxi[2, 1] * ddxdxi[2, 4];
            //  jac23[5, 3] = 2 * (dxdxi[0, 1] * ddxdxi[1, 4] + dxdxi[1, 1] * ddxdxi[0, 4]) + dxdxi[0, 2] * ddxdxi[1, 1] + dxdxi[1, 2] * ddxdxi[0, 1];
            //  jac23[5, 4] = 2 * (dxdxi[1, 1] * ddxdxi[2, 4] + dxdxi[2, 1] * ddxdxi[1, 4]) + dxdxi[1, 2] * ddxdxi[2, 1] + dxdxi[2, 2] * ddxdxi[1, 1];
            //  jac23[5, 5] = 2 * (dxdxi[0, 1] * ddxdxi[2, 4] + dxdxi[2, 1] * ddxdxi[0, 4]) + dxdxi[0, 2] * ddxdxi[2, 1] + dxdxi[2, 2] * ddxdxi[0, 1];
            //  jac23[6, 0] = dxdxi[0, 1] * ddxdxi[0, 2] + 2 * dxdxi[0, 2] * ddxdxi[0, 4];
            //  jac23[6, 1] = dxdxi[1, 1] * ddxdxi[1, 2] + 2 * dxdxi[1, 2] * ddxdxi[1, 4];
            //  jac23[6, 2] = dxdxi[2, 1] * ddxdxi[2, 2] + 2 * dxdxi[2, 2] * ddxdxi[2, 4];
            //  jac23[6, 3] = 2 * (dxdxi[0, 2] * ddxdxi[1, 4] + dxdxi[1, 2] * ddxdxi[0, 4]) + dxdxi[0, 1] * ddxdxi[1, 2] + dxdxi[1, 1] * ddxdxi[0, 2];
            //  jac23[6, 4] = 2 * (dxdxi[1, 2] * ddxdxi[2, 4] + dxdxi[2, 2] * ddxdxi[1, 4]) + dxdxi[1, 1] * ddxdxi[2, 2] + dxdxi[2, 1] * ddxdxi[1, 2];
            //  jac23[6, 5] = 2 * (dxdxi[0, 2] * ddxdxi[2, 4] + dxdxi[2, 2] * ddxdxi[0, 4]) + dxdxi[0, 1] * ddxdxi[2, 2] + dxdxi[2, 1] * ddxdxi[0, 2];
            //  jac23[7, 0] = dxdxi[0, 2] * ddxdxi[0, 0] + 2 * dxdxi[0, 0] * ddxdxi[0, 5];
            //  jac23[7, 1] = dxdxi[1, 2] * ddxdxi[1, 0] + 2 * dxdxi[1, 0] * ddxdxi[1, 5];
            //  jac23[7, 2] = dxdxi[2, 2] * ddxdxi[2, 0] + 2 * dxdxi[2, 0] * ddxdxi[2, 5];
            //  jac23[7, 3] = 2 * (dxdxi[0, 0] * ddxdxi[1, 5] + dxdxi[1, 0] * ddxdxi[0, 5]) + dxdxi[0, 2] * ddxdxi[1, 0] + dxdxi[1, 2] * ddxdxi[0, 0];
            //  jac23[7, 4] = 2 * (dxdxi[1, 0] * ddxdxi[2, 5] + dxdxi[2, 0] * ddxdxi[1, 5]) + dxdxi[1, 2] * ddxdxi[2, 0] + dxdxi[2, 2] * ddxdxi[1, 0];
            //  jac23[7, 5] = 2 * (dxdxi[0, 0] * ddxdxi[2, 5] + dxdxi[2, 0] * ddxdxi[0, 5]) + dxdxi[0, 2] * ddxdxi[2, 0] + dxdxi[2, 2] * ddxdxi[0, 0];
            //  jac23[8, 0] = dxdxi[0, 0] * ddxdxi[0, 2] + 2 * dxdxi[0, 2] * ddxdxi[0, 5];
            //  jac23[8, 1] = dxdxi[1, 0] * ddxdxi[1, 2] + 2 * dxdxi[1, 2] * ddxdxi[1, 5];
            //  jac23[8, 2] = dxdxi[2, 0] * ddxdxi[2, 2] + 2 * dxdxi[2, 2] * ddxdxi[2, 5];
            //  jac23[8, 3] = 2 * (dxdxi[0, 2] * ddxdxi[1, 5] + dxdxi[1, 2] * ddxdxi[0, 5]) + dxdxi[0, 0] * ddxdxi[1, 2] + dxdxi[1, 0] * ddxdxi[0, 2];
            //  jac23[8, 4] = 2 * (dxdxi[1, 2] * ddxdxi[2, 5] + dxdxi[2, 2] * ddxdxi[1, 5]) + dxdxi[1, 0] * ddxdxi[2, 2] + dxdxi[2, 0] * ddxdxi[1, 2];
            //  jac23[8, 5] = 2 * (dxdxi[0, 2] * ddxdxi[2, 5] + dxdxi[2, 2] * ddxdxi[0, 5]) + dxdxi[0, 0] * ddxdxi[2, 2] + dxdxi[2, 0] * ddxdxi[0, 2];
            //  jac23[9, 0] = dxdxi[0, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * ddxdxi[0, 5] + dxdxi[0, 2] * ddxdxi[0, 3];
            //  jac23[9, 1] = dxdxi[1, 0] * ddxdxi[1, 4] + dxdxi[1, 1] * ddxdxi[1, 5] + dxdxi[1, 2] * ddxdxi[1, 3];
            //  jac23[9, 2] = dxdxi[2, 0] * ddxdxi[2, 4] + dxdxi[2, 1] * ddxdxi[2, 5] + dxdxi[2, 2] * ddxdxi[2, 3];
            //  jac23[9, 3] = dxdxi[0, 0] * ddxdxi[1, 4] + dxdxi[0, 1] * ddxdxi[1, 5] + dxdxi[0, 2] * ddxdxi[1, 3] + dxdxi[1, 0] * ddxdxi[0, 4] + dxdxi[1, 1] * ddxdxi[0, 5] + dxdxi[1, 2] * ddxdxi[0, 3];
            //  jac23[9, 4] = dxdxi[1, 0] * ddxdxi[2, 4] + dxdxi[1, 1] * ddxdxi[2, 5] + dxdxi[1, 2] * ddxdxi[2, 3] + dxdxi[2, 0] * ddxdxi[1, 4] + dxdxi[2, 1] * ddxdxi[1, 5] + dxdxi[2, 2] * ddxdxi[1, 3];
            //  jac23[9, 5] = dxdxi[0, 0] * ddxdxi[2, 4] + dxdxi[0, 1] * ddxdxi[2, 5] + dxdxi[0, 2] * ddxdxi[2, 3] + dxdxi[2, 0] * ddxdxi[0, 4] + dxdxi[2, 1] * ddxdxi[0, 5] + dxdxi[2, 2] * ddxdxi[0, 3];
            //  DoubleMatrix term23 = MatrixFunctions.Product(ddNdX, jac23.Transpose());
            //  DoubleMatrix term3 = d3Ndxi - term13 - term23;

            //  DoubleMatrix jac13 = new DoubleMatrix(10, 10);
            //  jac13[0, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0];
            //  jac13[0, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
            //  jac13[0, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac13[0, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0];
            //  jac13[0, 4] = 3 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0];
            //  jac13[0, 5] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0];
            //  jac13[0, 6] = 3 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac13[0, 7] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 0];
            //  jac13[0, 8] = 3 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac13[0, 9] = 6 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 0];
            //  jac13[1, 0] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1];
            //  jac13[1, 1] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
            //  jac13[1, 2] = dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac13[1, 3] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1];
            //  jac13[1, 4] = 3 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1];
            //  jac13[1, 5] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1];
            //  jac13[1, 6] = 3 * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac13[1, 7] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 1];
            //  jac13[1, 8] = 3 * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac13[1, 9] = 6 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 1];
            //  jac13[2, 0] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2];
            //  jac13[2, 1] = dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac13[2, 2] = dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac13[2, 3] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2];
            //  jac13[2, 4] = 3 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac13[2, 5] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2];
            //  jac13[2, 6] = 3 * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac13[2, 7] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 2];
            //  jac13[2, 8] = 3 * dxdxi[0, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac13[2, 9] = 6 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 2];
            //  jac13[3, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1];
            //  jac13[3, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1];
            //  jac13[3, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac13[3, 3] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0];
            //  jac13[3, 4] = dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1];
            //  jac13[3, 5] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0];
            //  jac13[3, 6] = dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac13[3, 7] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 0];
            //  jac13[3, 8] = dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac13[3, 9] = 2 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 0];
            //  jac13[4, 0] = dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1];
            //  jac13[4, 1] = dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1];
            //  jac13[4, 2] = dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac13[4, 3] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1];
            //  jac13[4, 4] = dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1];
            //  jac13[4, 5] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 1];
            //  jac13[4, 6] = dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac13[4, 7] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 1];
            //  jac13[4, 8] = dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac13[4, 9] = 2 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 0];
            //  jac13[5, 0] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2];
            //  jac13[5, 1] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2];
            //  jac13[5, 2] = dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac13[5, 3] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 2] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1];
            //  jac13[5, 4] = dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 2];
            //  jac13[5, 5] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 2] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1];
            //  jac13[5, 6] = dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac13[5, 7] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 2] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 1];
            //  jac13[5, 8] = dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac13[5, 9] = 2 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 1];
            //  jac13[6, 0] = dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2];
            //  jac13[6, 1] = dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac13[6, 2] = dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac13[6, 3] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 2];
            //  jac13[6, 4] = dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 2];
            //  jac13[6, 5] = dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 2];
            //  jac13[6, 6] = dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac13[6, 7] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 2];
            //  jac13[6, 8] = dxdxi[0, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac13[6, 9] = 2 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 1];
            //  jac13[7, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2];
            //  jac13[7, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2];
            //  jac13[7, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac13[7, 3] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 2] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0];
            //  jac13[7, 4] = dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 2];
            //  jac13[7, 5] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 2] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0];
            //  jac13[7, 6] = dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac13[7, 7] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 2] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 0];
            //  jac13[7, 8] = dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac13[7, 9] = 2 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 0];
            //  jac13[8, 0] = dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2];
            //  jac13[8, 1] = dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac13[8, 2] = dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac13[8, 3] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 2];
            //  jac13[8, 4] = dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 2];
            //  jac13[8, 5] = dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 2];
            //  jac13[8, 6] = dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac13[8, 7] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 2];
            //  jac13[8, 8] = dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac13[8, 9] = 2 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 0];
            //  jac13[9, 0] = dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 2];
            //  jac13[9, 1] = dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 2];
            //  jac13[9, 2] = dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac13[9, 3] = dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 2] + dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 0] + dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 1];
            //  jac13[9, 4] = dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 2] + dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 1];
            //  jac13[9, 5] = dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 0] + dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 1];
            //  jac13[9, 6] = dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 1];
            //  jac13[9, 7] = dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 0] + dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 1];
            //  jac13[9, 8] = dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac13[9, 9] = dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 2] + dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 1] + dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 0] + dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 1] + dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 0];

            //  d3NdX = MatrixFunctions.Product(MatrixFunctions.Inverse(jac13), term3.Transpose()).Transpose();
            //  gpsijk.SetValue(DataInGausspoint.d3NdX, d3NdX/*.Transpose()*/);

            //  DoubleMatrix jac14 = new DoubleMatrix(15, 15);
            //  jac14[0, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0];
            //  jac14[0, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
            //  jac14[0, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[0, 3] = 4 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0];
            //  jac14[0, 4] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0];
            //  jac14[0, 5] = 4 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
            //  jac14[0, 6] = 4 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0];
            //  jac14[0, 7] = 6 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[0, 8] = 4 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[0, 9] = 4 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 0];
            //  jac14[0, 10] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[0, 11] = 4 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[0, 12] = 12 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 0];
            //  jac14[0, 13] = 12 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0];
            //  jac14[0, 14] = 12 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[1, 0] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1];
            //  jac14[1, 1] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
            //  jac14[1, 2] = dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[1, 3] = 4 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1];
            //  jac14[1, 4] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1];
            //  jac14[1, 5] = 4 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
            //  jac14[1, 6] = 4 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1];
            //  jac14[1, 7] = 6 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[1, 8] = 4 * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[1, 9] = 4 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 1];
            //  jac14[1, 10] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[1, 11] = 4 * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[1, 12] = 12 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 1];
            //  jac14[1, 13] = 12 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1];
            //  jac14[1, 14] = 12 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[2, 0] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2];
            //  jac14[2, 1] = dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac14[2, 2] = dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[2, 3] = 4 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2];
            //  jac14[2, 4] = 6 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac14[2, 5] = 4 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac14[2, 6] = 4 * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2];
            //  jac14[2, 7] = 6 * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[2, 8] = 4 * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[2, 9] = 4 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 2];
            //  jac14[2, 10] = 6 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[2, 11] = 4 * dxdxi[0, 2] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[2, 12] = 12 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 2];
            //  jac14[2, 13] = 12 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2];
            //  jac14[2, 14] = 12 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2];

            //  jac14[3, 0] = dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0];
            //  jac14[3, 1] = dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
            //  jac14[3, 2] = dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[3, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1];
            //  jac14[3, 4] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1];
            //  jac14[3, 5] = 3 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[1, 0] + dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
            //  jac14[3, 6] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0] + dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 1];
            //  jac14[3, 7] = 3 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] + 3 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac14[3, 8] = 3 * dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[2, 0] + dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[3, 9] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 0] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 1];
            //  jac14[3, 10] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac14[3, 11] = 3 * dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[2, 0] + dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[3, 12] = 6 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 1];
            //  jac14[3, 13] = 6 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0];
            //  jac14[3, 14] = 6 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 1] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 0];

            //  jac14[4, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1];
            //  jac14[4, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1];
            //  jac14[4, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[4, 3] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0];
            //  jac14[4, 4] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0];
            //  jac14[4, 5] = 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1];
            //  jac14[4, 6] = 2 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 1] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0];
            //  jac14[4, 7] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 4 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 1] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[4, 8] = 2 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac14[4, 9] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 0];
            //  jac14[4, 10] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[4, 11] = 2 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac14[4, 12] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 1];
            //  jac14[4, 13] = 2 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0] + 4 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 1] + 4 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0];
            //  jac14[4, 14] = 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 4 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 1] + 4 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 1];

            //  jac14[5, 0] = dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1];
            //  jac14[5, 1] = dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
            //  jac14[5, 2] = dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[5, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0];
            //  jac14[5, 4] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1];
            //  jac14[5, 5] = 3 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] + dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
            //  jac14[5, 6] = 3 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0];
            //  jac14[5, 7] = 3 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1] + 3 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac14[5, 8] = 3 * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1] + dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[5, 9] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 0];
            //  jac14[5, 10] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac14[5, 11] = 3 * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1] + dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[5, 12] = 6 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 0] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 1];
            //  jac14[5, 13] = 6 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0];
            //  jac14[5, 14] = 6 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 1] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 1];

            //  jac14[6, 0] = dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1];
            //  jac14[6, 1] = dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
            //  jac14[6, 2] = dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[6, 3] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 2];
            //  jac14[6, 4] = 3 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 2];
            //  jac14[6, 5] = 3 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[1, 1] + dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
            //  jac14[6, 6] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 2];
            //  jac14[6, 7] = 3 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 3 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac14[6, 8] = 3 * dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[2, 1] + dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[6, 9] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 2];
            //  jac14[6, 10] = 3 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac14[6, 11] = 3 * dxdxi[0, 1] * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[2, 1] + dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[6, 12] = 6 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 2];
            //  jac14[6, 13] = 6 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1];
            //  jac14[6, 14] = 6 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 2] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1];

            //  jac14[7, 0] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2];
            //  jac14[7, 1] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac14[7, 2] = dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[7, 3] = 2 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1];
            //  jac14[7, 4] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 2] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1];
            //  jac14[7, 5] = 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2];
            //  jac14[7, 6] = 2 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1];
            //  jac14[7, 7] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 4 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[7, 8] = 2 * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac14[7, 9] = 2 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 1];
            //  jac14[7, 10] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1];
            //  jac14[7, 11] = 2 * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac14[7, 12] = 2 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 1] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 1] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 2];
            //  jac14[7, 13] = 2 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 2] + 2 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1] + 4 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 2] + 4 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1];
            //  jac14[7, 14] = 2 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 4 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 2] + 4 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 2];

            //  jac14[8, 0] = dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2];
            //  jac14[8, 1] = dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac14[8, 2] = dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[8, 3] = 3 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1];
            //  jac14[8, 4] = 3 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 2];
            //  jac14[8, 5] = 3 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2] + dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac14[8, 6] = 3 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1];
            //  jac14[8, 7] = 3 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2] + 3 * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac14[8, 8] = 3 * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[8, 9] = 3 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 1];
            //  jac14[8, 10] = 3 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac14[8, 11] = 3 * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[8, 12] = 6 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 1] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 2];
            //  jac14[8, 13] = 6 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1];
            //  jac14[8, 14] = 6 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 2] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 2];

            //  jac14[9, 0] = dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0];
            //  jac14[9, 1] = dxdxi[1, 2] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
            //  jac14[9, 2] = dxdxi[2, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[9, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 2];
            //  jac14[9, 4] = 3 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 2];
            //  jac14[9, 5] = 3 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 0] * dxdxi[1, 0] + dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
            //  jac14[9, 6] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0] + dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 2];
            //  jac14[9, 7] = 3 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 3 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac14[9, 8] = 3 * dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 0] * dxdxi[2, 0] + dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[9, 9] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 0] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 2];
            //  jac14[9, 10] = 3 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac14[9, 11] = 3 * dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 0] * dxdxi[2, 0] + dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[9, 12] = 6 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 2];
            //  jac14[9, 13] = 6 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0];
            //  jac14[9, 14] = 6 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 2] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 0];

            //  jac14[10, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2];
            //  jac14[10, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac14[10, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[10, 3] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0];
            //  jac14[10, 4] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0];
            //  jac14[10, 5] = 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2];
            //  jac14[10, 6] = 2 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0];
            //  jac14[10, 7] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 4 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[10, 8] = 2 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac14[10, 9] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 0];
            //  jac14[10, 10] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0];
            //  jac14[10, 11] = 2 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac14[10, 12] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 2];
            //  jac14[10, 13] = 2 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 2] + 2 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0] + 4 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 2] + 4 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0];
            //  jac14[10, 14] = 2 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 4 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 2] + 4 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 2];

            //  jac14[11, 0] = dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2];
            //  jac14[11, 1] = dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac14[11, 2] = dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[11, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0];
            //  jac14[11, 4] = 3 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 2];
            //  jac14[11, 5] = 3 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2] + dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
            //  jac14[11, 6] = 3 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0];
            //  jac14[11, 7] = 3 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2] + 3 * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac14[11, 8] = 3 * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2] + dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[11, 9] = 3 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 0];
            //  jac14[11, 10] = 3 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac14[11, 11] = 3 * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2] + dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
            //  jac14[11, 12] = 6 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 0] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 2];
            //  jac14[11, 13] = 6 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0];
            //  jac14[11, 14] = 6 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 2] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2];

            //  jac14[12, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 2];
            //  jac14[12, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 2];
            //  jac14[12, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac14[12, 3] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 2] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 0];
            //  jac14[12, 4] = 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 2] + dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 1] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 2];
            //  jac14[12, 5] = 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 2] + dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1];
            //  jac14[12, 6] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 2] + dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 1] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 0];
            //  jac14[12, 7] = 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 1] + dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac14[12, 8] = 2 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac14[12, 9] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 2] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 0];
            //  jac14[12, 10] = 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 1] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 2];
            //  jac14[12, 11] = 2 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1];
            //  jac14[12, 12] = 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 2] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 1] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 2];
            //  jac14[12, 13] = 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[0, 0] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[0, 1] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[0, 2] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[0, 0] * dxdxi[2, 1] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[0, 0] * dxdxi[2, 2] + dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[0, 2] * dxdxi[2, 1] + dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[0, 1] * dxdxi[2, 2];
            //  jac14[12, 14] = 2 * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[1, 0] * dxdxi[0, 0] + 2 * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[1, 1] * dxdxi[0, 0] + 2 * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[1, 2] * dxdxi[0, 0] + 2 * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[1, 0] * dxdxi[0, 1] + 2 * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[1, 0] * dxdxi[0, 2] + dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[1, 2] * dxdxi[0, 1] + dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[1, 1] * dxdxi[0, 2];

            //  jac14[13, 0] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[0, 2];
            //  jac14[13, 1] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[1, 2];
            //  jac14[13, 2] = dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac14[13, 3] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 2] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 0] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 1];
            //  jac14[13, 4] = 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 2] + dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 0] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 2];
            //  jac14[13, 5] = 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[1, 2] + dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 0];
            //  jac14[13, 6] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 0] + 2 * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 1];
            //  jac14[13, 7] = 2 * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 0] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac14[13, 8] = 2 * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 0];
            //  jac14[13, 9] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 0] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 1];
            //  jac14[13, 10] = 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 0] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 2];
            //  jac14[13, 11] = 2 * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 0];
            //  jac14[13, 12] = 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 0] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 0] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 2];
            //  jac14[13, 13] = 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[0, 1] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[0, 0] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[0, 2] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[0, 1] * dxdxi[2, 0] + 2 * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[0, 1] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[0, 2] * dxdxi[2, 0] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[0, 0] * dxdxi[2, 2];
            //  jac14[13, 14] = 2 * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[1, 1] * dxdxi[0, 1] + 2 * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[1, 0] * dxdxi[0, 1] + 2 * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[1, 2] * dxdxi[0, 1] + 2 * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[1, 1] * dxdxi[0, 0] + 2 * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[1, 1] * dxdxi[0, 2] + dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[1, 2] * dxdxi[0, 0] + dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[1, 0] * dxdxi[0, 2];

            //  jac14[14, 0] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[0, 0];
            //  jac14[14, 1] = dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[1, 0];
            //  jac14[14, 2] = dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[2, 0];
            //  jac14[14, 3] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[1, 0] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[1, 1] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 2];
            //  jac14[14, 4] = 2 * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 0] + dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 1] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 0];
            //  jac14[14, 5] = 2 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[1, 0] + dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 0] + dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 1];
            //  jac14[14, 6] = dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[2, 0] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[2, 2];
            //  jac14[14, 7] = 2 * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 0] + dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 1] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 0];
            //  jac14[14, 8] = 2 * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[2, 0] + dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 0] + dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 1];
            //  jac14[14, 9] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[2, 0] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[2, 2];
            //  jac14[14, 10] = 2 * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[2, 2] * dxdxi[2, 0] + dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 1] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 0];
            //  jac14[14, 11] = 2 * dxdxi[0, 2] * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[2, 0] + dxdxi[0, 1] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 0] + dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 1];
            //  jac14[14, 12] = 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 1] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 0] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 1] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 0];
            //  jac14[14, 13] = 2 * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[0, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[1, 0] * dxdxi[0, 1] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[0, 0] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[1, 0] * dxdxi[0, 2] * dxdxi[2, 1] + 2 * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[0, 2] * dxdxi[2, 0] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[0, 0] * dxdxi[2, 1] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[0, 1] * dxdxi[2, 0];
            //  jac14[14, 14] = 2 * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[1, 2] * dxdxi[0, 2] + 2 * dxdxi[2, 2] * dxdxi[2, 0] * dxdxi[1, 1] * dxdxi[0, 2] + 2 * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[1, 0] * dxdxi[0, 2] + 2 * dxdxi[2, 2] * dxdxi[2, 0] * dxdxi[1, 2] * dxdxi[0, 1] + 2 * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[1, 2] * dxdxi[0, 0] + dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[1, 0] * dxdxi[0, 1] + dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[1, 1] * dxdxi[0, 0];

            //  DoubleMatrix d4xdxi = Jacobian14At(d4Ndxi);
            //  DoubleMatrix term14 = MatrixFunctions.Product(dNdx, d4xdxi);

            //  DoubleMatrix jac24 = new DoubleMatrix(15, 6);

            //  jac24[0, 0] = 3 * ddxdxi[0, 0] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * d3xdxi[0, 0];
            //  jac24[0, 1] = 3 * ddxdxi[1, 0] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * d3xdxi[1, 0];
            //  jac24[0, 2] = 3 * ddxdxi[2, 0] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * d3xdxi[2, 0];
            //  jac24[0, 3] = 4 * dxdxi[0, 0] * d3xdxi[1, 0] + 4 * dxdxi[1, 0] * d3xdxi[0, 0] + 6 * ddxdxi[0, 0] * ddxdxi[1, 0];
            //  jac24[0, 4] = 4 * dxdxi[1, 0] * d3xdxi[2, 0] + 4 * dxdxi[2, 0] * d3xdxi[1, 0] + 6 * ddxdxi[1, 0] * ddxdxi[2, 0];
            //  jac24[0, 5] = 4 * dxdxi[0, 0] * d3xdxi[2, 0] + 4 * dxdxi[2, 0] * d3xdxi[0, 0] + 6 * ddxdxi[0, 0] * ddxdxi[2, 0];

            //  jac24[1, 0] = 3 * ddxdxi[0, 1] * ddxdxi[0, 1] + 4 * dxdxi[0, 1] * d3xdxi[0, 1];
            //  jac24[1, 1] = 3 * ddxdxi[1, 1] * ddxdxi[1, 1] + 4 * dxdxi[1, 1] * d3xdxi[1, 1];
            //  jac24[1, 2] = 3 * ddxdxi[2, 1] * ddxdxi[2, 1] + 4 * dxdxi[2, 1] * d3xdxi[2, 1];
            //  jac24[1, 3] = 4 * dxdxi[0, 1] * d3xdxi[1, 1] + 4 * dxdxi[1, 1] * d3xdxi[0, 1] + 6 * ddxdxi[0, 1] * ddxdxi[1, 1];
            //  jac24[1, 4] = 4 * dxdxi[1, 1] * d3xdxi[2, 1] + 4 * dxdxi[2, 1] * d3xdxi[1, 1] + 6 * ddxdxi[1, 1] * ddxdxi[2, 1];
            //  jac24[1, 5] = 4 * dxdxi[0, 1] * d3xdxi[2, 1] + 4 * dxdxi[2, 1] * d3xdxi[0, 1] + 6 * ddxdxi[0, 1] * ddxdxi[2, 1];

            //  jac24[2, 0] = 3 * ddxdxi[0, 2] * ddxdxi[0, 2] + 4 * dxdxi[0, 2] * d3xdxi[0, 2];
            //  jac24[2, 1] = 3 * ddxdxi[1, 2] * ddxdxi[1, 2] + 4 * dxdxi[1, 2] * d3xdxi[1, 2];
            //  jac24[2, 2] = 3 * ddxdxi[2, 2] * ddxdxi[2, 2] + 4 * dxdxi[2, 2] * d3xdxi[2, 2];
            //  jac24[2, 3] = 4 * dxdxi[0, 2] * d3xdxi[1, 2] + 4 * dxdxi[1, 2] * d3xdxi[0, 2] + 6 * ddxdxi[0, 2] * ddxdxi[1, 2];
            //  jac24[2, 4] = 4 * dxdxi[1, 2] * d3xdxi[2, 2] + 4 * dxdxi[2, 2] * d3xdxi[1, 2] + 6 * ddxdxi[1, 2] * ddxdxi[2, 2];
            //  jac24[2, 5] = 4 * dxdxi[0, 2] * d3xdxi[2, 2] + 4 * dxdxi[2, 2] * d3xdxi[0, 2] + 6 * ddxdxi[0, 2] * ddxdxi[2, 2];

            //  jac24[3, 0] = dxdxi[0, 1] * d3xdxi[0, 0] + 3 * ddxdxi[0, 0] * ddxdxi[0, 3] + 3 * dxdxi[0, 0] * d3xdxi[0, 3];
            //  jac24[3, 1] = dxdxi[1, 1] * d3xdxi[1, 0] + 3 * ddxdxi[1, 0] * ddxdxi[1, 3] + 3 * dxdxi[1, 0] * d3xdxi[1, 3];
            //  jac24[3, 2] = dxdxi[2, 1] * d3xdxi[2, 0] + 3 * ddxdxi[2, 0] * ddxdxi[2, 3] + 3 * dxdxi[2, 0] * d3xdxi[2, 3];
            //  jac24[3, 3] = 3 * dxdxi[0, 0] * d3xdxi[1, 3] + 3 * dxdxi[1, 0] * d3xdxi[0, 3] + 3 * ddxdxi[0, 3] * ddxdxi[1, 0] + 3 * ddxdxi[0, 0] * ddxdxi[1, 3] + dxdxi[0, 1] * d3xdxi[1, 0] + dxdxi[1, 1] * d3xdxi[0, 0];
            //  jac24[3, 4] = 3 * dxdxi[1, 0] * d3xdxi[2, 3] + 3 * dxdxi[2, 0] * d3xdxi[1, 3] + 3 * ddxdxi[1, 3] * ddxdxi[2, 0] + 3 * ddxdxi[1, 0] * ddxdxi[2, 3] + dxdxi[1, 1] * d3xdxi[2, 0] + dxdxi[2, 1] * d3xdxi[1, 0];
            //  jac24[3, 5] = 3 * dxdxi[0, 0] * d3xdxi[2, 3] + 3 * dxdxi[2, 0] * d3xdxi[0, 3] + 3 * ddxdxi[0, 3] * ddxdxi[2, 0] + 3 * ddxdxi[0, 0] * ddxdxi[2, 3] + dxdxi[0, 1] * d3xdxi[2, 0] + dxdxi[2, 1] * d3xdxi[0, 0];

            //  jac24[4, 0] = ddxdxi[0, 0] * ddxdxi[0, 1] + 2 * ddxdxi[0, 3] * ddxdxi[0, 3] + 2 * dxdxi[0, 0] * d3xdxi[0, 4] + 2 * dxdxi[0, 1] * d3xdxi[0, 3];
            //  jac24[4, 1] = ddxdxi[1, 0] * ddxdxi[1, 1] + 2 * ddxdxi[1, 3] * ddxdxi[1, 3] + 2 * dxdxi[1, 0] * d3xdxi[1, 4] + 2 * dxdxi[1, 1] * d3xdxi[1, 3];
            //  jac24[4, 2] = ddxdxi[2, 0] * ddxdxi[2, 1] + 2 * ddxdxi[2, 3] * ddxdxi[2, 3] + 2 * dxdxi[2, 0] * d3xdxi[2, 4] + 2 * dxdxi[2, 1] * d3xdxi[2, 3];
            //  jac24[4, 3] = 2 * dxdxi[0, 0] * d3xdxi[1, 4] + 2 * dxdxi[0, 1] * d3xdxi[1, 3] + 2 * dxdxi[1, 0] * d3xdxi[0, 4] + 2 * dxdxi[1, 1] * d3xdxi[0, 3] + 4 * ddxdxi[0, 3] * ddxdxi[1, 3] + ddxdxi[0, 0] * ddxdxi[1, 1] + ddxdxi[0, 1] * ddxdxi[1, 0];
            //  jac24[4, 4] = 2 * dxdxi[1, 0] * d3xdxi[2, 4] + 2 * dxdxi[1, 1] * d3xdxi[2, 3] + 2 * dxdxi[2, 0] * d3xdxi[1, 4] + 2 * dxdxi[2, 1] * d3xdxi[1, 3] + 4 * ddxdxi[1, 3] * ddxdxi[2, 3] + ddxdxi[1, 0] * ddxdxi[2, 1] + ddxdxi[1, 1] * ddxdxi[2, 0];
            //  jac24[4, 5] = 2 * dxdxi[0, 0] * d3xdxi[2, 4] + 2 * dxdxi[0, 1] * d3xdxi[2, 3] + 2 * dxdxi[2, 0] * d3xdxi[0, 4] + 2 * dxdxi[2, 1] * d3xdxi[0, 3] + 4 * ddxdxi[0, 3] * ddxdxi[2, 3] + ddxdxi[0, 0] * ddxdxi[2, 1] + ddxdxi[0, 1] * ddxdxi[2, 0];

            //  jac24[5, 0] = dxdxi[0, 0] * d3xdxi[0, 1] + 3 * ddxdxi[0, 1] * ddxdxi[0, 3] + 3 * dxdxi[0, 1] * d3xdxi[0, 4];
            //  jac24[5, 1] = dxdxi[1, 0] * d3xdxi[1, 1] + 3 * ddxdxi[1, 1] * ddxdxi[1, 3] + 3 * dxdxi[1, 1] * d3xdxi[1, 4];
            //  jac24[5, 2] = dxdxi[2, 0] * d3xdxi[2, 1] + 3 * ddxdxi[2, 1] * ddxdxi[2, 3] + 3 * dxdxi[2, 1] * d3xdxi[2, 4];
            //  jac24[5, 3] = 3 * dxdxi[0, 1] * d3xdxi[1, 4] + 3 * dxdxi[1, 1] * d3xdxi[0, 4] + 3 * ddxdxi[0, 1] * ddxdxi[1, 3] + 3 * ddxdxi[0, 3] * ddxdxi[1, 1] + dxdxi[0, 0] * d3xdxi[1, 1] + dxdxi[1, 0] * d3xdxi[0, 1];
            //  jac24[5, 4] = 3 * dxdxi[1, 1] * d3xdxi[2, 4] + 3 * dxdxi[2, 1] * d3xdxi[1, 4] + 3 * ddxdxi[1, 1] * ddxdxi[2, 3] + 3 * ddxdxi[1, 3] * ddxdxi[2, 1] + dxdxi[1, 0] * d3xdxi[2, 1] + dxdxi[2, 0] * d3xdxi[1, 1];
            //  jac24[5, 5] = 3 * dxdxi[0, 1] * d3xdxi[2, 4] + 3 * dxdxi[2, 1] * d3xdxi[0, 4] + 3 * ddxdxi[0, 1] * ddxdxi[2, 3] + 3 * ddxdxi[0, 3] * ddxdxi[2, 1] + dxdxi[0, 0] * d3xdxi[2, 1] + dxdxi[2, 0] * d3xdxi[0, 1];

            //  jac24[6, 0] = dxdxi[0, 2] * d3xdxi[0, 1] + 3 * ddxdxi[0, 1] * ddxdxi[0, 4] + 3 * dxdxi[0, 1] * d3xdxi[0, 5];
            //  jac24[6, 1] = dxdxi[1, 2] * d3xdxi[1, 1] + 3 * ddxdxi[1, 1] * ddxdxi[1, 4] + 3 * dxdxi[1, 1] * d3xdxi[1, 5];
            //  jac24[6, 2] = dxdxi[2, 2] * d3xdxi[2, 1] + 3 * ddxdxi[2, 1] * ddxdxi[2, 4] + 3 * dxdxi[2, 1] * d3xdxi[2, 5];
            //  jac24[6, 3] = 3 * dxdxi[0, 1] * d3xdxi[1, 5] + 3 * dxdxi[1, 1] * d3xdxi[0, 5] + 3 * ddxdxi[0, 4] * ddxdxi[1, 1] + 3 * ddxdxi[0, 1] * ddxdxi[1, 4] + dxdxi[0, 2] * d3xdxi[1, 1] + dxdxi[1, 2] * d3xdxi[0, 1];
            //  jac24[6, 4] = 3 * dxdxi[1, 1] * d3xdxi[2, 5] + 3 * dxdxi[2, 1] * d3xdxi[1, 5] + 3 * ddxdxi[1, 4] * ddxdxi[2, 1] + 3 * ddxdxi[1, 1] * ddxdxi[2, 4] + dxdxi[1, 2] * d3xdxi[2, 1] + dxdxi[2, 2] * d3xdxi[1, 1];
            //  jac24[6, 5] = 3 * dxdxi[0, 1] * d3xdxi[2, 5] + 3 * dxdxi[2, 1] * d3xdxi[0, 5] + 3 * ddxdxi[0, 4] * ddxdxi[2, 1] + 3 * ddxdxi[0, 1] * ddxdxi[2, 4] + dxdxi[0, 2] * d3xdxi[2, 1] + dxdxi[2, 2] * d3xdxi[0, 1];

            //  jac24[7, 0] = ddxdxi[0, 1] * ddxdxi[0, 2] + 2 * ddxdxi[0, 4] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * d3xdxi[0, 6] + 2 * dxdxi[0, 2] * d3xdxi[0, 5];
            //  jac24[7, 1] = ddxdxi[1, 1] * ddxdxi[1, 2] + 2 * ddxdxi[1, 4] * ddxdxi[1, 4] + 2 * dxdxi[1, 1] * d3xdxi[1, 6] + 2 * dxdxi[1, 2] * d3xdxi[1, 5];
            //  jac24[7, 2] = ddxdxi[2, 1] * ddxdxi[2, 2] + 2 * ddxdxi[2, 4] * ddxdxi[2, 4] + 2 * dxdxi[2, 1] * d3xdxi[2, 6] + 2 * dxdxi[2, 2] * d3xdxi[2, 5];
            //  jac24[7, 3] = 2 * dxdxi[0, 1] * d3xdxi[1, 6] + 2 * dxdxi[0, 2] * d3xdxi[1, 5] + 2 * dxdxi[1, 1] * d3xdxi[0, 6] + 2 * dxdxi[1, 2] * d3xdxi[0, 5] + 4 * ddxdxi[0, 4] * ddxdxi[1, 4] + ddxdxi[0, 1] * ddxdxi[1, 2] + ddxdxi[0, 2] * ddxdxi[1, 1];
            //  jac24[7, 4] = 2 * dxdxi[1, 1] * d3xdxi[2, 6] + 2 * dxdxi[1, 2] * d3xdxi[2, 5] + 2 * dxdxi[2, 1] * d3xdxi[1, 6] + 2 * dxdxi[2, 2] * d3xdxi[1, 5] + 4 * ddxdxi[1, 4] * ddxdxi[2, 4] + ddxdxi[1, 1] * ddxdxi[2, 2] + ddxdxi[1, 2] * ddxdxi[2, 1];
            //  jac24[7, 5] = 2 * dxdxi[0, 1] * d3xdxi[2, 6] + 2 * dxdxi[0, 2] * d3xdxi[2, 5] + 2 * dxdxi[2, 1] * d3xdxi[0, 6] + 2 * dxdxi[2, 2] * d3xdxi[0, 5] + 4 * ddxdxi[0, 4] * ddxdxi[2, 4] + ddxdxi[0, 1] * ddxdxi[2, 2] + ddxdxi[0, 2] * ddxdxi[2, 1];

            //  jac24[8, 0] = dxdxi[0, 1] * d3xdxi[0, 2] + 3 * ddxdxi[0, 2] * ddxdxi[0, 4] + 3 * dxdxi[0, 2] * d3xdxi[0, 6];
            //  jac24[8, 1] = dxdxi[1, 1] * d3xdxi[1, 2] + 3 * ddxdxi[1, 2] * ddxdxi[1, 4] + 3 * dxdxi[1, 2] * d3xdxi[1, 6];
            //  jac24[8, 2] = dxdxi[2, 1] * d3xdxi[2, 2] + 3 * ddxdxi[2, 2] * ddxdxi[2, 4] + 3 * dxdxi[2, 2] * d3xdxi[2, 6];
            //  jac24[8, 3] = 3 * dxdxi[0, 2] * d3xdxi[1, 6] + 3 * dxdxi[1, 2] * d3xdxi[0, 6] + 3 * ddxdxi[0, 2] * ddxdxi[1, 4] + 3 * ddxdxi[0, 4] * ddxdxi[1, 2] + dxdxi[0, 1] * d3xdxi[1, 2] + dxdxi[1, 1] * d3xdxi[0, 2];
            //  jac24[8, 4] = 3 * dxdxi[1, 2] * d3xdxi[2, 6] + 3 * dxdxi[2, 2] * d3xdxi[1, 6] + 3 * ddxdxi[1, 2] * ddxdxi[2, 4] + 3 * ddxdxi[1, 4] * ddxdxi[2, 2] + dxdxi[1, 1] * d3xdxi[2, 2] + dxdxi[2, 1] * d3xdxi[1, 2];
            //  jac24[8, 5] = 3 * dxdxi[0, 2] * d3xdxi[2, 6] + 3 * dxdxi[2, 2] * d3xdxi[0, 6] + 3 * ddxdxi[0, 2] * ddxdxi[2, 4] + 3 * ddxdxi[0, 4] * ddxdxi[2, 2] + dxdxi[0, 1] * d3xdxi[2, 2] + dxdxi[2, 1] * d3xdxi[0, 2];

            //  jac24[9, 0] = dxdxi[0, 2] * d3xdxi[0, 0] + 3 * ddxdxi[0, 0] * ddxdxi[0, 5] + 3 * dxdxi[0, 0] * d3xdxi[0, 7];
            //  jac24[9, 1] = dxdxi[1, 2] * d3xdxi[1, 0] + 3 * ddxdxi[1, 0] * ddxdxi[1, 5] + 3 * dxdxi[1, 0] * d3xdxi[1, 7];
            //  jac24[9, 2] = dxdxi[2, 2] * d3xdxi[2, 0] + 3 * ddxdxi[2, 0] * ddxdxi[2, 5] + 3 * dxdxi[2, 0] * d3xdxi[2, 7];
            //  jac24[9, 3] = 3 * dxdxi[0, 0] * d3xdxi[1, 7] + 3 * dxdxi[1, 0] * d3xdxi[0, 7] + 3 * ddxdxi[0, 5] * ddxdxi[1, 0] + 3 * ddxdxi[0, 0] * ddxdxi[1, 5] + dxdxi[0, 2] * d3xdxi[1, 0] + dxdxi[1, 2] * d3xdxi[0, 0];
            //  jac24[9, 4] = 3 * dxdxi[1, 0] * d3xdxi[2, 7] + 3 * dxdxi[2, 0] * d3xdxi[1, 7] + 3 * ddxdxi[1, 5] * ddxdxi[2, 0] + 3 * ddxdxi[1, 0] * ddxdxi[2, 5] + dxdxi[1, 2] * d3xdxi[2, 0] + dxdxi[2, 2] * d3xdxi[1, 0];
            //  jac24[9, 5] = 3 * dxdxi[0, 0] * d3xdxi[2, 7] + 3 * dxdxi[2, 0] * d3xdxi[0, 7] + 3 * ddxdxi[0, 5] * ddxdxi[2, 0] + 3 * ddxdxi[0, 0] * ddxdxi[2, 5] + dxdxi[0, 2] * d3xdxi[2, 0] + dxdxi[2, 2] * d3xdxi[0, 0];

            //  jac24[10, 0] = ddxdxi[0, 0] * ddxdxi[0, 2] + 2 * ddxdxi[0, 5] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * d3xdxi[0, 8] + 2 * dxdxi[0, 2] * d3xdxi[0, 7];
            //  jac24[10, 1] = ddxdxi[1, 0] * ddxdxi[1, 2] + 2 * ddxdxi[1, 5] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * d3xdxi[1, 8] + 2 * dxdxi[1, 2] * d3xdxi[1, 7];
            //  jac24[10, 2] = ddxdxi[2, 0] * ddxdxi[2, 2] + 2 * ddxdxi[2, 5] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * d3xdxi[2, 8] + 2 * dxdxi[2, 2] * d3xdxi[2, 7];
            //  jac24[10, 3] = 2 * dxdxi[0, 0] * d3xdxi[1, 8] + 2 * dxdxi[0, 2] * d3xdxi[1, 7] + 2 * dxdxi[1, 0] * d3xdxi[0, 8] + 2 * dxdxi[1, 2] * d3xdxi[0, 7] + 4 * ddxdxi[0, 5] * ddxdxi[1, 5] + ddxdxi[0, 0] * ddxdxi[1, 2] + ddxdxi[0, 2] * ddxdxi[1, 0];
            //  jac24[10, 4] = 2 * dxdxi[1, 0] * d3xdxi[2, 8] + 2 * dxdxi[1, 2] * d3xdxi[2, 7] + 2 * dxdxi[2, 0] * d3xdxi[1, 8] + 2 * dxdxi[2, 2] * d3xdxi[1, 7] + 4 * ddxdxi[1, 5] * ddxdxi[2, 5] + ddxdxi[1, 0] * ddxdxi[2, 2] + ddxdxi[1, 2] * ddxdxi[2, 0];
            //  jac24[10, 5] = 2 * dxdxi[0, 0] * d3xdxi[2, 8] + 2 * dxdxi[0, 2] * d3xdxi[2, 7] + 2 * dxdxi[2, 0] * d3xdxi[0, 8] + 2 * dxdxi[2, 2] * d3xdxi[0, 7] + 4 * ddxdxi[0, 5] * ddxdxi[2, 5] + ddxdxi[0, 0] * ddxdxi[2, 2] + ddxdxi[0, 2] * ddxdxi[2, 0];

            //  jac24[11, 0] = dxdxi[0, 0] * d3xdxi[0, 2] + 3 * ddxdxi[0, 2] * ddxdxi[0, 5] + 3 * dxdxi[0, 2] * d3xdxi[0, 8];
            //  jac24[11, 1] = dxdxi[1, 0] * d3xdxi[1, 2] + 3 * ddxdxi[1, 2] * ddxdxi[1, 5] + 3 * dxdxi[1, 2] * d3xdxi[1, 8];
            //  jac24[11, 2] = dxdxi[2, 0] * d3xdxi[2, 2] + 3 * ddxdxi[2, 2] * ddxdxi[2, 5] + 3 * dxdxi[2, 2] * d3xdxi[2, 8];
            //  jac24[11, 3] = 3 * dxdxi[0, 2] * d3xdxi[1, 8] + 3 * dxdxi[1, 2] * d3xdxi[0, 8] + 3 * ddxdxi[0, 2] * ddxdxi[1, 5] + 3 * ddxdxi[0, 5] * ddxdxi[1, 2] + dxdxi[0, 0] * d3xdxi[1, 2] + dxdxi[1, 0] * d3xdxi[0, 2];
            //  jac24[11, 4] = 3 * dxdxi[1, 2] * d3xdxi[2, 8] + 3 * dxdxi[2, 2] * d3xdxi[1, 8] + 3 * ddxdxi[1, 2] * ddxdxi[2, 5] + 3 * ddxdxi[1, 5] * ddxdxi[2, 2] + dxdxi[1, 0] * d3xdxi[2, 2] + dxdxi[2, 0] * d3xdxi[1, 2];
            //  jac24[11, 5] = 3 * dxdxi[0, 2] * d3xdxi[2, 8] + 3 * dxdxi[2, 2] * d3xdxi[0, 8] + 3 * ddxdxi[0, 2] * ddxdxi[2, 5] + 3 * ddxdxi[0, 5] * ddxdxi[2, 2] + dxdxi[0, 0] * d3xdxi[2, 2] + dxdxi[2, 0] * d3xdxi[0, 2];

            //  jac24[12, 0] = ddxdxi[0, 0] * ddxdxi[0, 4] + 2 * ddxdxi[0, 3] * ddxdxi[0, 5] + dxdxi[0, 1] * d3xdxi[0, 7] + dxdxi[0, 2] * d3xdxi[0, 3] + 2 * dxdxi[0, 0] * d3xdxi[0, 9];
            //  jac24[12, 1] = ddxdxi[1, 0] * ddxdxi[1, 4] + 2 * ddxdxi[1, 3] * ddxdxi[1, 5] + dxdxi[1, 1] * d3xdxi[1, 7] + dxdxi[1, 2] * d3xdxi[1, 3] + 2 * dxdxi[1, 0] * d3xdxi[1, 9];
            //  jac24[12, 2] = ddxdxi[2, 0] * ddxdxi[2, 4] + 2 * ddxdxi[2, 3] * ddxdxi[2, 5] + dxdxi[2, 1] * d3xdxi[2, 7] + dxdxi[2, 2] * d3xdxi[2, 3] + 2 * dxdxi[2, 0] * d3xdxi[2, 9];
            //  jac24[12, 3] = 2 * dxdxi[1, 0] * d3xdxi[0, 9] + 2 * dxdxi[0, 0] * d3xdxi[1, 9] + ddxdxi[0, 0] * ddxdxi[1, 4] + 2 * ddxdxi[0, 3] * ddxdxi[1, 5] + 2 * ddxdxi[0, 5] * ddxdxi[1, 3] + ddxdxi[1, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * d3xdxi[1, 7] + dxdxi[0, 2] * d3xdxi[1, 3] + dxdxi[1, 1] * d3xdxi[0, 7] + dxdxi[1, 2] * d3xdxi[0, 3];
            //  jac24[12, 4] = 2 * dxdxi[2, 0] * d3xdxi[1, 9] + 2 * dxdxi[1, 0] * d3xdxi[2, 9] + ddxdxi[1, 0] * ddxdxi[2, 4] + 2 * ddxdxi[1, 3] * ddxdxi[2, 5] + 2 * ddxdxi[1, 5] * ddxdxi[2, 3] + ddxdxi[2, 0] * ddxdxi[1, 4] + dxdxi[1, 1] * d3xdxi[2, 7] + dxdxi[1, 2] * d3xdxi[2, 3] + dxdxi[2, 1] * d3xdxi[1, 7] + dxdxi[2, 2] * d3xdxi[1, 3];
            //  jac24[12, 5] = 2 * dxdxi[2, 0] * d3xdxi[0, 9] + 2 * dxdxi[0, 0] * d3xdxi[2, 9] + ddxdxi[0, 0] * ddxdxi[2, 4] + 2 * ddxdxi[0, 3] * ddxdxi[2, 5] + 2 * ddxdxi[0, 5] * ddxdxi[2, 3] + ddxdxi[2, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * d3xdxi[2, 7] + dxdxi[0, 2] * d3xdxi[2, 3] + dxdxi[2, 1] * d3xdxi[0, 7] + dxdxi[2, 2] * d3xdxi[0, 3];

            //  jac24[13, 0] = ddxdxi[0, 1] * ddxdxi[0, 5] + 2 * ddxdxi[0, 3] * ddxdxi[0, 4] + dxdxi[0, 0] * d3xdxi[0, 5] + dxdxi[0, 2] * d3xdxi[0, 4] + 2 * dxdxi[0, 1] * d3xdxi[0, 9];
            //  jac24[13, 1] = ddxdxi[1, 1] * ddxdxi[1, 5] + 2 * ddxdxi[1, 3] * ddxdxi[1, 4] + dxdxi[1, 0] * d3xdxi[1, 5] + dxdxi[1, 2] * d3xdxi[1, 4] + 2 * dxdxi[1, 1] * d3xdxi[1, 9];
            //  jac24[13, 2] = ddxdxi[2, 1] * ddxdxi[2, 5] + 2 * ddxdxi[2, 3] * ddxdxi[2, 4] + dxdxi[2, 0] * d3xdxi[2, 5] + dxdxi[2, 2] * d3xdxi[2, 4] + 2 * dxdxi[2, 1] * d3xdxi[2, 9];
            //  jac24[13, 3] = 2 * dxdxi[1, 1] * d3xdxi[0, 9] + 2 * dxdxi[0, 1] * d3xdxi[1, 9] + ddxdxi[0, 1] * ddxdxi[1, 5] + 2 * ddxdxi[0, 3] * ddxdxi[1, 4] + 2 * ddxdxi[0, 4] * ddxdxi[1, 3] + ddxdxi[1, 1] * ddxdxi[0, 5] + dxdxi[0, 0] * d3xdxi[1, 5] + dxdxi[0, 2] * d3xdxi[1, 4] + dxdxi[1, 0] * d3xdxi[0, 5] + dxdxi[1, 2] * d3xdxi[0, 4];
            //  jac24[13, 4] = 2 * dxdxi[2, 1] * d3xdxi[1, 9] + 2 * dxdxi[1, 1] * d3xdxi[2, 9] + ddxdxi[1, 1] * ddxdxi[2, 5] + 2 * ddxdxi[1, 3] * ddxdxi[2, 4] + 2 * ddxdxi[1, 4] * ddxdxi[2, 3] + ddxdxi[2, 1] * ddxdxi[1, 5] + dxdxi[1, 0] * d3xdxi[2, 5] + dxdxi[1, 2] * d3xdxi[2, 4] + dxdxi[2, 0] * d3xdxi[1, 5] + dxdxi[2, 2] * d3xdxi[1, 4];
            //  jac24[13, 5] = 2 * dxdxi[2, 1] * d3xdxi[0, 9] + 2 * dxdxi[0, 1] * d3xdxi[2, 9] + ddxdxi[0, 1] * ddxdxi[2, 5] + 2 * ddxdxi[0, 3] * ddxdxi[2, 4] + 2 * ddxdxi[0, 4] * ddxdxi[2, 3] + ddxdxi[2, 1] * ddxdxi[0, 5] + dxdxi[0, 0] * d3xdxi[2, 5] + dxdxi[0, 2] * d3xdxi[2, 4] + dxdxi[2, 0] * d3xdxi[0, 5] + dxdxi[2, 2] * d3xdxi[0, 4];

            //  jac24[14, 0] = ddxdxi[0, 2] * ddxdxi[0, 3] + 2 * ddxdxi[0, 4] * ddxdxi[0, 5] + dxdxi[0, 1] * d3xdxi[0, 8] + dxdxi[0, 0] * d3xdxi[0, 6] + 2 * dxdxi[0, 2] * d3xdxi[0, 9];
            //  jac24[14, 1] = ddxdxi[1, 2] * ddxdxi[1, 3] + 2 * ddxdxi[1, 4] * ddxdxi[1, 5] + dxdxi[1, 1] * d3xdxi[1, 8] + dxdxi[1, 0] * d3xdxi[1, 6] + 2 * dxdxi[1, 2] * d3xdxi[1, 9];
            //  jac24[14, 2] = ddxdxi[2, 2] * ddxdxi[2, 3] + 2 * ddxdxi[2, 4] * ddxdxi[2, 5] + dxdxi[2, 1] * d3xdxi[2, 8] + dxdxi[2, 0] * d3xdxi[2, 6] + 2 * dxdxi[2, 2] * d3xdxi[2, 9];

            //  jac24[14, 3] = 2 * dxdxi[1, 2] * d3xdxi[0, 9] + 2 * dxdxi[0, 2] * d3xdxi[1, 9] + ddxdxi[0, 2] * ddxdxi[1, 3] + 2 * ddxdxi[0, 4] * ddxdxi[1, 5] + 2 * ddxdxi[0, 5] * ddxdxi[1, 4] + ddxdxi[1, 2] * ddxdxi[0, 3] + dxdxi[0, 1] * d3xdxi[1, 8] + dxdxi[0, 0] * d3xdxi[1, 6] + dxdxi[1, 1] * d3xdxi[0, 8] + dxdxi[1, 0] * d3xdxi[0, 6];
            //  jac24[14, 4] = 2 * dxdxi[2, 2] * d3xdxi[1, 9] + 2 * dxdxi[1, 2] * d3xdxi[2, 9] + ddxdxi[1, 2] * ddxdxi[2, 3] + 2 * ddxdxi[1, 4] * ddxdxi[2, 5] + 2 * ddxdxi[1, 5] * ddxdxi[2, 4] + ddxdxi[2, 2] * ddxdxi[1, 3] + dxdxi[1, 1] * d3xdxi[2, 8] + dxdxi[1, 0] * d3xdxi[2, 6] + dxdxi[2, 1] * d3xdxi[1, 8] + dxdxi[2, 0] * d3xdxi[1, 6];
            //  jac24[14, 5] = 2 * dxdxi[2, 2] * d3xdxi[0, 9] + 2 * dxdxi[0, 2] * d3xdxi[2, 9] + ddxdxi[0, 2] * ddxdxi[2, 3] + 2 * ddxdxi[0, 4] * ddxdxi[2, 5] + 2 * ddxdxi[0, 5] * ddxdxi[2, 4] + ddxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[0, 1] * d3xdxi[2, 8] + dxdxi[0, 0] * d3xdxi[2, 6] + dxdxi[2, 1] * d3xdxi[0, 8] + dxdxi[2, 0] * d3xdxi[0, 6];

            //  DoubleMatrix term24 = MatrixFunctions.Product(ddNdX, jac24.Transpose());

            //  DoubleMatrix jac34 = new DoubleMatrix(15, 10);

            //  jac34[0, 0] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 0];
            //  jac34[0, 1] = 6 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 0];
            //  jac34[0, 2] = 6 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[2, 0];
            //  jac34[0, 3] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 0] + 12 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 0];
            //  jac34[0, 4] = 6 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 0] + 12 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[1, 0];
            //  jac34[0, 5] = 6 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[2, 0] + 12 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[1, 0];
            //  jac34[0, 6] = 6 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[1, 0] + 12 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[2, 0];
            //  jac34[0, 7] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[2, 0] + 12 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[0, 0];
            //  jac34[0, 8] = 6 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[0, 0] + 12 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[2, 0];
            //  jac34[0, 9] = 12 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[2, 0] + 12 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[1, 0] + 12 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[0, 0];

            //  jac34[1, 0] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 1];
            //  jac34[1, 1] = 6 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 1];
            //  jac34[1, 2] = 6 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[2, 1];
            //  jac34[1, 3] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 1] + 12 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 1];
            //  jac34[1, 4] = 6 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 1] + 12 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[1, 1];
            //  jac34[1, 5] = 6 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[2, 1] + 12 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[1, 1];
            //  jac34[1, 6] = 6 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[1, 1] + 12 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[2, 1];
            //  jac34[1, 7] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[2, 1] + 12 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[0, 1];
            //  jac34[1, 8] = 6 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[0, 1] + 12 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[2, 1];
            //  jac34[1, 9] = 12 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[2, 1] + 12 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[1, 1] + 12 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[0, 1];

            //  jac34[2, 0] = 6 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[0, 2];
            //  jac34[2, 1] = 6 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[1, 2];
            //  jac34[2, 2] = 6 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[2, 2];
            //  jac34[2, 3] = 6 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[1, 2] + 12 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[0, 2];
            //  jac34[2, 4] = 6 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[0, 2] + 12 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[1, 2];
            //  jac34[2, 5] = 6 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[2, 2] + 12 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[1, 2];
            //  jac34[2, 6] = 6 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[1, 2] + 12 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[2, 2];
            //  jac34[2, 7] = 6 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[2, 2] + 12 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[0, 2];
            //  jac34[2, 8] = 6 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[0, 2] + 12 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[2, 2];
            //  jac34[2, 9] = 12 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[2, 2] + 12 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[1, 2] + 12 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[0, 2];

            //  jac34[3, 0] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 3] + 3 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[0, 0];
            //  jac34[3, 1] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 3] + 3 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[1, 0];
            //  jac34[3, 2] = 3 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[2, 3] + 3 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[2, 0];
            //  jac34[3, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 3] + 3 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[1, 0] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 3] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 0];
            //  jac34[3, 4] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 3] + 3 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[0, 0] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[1, 3] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[1, 0] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[1, 0];
            //  jac34[3, 5] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[2, 3] + 3 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[2, 0] + 6 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[1, 3] + 3 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[1, 0] + 3 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[1, 0];
            //  jac34[3, 6] = 3 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[1, 3] + 3 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[1, 0] + 6 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[2, 3] + 3 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[2, 0] + 3 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[2, 0];
            //  jac34[3, 7] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[2, 3] + 3 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[2, 0] + 6 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[0, 3] + 3 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[0, 0] + 3 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[0, 0];
            //  jac34[3, 8] = 3 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[0, 3] + 3 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[0, 0] + 6 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[2, 3] + 3 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[2, 0] + 3 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[2, 0];
            //  jac34[3, 9] = 3 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[0, 0] + 3 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[0, 0] + 3 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[2, 0] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[2, 0] + 6 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[0, 3] + 6 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[1, 3] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[2, 3];

            //  jac34[4, 0] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[0, 3] + dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 1];
            //  jac34[4, 1] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[1, 3] + dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 1];
            //  jac34[4, 2] = dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[2, 3] + dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[2, 1];
            //  jac34[4, 3] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 1] + dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[1, 3] + 4 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 3] + 4 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 3];
            //  jac34[4, 4] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 1] + dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 0] + 2 * dxdxi[1, 0] * dxdxi[0, 0] * ddxdxi[1, 1] + 2 * dxdxi[1, 1] * dxdxi[0, 1] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[0, 3] + 4 * dxdxi[1, 0] * dxdxi[0, 1] * ddxdxi[1, 3] + 4 * dxdxi[1, 1] * dxdxi[0, 0] * ddxdxi[1, 3];
            //  jac34[4, 5] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[2, 1] + dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[1, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[2, 3] + 4 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[1, 3] + 4 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[1, 3];
            //  jac34[4, 6] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[1, 1] + dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[1, 0] + 2 * dxdxi[2, 0] * dxdxi[1, 0] * ddxdxi[2, 1] + 2 * dxdxi[2, 1] * dxdxi[1, 1] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[1, 3] + 4 * dxdxi[2, 0] * dxdxi[1, 1] * ddxdxi[2, 3] + 4 * dxdxi[2, 1] * dxdxi[1, 0] * ddxdxi[2, 3];
            //  jac34[4, 7] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[2, 1] + dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[0, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[2, 3] + 4 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[0, 3] + 4 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[0, 3];
            //  jac34[4, 8] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[0, 1] + dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[0, 0] + 2 * dxdxi[2, 0] * dxdxi[0, 0] * ddxdxi[2, 1] + 2 * dxdxi[2, 1] * dxdxi[0, 1] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[0, 3] + 4 * dxdxi[2, 0] * dxdxi[0, 1] * ddxdxi[2, 3] + 4 * dxdxi[2, 1] * dxdxi[0, 0] * ddxdxi[2, 3];
            //  jac34[4, 9] = 2 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[0, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[0, 0] + 4 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[0, 3] + 4 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[0, 3] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[1, 0] + 4 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[1, 3] + 4 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[1, 3] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[2, 0] + 4 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[2, 3] + 4 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[2, 3] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[2, 1];

            //  jac34[5, 0] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 3] + 3 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[0, 1];
            //  jac34[5, 1] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 3] + 3 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[1, 1];
            //  jac34[5, 2] = 3 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[2, 3] + 3 * dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[2, 1];
            //  jac34[5, 3] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 3] + 3 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[1, 1] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 3] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 1];
            //  jac34[5, 4] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 3] + 3 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[0, 1] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[1, 3] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[1, 1] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[1, 1];
            //  jac34[5, 5] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[2, 3] + 3 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[2, 1] + 6 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[1, 3] + 3 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[1, 1] + 3 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[1, 1];
            //  jac34[5, 6] = 3 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[1, 3] + 3 * dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[1, 1] + 6 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[2, 3] + 3 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[2, 1] + 3 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[2, 1];
            //  jac34[5, 7] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[2, 3] + 3 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[2, 1] + 6 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[0, 3] + 3 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[0, 1] + 3 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[0, 1];
            //  jac34[5, 8] = 3 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[0, 3] + 3 * dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[0, 1] + 6 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[2, 3] + 3 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[2, 1] + 3 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[2, 1];
            //  jac34[5, 9] = 3 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[0, 1] + 3 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[0, 1] + 3 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[2, 1] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[2, 1] + 6 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[0, 3] + 6 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[1, 3] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[2, 3];

            //  jac34[6, 0] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 4] + 3 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[0, 1];
            //  jac34[6, 1] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 4] + 3 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[1, 1];
            //  jac34[6, 2] = 3 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[2, 4] + 3 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[2, 1];
            //  jac34[6, 3] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 4] + 3 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[1, 1] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 4] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[0, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[0, 1];
            //  jac34[6, 4] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 4] + 3 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[0, 1] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[1, 4] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[1, 1] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[1, 1];
            //  jac34[6, 5] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[2, 4] + 3 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[2, 1] + 6 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[1, 4] + 3 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[1, 1] + 3 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[1, 1];
            //  jac34[6, 6] = 3 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[1, 4] + 3 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[1, 1] + 6 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[2, 4] + 3 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[2, 1] + 3 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[2, 1];
            //  jac34[6, 7] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[2, 4] + 3 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[2, 1] + 6 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[0, 4] + 3 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[0, 1] + 3 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[0, 1];
            //  jac34[6, 8] = 3 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[0, 4] + 3 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[0, 1] + 6 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[2, 4] + 3 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[2, 1] + 3 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[2, 1];
            //  jac34[6, 9] = 3 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[0, 1] + 3 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[0, 1] + 3 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[2, 1] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[2, 1] + 6 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[0, 4] + 6 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[1, 4] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[2, 4];

            //  jac34[7, 0] = dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[0, 1] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[0, 4] + dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 2];
            //  jac34[7, 1] = dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[1, 1] + 4 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[1, 4] + dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 2];
            //  jac34[7, 2] = dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[2, 1] + 4 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[2, 4] + dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[2, 2];
            //  jac34[7, 3] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 2] + dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[0, 1] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[1, 4] + 4 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[0, 4] + 4 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[0, 4];
            //  jac34[7, 4] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 2] + dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[0, 1] + 2 * dxdxi[1, 1] * dxdxi[0, 1] * ddxdxi[1, 2] + 2 * dxdxi[1, 2] * dxdxi[0, 2] * ddxdxi[1, 1] + 4 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[0, 4] + 4 * dxdxi[1, 1] * dxdxi[0, 2] * ddxdxi[1, 4] + 4 * dxdxi[1, 2] * dxdxi[0, 1] * ddxdxi[1, 4];
            //  jac34[7, 5] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[2, 2] + dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[1, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[1, 1] + 4 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[2, 4] + 4 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[1, 4] + 4 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[1, 4];
            //  jac34[7, 6] = dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[1, 2] + dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[1, 1] + 2 * dxdxi[2, 1] * dxdxi[1, 1] * ddxdxi[2, 2] + 2 * dxdxi[2, 2] * dxdxi[1, 2] * ddxdxi[2, 1] + 4 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[1, 4] + 4 * dxdxi[2, 1] * dxdxi[1, 2] * ddxdxi[2, 4] + 4 * dxdxi[2, 2] * dxdxi[1, 1] * ddxdxi[2, 4];
            //  jac34[7, 7] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[2, 2] + dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[0, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[0, 1] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[2, 4] + 4 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[0, 4] + 4 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[0, 4];
            //  jac34[7, 8] = dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[0, 2] + dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[0, 1] + 2 * dxdxi[2, 1] * dxdxi[0, 1] * ddxdxi[2, 2] + 2 * dxdxi[2, 2] * dxdxi[0, 2] * ddxdxi[2, 1] + 4 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[0, 4] + 4 * dxdxi[2, 1] * dxdxi[0, 2] * ddxdxi[2, 4] + 4 * dxdxi[2, 2] * dxdxi[0, 1] * ddxdxi[2, 4];
            //  jac34[7, 9] = 2 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[0, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[0, 1] + 4 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[0, 4] + 4 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[1, 1] + 4 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[1, 4] + 4 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[1, 4] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[2, 1] + 4 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[2, 4] + 4 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[2, 4] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[2, 2];

            //  jac34[8, 0] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[0, 4] + 3 * dxdxi[0, 2] * dxdxi[0, 1] * ddxdxi[0, 2];
            //  jac34[8, 1] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[1, 4] + 3 * dxdxi[1, 2] * dxdxi[1, 1] * ddxdxi[1, 2];
            //  jac34[8, 2] = 3 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[2, 4] + 3 * dxdxi[2, 2] * dxdxi[2, 1] * ddxdxi[2, 2];
            //  jac34[8, 3] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[1, 4] + 3 * dxdxi[0, 2] * dxdxi[0, 1] * ddxdxi[1, 2] + 6 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[0, 4] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[0, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[0, 2];
            //  jac34[8, 4] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[0, 4] + 3 * dxdxi[1, 2] * dxdxi[1, 1] * ddxdxi[0, 2] + 6 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[1, 4] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[1, 2] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[1, 2];
            //  jac34[8, 5] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[2, 4] + 3 * dxdxi[1, 2] * dxdxi[1, 1] * ddxdxi[2, 2] + 6 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[1, 4] + 3 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[1, 2] + 3 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[1, 2];
            //  jac34[8, 6] = 3 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[1, 4] + 3 * dxdxi[2, 2] * dxdxi[2, 1] * ddxdxi[1, 2] + 6 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[2, 4] + 3 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[2, 2] + 3 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[2, 2];
            //  jac34[8, 7] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[2, 4] + 3 * dxdxi[0, 2] * dxdxi[0, 1] * ddxdxi[2, 2] + 6 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[0, 4] + 3 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[0, 2] + 3 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[0, 2];
            //  jac34[8, 8] = 3 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[0, 4] + 3 * dxdxi[2, 2] * dxdxi[2, 1] * ddxdxi[0, 2] + 6 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[2, 4] + 3 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[2, 2] + 3 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[2, 2];
            //  jac34[8, 9] = 3 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[0, 2] + 3 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[0, 2] + 3 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[1, 2] + 3 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[1, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[2, 2] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[2, 2] + 6 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[0, 4] + 6 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[1, 4] + 6 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[2, 4];

            //  jac34[9, 0] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 5] + 3 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[0, 0];
            //  jac34[9, 1] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 5] + 3 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[1, 0];
            //  jac34[9, 2] = 3 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[2, 5] + 3 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[2, 0];
            //  jac34[9, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 5] + 3 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[1, 0] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 5] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[0, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[0, 0];
            //  jac34[9, 4] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 5] + 3 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[0, 0] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[1, 5] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[1, 0] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[1, 0];
            //  jac34[9, 5] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[2, 5] + 3 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[2, 0] + 6 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[1, 5] + 3 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[1, 0] + 3 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[1, 0];
            //  jac34[9, 6] = 3 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[1, 5] + 3 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[1, 0] + 6 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[2, 5] + 3 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[2, 0] + 3 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[2, 0];
            //  jac34[9, 7] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[2, 5] + 3 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[2, 0] + 6 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[0, 5] + 3 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[0, 0] + 3 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[0, 0];
            //  jac34[9, 8] = 3 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[0, 5] + 3 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[0, 0] + 6 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[2, 5] + 3 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[2, 0] + 3 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[2, 0];
            //  jac34[9, 9] = 3 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[0, 0] + 3 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[0, 0] + 3 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[2, 0] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[2, 0] + 6 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[0, 5] + 6 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[1, 5] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[2, 5];

            //  jac34[10, 0] = dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[0, 5] + dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 2];
            //  jac34[10, 1] = dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[1, 5] + dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 2];
            //  jac34[10, 2] = dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[2, 5] + dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[2, 2];
            //  jac34[10, 3] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 2] + dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[1, 5] + 4 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[0, 5] + 4 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[0, 5];
            //  jac34[10, 4] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 2] + dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[0, 0] + 2 * dxdxi[1, 0] * dxdxi[0, 0] * ddxdxi[1, 2] + 2 * dxdxi[1, 2] * dxdxi[0, 2] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[0, 5] + 4 * dxdxi[1, 0] * dxdxi[0, 2] * ddxdxi[1, 5] + 4 * dxdxi[1, 2] * dxdxi[0, 0] * ddxdxi[1, 5];
            //  jac34[10, 5] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[2, 2] + dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[1, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[2, 5] + 4 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[1, 5] + 4 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[1, 5];
            //  jac34[10, 6] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[1, 2] + dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[1, 0] + 2 * dxdxi[2, 0] * dxdxi[1, 0] * ddxdxi[2, 2] + 2 * dxdxi[2, 2] * dxdxi[1, 2] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[1, 5] + 4 * dxdxi[2, 0] * dxdxi[1, 2] * ddxdxi[2, 5] + 4 * dxdxi[2, 2] * dxdxi[1, 0] * ddxdxi[2, 5];
            //  jac34[10, 7] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[2, 2] + dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[0, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[2, 5] + 4 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[0, 5] + 4 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[0, 5];
            //  jac34[10, 8] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[0, 2] + dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[0, 0] + 2 * dxdxi[2, 0] * dxdxi[0, 0] * ddxdxi[2, 2] + 2 * dxdxi[2, 2] * dxdxi[0, 2] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[0, 5] + 4 * dxdxi[2, 0] * dxdxi[0, 2] * ddxdxi[2, 5] + 4 * dxdxi[2, 2] * dxdxi[0, 0] * ddxdxi[2, 5];
            //  jac34[10, 9] = 2 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[0, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[0, 0] + 4 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[0, 5] + 4 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[1, 0] + 4 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[1, 5] + 4 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[1, 5] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[2, 0] + 4 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[2, 5] + 4 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[2, 5] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[2, 2];

            //  jac34[11, 0] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[0, 5] + 3 * dxdxi[0, 2] * dxdxi[0, 0] * ddxdxi[0, 2];
            //  jac34[11, 1] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[1, 5] + 3 * dxdxi[1, 2] * dxdxi[1, 0] * ddxdxi[1, 2];
            //  jac34[11, 2] = 3 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[2, 5] + 3 * dxdxi[2, 2] * dxdxi[2, 0] * ddxdxi[2, 2];
            //  jac34[11, 3] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[1, 5] + 3 * dxdxi[0, 2] * dxdxi[0, 0] * ddxdxi[1, 2] + 6 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[0, 5] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[0, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[0, 2];
            //  jac34[11, 4] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[0, 5] + 3 * dxdxi[1, 2] * dxdxi[1, 0] * ddxdxi[0, 2] + 6 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[1, 5] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[1, 2] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[1, 2];
            //  jac34[11, 5] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[2, 5] + 3 * dxdxi[1, 2] * dxdxi[1, 0] * ddxdxi[2, 2] + 6 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[1, 5] + 3 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[1, 2] + 3 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[1, 2];
            //  jac34[11, 6] = 3 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[1, 5] + 3 * dxdxi[2, 2] * dxdxi[2, 0] * ddxdxi[1, 2] + 6 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[2, 5] + 3 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[2, 2] + 3 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[2, 2];
            //  jac34[11, 7] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[2, 5] + 3 * dxdxi[0, 2] * dxdxi[0, 0] * ddxdxi[2, 2] + 6 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[0, 5] + 3 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[0, 2] + 3 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[0, 2];
            //  jac34[11, 8] = 3 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[0, 5] + 3 * dxdxi[2, 2] * dxdxi[2, 0] * ddxdxi[0, 2] + 6 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[2, 5] + 3 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[2, 2] + 3 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[2, 2];
            //  jac34[11, 9] = 3 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[0, 2] + 3 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[0, 2] + 3 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[1, 2] + 3 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[1, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[2, 2] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[2, 2] + 6 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[0, 5] + 6 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[1, 5] + 6 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[2, 5];

            //  jac34[12, 0] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 4] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[0, 3] + dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[0, 0];
            //  jac34[12, 1] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[1, 3] + dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[1, 0];
            //  jac34[12, 2] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[2, 4] + 2 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[2, 3] + dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[2, 0];
            //  jac34[12, 3] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 4] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[1, 5] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[1, 3] + dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 5] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[0, 3] + 2 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[0, 3] + dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[0, 0] + dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[0, 0];
            //  jac34[12, 4] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 4] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[0, 5] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[0, 3] + dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[0, 0] + 2 * dxdxi[1, 0] * dxdxi[0, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 1] * dxdxi[0, 0] * ddxdxi[1, 5] + 2 * dxdxi[1, 2] * dxdxi[0, 0] * ddxdxi[1, 3] + 2 * dxdxi[1, 0] * dxdxi[0, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[0, 2] * ddxdxi[1, 3] + dxdxi[1, 1] * dxdxi[0, 2] * ddxdxi[1, 0] + dxdxi[1, 2] * dxdxi[0, 1] * ddxdxi[1, 0];
            //  jac34[12, 5] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[2, 3] + dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[1, 5] + 2 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[1, 3] + 2 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[1, 3] + dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[1, 0] + dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[1, 0];
            //  jac34[12, 6] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[1, 3] + dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[1, 0] + 2 * dxdxi[2, 0] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[2, 1] * dxdxi[1, 0] * ddxdxi[2, 5] + 2 * dxdxi[2, 2] * dxdxi[1, 0] * ddxdxi[2, 3] + 2 * dxdxi[2, 0] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[1, 2] * ddxdxi[2, 3] + dxdxi[2, 1] * dxdxi[1, 2] * ddxdxi[2, 0] + dxdxi[2, 2] * dxdxi[1, 1] * ddxdxi[2, 0];
            //  jac34[12, 7] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[2, 4] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[2, 5] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[2, 3] + dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[0, 5] + 2 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[0, 3] + 2 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[0, 0] + dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[0, 0];
            //  jac34[12, 8] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[0, 4] + 2 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[0, 0] + 2 * dxdxi[2, 0] * dxdxi[0, 0] * ddxdxi[2, 4] + 2 * dxdxi[2, 1] * dxdxi[0, 0] * ddxdxi[2, 5] + 2 * dxdxi[2, 2] * dxdxi[0, 0] * ddxdxi[2, 3] + 2 * dxdxi[2, 0] * dxdxi[0, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[0, 2] * ddxdxi[2, 3] + dxdxi[2, 1] * dxdxi[0, 2] * ddxdxi[2, 0] + dxdxi[2, 2] * dxdxi[0, 1] * ddxdxi[2, 0];
            //  jac34[12, 9] = 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[2, 3] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[2, 5] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[2, 3] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[1, 5] + 2 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[1, 3] + 2 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[1, 3] + 2 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[0, 4] + 2 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[0, 5] + 2 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[0, 3] + 2 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[2, 0] + dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[2, 0] + dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[1, 0] + dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[1, 0] + dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[0, 0] + dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[0, 0];

            //  jac34[13, 0] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[0, 3] + dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[0, 1];
            //  jac34[13, 1] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[1, 3] + dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[1, 1];
            //  jac34[13, 2] = dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[2, 4] + 2 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[2, 3] + dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[2, 1];
            //  jac34[13, 3] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 5] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[1, 4] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[1, 3] + dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 4] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[0, 3] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[0, 3] + dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[0, 1] + dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[0, 1];
            //  jac34[13, 4] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 5] + 2 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[0, 4] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[0, 3] + dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[0, 1] + 2 * dxdxi[1, 1] * dxdxi[0, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[0, 1] * ddxdxi[1, 4] + 2 * dxdxi[1, 2] * dxdxi[0, 1] * ddxdxi[1, 3] + 2 * dxdxi[1, 1] * dxdxi[0, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 1] * dxdxi[0, 2] * ddxdxi[1, 3] + dxdxi[1, 0] * dxdxi[0, 2] * ddxdxi[1, 1] + dxdxi[1, 2] * dxdxi[0, 0] * ddxdxi[1, 1];
            //  jac34[13, 5] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[2, 3] + dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[1, 4] + 2 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[1, 3] + 2 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[1, 3] + dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[1, 1] + dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[1, 1];
            //  jac34[13, 6] = dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[1, 3] + dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[1, 1] + 2 * dxdxi[2, 1] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[1, 1] * ddxdxi[2, 4] + 2 * dxdxi[2, 2] * dxdxi[1, 1] * ddxdxi[2, 3] + 2 * dxdxi[2, 1] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[2, 1] * dxdxi[1, 2] * ddxdxi[2, 3] + dxdxi[2, 0] * dxdxi[1, 2] * ddxdxi[2, 1] + dxdxi[2, 2] * dxdxi[1, 0] * ddxdxi[2, 1];
            //  jac34[13, 7] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[2, 5] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[2, 4] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[2, 3] + dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[0, 4] + 2 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[0, 3] + 2 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[0, 1] + dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[0, 1];
            //  jac34[13, 8] = dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[0, 4] + 2 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[0, 1] + 2 * dxdxi[2, 1] * dxdxi[0, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[0, 1] * ddxdxi[2, 4] + 2 * dxdxi[2, 2] * dxdxi[0, 1] * ddxdxi[2, 3] + 2 * dxdxi[2, 1] * dxdxi[0, 0] * ddxdxi[2, 4] + 2 * dxdxi[2, 1] * dxdxi[0, 2] * ddxdxi[2, 3] + dxdxi[2, 0] * dxdxi[0, 2] * ddxdxi[2, 1] + dxdxi[2, 2] * dxdxi[0, 0] * ddxdxi[2, 1];
            //  jac34[13, 9] = 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[2, 3] + 2 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[2, 4] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[2, 3] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[1, 4] + 2 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[1, 3] + 2 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[1, 3] + 2 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[0, 4] + 2 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[0, 3] + 2 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[0, 4] + 2 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[2, 1] + dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[2, 1] + dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[1, 1] + dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[1, 1] + dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[0, 1] + dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[0, 1];

            //  jac34[14, 0] = dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[0, 3] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[0, 2];
            //  jac34[14, 1] = dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[1, 3] + 2 * dxdxi[1, 2] * dxdxi[1, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 2] * dxdxi[1, 0] * ddxdxi[1, 4] + dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[1, 2];
            //  jac34[14, 2] = dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[2, 3] + 2 * dxdxi[2, 2] * dxdxi[2, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 2] * dxdxi[2, 0] * ddxdxi[2, 4] + dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[2, 2];
            //  jac34[14, 3] = dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[1, 3] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * ddxdxi[1, 5] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * ddxdxi[1, 4] + dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[0, 3] + 2 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[0, 4] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 2] + dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 2];
            //  jac34[14, 4] = dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[0, 3] + 2 * dxdxi[1, 2] * dxdxi[1, 1] * ddxdxi[0, 5] + 2 * dxdxi[1, 2] * dxdxi[1, 0] * ddxdxi[0, 4] + dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[0, 2] + 2 * dxdxi[1, 2] * dxdxi[0, 2] * ddxdxi[1, 3] + 2 * dxdxi[1, 1] * dxdxi[0, 2] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[0, 2] * ddxdxi[1, 4] + 2 * dxdxi[1, 2] * dxdxi[0, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 2] * dxdxi[0, 0] * ddxdxi[1, 4] + dxdxi[1, 1] * dxdxi[0, 0] * ddxdxi[1, 2] + dxdxi[1, 0] * dxdxi[0, 1] * ddxdxi[1, 2];
            //  jac34[14, 5] = dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[2, 3] + 2 * dxdxi[1, 2] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[1, 2] * dxdxi[1, 0] * ddxdxi[2, 4] + dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[1, 3] + 2 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[1, 4] + 2 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[1, 4] + dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[1, 2] + dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[1, 2];
            //  jac34[14, 6] = dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[1, 3] + 2 * dxdxi[2, 2] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[2, 2] * dxdxi[2, 0] * ddxdxi[1, 4] + dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[1, 2] + 2 * dxdxi[2, 2] * dxdxi[1, 2] * ddxdxi[2, 3] + 2 * dxdxi[2, 1] * dxdxi[1, 2] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[1, 2] * ddxdxi[2, 4] + 2 * dxdxi[2, 2] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 2] * dxdxi[1, 0] * ddxdxi[2, 4] + dxdxi[2, 1] * dxdxi[1, 0] * ddxdxi[2, 2] + dxdxi[2, 0] * dxdxi[1, 1] * ddxdxi[2, 2];
            //  jac34[14, 7] = dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[2, 3] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * ddxdxi[2, 5] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * ddxdxi[2, 4] + dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[0, 3] + 2 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[0, 4] + 2 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[0, 2] + dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[0, 2];
            //  jac34[14, 8] = dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[0, 3] + 2 * dxdxi[2, 2] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[2, 2] * dxdxi[2, 0] * ddxdxi[0, 4] + dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[0, 2] + 2 * dxdxi[2, 2] * dxdxi[0, 2] * ddxdxi[2, 3] + 2 * dxdxi[2, 1] * dxdxi[0, 2] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[0, 2] * ddxdxi[2, 4] + 2 * dxdxi[2, 2] * dxdxi[0, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 2] * dxdxi[0, 0] * ddxdxi[2, 4] + dxdxi[2, 1] * dxdxi[0, 0] * ddxdxi[2, 2] + dxdxi[2, 0] * dxdxi[0, 1] * ddxdxi[2, 2];
            //  jac34[14, 9] = 2 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[2, 3] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[2, 5] + 2 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[2, 4] + 2 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[1, 3] + 2 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[1, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[1, 4] + 2 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[0, 3] + 2 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[0, 5] + 2 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[0, 4] + 2 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[2, 2] + dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[2, 2] + dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[1, 2] + dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[1, 2] + dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[0, 2] + dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[0, 2];

            //  DoubleMatrix term34 = MatrixFunctions.Product(d3NdX, jac34.Transpose());

            //  DoubleMatrix term4 = d4Ndxi - term14 - term24 - term34;
            //  d4NdX = MatrixFunctions.Product(MatrixFunctions.Inverse(jac14), term4.Transpose()).Transpose();

            //  gpsijk.SetValue(DataInGausspoint.d4NdX, d4NdX);
            //}
            if (isNeed_d3NdX)
            {
              d3xdxi = Jacobian13At(d3Ndxi);
              DoubleMatrix term13 = MatrixFunctions.Product(dNdx, d3xdxi);
              DoubleMatrix jac23 = new DoubleMatrix(10, 6);

              jac23[0, 0] = 3 * dxdxi[0, 0] * ddxdxi[0, 0];
              jac23[0, 1] = 3 * dxdxi[1, 0] * ddxdxi[1, 0];
              jac23[0, 2] = 3 * dxdxi[2, 0] * ddxdxi[2, 0];
              jac23[0, 3] = 3 * (dxdxi[0, 0] * ddxdxi[1, 0] + dxdxi[1, 0] * ddxdxi[0, 0]);
              jac23[0, 4] = 3 * (dxdxi[1, 0] * ddxdxi[2, 0] + dxdxi[2, 0] * ddxdxi[1, 0]);
              jac23[0, 5] = 3 * (dxdxi[0, 0] * ddxdxi[2, 0] + dxdxi[2, 0] * ddxdxi[0, 0]);
              jac23[1, 0] = 3 * dxdxi[0, 1] * ddxdxi[0, 1];
              jac23[1, 1] = 3 * dxdxi[1, 1] * ddxdxi[1, 1];
              jac23[1, 2] = 3 * dxdxi[2, 1] * ddxdxi[2, 1];
              jac23[1, 3] = 3 * (dxdxi[0, 1] * ddxdxi[1, 1] + dxdxi[1, 1] * ddxdxi[0, 1]);
              jac23[1, 4] = 3 * (dxdxi[1, 1] * ddxdxi[2, 1] + dxdxi[2, 1] * ddxdxi[1, 1]);
              jac23[1, 5] = 3 * (dxdxi[0, 1] * ddxdxi[2, 1] + dxdxi[2, 1] * ddxdxi[0, 1]);
              jac23[2, 0] = 3 * dxdxi[0, 2] * ddxdxi[0, 2];
              jac23[2, 1] = 3 * dxdxi[1, 2] * ddxdxi[1, 2];
              jac23[2, 2] = 3 * dxdxi[2, 2] * ddxdxi[2, 2];
              jac23[2, 3] = 3 * (dxdxi[0, 2] * ddxdxi[1, 2] + dxdxi[1, 2] * ddxdxi[0, 2]);
              jac23[2, 4] = 3 * (dxdxi[1, 2] * ddxdxi[2, 2] + dxdxi[2, 2] * ddxdxi[1, 2]);
              jac23[2, 5] = 3 * (dxdxi[0, 2] * ddxdxi[2, 2] + dxdxi[2, 2] * ddxdxi[0, 2]);
              jac23[3, 0] = dxdxi[0, 1] * ddxdxi[0, 0] + 2 * dxdxi[0, 0] * ddxdxi[0, 3];
              jac23[3, 1] = dxdxi[1, 1] * ddxdxi[1, 0] + 2 * dxdxi[1, 0] * ddxdxi[1, 3];
              jac23[3, 2] = dxdxi[2, 1] * ddxdxi[2, 0] + 2 * dxdxi[2, 0] * ddxdxi[2, 3];
              jac23[3, 3] = 2 * (dxdxi[0, 0] * ddxdxi[1, 3] + dxdxi[1, 0] * ddxdxi[0, 3]) + dxdxi[0, 1] * ddxdxi[1, 0] + dxdxi[1, 1] * ddxdxi[0, 0];
              jac23[3, 4] = 2 * (dxdxi[1, 0] * ddxdxi[2, 3] + dxdxi[2, 0] * ddxdxi[1, 3]) + dxdxi[1, 1] * ddxdxi[2, 0] + dxdxi[2, 1] * ddxdxi[1, 0];
              jac23[3, 5] = 2 * (dxdxi[0, 0] * ddxdxi[2, 3] + dxdxi[2, 0] * ddxdxi[0, 3]) + dxdxi[0, 1] * ddxdxi[2, 0] + dxdxi[2, 1] * ddxdxi[0, 0];
              jac23[4, 0] = dxdxi[0, 0] * ddxdxi[0, 1] + 2 * dxdxi[0, 1] * ddxdxi[0, 3];
              jac23[4, 1] = dxdxi[1, 0] * ddxdxi[1, 1] + 2 * dxdxi[1, 1] * ddxdxi[1, 3];
              jac23[4, 2] = dxdxi[2, 0] * ddxdxi[2, 1] + 2 * dxdxi[2, 1] * ddxdxi[2, 3];
              jac23[4, 3] = 2 * (dxdxi[0, 1] * ddxdxi[1, 3] + dxdxi[1, 1] * ddxdxi[0, 3]) + dxdxi[0, 0] * ddxdxi[1, 1] + dxdxi[1, 0] * ddxdxi[0, 1];
              jac23[4, 4] = 2 * (dxdxi[1, 1] * ddxdxi[2, 3] + dxdxi[2, 1] * ddxdxi[1, 3]) + dxdxi[1, 0] * ddxdxi[2, 1] + dxdxi[2, 0] * ddxdxi[1, 1];
              jac23[4, 5] = 2 * (dxdxi[0, 1] * ddxdxi[2, 3] + dxdxi[2, 1] * ddxdxi[0, 3]) + dxdxi[0, 0] * ddxdxi[2, 1] + dxdxi[2, 0] * ddxdxi[0, 1];
              jac23[5, 0] = dxdxi[0, 2] * ddxdxi[0, 1] + 2 * dxdxi[0, 1] * ddxdxi[0, 4];
              jac23[5, 1] = dxdxi[1, 2] * ddxdxi[1, 1] + 2 * dxdxi[1, 1] * ddxdxi[1, 4];
              jac23[5, 2] = dxdxi[2, 2] * ddxdxi[2, 1] + 2 * dxdxi[2, 1] * ddxdxi[2, 4];
              jac23[5, 3] = 2 * (dxdxi[0, 1] * ddxdxi[1, 4] + dxdxi[1, 1] * ddxdxi[0, 4]) + dxdxi[0, 2] * ddxdxi[1, 1] + dxdxi[1, 2] * ddxdxi[0, 1];
              jac23[5, 4] = 2 * (dxdxi[1, 1] * ddxdxi[2, 4] + dxdxi[2, 1] * ddxdxi[1, 4]) + dxdxi[1, 2] * ddxdxi[2, 1] + dxdxi[2, 2] * ddxdxi[1, 1];
              jac23[5, 5] = 2 * (dxdxi[0, 1] * ddxdxi[2, 4] + dxdxi[2, 1] * ddxdxi[0, 4]) + dxdxi[0, 2] * ddxdxi[2, 1] + dxdxi[2, 2] * ddxdxi[0, 1];
              jac23[6, 0] = dxdxi[0, 1] * ddxdxi[0, 2] + 2 * dxdxi[0, 2] * ddxdxi[0, 4];
              jac23[6, 1] = dxdxi[1, 1] * ddxdxi[1, 2] + 2 * dxdxi[1, 2] * ddxdxi[1, 4];
              jac23[6, 2] = dxdxi[2, 1] * ddxdxi[2, 2] + 2 * dxdxi[2, 2] * ddxdxi[2, 4];
              jac23[6, 3] = 2 * (dxdxi[0, 2] * ddxdxi[1, 4] + dxdxi[1, 2] * ddxdxi[0, 4]) + dxdxi[0, 1] * ddxdxi[1, 2] + dxdxi[1, 1] * ddxdxi[0, 2];
              jac23[6, 4] = 2 * (dxdxi[1, 2] * ddxdxi[2, 4] + dxdxi[2, 2] * ddxdxi[1, 4]) + dxdxi[1, 1] * ddxdxi[2, 2] + dxdxi[2, 1] * ddxdxi[1, 2];
              jac23[6, 5] = 2 * (dxdxi[0, 2] * ddxdxi[2, 4] + dxdxi[2, 2] * ddxdxi[0, 4]) + dxdxi[0, 1] * ddxdxi[2, 2] + dxdxi[2, 1] * ddxdxi[0, 2];
              jac23[7, 0] = dxdxi[0, 2] * ddxdxi[0, 0] + 2 * dxdxi[0, 0] * ddxdxi[0, 5];
              jac23[7, 1] = dxdxi[1, 2] * ddxdxi[1, 0] + 2 * dxdxi[1, 0] * ddxdxi[1, 5];
              jac23[7, 2] = dxdxi[2, 2] * ddxdxi[2, 0] + 2 * dxdxi[2, 0] * ddxdxi[2, 5];
              jac23[7, 3] = 2 * (dxdxi[0, 0] * ddxdxi[1, 5] + dxdxi[1, 0] * ddxdxi[0, 5]) + dxdxi[0, 2] * ddxdxi[1, 0] + dxdxi[1, 2] * ddxdxi[0, 0];
              jac23[7, 4] = 2 * (dxdxi[1, 0] * ddxdxi[2, 5] + dxdxi[2, 0] * ddxdxi[1, 5]) + dxdxi[1, 2] * ddxdxi[2, 0] + dxdxi[2, 2] * ddxdxi[1, 0];
              jac23[7, 5] = 2 * (dxdxi[0, 0] * ddxdxi[2, 5] + dxdxi[2, 0] * ddxdxi[0, 5]) + dxdxi[0, 2] * ddxdxi[2, 0] + dxdxi[2, 2] * ddxdxi[0, 0];
              jac23[8, 0] = dxdxi[0, 0] * ddxdxi[0, 2] + 2 * dxdxi[0, 2] * ddxdxi[0, 5];
              jac23[8, 1] = dxdxi[1, 0] * ddxdxi[1, 2] + 2 * dxdxi[1, 2] * ddxdxi[1, 5];
              jac23[8, 2] = dxdxi[2, 0] * ddxdxi[2, 2] + 2 * dxdxi[2, 2] * ddxdxi[2, 5];
              jac23[8, 3] = 2 * (dxdxi[0, 2] * ddxdxi[1, 5] + dxdxi[1, 2] * ddxdxi[0, 5]) + dxdxi[0, 0] * ddxdxi[1, 2] + dxdxi[1, 0] * ddxdxi[0, 2];
              jac23[8, 4] = 2 * (dxdxi[1, 2] * ddxdxi[2, 5] + dxdxi[2, 2] * ddxdxi[1, 5]) + dxdxi[1, 0] * ddxdxi[2, 2] + dxdxi[2, 0] * ddxdxi[1, 2];
              jac23[8, 5] = 2 * (dxdxi[0, 2] * ddxdxi[2, 5] + dxdxi[2, 2] * ddxdxi[0, 5]) + dxdxi[0, 0] * ddxdxi[2, 2] + dxdxi[2, 0] * ddxdxi[0, 2];
              jac23[9, 0] = dxdxi[0, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * ddxdxi[0, 5] + dxdxi[0, 2] * ddxdxi[0, 3];
              jac23[9, 1] = dxdxi[1, 0] * ddxdxi[1, 4] + dxdxi[1, 1] * ddxdxi[1, 5] + dxdxi[1, 2] * ddxdxi[1, 3];
              jac23[9, 2] = dxdxi[2, 0] * ddxdxi[2, 4] + dxdxi[2, 1] * ddxdxi[2, 5] + dxdxi[2, 2] * ddxdxi[2, 3];
              jac23[9, 3] = dxdxi[0, 0] * ddxdxi[1, 4] + dxdxi[0, 1] * ddxdxi[1, 5] + dxdxi[0, 2] * ddxdxi[1, 3] + dxdxi[1, 0] * ddxdxi[0, 4] + dxdxi[1, 1] * ddxdxi[0, 5] + dxdxi[1, 2] * ddxdxi[0, 3];
              jac23[9, 4] = dxdxi[1, 0] * ddxdxi[2, 4] + dxdxi[1, 1] * ddxdxi[2, 5] + dxdxi[1, 2] * ddxdxi[2, 3] + dxdxi[2, 0] * ddxdxi[1, 4] + dxdxi[2, 1] * ddxdxi[1, 5] + dxdxi[2, 2] * ddxdxi[1, 3];
              jac23[9, 5] = dxdxi[0, 0] * ddxdxi[2, 4] + dxdxi[0, 1] * ddxdxi[2, 5] + dxdxi[0, 2] * ddxdxi[2, 3] + dxdxi[2, 0] * ddxdxi[0, 4] + dxdxi[2, 1] * ddxdxi[0, 5] + dxdxi[2, 2] * ddxdxi[0, 3];
              DoubleMatrix term23 = MatrixFunctions.Product(ddNdX, jac23.Transpose());
              DoubleMatrix term3 = d3Ndxi - term13 - term23;

              DoubleMatrix jac13 = new DoubleMatrix(10, 10);
              jac13[0, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0];
              jac13[0, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
              jac13[0, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac13[0, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0];
              jac13[0, 4] = 3 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0];
              jac13[0, 5] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0];
              jac13[0, 6] = 3 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac13[0, 7] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 0];
              jac13[0, 8] = 3 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac13[0, 9] = 6 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 0];
              jac13[1, 0] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1];
              jac13[1, 1] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
              jac13[1, 2] = dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac13[1, 3] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1];
              jac13[1, 4] = 3 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1];
              jac13[1, 5] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1];
              jac13[1, 6] = 3 * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac13[1, 7] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 1];
              jac13[1, 8] = 3 * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac13[1, 9] = 6 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 1];
              jac13[2, 0] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2];
              jac13[2, 1] = dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
              jac13[2, 2] = dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac13[2, 3] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2];
              jac13[2, 4] = 3 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2];
              jac13[2, 5] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2];
              jac13[2, 6] = 3 * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac13[2, 7] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 2];
              jac13[2, 8] = 3 * dxdxi[0, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac13[2, 9] = 6 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 2];
              jac13[3, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1];
              jac13[3, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1];
              jac13[3, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1];
              jac13[3, 3] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0];
              jac13[3, 4] = dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1];
              jac13[3, 5] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0];
              jac13[3, 6] = dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 1];
              jac13[3, 7] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 0];
              jac13[3, 8] = dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 1];
              jac13[3, 9] = 2 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 0];
              jac13[4, 0] = dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1];
              jac13[4, 1] = dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1];
              jac13[4, 2] = dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1];
              jac13[4, 3] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1];
              jac13[4, 4] = dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1];
              jac13[4, 5] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 1];
              jac13[4, 6] = dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 1];
              jac13[4, 7] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 1];
              jac13[4, 8] = dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 1];
              jac13[4, 9] = 2 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 0];
              jac13[5, 0] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2];
              jac13[5, 1] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2];
              jac13[5, 2] = dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 2];
              jac13[5, 3] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 2] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1];
              jac13[5, 4] = dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 2];
              jac13[5, 5] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 2] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1];
              jac13[5, 6] = dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 2];
              jac13[5, 7] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 2] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 1];
              jac13[5, 8] = dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 2];
              jac13[5, 9] = 2 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 1];
              jac13[6, 0] = dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2];
              jac13[6, 1] = dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2];
              jac13[6, 2] = dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2];
              jac13[6, 3] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 2];
              jac13[6, 4] = dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 2];
              jac13[6, 5] = dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 2];
              jac13[6, 6] = dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 2];
              jac13[6, 7] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 2];
              jac13[6, 8] = dxdxi[0, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 2];
              jac13[6, 9] = 2 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 1];
              jac13[7, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2];
              jac13[7, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2];
              jac13[7, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 2];
              jac13[7, 3] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 2] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0];
              jac13[7, 4] = dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 2];
              jac13[7, 5] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 2] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0];
              jac13[7, 6] = dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 2];
              jac13[7, 7] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 2] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 0];
              jac13[7, 8] = dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 2];
              jac13[7, 9] = 2 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 0];
              jac13[8, 0] = dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2];
              jac13[8, 1] = dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2];
              jac13[8, 2] = dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2];
              jac13[8, 3] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 2];
              jac13[8, 4] = dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 2];
              jac13[8, 5] = dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 2];
              jac13[8, 6] = dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 2];
              jac13[8, 7] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 2];
              jac13[8, 8] = dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 2];
              jac13[8, 9] = 2 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 0];
              jac13[9, 0] = dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 2];
              jac13[9, 1] = dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 2];
              jac13[9, 2] = dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 2];
              jac13[9, 3] = dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 2] + dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 0] + dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 1];
              jac13[9, 4] = dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 2] + dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 1];
              jac13[9, 5] = dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 0] + dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 1];
              jac13[9, 6] = dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 1];
              jac13[9, 7] = dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 0] + dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 1];
              jac13[9, 8] = dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 1];
              jac13[9, 9] = dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 2] + dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 1] + dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 0] + dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 1] + dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 0];

              d3NdX = MatrixFunctions.Product(MatrixFunctions.Inverse(jac13), term3.Transpose()).Transpose();
              gpsijk.SetValue(DataInGausspoint.d3NdX, d3NdX/*.Transpose()*/);
            }

            if (isNeed_d4NdX)
            {
              DoubleMatrix jac14 = new DoubleMatrix(15, 15);
              jac14[0, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0];
              jac14[0, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
              jac14[0, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[0, 3] = 4 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0];
              jac14[0, 4] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0];
              jac14[0, 5] = 4 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
              jac14[0, 6] = 4 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0];
              jac14[0, 7] = 6 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[0, 8] = 4 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[0, 9] = 4 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 0];
              jac14[0, 10] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[0, 11] = 4 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[0, 12] = 12 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 0];
              jac14[0, 13] = 12 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0];
              jac14[0, 14] = 12 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[1, 0] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1];
              jac14[1, 1] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
              jac14[1, 2] = dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[1, 3] = 4 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1];
              jac14[1, 4] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1];
              jac14[1, 5] = 4 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
              jac14[1, 6] = 4 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1];
              jac14[1, 7] = 6 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[1, 8] = 4 * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[1, 9] = 4 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 1];
              jac14[1, 10] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[1, 11] = 4 * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[1, 12] = 12 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 1];
              jac14[1, 13] = 12 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1];
              jac14[1, 14] = 12 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[2, 0] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2];
              jac14[2, 1] = dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
              jac14[2, 2] = dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[2, 3] = 4 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2];
              jac14[2, 4] = 6 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2];
              jac14[2, 5] = 4 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
              jac14[2, 6] = 4 * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2];
              jac14[2, 7] = 6 * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[2, 8] = 4 * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[2, 9] = 4 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 2];
              jac14[2, 10] = 6 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[2, 11] = 4 * dxdxi[0, 2] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[2, 12] = 12 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 2];
              jac14[2, 13] = 12 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2];
              jac14[2, 14] = 12 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2];

              jac14[3, 0] = dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0];
              jac14[3, 1] = dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
              jac14[3, 2] = dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[3, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1];
              jac14[3, 4] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1];
              jac14[3, 5] = 3 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[1, 0] + dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
              jac14[3, 6] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0] + dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 1];
              jac14[3, 7] = 3 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] + 3 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 1];
              jac14[3, 8] = 3 * dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[2, 0] + dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[3, 9] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 0] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 1];
              jac14[3, 10] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 1];
              jac14[3, 11] = 3 * dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[2, 0] + dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[3, 12] = 6 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 1];
              jac14[3, 13] = 6 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0];
              jac14[3, 14] = 6 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 1] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 0];

              jac14[4, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1];
              jac14[4, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1];
              jac14[4, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[4, 3] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0];
              jac14[4, 4] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0];
              jac14[4, 5] = 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1];
              jac14[4, 6] = 2 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 1] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0];
              jac14[4, 7] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 4 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 1] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[4, 8] = 2 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1];
              jac14[4, 9] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 0];
              jac14[4, 10] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[4, 11] = 2 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1];
              jac14[4, 12] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 1];
              jac14[4, 13] = 2 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0] + 4 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 1] + 4 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0];
              jac14[4, 14] = 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 1] + 4 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 1] + 4 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 1];

              jac14[5, 0] = dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1];
              jac14[5, 1] = dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
              jac14[5, 2] = dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[5, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0];
              jac14[5, 4] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1];
              jac14[5, 5] = 3 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] + dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
              jac14[5, 6] = 3 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0];
              jac14[5, 7] = 3 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1] + 3 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 1];
              jac14[5, 8] = 3 * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1] + dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[5, 9] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 0];
              jac14[5, 10] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 1];
              jac14[5, 11] = 3 * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 1] + dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[5, 12] = 6 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 0] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 1];
              jac14[5, 13] = 6 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0];
              jac14[5, 14] = 6 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 1] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 1];

              jac14[6, 0] = dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1];
              jac14[6, 1] = dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
              jac14[6, 2] = dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[6, 3] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 2];
              jac14[6, 4] = 3 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 2];
              jac14[6, 5] = 3 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[1, 1] + dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
              jac14[6, 6] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 2];
              jac14[6, 7] = 3 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 3 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 2];
              jac14[6, 8] = 3 * dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[2, 1] + dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[6, 9] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 2];
              jac14[6, 10] = 3 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 2];
              jac14[6, 11] = 3 * dxdxi[0, 1] * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[2, 1] + dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[6, 12] = 6 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 2];
              jac14[6, 13] = 6 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 1];
              jac14[6, 14] = 6 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 2] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 1];

              jac14[7, 0] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2];
              jac14[7, 1] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2];
              jac14[7, 2] = dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[7, 3] = 2 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1];
              jac14[7, 4] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 2] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1];
              jac14[7, 5] = 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2];
              jac14[7, 6] = 2 * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1];
              jac14[7, 7] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 4 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[7, 8] = 2 * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 2];
              jac14[7, 9] = 2 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 1];
              jac14[7, 10] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1];
              jac14[7, 11] = 2 * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 2];
              jac14[7, 12] = 2 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 1] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 1] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 2];
              jac14[7, 13] = 2 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 2] + 2 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1] + 4 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 2] + 4 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1];
              jac14[7, 14] = 2 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 2] + 4 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 2] + 4 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 2];

              jac14[8, 0] = dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2];
              jac14[8, 1] = dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
              jac14[8, 2] = dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[8, 3] = 3 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1];
              jac14[8, 4] = 3 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 2];
              jac14[8, 5] = 3 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2] + dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
              jac14[8, 6] = 3 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1];
              jac14[8, 7] = 3 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2] + 3 * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 2];
              jac14[8, 8] = 3 * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[8, 9] = 3 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 1];
              jac14[8, 10] = 3 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 2];
              jac14[8, 11] = 3 * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[8, 12] = 6 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 1] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 2];
              jac14[8, 13] = 6 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1];
              jac14[8, 14] = 6 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 2] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 2];

              jac14[9, 0] = dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0];
              jac14[9, 1] = dxdxi[1, 2] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
              jac14[9, 2] = dxdxi[2, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[9, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 2];
              jac14[9, 4] = 3 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 2];
              jac14[9, 5] = 3 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 0] * dxdxi[1, 0] + dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
              jac14[9, 6] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0] + dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 2];
              jac14[9, 7] = 3 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 3 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 2];
              jac14[9, 8] = 3 * dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 0] * dxdxi[2, 0] + dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[9, 9] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 0] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 2];
              jac14[9, 10] = 3 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 2];
              jac14[9, 11] = 3 * dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 0] * dxdxi[2, 0] + dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[9, 12] = 6 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 2];
              jac14[9, 13] = 6 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 0];
              jac14[9, 14] = 6 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 2] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 0];

              jac14[10, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2];
              jac14[10, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2];
              jac14[10, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[10, 3] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0];
              jac14[10, 4] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0];
              jac14[10, 5] = 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2];
              jac14[10, 6] = 2 * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0];
              jac14[10, 7] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 4 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[10, 8] = 2 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 2];
              jac14[10, 9] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 0];
              jac14[10, 10] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0];
              jac14[10, 11] = 2 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 2];
              jac14[10, 12] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 2];
              jac14[10, 13] = 2 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 2] + 2 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0] + 4 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 2] + 4 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0];
              jac14[10, 14] = 2 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 4 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 2] + 4 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 2];

              jac14[11, 0] = dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2];
              jac14[11, 1] = dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
              jac14[11, 2] = dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[11, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0];
              jac14[11, 4] = 3 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 2];
              jac14[11, 5] = 3 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2] + dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2];
              jac14[11, 6] = 3 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0];
              jac14[11, 7] = 3 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2] + 3 * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 2];
              jac14[11, 8] = 3 * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2] + dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[11, 9] = 3 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 0];
              jac14[11, 10] = 3 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 2];
              jac14[11, 11] = 3 * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[2, 2] + dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 2];
              jac14[11, 12] = 6 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 0] + 3 * dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 2];
              jac14[11, 13] = 6 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 0];
              jac14[11, 14] = 6 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 2] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2];

              jac14[12, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 2];
              jac14[12, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 2];
              jac14[12, 2] = dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 2];
              jac14[12, 3] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 2] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 0];
              jac14[12, 4] = 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 2] + dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 1] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 2];
              jac14[12, 5] = 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 2] + dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1];
              jac14[12, 6] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 2] + dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 1] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 0];
              jac14[12, 7] = 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 1] + dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 2];
              jac14[12, 8] = 2 * dxdxi[1, 0] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1];
              jac14[12, 9] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 2] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 0];
              jac14[12, 10] = 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 1] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 2];
              jac14[12, 11] = 2 * dxdxi[0, 0] * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[2, 1];
              jac14[12, 12] = 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 2] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 1] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 2];
              jac14[12, 13] = 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[0, 0] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[0, 1] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[0, 2] * dxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[0, 0] * dxdxi[2, 1] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[0, 0] * dxdxi[2, 2] + dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[0, 2] * dxdxi[2, 1] + dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[0, 1] * dxdxi[2, 2];
              jac14[12, 14] = 2 * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[1, 0] * dxdxi[0, 0] + 2 * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[1, 1] * dxdxi[0, 0] + 2 * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[1, 2] * dxdxi[0, 0] + 2 * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[1, 0] * dxdxi[0, 1] + 2 * dxdxi[2, 0] * dxdxi[2, 1] * dxdxi[1, 0] * dxdxi[0, 2] + dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[1, 2] * dxdxi[0, 1] + dxdxi[2, 0] * dxdxi[2, 0] * dxdxi[1, 1] * dxdxi[0, 2];

              jac14[13, 0] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[0, 2];
              jac14[13, 1] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[1, 2];
              jac14[13, 2] = dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[2, 2];
              jac14[13, 3] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 2] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 0] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 1];
              jac14[13, 4] = 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 2] + dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 0] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 2];
              jac14[13, 5] = 2 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[1, 2] + dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2] + dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 0];
              jac14[13, 6] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 0] + 2 * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 1];
              jac14[13, 7] = 2 * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 0] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[2, 0] * dxdxi[2, 2];
              jac14[13, 8] = 2 * dxdxi[1, 1] * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[1, 0] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 0];
              jac14[13, 9] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 0] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 1];
              jac14[13, 10] = 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 0] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[2, 0] * dxdxi[2, 2];
              jac14[13, 11] = 2 * dxdxi[0, 1] * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[2, 2] + dxdxi[0, 0] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 2] + dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[2, 0];
              jac14[13, 12] = 2 * dxdxi[0, 0] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 0] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 2] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 0] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 2];
              jac14[13, 13] = 2 * dxdxi[1, 0] * dxdxi[1, 2] * dxdxi[0, 1] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[0, 0] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[0, 2] * dxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * dxdxi[0, 1] * dxdxi[2, 0] + 2 * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[0, 1] * dxdxi[2, 2] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[0, 2] * dxdxi[2, 0] + dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[0, 0] * dxdxi[2, 2];
              jac14[13, 14] = 2 * dxdxi[2, 0] * dxdxi[2, 2] * dxdxi[1, 1] * dxdxi[0, 1] + 2 * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[1, 0] * dxdxi[0, 1] + 2 * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[1, 2] * dxdxi[0, 1] + 2 * dxdxi[2, 1] * dxdxi[2, 2] * dxdxi[1, 1] * dxdxi[0, 0] + 2 * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[1, 1] * dxdxi[0, 2] + dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[1, 2] * dxdxi[0, 0] + dxdxi[2, 1] * dxdxi[2, 1] * dxdxi[1, 0] * dxdxi[0, 2];

              jac14[14, 0] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[0, 0];
              jac14[14, 1] = dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[1, 0];
              jac14[14, 2] = dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[2, 0];
              jac14[14, 3] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[1, 0] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[1, 1] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 2];
              jac14[14, 4] = 2 * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 0] + dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 1] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[1, 0];
              jac14[14, 5] = 2 * dxdxi[0, 2] * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[1, 0] + dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 0] + dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 1];
              jac14[14, 6] = dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[2, 0] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[1, 0] * dxdxi[2, 1] + 2 * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[2, 2];
              jac14[14, 7] = 2 * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 0] + dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 1] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[2, 1] * dxdxi[2, 0];
              jac14[14, 8] = 2 * dxdxi[1, 2] * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[2, 0] + dxdxi[1, 1] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 0] + dxdxi[1, 0] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 1];
              jac14[14, 9] = dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[2, 0] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[2, 1] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[2, 2];
              jac14[14, 10] = 2 * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[2, 2] * dxdxi[2, 0] + dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 1] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[2, 1] * dxdxi[2, 0];
              jac14[14, 11] = 2 * dxdxi[0, 2] * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[2, 0] + dxdxi[0, 1] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 0] + dxdxi[0, 0] * dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[2, 1];
              jac14[14, 12] = 2 * dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * dxdxi[1, 2] * dxdxi[2, 1] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * dxdxi[1, 2] * dxdxi[2, 0] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 0] * dxdxi[2, 1] + dxdxi[0, 2] * dxdxi[0, 2] * dxdxi[1, 1] * dxdxi[2, 0];
              jac14[14, 13] = 2 * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[0, 2] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[1, 0] * dxdxi[0, 1] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[0, 0] * dxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[1, 0] * dxdxi[0, 2] * dxdxi[2, 1] + 2 * dxdxi[1, 2] * dxdxi[1, 1] * dxdxi[0, 2] * dxdxi[2, 0] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[0, 0] * dxdxi[2, 1] + dxdxi[1, 2] * dxdxi[1, 2] * dxdxi[0, 1] * dxdxi[2, 0];
              jac14[14, 14] = 2 * dxdxi[2, 1] * dxdxi[2, 0] * dxdxi[1, 2] * dxdxi[0, 2] + 2 * dxdxi[2, 2] * dxdxi[2, 0] * dxdxi[1, 1] * dxdxi[0, 2] + 2 * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[1, 0] * dxdxi[0, 2] + 2 * dxdxi[2, 2] * dxdxi[2, 0] * dxdxi[1, 2] * dxdxi[0, 1] + 2 * dxdxi[2, 2] * dxdxi[2, 1] * dxdxi[1, 2] * dxdxi[0, 0] + dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[1, 0] * dxdxi[0, 1] + dxdxi[2, 2] * dxdxi[2, 2] * dxdxi[1, 1] * dxdxi[0, 0];

              DoubleMatrix d4xdxi = Jacobian14At(d4Ndxi);
              DoubleMatrix term14 = MatrixFunctions.Product(dNdx, d4xdxi);

              DoubleMatrix jac24 = new DoubleMatrix(15, 6);

              jac24[0, 0] = 3 * ddxdxi[0, 0] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * d3xdxi[0, 0];
              jac24[0, 1] = 3 * ddxdxi[1, 0] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * d3xdxi[1, 0];
              jac24[0, 2] = 3 * ddxdxi[2, 0] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * d3xdxi[2, 0];
              jac24[0, 3] = 4 * dxdxi[0, 0] * d3xdxi[1, 0] + 4 * dxdxi[1, 0] * d3xdxi[0, 0] + 6 * ddxdxi[0, 0] * ddxdxi[1, 0];
              jac24[0, 4] = 4 * dxdxi[1, 0] * d3xdxi[2, 0] + 4 * dxdxi[2, 0] * d3xdxi[1, 0] + 6 * ddxdxi[1, 0] * ddxdxi[2, 0];
              jac24[0, 5] = 4 * dxdxi[0, 0] * d3xdxi[2, 0] + 4 * dxdxi[2, 0] * d3xdxi[0, 0] + 6 * ddxdxi[0, 0] * ddxdxi[2, 0];

              jac24[1, 0] = 3 * ddxdxi[0, 1] * ddxdxi[0, 1] + 4 * dxdxi[0, 1] * d3xdxi[0, 1];
              jac24[1, 1] = 3 * ddxdxi[1, 1] * ddxdxi[1, 1] + 4 * dxdxi[1, 1] * d3xdxi[1, 1];
              jac24[1, 2] = 3 * ddxdxi[2, 1] * ddxdxi[2, 1] + 4 * dxdxi[2, 1] * d3xdxi[2, 1];
              jac24[1, 3] = 4 * dxdxi[0, 1] * d3xdxi[1, 1] + 4 * dxdxi[1, 1] * d3xdxi[0, 1] + 6 * ddxdxi[0, 1] * ddxdxi[1, 1];
              jac24[1, 4] = 4 * dxdxi[1, 1] * d3xdxi[2, 1] + 4 * dxdxi[2, 1] * d3xdxi[1, 1] + 6 * ddxdxi[1, 1] * ddxdxi[2, 1];
              jac24[1, 5] = 4 * dxdxi[0, 1] * d3xdxi[2, 1] + 4 * dxdxi[2, 1] * d3xdxi[0, 1] + 6 * ddxdxi[0, 1] * ddxdxi[2, 1];

              jac24[2, 0] = 3 * ddxdxi[0, 2] * ddxdxi[0, 2] + 4 * dxdxi[0, 2] * d3xdxi[0, 2];
              jac24[2, 1] = 3 * ddxdxi[1, 2] * ddxdxi[1, 2] + 4 * dxdxi[1, 2] * d3xdxi[1, 2];
              jac24[2, 2] = 3 * ddxdxi[2, 2] * ddxdxi[2, 2] + 4 * dxdxi[2, 2] * d3xdxi[2, 2];
              jac24[2, 3] = 4 * dxdxi[0, 2] * d3xdxi[1, 2] + 4 * dxdxi[1, 2] * d3xdxi[0, 2] + 6 * ddxdxi[0, 2] * ddxdxi[1, 2];
              jac24[2, 4] = 4 * dxdxi[1, 2] * d3xdxi[2, 2] + 4 * dxdxi[2, 2] * d3xdxi[1, 2] + 6 * ddxdxi[1, 2] * ddxdxi[2, 2];
              jac24[2, 5] = 4 * dxdxi[0, 2] * d3xdxi[2, 2] + 4 * dxdxi[2, 2] * d3xdxi[0, 2] + 6 * ddxdxi[0, 2] * ddxdxi[2, 2];

              jac24[3, 0] = dxdxi[0, 1] * d3xdxi[0, 0] + 3 * ddxdxi[0, 0] * ddxdxi[0, 3] + 3 * dxdxi[0, 0] * d3xdxi[0, 3];
              jac24[3, 1] = dxdxi[1, 1] * d3xdxi[1, 0] + 3 * ddxdxi[1, 0] * ddxdxi[1, 3] + 3 * dxdxi[1, 0] * d3xdxi[1, 3];
              jac24[3, 2] = dxdxi[2, 1] * d3xdxi[2, 0] + 3 * ddxdxi[2, 0] * ddxdxi[2, 3] + 3 * dxdxi[2, 0] * d3xdxi[2, 3];
              jac24[3, 3] = 3 * dxdxi[0, 0] * d3xdxi[1, 3] + 3 * dxdxi[1, 0] * d3xdxi[0, 3] + 3 * ddxdxi[0, 3] * ddxdxi[1, 0] + 3 * ddxdxi[0, 0] * ddxdxi[1, 3] + dxdxi[0, 1] * d3xdxi[1, 0] + dxdxi[1, 1] * d3xdxi[0, 0];
              jac24[3, 4] = 3 * dxdxi[1, 0] * d3xdxi[2, 3] + 3 * dxdxi[2, 0] * d3xdxi[1, 3] + 3 * ddxdxi[1, 3] * ddxdxi[2, 0] + 3 * ddxdxi[1, 0] * ddxdxi[2, 3] + dxdxi[1, 1] * d3xdxi[2, 0] + dxdxi[2, 1] * d3xdxi[1, 0];
              jac24[3, 5] = 3 * dxdxi[0, 0] * d3xdxi[2, 3] + 3 * dxdxi[2, 0] * d3xdxi[0, 3] + 3 * ddxdxi[0, 3] * ddxdxi[2, 0] + 3 * ddxdxi[0, 0] * ddxdxi[2, 3] + dxdxi[0, 1] * d3xdxi[2, 0] + dxdxi[2, 1] * d3xdxi[0, 0];

              jac24[4, 0] = ddxdxi[0, 0] * ddxdxi[0, 1] + 2 * ddxdxi[0, 3] * ddxdxi[0, 3] + 2 * dxdxi[0, 0] * d3xdxi[0, 4] + 2 * dxdxi[0, 1] * d3xdxi[0, 3];
              jac24[4, 1] = ddxdxi[1, 0] * ddxdxi[1, 1] + 2 * ddxdxi[1, 3] * ddxdxi[1, 3] + 2 * dxdxi[1, 0] * d3xdxi[1, 4] + 2 * dxdxi[1, 1] * d3xdxi[1, 3];
              jac24[4, 2] = ddxdxi[2, 0] * ddxdxi[2, 1] + 2 * ddxdxi[2, 3] * ddxdxi[2, 3] + 2 * dxdxi[2, 0] * d3xdxi[2, 4] + 2 * dxdxi[2, 1] * d3xdxi[2, 3];
              jac24[4, 3] = 2 * dxdxi[0, 0] * d3xdxi[1, 4] + 2 * dxdxi[0, 1] * d3xdxi[1, 3] + 2 * dxdxi[1, 0] * d3xdxi[0, 4] + 2 * dxdxi[1, 1] * d3xdxi[0, 3] + 4 * ddxdxi[0, 3] * ddxdxi[1, 3] + ddxdxi[0, 0] * ddxdxi[1, 1] + ddxdxi[0, 1] * ddxdxi[1, 0];
              jac24[4, 4] = 2 * dxdxi[1, 0] * d3xdxi[2, 4] + 2 * dxdxi[1, 1] * d3xdxi[2, 3] + 2 * dxdxi[2, 0] * d3xdxi[1, 4] + 2 * dxdxi[2, 1] * d3xdxi[1, 3] + 4 * ddxdxi[1, 3] * ddxdxi[2, 3] + ddxdxi[1, 0] * ddxdxi[2, 1] + ddxdxi[1, 1] * ddxdxi[2, 0];
              jac24[4, 5] = 2 * dxdxi[0, 0] * d3xdxi[2, 4] + 2 * dxdxi[0, 1] * d3xdxi[2, 3] + 2 * dxdxi[2, 0] * d3xdxi[0, 4] + 2 * dxdxi[2, 1] * d3xdxi[0, 3] + 4 * ddxdxi[0, 3] * ddxdxi[2, 3] + ddxdxi[0, 0] * ddxdxi[2, 1] + ddxdxi[0, 1] * ddxdxi[2, 0];

              jac24[5, 0] = dxdxi[0, 0] * d3xdxi[0, 1] + 3 * ddxdxi[0, 1] * ddxdxi[0, 3] + 3 * dxdxi[0, 1] * d3xdxi[0, 4];
              jac24[5, 1] = dxdxi[1, 0] * d3xdxi[1, 1] + 3 * ddxdxi[1, 1] * ddxdxi[1, 3] + 3 * dxdxi[1, 1] * d3xdxi[1, 4];
              jac24[5, 2] = dxdxi[2, 0] * d3xdxi[2, 1] + 3 * ddxdxi[2, 1] * ddxdxi[2, 3] + 3 * dxdxi[2, 1] * d3xdxi[2, 4];
              jac24[5, 3] = 3 * dxdxi[0, 1] * d3xdxi[1, 4] + 3 * dxdxi[1, 1] * d3xdxi[0, 4] + 3 * ddxdxi[0, 1] * ddxdxi[1, 3] + 3 * ddxdxi[0, 3] * ddxdxi[1, 1] + dxdxi[0, 0] * d3xdxi[1, 1] + dxdxi[1, 0] * d3xdxi[0, 1];
              jac24[5, 4] = 3 * dxdxi[1, 1] * d3xdxi[2, 4] + 3 * dxdxi[2, 1] * d3xdxi[1, 4] + 3 * ddxdxi[1, 1] * ddxdxi[2, 3] + 3 * ddxdxi[1, 3] * ddxdxi[2, 1] + dxdxi[1, 0] * d3xdxi[2, 1] + dxdxi[2, 0] * d3xdxi[1, 1];
              jac24[5, 5] = 3 * dxdxi[0, 1] * d3xdxi[2, 4] + 3 * dxdxi[2, 1] * d3xdxi[0, 4] + 3 * ddxdxi[0, 1] * ddxdxi[2, 3] + 3 * ddxdxi[0, 3] * ddxdxi[2, 1] + dxdxi[0, 0] * d3xdxi[2, 1] + dxdxi[2, 0] * d3xdxi[0, 1];

              jac24[6, 0] = dxdxi[0, 2] * d3xdxi[0, 1] + 3 * ddxdxi[0, 1] * ddxdxi[0, 4] + 3 * dxdxi[0, 1] * d3xdxi[0, 5];
              jac24[6, 1] = dxdxi[1, 2] * d3xdxi[1, 1] + 3 * ddxdxi[1, 1] * ddxdxi[1, 4] + 3 * dxdxi[1, 1] * d3xdxi[1, 5];
              jac24[6, 2] = dxdxi[2, 2] * d3xdxi[2, 1] + 3 * ddxdxi[2, 1] * ddxdxi[2, 4] + 3 * dxdxi[2, 1] * d3xdxi[2, 5];
              jac24[6, 3] = 3 * dxdxi[0, 1] * d3xdxi[1, 5] + 3 * dxdxi[1, 1] * d3xdxi[0, 5] + 3 * ddxdxi[0, 4] * ddxdxi[1, 1] + 3 * ddxdxi[0, 1] * ddxdxi[1, 4] + dxdxi[0, 2] * d3xdxi[1, 1] + dxdxi[1, 2] * d3xdxi[0, 1];
              jac24[6, 4] = 3 * dxdxi[1, 1] * d3xdxi[2, 5] + 3 * dxdxi[2, 1] * d3xdxi[1, 5] + 3 * ddxdxi[1, 4] * ddxdxi[2, 1] + 3 * ddxdxi[1, 1] * ddxdxi[2, 4] + dxdxi[1, 2] * d3xdxi[2, 1] + dxdxi[2, 2] * d3xdxi[1, 1];
              jac24[6, 5] = 3 * dxdxi[0, 1] * d3xdxi[2, 5] + 3 * dxdxi[2, 1] * d3xdxi[0, 5] + 3 * ddxdxi[0, 4] * ddxdxi[2, 1] + 3 * ddxdxi[0, 1] * ddxdxi[2, 4] + dxdxi[0, 2] * d3xdxi[2, 1] + dxdxi[2, 2] * d3xdxi[0, 1];

              jac24[7, 0] = ddxdxi[0, 1] * ddxdxi[0, 2] + 2 * ddxdxi[0, 4] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * d3xdxi[0, 6] + 2 * dxdxi[0, 2] * d3xdxi[0, 5];
              jac24[7, 1] = ddxdxi[1, 1] * ddxdxi[1, 2] + 2 * ddxdxi[1, 4] * ddxdxi[1, 4] + 2 * dxdxi[1, 1] * d3xdxi[1, 6] + 2 * dxdxi[1, 2] * d3xdxi[1, 5];
              jac24[7, 2] = ddxdxi[2, 1] * ddxdxi[2, 2] + 2 * ddxdxi[2, 4] * ddxdxi[2, 4] + 2 * dxdxi[2, 1] * d3xdxi[2, 6] + 2 * dxdxi[2, 2] * d3xdxi[2, 5];
              jac24[7, 3] = 2 * dxdxi[0, 1] * d3xdxi[1, 6] + 2 * dxdxi[0, 2] * d3xdxi[1, 5] + 2 * dxdxi[1, 1] * d3xdxi[0, 6] + 2 * dxdxi[1, 2] * d3xdxi[0, 5] + 4 * ddxdxi[0, 4] * ddxdxi[1, 4] + ddxdxi[0, 1] * ddxdxi[1, 2] + ddxdxi[0, 2] * ddxdxi[1, 1];
              jac24[7, 4] = 2 * dxdxi[1, 1] * d3xdxi[2, 6] + 2 * dxdxi[1, 2] * d3xdxi[2, 5] + 2 * dxdxi[2, 1] * d3xdxi[1, 6] + 2 * dxdxi[2, 2] * d3xdxi[1, 5] + 4 * ddxdxi[1, 4] * ddxdxi[2, 4] + ddxdxi[1, 1] * ddxdxi[2, 2] + ddxdxi[1, 2] * ddxdxi[2, 1];
              jac24[7, 5] = 2 * dxdxi[0, 1] * d3xdxi[2, 6] + 2 * dxdxi[0, 2] * d3xdxi[2, 5] + 2 * dxdxi[2, 1] * d3xdxi[0, 6] + 2 * dxdxi[2, 2] * d3xdxi[0, 5] + 4 * ddxdxi[0, 4] * ddxdxi[2, 4] + ddxdxi[0, 1] * ddxdxi[2, 2] + ddxdxi[0, 2] * ddxdxi[2, 1];

              jac24[8, 0] = dxdxi[0, 1] * d3xdxi[0, 2] + 3 * ddxdxi[0, 2] * ddxdxi[0, 4] + 3 * dxdxi[0, 2] * d3xdxi[0, 6];
              jac24[8, 1] = dxdxi[1, 1] * d3xdxi[1, 2] + 3 * ddxdxi[1, 2] * ddxdxi[1, 4] + 3 * dxdxi[1, 2] * d3xdxi[1, 6];
              jac24[8, 2] = dxdxi[2, 1] * d3xdxi[2, 2] + 3 * ddxdxi[2, 2] * ddxdxi[2, 4] + 3 * dxdxi[2, 2] * d3xdxi[2, 6];
              jac24[8, 3] = 3 * dxdxi[0, 2] * d3xdxi[1, 6] + 3 * dxdxi[1, 2] * d3xdxi[0, 6] + 3 * ddxdxi[0, 2] * ddxdxi[1, 4] + 3 * ddxdxi[0, 4] * ddxdxi[1, 2] + dxdxi[0, 1] * d3xdxi[1, 2] + dxdxi[1, 1] * d3xdxi[0, 2];
              jac24[8, 4] = 3 * dxdxi[1, 2] * d3xdxi[2, 6] + 3 * dxdxi[2, 2] * d3xdxi[1, 6] + 3 * ddxdxi[1, 2] * ddxdxi[2, 4] + 3 * ddxdxi[1, 4] * ddxdxi[2, 2] + dxdxi[1, 1] * d3xdxi[2, 2] + dxdxi[2, 1] * d3xdxi[1, 2];
              jac24[8, 5] = 3 * dxdxi[0, 2] * d3xdxi[2, 6] + 3 * dxdxi[2, 2] * d3xdxi[0, 6] + 3 * ddxdxi[0, 2] * ddxdxi[2, 4] + 3 * ddxdxi[0, 4] * ddxdxi[2, 2] + dxdxi[0, 1] * d3xdxi[2, 2] + dxdxi[2, 1] * d3xdxi[0, 2];

              jac24[9, 0] = dxdxi[0, 2] * d3xdxi[0, 0] + 3 * ddxdxi[0, 0] * ddxdxi[0, 5] + 3 * dxdxi[0, 0] * d3xdxi[0, 7];
              jac24[9, 1] = dxdxi[1, 2] * d3xdxi[1, 0] + 3 * ddxdxi[1, 0] * ddxdxi[1, 5] + 3 * dxdxi[1, 0] * d3xdxi[1, 7];
              jac24[9, 2] = dxdxi[2, 2] * d3xdxi[2, 0] + 3 * ddxdxi[2, 0] * ddxdxi[2, 5] + 3 * dxdxi[2, 0] * d3xdxi[2, 7];
              jac24[9, 3] = 3 * dxdxi[0, 0] * d3xdxi[1, 7] + 3 * dxdxi[1, 0] * d3xdxi[0, 7] + 3 * ddxdxi[0, 5] * ddxdxi[1, 0] + 3 * ddxdxi[0, 0] * ddxdxi[1, 5] + dxdxi[0, 2] * d3xdxi[1, 0] + dxdxi[1, 2] * d3xdxi[0, 0];
              jac24[9, 4] = 3 * dxdxi[1, 0] * d3xdxi[2, 7] + 3 * dxdxi[2, 0] * d3xdxi[1, 7] + 3 * ddxdxi[1, 5] * ddxdxi[2, 0] + 3 * ddxdxi[1, 0] * ddxdxi[2, 5] + dxdxi[1, 2] * d3xdxi[2, 0] + dxdxi[2, 2] * d3xdxi[1, 0];
              jac24[9, 5] = 3 * dxdxi[0, 0] * d3xdxi[2, 7] + 3 * dxdxi[2, 0] * d3xdxi[0, 7] + 3 * ddxdxi[0, 5] * ddxdxi[2, 0] + 3 * ddxdxi[0, 0] * ddxdxi[2, 5] + dxdxi[0, 2] * d3xdxi[2, 0] + dxdxi[2, 2] * d3xdxi[0, 0];

              jac24[10, 0] = ddxdxi[0, 0] * ddxdxi[0, 2] + 2 * ddxdxi[0, 5] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * d3xdxi[0, 8] + 2 * dxdxi[0, 2] * d3xdxi[0, 7];
              jac24[10, 1] = ddxdxi[1, 0] * ddxdxi[1, 2] + 2 * ddxdxi[1, 5] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * d3xdxi[1, 8] + 2 * dxdxi[1, 2] * d3xdxi[1, 7];
              jac24[10, 2] = ddxdxi[2, 0] * ddxdxi[2, 2] + 2 * ddxdxi[2, 5] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * d3xdxi[2, 8] + 2 * dxdxi[2, 2] * d3xdxi[2, 7];
              jac24[10, 3] = 2 * dxdxi[0, 0] * d3xdxi[1, 8] + 2 * dxdxi[0, 2] * d3xdxi[1, 7] + 2 * dxdxi[1, 0] * d3xdxi[0, 8] + 2 * dxdxi[1, 2] * d3xdxi[0, 7] + 4 * ddxdxi[0, 5] * ddxdxi[1, 5] + ddxdxi[0, 0] * ddxdxi[1, 2] + ddxdxi[0, 2] * ddxdxi[1, 0];
              jac24[10, 4] = 2 * dxdxi[1, 0] * d3xdxi[2, 8] + 2 * dxdxi[1, 2] * d3xdxi[2, 7] + 2 * dxdxi[2, 0] * d3xdxi[1, 8] + 2 * dxdxi[2, 2] * d3xdxi[1, 7] + 4 * ddxdxi[1, 5] * ddxdxi[2, 5] + ddxdxi[1, 0] * ddxdxi[2, 2] + ddxdxi[1, 2] * ddxdxi[2, 0];
              jac24[10, 5] = 2 * dxdxi[0, 0] * d3xdxi[2, 8] + 2 * dxdxi[0, 2] * d3xdxi[2, 7] + 2 * dxdxi[2, 0] * d3xdxi[0, 8] + 2 * dxdxi[2, 2] * d3xdxi[0, 7] + 4 * ddxdxi[0, 5] * ddxdxi[2, 5] + ddxdxi[0, 0] * ddxdxi[2, 2] + ddxdxi[0, 2] * ddxdxi[2, 0];

              jac24[11, 0] = dxdxi[0, 0] * d3xdxi[0, 2] + 3 * ddxdxi[0, 2] * ddxdxi[0, 5] + 3 * dxdxi[0, 2] * d3xdxi[0, 8];
              jac24[11, 1] = dxdxi[1, 0] * d3xdxi[1, 2] + 3 * ddxdxi[1, 2] * ddxdxi[1, 5] + 3 * dxdxi[1, 2] * d3xdxi[1, 8];
              jac24[11, 2] = dxdxi[2, 0] * d3xdxi[2, 2] + 3 * ddxdxi[2, 2] * ddxdxi[2, 5] + 3 * dxdxi[2, 2] * d3xdxi[2, 8];
              jac24[11, 3] = 3 * dxdxi[0, 2] * d3xdxi[1, 8] + 3 * dxdxi[1, 2] * d3xdxi[0, 8] + 3 * ddxdxi[0, 2] * ddxdxi[1, 5] + 3 * ddxdxi[0, 5] * ddxdxi[1, 2] + dxdxi[0, 0] * d3xdxi[1, 2] + dxdxi[1, 0] * d3xdxi[0, 2];
              jac24[11, 4] = 3 * dxdxi[1, 2] * d3xdxi[2, 8] + 3 * dxdxi[2, 2] * d3xdxi[1, 8] + 3 * ddxdxi[1, 2] * ddxdxi[2, 5] + 3 * ddxdxi[1, 5] * ddxdxi[2, 2] + dxdxi[1, 0] * d3xdxi[2, 2] + dxdxi[2, 0] * d3xdxi[1, 2];
              jac24[11, 5] = 3 * dxdxi[0, 2] * d3xdxi[2, 8] + 3 * dxdxi[2, 2] * d3xdxi[0, 8] + 3 * ddxdxi[0, 2] * ddxdxi[2, 5] + 3 * ddxdxi[0, 5] * ddxdxi[2, 2] + dxdxi[0, 0] * d3xdxi[2, 2] + dxdxi[2, 0] * d3xdxi[0, 2];

              jac24[12, 0] = ddxdxi[0, 0] * ddxdxi[0, 4] + 2 * ddxdxi[0, 3] * ddxdxi[0, 5] + dxdxi[0, 1] * d3xdxi[0, 7] + dxdxi[0, 2] * d3xdxi[0, 3] + 2 * dxdxi[0, 0] * d3xdxi[0, 9];
              jac24[12, 1] = ddxdxi[1, 0] * ddxdxi[1, 4] + 2 * ddxdxi[1, 3] * ddxdxi[1, 5] + dxdxi[1, 1] * d3xdxi[1, 7] + dxdxi[1, 2] * d3xdxi[1, 3] + 2 * dxdxi[1, 0] * d3xdxi[1, 9];
              jac24[12, 2] = ddxdxi[2, 0] * ddxdxi[2, 4] + 2 * ddxdxi[2, 3] * ddxdxi[2, 5] + dxdxi[2, 1] * d3xdxi[2, 7] + dxdxi[2, 2] * d3xdxi[2, 3] + 2 * dxdxi[2, 0] * d3xdxi[2, 9];
              jac24[12, 3] = 2 * dxdxi[1, 0] * d3xdxi[0, 9] + 2 * dxdxi[0, 0] * d3xdxi[1, 9] + ddxdxi[0, 0] * ddxdxi[1, 4] + 2 * ddxdxi[0, 3] * ddxdxi[1, 5] + 2 * ddxdxi[0, 5] * ddxdxi[1, 3] + ddxdxi[1, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * d3xdxi[1, 7] + dxdxi[0, 2] * d3xdxi[1, 3] + dxdxi[1, 1] * d3xdxi[0, 7] + dxdxi[1, 2] * d3xdxi[0, 3];
              jac24[12, 4] = 2 * dxdxi[2, 0] * d3xdxi[1, 9] + 2 * dxdxi[1, 0] * d3xdxi[2, 9] + ddxdxi[1, 0] * ddxdxi[2, 4] + 2 * ddxdxi[1, 3] * ddxdxi[2, 5] + 2 * ddxdxi[1, 5] * ddxdxi[2, 3] + ddxdxi[2, 0] * ddxdxi[1, 4] + dxdxi[1, 1] * d3xdxi[2, 7] + dxdxi[1, 2] * d3xdxi[2, 3] + dxdxi[2, 1] * d3xdxi[1, 7] + dxdxi[2, 2] * d3xdxi[1, 3];
              jac24[12, 5] = 2 * dxdxi[2, 0] * d3xdxi[0, 9] + 2 * dxdxi[0, 0] * d3xdxi[2, 9] + ddxdxi[0, 0] * ddxdxi[2, 4] + 2 * ddxdxi[0, 3] * ddxdxi[2, 5] + 2 * ddxdxi[0, 5] * ddxdxi[2, 3] + ddxdxi[2, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * d3xdxi[2, 7] + dxdxi[0, 2] * d3xdxi[2, 3] + dxdxi[2, 1] * d3xdxi[0, 7] + dxdxi[2, 2] * d3xdxi[0, 3];

              jac24[13, 0] = ddxdxi[0, 1] * ddxdxi[0, 5] + 2 * ddxdxi[0, 3] * ddxdxi[0, 4] + dxdxi[0, 0] * d3xdxi[0, 5] + dxdxi[0, 2] * d3xdxi[0, 4] + 2 * dxdxi[0, 1] * d3xdxi[0, 9];
              jac24[13, 1] = ddxdxi[1, 1] * ddxdxi[1, 5] + 2 * ddxdxi[1, 3] * ddxdxi[1, 4] + dxdxi[1, 0] * d3xdxi[1, 5] + dxdxi[1, 2] * d3xdxi[1, 4] + 2 * dxdxi[1, 1] * d3xdxi[1, 9];
              jac24[13, 2] = ddxdxi[2, 1] * ddxdxi[2, 5] + 2 * ddxdxi[2, 3] * ddxdxi[2, 4] + dxdxi[2, 0] * d3xdxi[2, 5] + dxdxi[2, 2] * d3xdxi[2, 4] + 2 * dxdxi[2, 1] * d3xdxi[2, 9];
              jac24[13, 3] = 2 * dxdxi[1, 1] * d3xdxi[0, 9] + 2 * dxdxi[0, 1] * d3xdxi[1, 9] + ddxdxi[0, 1] * ddxdxi[1, 5] + 2 * ddxdxi[0, 3] * ddxdxi[1, 4] + 2 * ddxdxi[0, 4] * ddxdxi[1, 3] + ddxdxi[1, 1] * ddxdxi[0, 5] + dxdxi[0, 0] * d3xdxi[1, 5] + dxdxi[0, 2] * d3xdxi[1, 4] + dxdxi[1, 0] * d3xdxi[0, 5] + dxdxi[1, 2] * d3xdxi[0, 4];
              jac24[13, 4] = 2 * dxdxi[2, 1] * d3xdxi[1, 9] + 2 * dxdxi[1, 1] * d3xdxi[2, 9] + ddxdxi[1, 1] * ddxdxi[2, 5] + 2 * ddxdxi[1, 3] * ddxdxi[2, 4] + 2 * ddxdxi[1, 4] * ddxdxi[2, 3] + ddxdxi[2, 1] * ddxdxi[1, 5] + dxdxi[1, 0] * d3xdxi[2, 5] + dxdxi[1, 2] * d3xdxi[2, 4] + dxdxi[2, 0] * d3xdxi[1, 5] + dxdxi[2, 2] * d3xdxi[1, 4];
              jac24[13, 5] = 2 * dxdxi[2, 1] * d3xdxi[0, 9] + 2 * dxdxi[0, 1] * d3xdxi[2, 9] + ddxdxi[0, 1] * ddxdxi[2, 5] + 2 * ddxdxi[0, 3] * ddxdxi[2, 4] + 2 * ddxdxi[0, 4] * ddxdxi[2, 3] + ddxdxi[2, 1] * ddxdxi[0, 5] + dxdxi[0, 0] * d3xdxi[2, 5] + dxdxi[0, 2] * d3xdxi[2, 4] + dxdxi[2, 0] * d3xdxi[0, 5] + dxdxi[2, 2] * d3xdxi[0, 4];

              jac24[14, 0] = ddxdxi[0, 2] * ddxdxi[0, 3] + 2 * ddxdxi[0, 4] * ddxdxi[0, 5] + dxdxi[0, 1] * d3xdxi[0, 8] + dxdxi[0, 0] * d3xdxi[0, 6] + 2 * dxdxi[0, 2] * d3xdxi[0, 9];
              jac24[14, 1] = ddxdxi[1, 2] * ddxdxi[1, 3] + 2 * ddxdxi[1, 4] * ddxdxi[1, 5] + dxdxi[1, 1] * d3xdxi[1, 8] + dxdxi[1, 0] * d3xdxi[1, 6] + 2 * dxdxi[1, 2] * d3xdxi[1, 9];
              jac24[14, 2] = ddxdxi[2, 2] * ddxdxi[2, 3] + 2 * ddxdxi[2, 4] * ddxdxi[2, 5] + dxdxi[2, 1] * d3xdxi[2, 8] + dxdxi[2, 0] * d3xdxi[2, 6] + 2 * dxdxi[2, 2] * d3xdxi[2, 9];

              jac24[14, 3] = 2 * dxdxi[1, 2] * d3xdxi[0, 9] + 2 * dxdxi[0, 2] * d3xdxi[1, 9] + ddxdxi[0, 2] * ddxdxi[1, 3] + 2 * ddxdxi[0, 4] * ddxdxi[1, 5] + 2 * ddxdxi[0, 5] * ddxdxi[1, 4] + ddxdxi[1, 2] * ddxdxi[0, 3] + dxdxi[0, 1] * d3xdxi[1, 8] + dxdxi[0, 0] * d3xdxi[1, 6] + dxdxi[1, 1] * d3xdxi[0, 8] + dxdxi[1, 0] * d3xdxi[0, 6];
              jac24[14, 4] = 2 * dxdxi[2, 2] * d3xdxi[1, 9] + 2 * dxdxi[1, 2] * d3xdxi[2, 9] + ddxdxi[1, 2] * ddxdxi[2, 3] + 2 * ddxdxi[1, 4] * ddxdxi[2, 5] + 2 * ddxdxi[1, 5] * ddxdxi[2, 4] + ddxdxi[2, 2] * ddxdxi[1, 3] + dxdxi[1, 1] * d3xdxi[2, 8] + dxdxi[1, 0] * d3xdxi[2, 6] + dxdxi[2, 1] * d3xdxi[1, 8] + dxdxi[2, 0] * d3xdxi[1, 6];
              jac24[14, 5] = 2 * dxdxi[2, 2] * d3xdxi[0, 9] + 2 * dxdxi[0, 2] * d3xdxi[2, 9] + ddxdxi[0, 2] * ddxdxi[2, 3] + 2 * ddxdxi[0, 4] * ddxdxi[2, 5] + 2 * ddxdxi[0, 5] * ddxdxi[2, 4] + ddxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[0, 1] * d3xdxi[2, 8] + dxdxi[0, 0] * d3xdxi[2, 6] + dxdxi[2, 1] * d3xdxi[0, 8] + dxdxi[2, 0] * d3xdxi[0, 6];

              DoubleMatrix term24 = MatrixFunctions.Product(ddNdX, jac24.Transpose());

              DoubleMatrix jac34 = new DoubleMatrix(15, 10);

              jac34[0, 0] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 0];
              jac34[0, 1] = 6 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 0];
              jac34[0, 2] = 6 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[2, 0];
              jac34[0, 3] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 0] + 12 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 0];
              jac34[0, 4] = 6 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 0] + 12 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[1, 0];
              jac34[0, 5] = 6 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[2, 0] + 12 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[1, 0];
              jac34[0, 6] = 6 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[1, 0] + 12 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[2, 0];
              jac34[0, 7] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[2, 0] + 12 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[0, 0];
              jac34[0, 8] = 6 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[0, 0] + 12 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[2, 0];
              jac34[0, 9] = 12 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[2, 0] + 12 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[1, 0] + 12 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[0, 0];

              jac34[1, 0] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 1];
              jac34[1, 1] = 6 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 1];
              jac34[1, 2] = 6 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[2, 1];
              jac34[1, 3] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 1] + 12 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 1];
              jac34[1, 4] = 6 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 1] + 12 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[1, 1];
              jac34[1, 5] = 6 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[2, 1] + 12 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[1, 1];
              jac34[1, 6] = 6 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[1, 1] + 12 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[2, 1];
              jac34[1, 7] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[2, 1] + 12 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[0, 1];
              jac34[1, 8] = 6 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[0, 1] + 12 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[2, 1];
              jac34[1, 9] = 12 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[2, 1] + 12 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[1, 1] + 12 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[0, 1];

              jac34[2, 0] = 6 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[0, 2];
              jac34[2, 1] = 6 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[1, 2];
              jac34[2, 2] = 6 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[2, 2];
              jac34[2, 3] = 6 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[1, 2] + 12 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[0, 2];
              jac34[2, 4] = 6 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[0, 2] + 12 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[1, 2];
              jac34[2, 5] = 6 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[2, 2] + 12 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[1, 2];
              jac34[2, 6] = 6 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[1, 2] + 12 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[2, 2];
              jac34[2, 7] = 6 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[2, 2] + 12 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[0, 2];
              jac34[2, 8] = 6 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[0, 2] + 12 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[2, 2];
              jac34[2, 9] = 12 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[2, 2] + 12 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[1, 2] + 12 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[0, 2];

              jac34[3, 0] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 3] + 3 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[0, 0];
              jac34[3, 1] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 3] + 3 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[1, 0];
              jac34[3, 2] = 3 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[2, 3] + 3 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[2, 0];
              jac34[3, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 3] + 3 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[1, 0] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 3] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 0];
              jac34[3, 4] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 3] + 3 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[0, 0] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[1, 3] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[1, 0] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[1, 0];
              jac34[3, 5] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[2, 3] + 3 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[2, 0] + 6 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[1, 3] + 3 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[1, 0] + 3 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[1, 0];
              jac34[3, 6] = 3 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[1, 3] + 3 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[1, 0] + 6 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[2, 3] + 3 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[2, 0] + 3 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[2, 0];
              jac34[3, 7] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[2, 3] + 3 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[2, 0] + 6 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[0, 3] + 3 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[0, 0] + 3 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[0, 0];
              jac34[3, 8] = 3 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[0, 3] + 3 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[0, 0] + 6 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[2, 3] + 3 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[2, 0] + 3 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[2, 0];
              jac34[3, 9] = 3 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[0, 0] + 3 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[0, 0] + 3 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[2, 0] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[2, 0] + 6 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[0, 3] + 6 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[1, 3] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[2, 3];

              jac34[4, 0] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[0, 3] + dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 1];
              jac34[4, 1] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[1, 3] + dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 1];
              jac34[4, 2] = dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[2, 3] + dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[2, 1];
              jac34[4, 3] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 1] + dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[1, 3] + 4 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 3] + 4 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 3];
              jac34[4, 4] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 1] + dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 0] + 2 * dxdxi[1, 0] * dxdxi[0, 0] * ddxdxi[1, 1] + 2 * dxdxi[1, 1] * dxdxi[0, 1] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[0, 3] + 4 * dxdxi[1, 0] * dxdxi[0, 1] * ddxdxi[1, 3] + 4 * dxdxi[1, 1] * dxdxi[0, 0] * ddxdxi[1, 3];
              jac34[4, 5] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[2, 1] + dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[1, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[2, 3] + 4 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[1, 3] + 4 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[1, 3];
              jac34[4, 6] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[1, 1] + dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[1, 0] + 2 * dxdxi[2, 0] * dxdxi[1, 0] * ddxdxi[2, 1] + 2 * dxdxi[2, 1] * dxdxi[1, 1] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[1, 3] + 4 * dxdxi[2, 0] * dxdxi[1, 1] * ddxdxi[2, 3] + 4 * dxdxi[2, 1] * dxdxi[1, 0] * ddxdxi[2, 3];
              jac34[4, 7] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[2, 1] + dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[0, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[2, 3] + 4 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[0, 3] + 4 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[0, 3];
              jac34[4, 8] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[0, 1] + dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[0, 0] + 2 * dxdxi[2, 0] * dxdxi[0, 0] * ddxdxi[2, 1] + 2 * dxdxi[2, 1] * dxdxi[0, 1] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[0, 3] + 4 * dxdxi[2, 0] * dxdxi[0, 1] * ddxdxi[2, 3] + 4 * dxdxi[2, 1] * dxdxi[0, 0] * ddxdxi[2, 3];
              jac34[4, 9] = 2 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[0, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[0, 0] + 4 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[0, 3] + 4 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[0, 3] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[1, 0] + 4 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[1, 3] + 4 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[1, 3] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[2, 0] + 4 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[2, 3] + 4 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[2, 3] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[2, 1];

              jac34[5, 0] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 3] + 3 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[0, 1];
              jac34[5, 1] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 3] + 3 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[1, 1];
              jac34[5, 2] = 3 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[2, 3] + 3 * dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[2, 1];
              jac34[5, 3] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 3] + 3 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[1, 1] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 3] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 1];
              jac34[5, 4] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 3] + 3 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[0, 1] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[1, 3] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[1, 1] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[1, 1];
              jac34[5, 5] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[2, 3] + 3 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[2, 1] + 6 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[1, 3] + 3 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[1, 1] + 3 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[1, 1];
              jac34[5, 6] = 3 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[1, 3] + 3 * dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[1, 1] + 6 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[2, 3] + 3 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[2, 1] + 3 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[2, 1];
              jac34[5, 7] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[2, 3] + 3 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[2, 1] + 6 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[0, 3] + 3 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[0, 1] + 3 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[0, 1];
              jac34[5, 8] = 3 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[0, 3] + 3 * dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[0, 1] + 6 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[2, 3] + 3 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[2, 1] + 3 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[2, 1];
              jac34[5, 9] = 3 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[0, 1] + 3 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[0, 1] + 3 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[2, 1] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[2, 1] + 6 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[0, 3] + 6 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[1, 3] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[2, 3];

              jac34[6, 0] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 4] + 3 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[0, 1];
              jac34[6, 1] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 4] + 3 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[1, 1];
              jac34[6, 2] = 3 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[2, 4] + 3 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[2, 1];
              jac34[6, 3] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 4] + 3 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[1, 1] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 4] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[0, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[0, 1];
              jac34[6, 4] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 4] + 3 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[0, 1] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[1, 4] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[1, 1] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[1, 1];
              jac34[6, 5] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[2, 4] + 3 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[2, 1] + 6 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[1, 4] + 3 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[1, 1] + 3 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[1, 1];
              jac34[6, 6] = 3 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[1, 4] + 3 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[1, 1] + 6 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[2, 4] + 3 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[2, 1] + 3 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[2, 1];
              jac34[6, 7] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[2, 4] + 3 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[2, 1] + 6 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[0, 4] + 3 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[0, 1] + 3 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[0, 1];
              jac34[6, 8] = 3 * dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[0, 4] + 3 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[0, 1] + 6 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[2, 4] + 3 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[2, 1] + 3 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[2, 1];
              jac34[6, 9] = 3 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[0, 1] + 3 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[0, 1] + 3 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[2, 1] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[2, 1] + 6 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[0, 4] + 6 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[1, 4] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[2, 4];

              jac34[7, 0] = dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[0, 1] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[0, 4] + dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 2];
              jac34[7, 1] = dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[1, 1] + 4 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[1, 4] + dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 2];
              jac34[7, 2] = dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[2, 1] + 4 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[2, 4] + dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[2, 2];
              jac34[7, 3] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 2] + dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[0, 1] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[1, 4] + 4 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[0, 4] + 4 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[0, 4];
              jac34[7, 4] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 2] + dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[0, 1] + 2 * dxdxi[1, 1] * dxdxi[0, 1] * ddxdxi[1, 2] + 2 * dxdxi[1, 2] * dxdxi[0, 2] * ddxdxi[1, 1] + 4 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[0, 4] + 4 * dxdxi[1, 1] * dxdxi[0, 2] * ddxdxi[1, 4] + 4 * dxdxi[1, 2] * dxdxi[0, 1] * ddxdxi[1, 4];
              jac34[7, 5] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[2, 2] + dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[1, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[1, 1] + 4 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[2, 4] + 4 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[1, 4] + 4 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[1, 4];
              jac34[7, 6] = dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[1, 2] + dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[1, 1] + 2 * dxdxi[2, 1] * dxdxi[1, 1] * ddxdxi[2, 2] + 2 * dxdxi[2, 2] * dxdxi[1, 2] * ddxdxi[2, 1] + 4 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[1, 4] + 4 * dxdxi[2, 1] * dxdxi[1, 2] * ddxdxi[2, 4] + 4 * dxdxi[2, 2] * dxdxi[1, 1] * ddxdxi[2, 4];
              jac34[7, 7] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[2, 2] + dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[0, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[0, 1] + 4 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[2, 4] + 4 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[0, 4] + 4 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[0, 4];
              jac34[7, 8] = dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[0, 2] + dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[0, 1] + 2 * dxdxi[2, 1] * dxdxi[0, 1] * ddxdxi[2, 2] + 2 * dxdxi[2, 2] * dxdxi[0, 2] * ddxdxi[2, 1] + 4 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[0, 4] + 4 * dxdxi[2, 1] * dxdxi[0, 2] * ddxdxi[2, 4] + 4 * dxdxi[2, 2] * dxdxi[0, 1] * ddxdxi[2, 4];
              jac34[7, 9] = 2 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[0, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[0, 1] + 4 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[0, 4] + 4 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[1, 1] + 4 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[1, 4] + 4 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[1, 4] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[2, 1] + 4 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[2, 4] + 4 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[2, 4] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[2, 2];

              jac34[8, 0] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[0, 4] + 3 * dxdxi[0, 2] * dxdxi[0, 1] * ddxdxi[0, 2];
              jac34[8, 1] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[1, 4] + 3 * dxdxi[1, 2] * dxdxi[1, 1] * ddxdxi[1, 2];
              jac34[8, 2] = 3 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[2, 4] + 3 * dxdxi[2, 2] * dxdxi[2, 1] * ddxdxi[2, 2];
              jac34[8, 3] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[1, 4] + 3 * dxdxi[0, 2] * dxdxi[0, 1] * ddxdxi[1, 2] + 6 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[0, 4] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[0, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[0, 2];
              jac34[8, 4] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[0, 4] + 3 * dxdxi[1, 2] * dxdxi[1, 1] * ddxdxi[0, 2] + 6 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[1, 4] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[1, 2] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[1, 2];
              jac34[8, 5] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[2, 4] + 3 * dxdxi[1, 2] * dxdxi[1, 1] * ddxdxi[2, 2] + 6 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[1, 4] + 3 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[1, 2] + 3 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[1, 2];
              jac34[8, 6] = 3 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[1, 4] + 3 * dxdxi[2, 2] * dxdxi[2, 1] * ddxdxi[1, 2] + 6 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[2, 4] + 3 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[2, 2] + 3 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[2, 2];
              jac34[8, 7] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[2, 4] + 3 * dxdxi[0, 2] * dxdxi[0, 1] * ddxdxi[2, 2] + 6 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[0, 4] + 3 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[0, 2] + 3 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[0, 2];
              jac34[8, 8] = 3 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[0, 4] + 3 * dxdxi[2, 2] * dxdxi[2, 1] * ddxdxi[0, 2] + 6 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[2, 4] + 3 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[2, 2] + 3 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[2, 2];
              jac34[8, 9] = 3 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[0, 2] + 3 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[0, 2] + 3 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[1, 2] + 3 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[1, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[2, 2] + 3 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[2, 2] + 6 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[0, 4] + 6 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[1, 4] + 6 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[2, 4];

              jac34[9, 0] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 5] + 3 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[0, 0];
              jac34[9, 1] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 5] + 3 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[1, 0];
              jac34[9, 2] = 3 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[2, 5] + 3 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[2, 0];
              jac34[9, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 5] + 3 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[1, 0] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 5] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[0, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[0, 0];
              jac34[9, 4] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 5] + 3 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[0, 0] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[1, 5] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[1, 0] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[1, 0];
              jac34[9, 5] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[2, 5] + 3 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[2, 0] + 6 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[1, 5] + 3 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[1, 0] + 3 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[1, 0];
              jac34[9, 6] = 3 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[1, 5] + 3 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[1, 0] + 6 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[2, 5] + 3 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[2, 0] + 3 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[2, 0];
              jac34[9, 7] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[2, 5] + 3 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[2, 0] + 6 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[0, 5] + 3 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[0, 0] + 3 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[0, 0];
              jac34[9, 8] = 3 * dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[0, 5] + 3 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[0, 0] + 6 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[2, 5] + 3 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[2, 0] + 3 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[2, 0];
              jac34[9, 9] = 3 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[0, 0] + 3 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[0, 0] + 3 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[2, 0] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[2, 0] + 6 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[0, 5] + 6 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[1, 5] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[2, 5];

              jac34[10, 0] = dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[0, 5] + dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 2];
              jac34[10, 1] = dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[1, 5] + dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 2];
              jac34[10, 2] = dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[2, 5] + dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[2, 2];
              jac34[10, 3] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 2] + dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[1, 5] + 4 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[0, 5] + 4 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[0, 5];
              jac34[10, 4] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 2] + dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[0, 0] + 2 * dxdxi[1, 0] * dxdxi[0, 0] * ddxdxi[1, 2] + 2 * dxdxi[1, 2] * dxdxi[0, 2] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[0, 5] + 4 * dxdxi[1, 0] * dxdxi[0, 2] * ddxdxi[1, 5] + 4 * dxdxi[1, 2] * dxdxi[0, 0] * ddxdxi[1, 5];
              jac34[10, 5] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[2, 2] + dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[1, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[2, 5] + 4 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[1, 5] + 4 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[1, 5];
              jac34[10, 6] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[1, 2] + dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[1, 0] + 2 * dxdxi[2, 0] * dxdxi[1, 0] * ddxdxi[2, 2] + 2 * dxdxi[2, 2] * dxdxi[1, 2] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[1, 5] + 4 * dxdxi[2, 0] * dxdxi[1, 2] * ddxdxi[2, 5] + 4 * dxdxi[2, 2] * dxdxi[1, 0] * ddxdxi[2, 5];
              jac34[10, 7] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[2, 2] + dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[0, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[2, 5] + 4 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[0, 5] + 4 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[0, 5];
              jac34[10, 8] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[0, 2] + dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[0, 0] + 2 * dxdxi[2, 0] * dxdxi[0, 0] * ddxdxi[2, 2] + 2 * dxdxi[2, 2] * dxdxi[0, 2] * ddxdxi[2, 0] + 4 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[0, 5] + 4 * dxdxi[2, 0] * dxdxi[0, 2] * ddxdxi[2, 5] + 4 * dxdxi[2, 2] * dxdxi[0, 0] * ddxdxi[2, 5];
              jac34[10, 9] = 2 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[0, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[0, 0] + 4 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[0, 5] + 4 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[1, 0] + 4 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[1, 5] + 4 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[1, 5] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[2, 0] + 4 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[2, 5] + 4 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[2, 5] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[2, 2];

              jac34[11, 0] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[0, 5] + 3 * dxdxi[0, 2] * dxdxi[0, 0] * ddxdxi[0, 2];
              jac34[11, 1] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[1, 5] + 3 * dxdxi[1, 2] * dxdxi[1, 0] * ddxdxi[1, 2];
              jac34[11, 2] = 3 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[2, 5] + 3 * dxdxi[2, 2] * dxdxi[2, 0] * ddxdxi[2, 2];
              jac34[11, 3] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[1, 5] + 3 * dxdxi[0, 2] * dxdxi[0, 0] * ddxdxi[1, 2] + 6 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[0, 5] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[0, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[0, 2];
              jac34[11, 4] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[0, 5] + 3 * dxdxi[1, 2] * dxdxi[1, 0] * ddxdxi[0, 2] + 6 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[1, 5] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[1, 2] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[1, 2];
              jac34[11, 5] = 3 * dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[2, 5] + 3 * dxdxi[1, 2] * dxdxi[1, 0] * ddxdxi[2, 2] + 6 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[1, 5] + 3 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[1, 2] + 3 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[1, 2];
              jac34[11, 6] = 3 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[1, 5] + 3 * dxdxi[2, 2] * dxdxi[2, 0] * ddxdxi[1, 2] + 6 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[2, 5] + 3 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[2, 2] + 3 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[2, 2];
              jac34[11, 7] = 3 * dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[2, 5] + 3 * dxdxi[0, 2] * dxdxi[0, 0] * ddxdxi[2, 2] + 6 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[0, 5] + 3 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[0, 2] + 3 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[0, 2];
              jac34[11, 8] = 3 * dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[0, 5] + 3 * dxdxi[2, 2] * dxdxi[2, 0] * ddxdxi[0, 2] + 6 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[2, 5] + 3 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[2, 2] + 3 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[2, 2];
              jac34[11, 9] = 3 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[0, 2] + 3 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[0, 2] + 3 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[1, 2] + 3 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[1, 2] + 3 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[2, 2] + 3 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[2, 2] + 6 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[0, 5] + 6 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[1, 5] + 6 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[2, 5];

              jac34[12, 0] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 4] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[0, 3] + dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[0, 0];
              jac34[12, 1] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[1, 3] + dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[1, 0];
              jac34[12, 2] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[2, 4] + 2 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[2, 3] + dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[2, 0];
              jac34[12, 3] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 4] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[1, 5] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[1, 3] + dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 5] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[0, 3] + 2 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[0, 3] + dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[0, 0] + dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[0, 0];
              jac34[12, 4] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 4] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[0, 5] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[0, 3] + dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[0, 0] + 2 * dxdxi[1, 0] * dxdxi[0, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 1] * dxdxi[0, 0] * ddxdxi[1, 5] + 2 * dxdxi[1, 2] * dxdxi[0, 0] * ddxdxi[1, 3] + 2 * dxdxi[1, 0] * dxdxi[0, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[0, 2] * ddxdxi[1, 3] + dxdxi[1, 1] * dxdxi[0, 2] * ddxdxi[1, 0] + dxdxi[1, 2] * dxdxi[0, 1] * ddxdxi[1, 0];
              jac34[12, 5] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[2, 3] + dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[2, 0] + 2 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[1, 5] + 2 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[1, 3] + 2 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[1, 3] + dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[1, 0] + dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[1, 0];
              jac34[12, 6] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[1, 3] + dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[1, 0] + 2 * dxdxi[2, 0] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[2, 1] * dxdxi[1, 0] * ddxdxi[2, 5] + 2 * dxdxi[2, 2] * dxdxi[1, 0] * ddxdxi[2, 3] + 2 * dxdxi[2, 0] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[1, 2] * ddxdxi[2, 3] + dxdxi[2, 1] * dxdxi[1, 2] * ddxdxi[2, 0] + dxdxi[2, 2] * dxdxi[1, 1] * ddxdxi[2, 0];
              jac34[12, 7] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[2, 4] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[2, 5] + 2 * dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[2, 3] + dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[2, 0] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[0, 5] + 2 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[0, 3] + 2 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[0, 0] + dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[0, 0];
              jac34[12, 8] = dxdxi[2, 0] * dxdxi[2, 0] * ddxdxi[0, 4] + 2 * dxdxi[2, 0] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[0, 0] + 2 * dxdxi[2, 0] * dxdxi[0, 0] * ddxdxi[2, 4] + 2 * dxdxi[2, 1] * dxdxi[0, 0] * ddxdxi[2, 5] + 2 * dxdxi[2, 2] * dxdxi[0, 0] * ddxdxi[2, 3] + 2 * dxdxi[2, 0] * dxdxi[0, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[0, 2] * ddxdxi[2, 3] + dxdxi[2, 1] * dxdxi[0, 2] * ddxdxi[2, 0] + dxdxi[2, 2] * dxdxi[0, 1] * ddxdxi[2, 0];
              jac34[12, 9] = 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[2, 3] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[2, 5] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[2, 3] + 2 * dxdxi[0, 0] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[1, 5] + 2 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[1, 3] + 2 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[1, 3] + 2 * dxdxi[1, 0] * dxdxi[2, 0] * ddxdxi[0, 4] + 2 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[0, 5] + 2 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[0, 3] + 2 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[2, 0] + dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[2, 0] + dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[1, 0] + dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[1, 0] + dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[0, 0] + dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[0, 0];

              jac34[13, 0] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[0, 3] + dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[0, 1];
              jac34[13, 1] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[1, 3] + dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[1, 1];
              jac34[13, 2] = dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[2, 4] + 2 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[2, 3] + dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[2, 1];
              jac34[13, 3] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 5] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[1, 4] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[1, 3] + dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 4] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[0, 3] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[0, 3] + dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[0, 1] + dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[0, 1];
              jac34[13, 4] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 5] + 2 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[0, 4] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[0, 3] + dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[0, 1] + 2 * dxdxi[1, 1] * dxdxi[0, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[0, 1] * ddxdxi[1, 4] + 2 * dxdxi[1, 2] * dxdxi[0, 1] * ddxdxi[1, 3] + 2 * dxdxi[1, 1] * dxdxi[0, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 1] * dxdxi[0, 2] * ddxdxi[1, 3] + dxdxi[1, 0] * dxdxi[0, 2] * ddxdxi[1, 1] + dxdxi[1, 2] * dxdxi[0, 0] * ddxdxi[1, 1];
              jac34[13, 5] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[1, 1] * dxdxi[1, 2] * ddxdxi[2, 3] + dxdxi[1, 0] * dxdxi[1, 2] * ddxdxi[2, 1] + 2 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[1, 4] + 2 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[1, 3] + 2 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[1, 3] + dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[1, 1] + dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[1, 1];
              jac34[13, 6] = dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[1, 3] + dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[1, 1] + 2 * dxdxi[2, 1] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[1, 1] * ddxdxi[2, 4] + 2 * dxdxi[2, 2] * dxdxi[1, 1] * ddxdxi[2, 3] + 2 * dxdxi[2, 1] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[2, 1] * dxdxi[1, 2] * ddxdxi[2, 3] + dxdxi[2, 0] * dxdxi[1, 2] * ddxdxi[2, 1] + dxdxi[2, 2] * dxdxi[1, 0] * ddxdxi[2, 1];
              jac34[13, 7] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[2, 5] + 2 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[2, 4] + 2 * dxdxi[0, 1] * dxdxi[0, 2] * ddxdxi[2, 3] + dxdxi[0, 0] * dxdxi[0, 2] * ddxdxi[2, 1] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[0, 4] + 2 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[0, 3] + 2 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[0, 4] + 2 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[0, 1] + dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[0, 1];
              jac34[13, 8] = dxdxi[2, 1] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[0, 4] + 2 * dxdxi[2, 1] * dxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[2, 0] * dxdxi[2, 2] * ddxdxi[0, 1] + 2 * dxdxi[2, 1] * dxdxi[0, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[0, 1] * ddxdxi[2, 4] + 2 * dxdxi[2, 2] * dxdxi[0, 1] * ddxdxi[2, 3] + 2 * dxdxi[2, 1] * dxdxi[0, 0] * ddxdxi[2, 4] + 2 * dxdxi[2, 1] * dxdxi[0, 2] * ddxdxi[2, 3] + dxdxi[2, 0] * dxdxi[0, 2] * ddxdxi[2, 1] + dxdxi[2, 2] * dxdxi[0, 0] * ddxdxi[2, 1];
              jac34[13, 9] = 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[2, 3] + 2 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[2, 4] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[2, 3] + 2 * dxdxi[0, 1] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[1, 4] + 2 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[1, 3] + 2 * dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[1, 3] + 2 * dxdxi[1, 1] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[0, 4] + 2 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[0, 3] + 2 * dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[0, 4] + 2 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[0, 3] + dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[2, 1] + dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[2, 1] + dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[1, 1] + dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[1, 1] + dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[0, 1] + dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[0, 1];

              jac34[14, 0] = dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[0, 3] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[0, 2];
              jac34[14, 1] = dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[1, 3] + 2 * dxdxi[1, 2] * dxdxi[1, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 2] * dxdxi[1, 0] * ddxdxi[1, 4] + dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[1, 2];
              jac34[14, 2] = dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[2, 3] + 2 * dxdxi[2, 2] * dxdxi[2, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 2] * dxdxi[2, 0] * ddxdxi[2, 4] + dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[2, 2];
              jac34[14, 3] = dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[1, 3] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * ddxdxi[1, 5] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * ddxdxi[1, 4] + dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[1, 2] + 2 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[0, 3] + 2 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[0, 4] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 2] + dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 2];
              jac34[14, 4] = dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[0, 3] + 2 * dxdxi[1, 2] * dxdxi[1, 1] * ddxdxi[0, 5] + 2 * dxdxi[1, 2] * dxdxi[1, 0] * ddxdxi[0, 4] + dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[0, 2] + 2 * dxdxi[1, 2] * dxdxi[0, 2] * ddxdxi[1, 3] + 2 * dxdxi[1, 1] * dxdxi[0, 2] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[0, 2] * ddxdxi[1, 4] + 2 * dxdxi[1, 2] * dxdxi[0, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 2] * dxdxi[0, 0] * ddxdxi[1, 4] + dxdxi[1, 1] * dxdxi[0, 0] * ddxdxi[1, 2] + dxdxi[1, 0] * dxdxi[0, 1] * ddxdxi[1, 2];
              jac34[14, 5] = dxdxi[1, 2] * dxdxi[1, 2] * ddxdxi[2, 3] + 2 * dxdxi[1, 2] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[1, 2] * dxdxi[1, 0] * ddxdxi[2, 4] + dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[2, 2] + 2 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[1, 3] + 2 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[1, 5] + 2 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[1, 4] + 2 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[1, 4] + dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[1, 2] + dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[1, 2];
              jac34[14, 6] = dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[1, 3] + 2 * dxdxi[2, 2] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[2, 2] * dxdxi[2, 0] * ddxdxi[1, 4] + dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[1, 2] + 2 * dxdxi[2, 2] * dxdxi[1, 2] * ddxdxi[2, 3] + 2 * dxdxi[2, 1] * dxdxi[1, 2] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[1, 2] * ddxdxi[2, 4] + 2 * dxdxi[2, 2] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 2] * dxdxi[1, 0] * ddxdxi[2, 4] + dxdxi[2, 1] * dxdxi[1, 0] * ddxdxi[2, 2] + dxdxi[2, 0] * dxdxi[1, 1] * ddxdxi[2, 2];
              jac34[14, 7] = dxdxi[0, 2] * dxdxi[0, 2] * ddxdxi[2, 3] + 2 * dxdxi[0, 2] * dxdxi[0, 1] * ddxdxi[2, 5] + 2 * dxdxi[0, 2] * dxdxi[0, 0] * ddxdxi[2, 4] + dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[2, 2] + 2 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[0, 3] + 2 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[0, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[0, 4] + 2 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[0, 2] + dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[0, 2];
              jac34[14, 8] = dxdxi[2, 2] * dxdxi[2, 2] * ddxdxi[0, 3] + 2 * dxdxi[2, 2] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[2, 2] * dxdxi[2, 0] * ddxdxi[0, 4] + dxdxi[2, 1] * dxdxi[2, 0] * ddxdxi[0, 2] + 2 * dxdxi[2, 2] * dxdxi[0, 2] * ddxdxi[2, 3] + 2 * dxdxi[2, 1] * dxdxi[0, 2] * ddxdxi[2, 5] + 2 * dxdxi[2, 0] * dxdxi[0, 2] * ddxdxi[2, 4] + 2 * dxdxi[2, 2] * dxdxi[0, 1] * ddxdxi[2, 5] + 2 * dxdxi[2, 2] * dxdxi[0, 0] * ddxdxi[2, 4] + dxdxi[2, 1] * dxdxi[0, 0] * ddxdxi[2, 2] + dxdxi[2, 0] * dxdxi[0, 1] * ddxdxi[2, 2];
              jac34[14, 9] = 2 * dxdxi[0, 2] * dxdxi[1, 2] * ddxdxi[2, 3] + 2 * dxdxi[0, 2] * dxdxi[1, 1] * ddxdxi[2, 5] + 2 * dxdxi[0, 2] * dxdxi[1, 0] * ddxdxi[2, 4] + 2 * dxdxi[0, 1] * dxdxi[1, 2] * ddxdxi[2, 5] + 2 * dxdxi[0, 0] * dxdxi[1, 2] * ddxdxi[2, 4] + 2 * dxdxi[0, 2] * dxdxi[2, 2] * ddxdxi[1, 3] + 2 * dxdxi[0, 1] * dxdxi[2, 2] * ddxdxi[1, 5] + 2 * dxdxi[0, 0] * dxdxi[2, 2] * ddxdxi[1, 4] + 2 * dxdxi[0, 2] * dxdxi[2, 1] * ddxdxi[1, 5] + 2 * dxdxi[0, 2] * dxdxi[2, 0] * ddxdxi[1, 4] + 2 * dxdxi[1, 2] * dxdxi[2, 2] * ddxdxi[0, 3] + 2 * dxdxi[1, 1] * dxdxi[2, 2] * ddxdxi[0, 5] + 2 * dxdxi[1, 0] * dxdxi[2, 2] * ddxdxi[0, 4] + 2 * dxdxi[1, 2] * dxdxi[2, 1] * ddxdxi[0, 5] + 2 * dxdxi[1, 2] * dxdxi[2, 0] * ddxdxi[0, 4] + dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[2, 2] + dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[2, 2] + dxdxi[0, 1] * dxdxi[2, 0] * ddxdxi[1, 2] + dxdxi[0, 0] * dxdxi[2, 1] * ddxdxi[1, 2] + dxdxi[1, 1] * dxdxi[2, 0] * ddxdxi[0, 2] + dxdxi[1, 0] * dxdxi[2, 1] * ddxdxi[0, 2];

              DoubleMatrix term34 = MatrixFunctions.Product(d3NdX, jac34.Transpose());

              DoubleMatrix term4 = d4Ndxi - term14 - term24 - term34;
              d4NdX = MatrixFunctions.Product(MatrixFunctions.Inverse(jac14), term4.Transpose()).Transpose();

              gpsijk.SetValue(DataInGausspoint.d4NdX, d4NdX);
            }
            MaterialProperty Gc = Material.GetProperty(MaterialPropertyName.CriticalEnergyReleaseRate);
            if (Gc != null)
            {
              gpsijk.SetValue(DataInGausspoint.currentPhase, 0.0);
              gpsijk.SetValue(DataInGausspoint.lastPhase, 0.0);
              gpsijk.SetValue(DataInGausspoint.previousPhase, 0.0);
              gpsijk.SetValue(DataInGausspoint.xiEpsilon, 0.0);
              gpsijk.SetValue(DataInGausspoint.currentStress, new DoubleVector(6));
              gpsijk.SetValue(DataInGausspoint.currentStrain, new DoubleVector(6));
              gpsijk.SetValue(DataInGausspoint.lastStress, new DoubleVector(6));
              gpsijk.SetValue(DataInGausspoint.lastStrain, new DoubleVector(6));
            }

            if (Material is PlasticityMaterial)
            {
              gpsijk.SetValue(DataInGausspoint.lastAlpha, 0.0);
              gpsijk.SetValue(DataInGausspoint.lastBackStress, new DoubleVector(6));
              gpsijk.SetValue(DataInGausspoint.lastPlasticStrain, new DoubleVector(6));
              gpsijk.SetValue(DataInGausspoint.lastStress, new DoubleVector(6));
            }
            //else if (Material is FGMOneGradedDirectionMaterial)
            //{
            //  double[] point = null;
            //  int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
            //  if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
            //  {
            //    point = PointAt(xi1, xi2, xi3);
            //  }
            //  else
            //  {
            //    point = new double[] { xi1, xi2, xi3 };
            //  }
            //  MaterialProperty matE = Material.GetProperty(MaterialPropertyName.YoungModulus);
            //  MaterialProperty matNu = Material.GetProperty(MaterialPropertyName.PoissonRatio);
            //  MaterialProperty matK = Material.GetProperty(MaterialPropertyName.IsotropicThermalConductivity);
            //  gpsijk.SetValue(DataInGausspoint.EModulus, (matE != null) ? matE.GetValueProperty(point[direction]) : 0.0);
            //  gpsijk.SetValue(DataInGausspoint.nu, (matNu != null) ? matNu.GetValueProperty(point[direction]) : 0.0);
            //  gpsijk.SetValue(DataInGausspoint.ThermalConductivity, (matK != null) ? matK.GetValueProperty(point[direction]) : 0.0);

            //  //gps[i, j, k].EModulus = (matE != null) ? matE.GetValueProperty(point[direction]) : 0;
            //  //gps[i, j, k].nu = (matNu != null) ? matNu.GetValueProperty(point[direction]) : 0;
            //  //gps[i, j, k].ThermalConductivity = (matK != null) ? matK.GetValueProperty(point[direction]) : 0;
            //  MaterialProperty matEpsilon = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion);
            //  gpsijk.SetValue(DataInGausspoint.CoefficientsThermalExpansion, (matEpsilon != null) ? matEpsilon.GetValueProperty(point[direction]) : 0.0);

            //  //gps[i, j, k].coefficientsthermalexpansion = (matEpsilon != null) ? matEpsilon.GetValueProperty(point[direction]) : 0;
            //}
            else if (Material is FGMUserDefinedGradedMaterial)
            {
              MaterialProperty matE = Material.GetProperty(MaterialPropertyName.YoungModulus);
              MaterialProperty matNu = Material.GetProperty(MaterialPropertyName.PoissonRatio);
              gpsijk.SetValue(DataInGausspoint.EModulus, (matE != null) ? ComputeParameterProperty(MaterialPropertyName.YoungModulus, xi1, xi2, xi3) : 0.0);
              gpsijk.SetValue(DataInGausspoint.nu, (matNu != null) ? ComputeParameterProperty(MaterialPropertyName.PoissonRatio, xi1, xi2, xi3) : 0.0);

              MaterialProperty matRho = Material.GetProperty(MaterialPropertyName.Density);
              if (matRho != null)
                gpsijk.SetValue(DataInGausspoint.Density, ComputeParameterProperty(MaterialPropertyName.Density, xi1, xi2, xi3));
              MaterialProperty matK = Material.GetProperty(MaterialPropertyName.IsotropicThermalConductivity);
              if (matK != null)
                gpsijk.SetValue(DataInGausspoint.ThermalConductivity, ComputeParameterProperty(MaterialPropertyName.IsotropicThermalConductivity, xi1, xi2, xi3));

              MaterialProperty matEpsilon = Material.GetProperty(MaterialPropertyName.CoefficientThermalExpansion);
              if (matEpsilon != null)
                gpsijk.SetValue(DataInGausspoint.CoefficientsThermalExpansion, ComputeParameterProperty(MaterialPropertyName.CoefficientThermalExpansion, xi1, xi2, xi3));

              MaterialProperty matNonlocal = Material.GetProperty(MaterialPropertyName.NonLocalParameter);
              if (matNonlocal != null)
                gpsijk.SetValue(DataInGausspoint.NonLocalParameter, ComputeParameterProperty(MaterialPropertyName.NonLocalParameter, xi1, xi2, xi3));

              MaterialProperty matGeneralNonlocal1 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter1);
              MaterialProperty matGeneralNonlocal2 = Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter2);
              if ((matGeneralNonlocal1 != null) && (matGeneralNonlocal2 != null))
              {
                gpsijk.SetValue(DataInGausspoint.GeneralNonLocalParameter1, ComputeParameterProperty(MaterialPropertyName.GeneralNonLocalParameter1, xi1, xi2, xi3));
                gpsijk.SetValue(DataInGausspoint.GeneralNonLocalParameter2, ComputeParameterProperty(MaterialPropertyName.GeneralNonLocalParameter2, xi1, xi2, xi3));
              }
            }
          }
        }
      }
    }

    protected DoubleMatrix Jacobian14At(DoubleMatrix d4Ndxi)
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      int nen = patch.GetCountLocalBasisFunctions();
      int d = patch.GetCountDimension();
      DoubleMatrix d4xdxi = new DoubleMatrix(d, 15);//15 là xxxx,yyyy,zzzz,xxxy,xxyy,xyyy,yyyz,yyzz,yzzz,xxxz,xxzz,xzzz,xxyz,xyyz,xyzz //d==3
      var cps = patch.GetVolume().ControlPoints;
      for (int i = 0; i < d; i++)//d==3
      {
        for (int k = 0; k < nen; k++)
        {
          int ien = patch.GetIEN(id, k);
          var controlPointCoord = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(ien, 2)].GetCoordinate();
          for (int j = 0; j < 15; j++)
          {
            d4xdxi[i, j] += d4Ndxi[k, j] * controlPointCoord[i];
          }
        }
      }
      return d4xdxi;
    }
    protected DoubleMatrix Jacobian13At(DoubleMatrix d3Ndxi)
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      int nen = patch.GetCountLocalBasisFunctions();
      int d = patch.GetCountDimension();
      DoubleMatrix d3xdxi = new DoubleMatrix(d, 10);//10 là xxx,yyy,zzz,xxy,xyy,yyz,yzz,xxz,xzz,xyz //d==3
      var cps = patch.GetVolume().ControlPoints;
      for (int i = 0; i < d; i++)//d==3
      {
        for (int k = 0; k < nen; k++)
        {
          int ien = patch.GetIEN(id, k);
          var controlPointCoord = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(ien, 2)].GetCoordinate();
          for (int j = 0; j < 10; j++)
          {
            d3xdxi[i, j] += d3Ndxi[k, j] * controlPointCoord[i];
          }
        }
      }
      return d3xdxi;
    }

    protected DoubleMatrix Jacobian2At(DoubleMatrix ddNdxi)
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      int nen = patch.GetCountLocalBasisFunctions();
      int d = patch.GetCountDimension();
      DoubleMatrix ddxdxi = new DoubleMatrix(d, 6);//6 là xx,yy,zz,xy,yz,xz //d==3
      var cps = patch.GetVolume().ControlPoints;
      for (int i = 0; i < d; i++)//d==3
      {
        for (int k = 0; k < nen; k++)
        {
          int ien = patch.GetIEN(id, k);
          var controlPointCoord = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(ien, 2)].GetCoordinate();
          for (int j = 0; j < 6; j++)
          {
            ddxdxi[i, j] += ddNdxi[k, j] * controlPointCoord[i];
          }
        }
      }
      return ddxdxi;
    }
    public double ComputeParameterProperty(MaterialPropertyName name, double xi, double eta, double zeta)
    {
      double v;
      if (Material is FGMUserDefinedGradedMaterial)
      {
        double[] point = null;
        if (((FGMUserDefinedGradedMaterial)Material).GetIsGlobal())
        {
          point = PointAt(xi, eta, zeta);
        }
        else
        {
          point = new double[] { xi, eta, zeta };
        }
        if (Material is FGMOneGradedDirectionMaterial)
        {
          int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
          v = Material.GetProperty(name).GetValueProperty(point[direction]);
        }
        else
          v = Material.GetProperty(name).GetValueProperty(point);
      }
      else
      {
        v = Material.GetProperty(name).GetValueProperty();
      }
      return v;
    }

    public override void ComputeMaterialPropertyValueAtGaussPoint(MaterialPropertyName name)
    {
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {
            if (Material.GetProperty(name) != null)
            {
              double xi1 = gps[i, j, k].location[0];
              double xi2 = gps[i, j, k].location[1];
              double xi3 = gps[i, j, k].location[2];
              //if (Material is FGMOneGradedDirectionMaterial)
              //{

              //  double[] point = null;
              //  int direction = ((FGMOneGradedDirectionMaterial)Material).GetDirectionGraded();
              //  if (((FGMOneGradedDirectionMaterial)Material).GetIsGlobal())
              //  {
              //    point = PointAt(xi1, xi2, xi3);
              //  }
              //  else
              //  {
              //    point = new double[] { xi1, xi2, xi3 };
              //  }
              //  gps[i, j, k].SetValue(DataInGausspoint.MaterialPropertyValue, Material.GetProperty(name).GetValueProperty(point[direction]));
              //  //gps[i, j, k].materialPropertyValue = Material.GetProperty(name).GetValueProperty(point[direction]);

              //}
              //else
              //{
              //  gps[i, j, k].SetValue(DataInGausspoint.MaterialPropertyValue, Material.GetProperty(name).GetValueProperty());
              //  //gps[i, j, k].materialPropertyValue = Material.GetProperty(name).GetValueProperty();
              //}
              gps[i, j, k].SetValue(DataInGausspoint.MaterialPropertyValue, ComputeParameterProperty(name, xi1, xi2, xi3));
            }
          }
        }
      }
    }
    public override void ComputeDrawValueAtGaussPoint(DataInGausspoint name)
    {
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {
            GaussPoints gpsij = gps[i, j, k];
            gps[i, j, k].SetValue(DataInGausspoint.DrawValue, gpsij.GetValue(name));
            //switch (name)
            //{
            //  case GaussPoints.Variable.xiEpsilon:
            //    gps[i, j, k].DrawValue = gps[i, j, k].xiEpsilon;
            //    break;
            //  case GaussPoints.Variable.currentPhase:
            //    gps[i, j, k].DrawValue = gps[i, j, k].currentPhase;
            //    break;
            //}
          }
        }
      }
    }
    public override int[] GetTArrayGlobal()
    {
      return volume.GetTArrayGlobal();
    }

    public double[] PointAt(double xi, double eta, double zeta)
    {
      return ((NURBSVolume)(patch.GetGeometry(0))).PointAt(xi, eta, zeta);
    }

    public override bool IsInsideRegion(IRegion loc)
    {
      var paramU = GetParameterTwoEndElement(0);
      var paramV = GetParameterTwoEndElement(1);
      var paramW = GetParameterTwoEndElement(2);
      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
          for (int k = 0; k < 2; k++)
          {
            var p = PointAt(paramU[i], paramV[j], paramW[k]);
            if (loc.IsContain(p[0], p[1], p[2]))
              return true;
          }
      return false;
    }

    public override bool IsCompleteInsideRegion(IRegion loc)
    {
      var paramU = GetParameterTwoEndElement(0);
      var paramV = GetParameterTwoEndElement(1);
      var paramW = GetParameterTwoEndElement(2);
      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
          for (int k = 0; k < 2; k++)
          {
            var p = PointAt(paramU[i], paramV[j], paramW[k]);
            if (!loc.IsContain(p[0], p[1], p[2]))
              return false;
          }
      return true;
    }

    public Volume GetVolume()
    {
      return volume;
    }

    /// <summary>
    /// Compute derivatives wrt local coords
    /// dNdxi
    /// </summary>
    /// <param name="xi1"></param>
    /// <param name="xi2"></param>
    /// <param name="xi3"></param>
    /// <returns></returns>
    public override DoubleMatrix GradBasisFunction(params double[] xi)
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      ///////// sua file goc Tien
      int d = patch.GetCountDimension();
      //////////
      NURBSVolume volume = patch.GetVolume();
      var basis = (TrivariateNURBSBasisFunction)volume.Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int nen = patch.GetCountLocalBasisFunctions(0);// (p + 1) * (q + 1); // number of local basis functions
      int span1 = basis.FindSpan(xi[0], 0);
      int span2 = basis.FindSpan(xi[1], 1);
      int span3 = basis.FindSpan(xi[2], 2);
      DoubleMatrix dNdxi = new DoubleMatrix(nen, d);
      var gradBasis = basis.GetDerivativeTrivariateBasisFunctions(xi[0], xi[1], xi[2], 1);
      for (int i = 0; i < nen; i++)
      {
        int ien = patch.GetIEN(0, id, i);
        double[,,] v = gradBasis[patch.GetINC(0, ien, 0) - span1 + p, patch.GetINC(0, ien, 1) - span2 + q, patch.GetINC(0, ien, 2) - span3 + r];
        dNdxi[i, 0] = v[1, 0, 0];
        dNdxi[i, 1] = v[0, 1, 0];
        dNdxi[i, 2] = v[0, 0, 1];
      }
      return dNdxi;
    }
    public override DoubleVector ValueBasisFunction(params double[] xi)
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      NURBSVolume volume = patch.GetVolume();
      var basis = (TrivariateNURBSBasisFunction)volume.Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int nen = patch.GetCountLocalBasisFunctions(0);// (p + 1) * (q + 1); // number of local basis functions
      int span1 = basis.FindSpan(xi[0], 0);
      int span2 = basis.FindSpan(xi[1], 1);
      int span3 = basis.FindSpan(xi[2], 2);
      DoubleVector Ni = new DoubleVector(nen);
      var valueBasis = basis.GetValueTrivariateBasisFunctions(xi[0], xi[1], xi[2]);
      for (int i = 0; i < nen; i++)
      {
        int ien = patch.GetIEN(0, id, i);
        Ni[i] = valueBasis[patch.GetINC(0, ien, 0) - span1 + p, patch.GetINC(0, ien, 1) - span2 + q, patch.GetINC(0, ien, 2) - span3 + r];
      }
      return Ni;
    }

    /// <summary>
    /// Compute the jacobian matrix
    /// dxdxi
    /// </summary>
    /// <param name="dNdxi">Derivatives wrt local coords dNdxi</param>
    /// <returns></returns>
    public override DoubleMatrix JacobianAt(DoubleMatrix dNdxi)
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      int nen = patch.GetCountLocalBasisFunctions();
      int d = patch.GetCountDimension();//3
      DoubleMatrix dxdxi = new DoubleMatrix(d, d);
      var cps = patch.GetVolume().ControlPoints;
      for (int j = 0; j < d; j++)
      {
        for (int k = 0; k < nen; k++)
        {
          int ien = patch.GetIEN(id, k);
          var controlPointCoord = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(ien, 2)].GetCoordinate();
          for (int i = 0; i < d; i++)
          {
            dxdxi[i, j] += dNdxi[k, j] * controlPointCoord[i];
          }
        }
      }
      return dxdxi;
    }

    public GaussPoints GetGaussPoint(int i, int j, int k)
    { return gps[i, j, k]; }

    public double ComupteVolumeOfElement()
    {
      double volume = 0;
      var d = GetPatch().GetCountDimension();
      int id = GetID();
      AbstractPatch3D patch = (AbstractPatch3D)GetPatch();
      NURBSVolume vol = patch.GetVolume();
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
      int numberOfGaussPointOnEachDirection = GetNumberOfGaussPointOnEachDirection();
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
            DoubleMatrix dNdxi = (DoubleMatrix)GetGaussPoint(i, j, k).GetValue(DataInGausspoint.dNdxi);//.dNdxi;
            DoubleMatrix J = JacobianAt(dNdxi);
            //double detJ = (double)gps[i, j, k].GetValue(DataInGausspoint.detJ);
            volume += GetGaussPoint(i, j, k).weight * detJbar * Math.Abs(MatrixFunctions.Determinant(J));
          }
        }
      }
      return Math.Abs(volume);
    }
    public override DoubleVector GetDisplacementLocal()
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      var cps = patch.GetVolume().ControlPoints;
      int d = patch.GetCountDimension();
      int nen = patch.GetCountLocalBasisFunctions();
      DoubleVector U = new DoubleVector(d * nen);
      for (int i = 0; i < nen; i++)
      {
        var ien = patch.GetIEN(id, i);
        ControlPoint cp = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(0, ien, 2)];
        U[d * i] = cp.GetResult(Result.UX);
        U[d * i + 1] = cp.GetResult(Result.UY);
        U[d * i + 2] = cp.GetResult(Result.UZ);
      }
      return U;
    }

    public override ControlPoint[] GetControlPointsLocal()
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      var cps = patch.GetVolume().ControlPoints;
      int nen = patch.GetCountLocalBasisFunctions();
      ControlPoint[] cpLocal = new ControlPoint[nen];
      for (int i = 0; i < nen; i++)
      {
        var ien = patch.GetIEN(id, i);
        cpLocal[i] = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(0, ien, 2)];
      }
      return cpLocal;
    }
    public override DoubleVector GetEachVariableLocal(Result variable)
    {
      AbstractPatch3D patch = (AbstractPatch3D)this.patch;
      var cps = patch.GetVolume().ControlPoints;
      int nen = patch.GetCountLocalBasisFunctions();
      DoubleVector U = new DoubleVector(nen);
      for (int i = 0; i < nen; i++)
      {
        var ien = patch.GetIEN(id, i);
        U[i] = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(0, ien, 2)].GetResult(variable);
      }
      return U;
    }


  }
}
