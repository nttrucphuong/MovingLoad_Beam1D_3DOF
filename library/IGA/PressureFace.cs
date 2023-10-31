using System;
using DEMSoft.NURBS;
using DEMSoft.Function;
using CenterSpace.NMath.Core;
using DEMSoft.EngineeringData;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Pressure load on face of 3D patch
  /// </summary>
  public class PressureFace : Pressure
  {
    private FunctionR2ToR[] press;

    /// <summary>
    /// Constructor class
    /// </summary>
    /// <param name="face">face of element which be applied load</param>
    /// <param name="isInGlobal">Direction coresponding on global coordinate (x-y) or local coordinate (n-t)</param>
    /// <param name="p">global coordinate (x-y) : p[0] = px, p[1] = py,  local coordinate (n-t) : p[0] = pn, p[1] = pt</param>
    public PressureFace(Face face, bool isInGlobal, params FunctionR2ToR[] p)
             : base(face, isInGlobal, true)
    { press = p; }

    /// <summary>
    /// Constructor class
    /// </summary>
    /// <param name="face">face of element which be applied load</param>
    /// <param name="isInGlobal">Direction coresponding on global coordinate (x-y) or local coordinate (n-t)</param>
    /// <param name="isNaturalCoordinate">Value coresponding on natural coordinate or physical coordinate</param>
    /// <param name="p">global coordinate (x-y) : p[0] = px, p[1] = py,  local coordinate (n-t) : p[0] = pn, p[1] = pt</param>
    public PressureFace(Face face, bool isInGlobal, bool isNaturalCoordinate, params FunctionR2ToR[] p)
             : base(face, isInGlobal, isNaturalCoordinate)
    { press = p; }

    /// <summary>
    /// Compute local load vector
    /// </summary>
    /// <returns></returns>
    public override DoubleVector ComputeLocalLoadVector(double time = 1)
    {
      int d = GetMeshPart().CountDimension();
      Face face = (Face)GetMeshPart();
      int p = face.GetDegree(0);
      int q = face.GetDegree(1);
      int nen = (p + 1) * (q + 1);
      int n = d * nen;
      DoubleVector re = new DoubleVector(n);

      var elem = face.GetVolumeBeAttached().GetElement();
      var patch = elem.GetPatch();
      var vol = (NURBSVolume)(patch.GetGeometry(0));
      var paraEndPatchU = face.GetVolumeBeAttached().GetElement().GetParameterTwoEndElement(0);
      var paraEndPatchV = face.GetVolumeBeAttached().GetElement().GetParameterTwoEndElement(1);
      var paraEndPatchW = face.GetVolumeBeAttached().GetElement().GetParameterTwoEndElement(2);
      double xi = 0;
      double eta = 0;
      double zeta = 0;
      //double[,] paraEndFace = face.GetParametricEndEdge();
      int numGaussPoint = Math.Max(p, q) + 1;
      for (int k = 0; k < numGaussPoint; k++)
      {
        double psi1 = GaussPoints.GetPoint(numGaussPoint, k);
        double w1 = GaussPoints.GetWeight(numGaussPoint, k);
        for (int kk = 0; kk < numGaussPoint; kk++)
        {
          double psi2 = GaussPoints.GetPoint(numGaussPoint, kk);
          double w2 = GaussPoints.GetWeight(numGaussPoint, kk);
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

          double[] paraLoad = { 0, 0 };

          if (isNaturalCoordinate)
          {
            switch (face.GetIndexCoordinate())
            {
              case 0:
                paraLoad[0] = xi;
                paraLoad[1] = eta;
                break;
              case 1:
                paraLoad[0] = xi;
                paraLoad[1] = zeta;
                break;
              case 2:
                paraLoad[0] = eta;
                paraLoad[1] = zeta;
                break;
            }
          }
          //else
          //{
          //    double[] pp = vol.PointAt(xi, eta, zeta);
          //    switch (face.GetIndexCoordinate())
          //    {
          //        case 0:
          //            paraLoad[0] = pp[0];
          //            paraLoad[1] = pp[1];
          //            break;
          //        case 1:
          //            paraLoad[0] = pp[0];
          //            paraLoad[1] = pp[2];
          //            break;
          //        case 2:
          //            paraLoad[0] = pp[1];
          //            paraLoad[1] = pp[2];
          //            break;
          //    }
          //}
          //DoubleMatrix ddNdX = (DoubleMatrix)gps[i, j].GetValue(DataInGausspoint.ddNdX);//gps[i, j].ddNdX;

          double[,] Nij = face.GetTrivariateBasisFunctionOnFace(xi, eta, zeta);
          //double[,] ddNdX=face.get 
          //DoubleVector N = new DoubleVector(n);

          //for (int j = 0; j <= q; j++)
          //  for (int i = 0; i <= p; i++)
          //  {
          //    N[(j * (p + 1) + i) * d + 0] = Nij[i, j];
          //    N[(j * (p + 1) + i) * d + 1] = Nij[i, j];
          //    N[(j * (p + 1) + i) * d + 2] = Nij[i, j];
          //  }
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

          DoubleMatrix gradNE = elem.GradBasisFunction(xi, eta, zeta);
          DoubleMatrix J = elem.JacobianAt(gradNE);

          //DoubleMatrix dNij = new DoubleMatrix(face.GetDerivativeTrivariateBasisFunctionOnFace(xi, eta, zeta));
          //DoubleMatrix J = JacobianAt(dNij);
          double nuy = 0;
          double nuy1 = 0;
          double nuy2 = 0;
          MaterialProperty matNonlocal = elem.Material.GetProperty(MaterialPropertyName.NonLocalParameter);
          MaterialProperty matGeneralNonlocal1 = elem.Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter1);
          MaterialProperty matGeneralNonlocal2 = elem.Material.GetProperty(MaterialPropertyName.GeneralNonLocalParameter2);
          if (matNonlocal != null)
          {
            nuy = elem.ComputeParameterProperty(MaterialPropertyName.NonLocalParameter, xi, eta, zeta);
          }
          if (matGeneralNonlocal1 != null && matGeneralNonlocal2 != null)
          {
            nuy1 = elem.ComputeParameterProperty(MaterialPropertyName.GeneralNonLocalParameter1, xi, eta, zeta);
            nuy2 = elem.ComputeParameterProperty(MaterialPropertyName.GeneralNonLocalParameter2, xi, eta, zeta);
          }
          //double nuy = elem.Material.GetProperty(EngineeringData.MaterialPropertyName.NonlocalParameter).GetValueProperty();

          //Jacobian of face mapping
          double normJ = Math.Abs(normVectorJacobian(J));
          DoubleVector nf = GetNormalFace(J);

          DoubleVector normalFace = nf / nf.TwoNorm();
          DoubleMatrix ddNdX = null;
          DoubleMatrix d4NdX = null;
          if (matNonlocal != null)
          {
            ddNdX = face.GetSecondDerivativePhysicalCoordinationOnFace(xi, eta, zeta);//[nen,3]
          }
          if (matGeneralNonlocal1 != null && matGeneralNonlocal2 != null)
          {
            ddNdX = face.GetSecondDerivativePhysicalCoordinationOnFace(xi, eta, zeta);//[nen,3]
            d4NdX = face.GetFourthDerivativePhysicalCoordinationOnFace(xi, eta, zeta);//[nen,3]
          }

          double coef = w1 * w2 * normJ * J2;
          for (int j = 0; j <= q; j++)
            for (int i = 0; i <= p; i++)
            {
              double Nijij = Nij[i, j];
              int idx = (j * (p + 1) + i);
              double NijGrad2 = Nijij;//(Nijij - nuy * LaplaceNij);
              if (matNonlocal != null)
              {
                double LaplaceNij = (ddNdX[idx, 0] + ddNdX[idx, 1]);
                NijGrad2 = (Nijij - nuy * LaplaceNij);
              }
              if (matGeneralNonlocal1 != null && matGeneralNonlocal2 != null)
              {
                double LaplaceNij = (ddNdX[idx, 0] + ddNdX[idx, 1]);
                double laplaceOfLaplaceNij = (d4NdX[idx, 0] + d4NdX[idx, 1] + 2 * d4NdX[idx, 3]);
                NijGrad2 = Nijij - (nuy1 + nuy2) * LaplaceNij + nuy1 * nuy2 * laplaceOfLaplaceNij;
              }
              if (isInGlobal)//Load distribute on x-y direction
              {
                for (int t = 0; t < d; t++)
                {
                  if (press[t] != null)
                  {
                    double pressValue = press[t].ValueAt(paraLoad[0], paraLoad[1]);
                    re[idx * d + t] += coef * NijGrad2 * pressValue;//x y
                  }
                }
              }
              else//n-t direction
              {
                for (int jj = 0; jj < d; jj++)
                {
                  for (int t = 0; t < d; t++)
                  {
                    if (press[t] != null)
                    {
                      switch (t)
                      {
                        case 0://Normal
                               //re[(j * (p + 1) + i) * d + jj] += w1 * w2 * normJ * J2 * Nij[i, j] * press[t].ValueAt(paraLoad[0], paraLoad[1]) * normalFace[jj];//tao, eta
                          double pressValue = press[t].ValueAt(paraLoad[0], paraLoad[1]);
                          re[idx * d + jj] += coef * NijGrad2 * pressValue * normalFace[jj];//tao, eta
                          break;
                      }
                    }
                    break;
                  }
                }
              }
            }
        }
      }
      return re;
    }

    private double normVectorJacobian(DoubleMatrix J)
    {
      Face face = (Face)GetMeshPart();

      double ee = 0;
      double gg = 0;
      double ff = 0;

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
      return Math.Sqrt(ee * gg - ff * ff);
    }

    private DoubleVector GetNormalFace(DoubleMatrix J)
    {
      DoubleVector vec = new DoubleVector(3);
      Face face = (Face)GetMeshPart();
      switch (face.GetIndexCoordinate())
      {
        case 0:
          if (face.GetIndexFrontBack() == 0)
          {
            vec[0] = J[1, 1] * J[2, 0] - J[2, 1] * J[1, 0];
            vec[1] = J[2, 1] * J[0, 0] - J[0, 1] * J[2, 0];
            vec[2] = J[0, 1] * J[1, 0] - J[1, 1] * J[0, 0];
          }
          else
          {
            vec[0] = -J[1, 1] * J[2, 0] + J[2, 1] * J[1, 0];
            vec[1] = -J[2, 1] * J[0, 0] + J[0, 1] * J[2, 0];
            vec[2] = -J[0, 1] * J[1, 0] + J[1, 1] * J[0, 0];
          }
          break;
        case 1:
          if (face.GetIndexFrontBack() == 0)
          {
            vec[0] = J[1, 0] * J[2, 2] - J[2, 0] * J[1, 2];
            vec[1] = J[2, 0] * J[0, 2] - J[0, 0] * J[2, 2];
            vec[2] = J[0, 0] * J[1, 2] - J[1, 0] * J[0, 2];
          }
          else
          {
            vec[0] = -J[1, 0] * J[2, 2] + J[2, 0] * J[1, 2];
            vec[1] = -J[2, 0] * J[0, 2] + J[0, 0] * J[2, 2];
            vec[2] = -J[0, 0] * J[1, 2] + J[1, 0] * J[0, 2];
          }
          break;
        case 2:
          if (face.GetIndexFrontBack() == 0)
          {
            vec[0] = -J[1, 1] * J[2, 2] + J[2, 1] * J[1, 2];
            vec[1] = -J[2, 1] * J[0, 2] + J[0, 1] * J[2, 2];
            vec[2] = -J[0, 1] * J[1, 2] + J[1, 1] * J[0, 2];
          }
          else
          {
            vec[0] = J[1, 1] * J[2, 2] - J[2, 1] * J[1, 2];
            vec[1] = J[2, 1] * J[0, 2] - J[0, 1] * J[2, 2];
            vec[2] = J[0, 1] * J[1, 2] - J[1, 1] * J[0, 2];
          }
          break;
      }
      return vec;
    }

    //private DoubleMatrix JacobianAt(DoubleMatrix dNij)
    //{
    //    int d = GetMeshPart().GetNumberOfFields() - 1;
    //    Face face = (Face)GetMeshPart();
    //    int p1 = face.GetDegree(0);
    //    int p2 = face.GetDegree(1);
    //    int nen = (p1 + 1) * (p2 + 1);
    //    DoubleMatrix grad = new DoubleMatrix(d + 1, d);
    //    var cpsOnFace = face.GetControlPointsOnFace();
    //    var inc = face.CreateIndexCoordinate();
    //    for (int j = 0; j < d; j++)
    //        for (int k = 0; k < nen; k++)
    //        {
    //            var controlPointCoord = cpsOnFace[inc[k, 0], inc[k, 1]].GetCoordinate();
    //            for (int i = 0; i <= d; i++)
    //            {
    //                grad[i, j] += dNij[k, j] * controlPointCoord[i];
    //            }
    //        }
    //    return grad;
    //}
  }
}
