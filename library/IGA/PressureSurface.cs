using DEMSoft.NURBS;
using DEMSoft.Function;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Pressure load on face of 3D patch
  /// </summary>
  public class PressureSurface : Pressure
  {
    private FunctionR2ToR press;

    /// <summary>
    /// Constructor class
    /// </summary>
    /// <param name="face">face of element which be applied load</param>
    /// <param name="isInGlobal">Coresponding on global coordinate (x-y) or local coordinate (n-t)</param>
    /// <param name="p">global coordinate (x-y) : p[0] = px, p[1] = py,  local coordinate (n-t) : p[0] = pn, p[1] = pt</param>
    public PressureSurface(Face face, bool isInGlobal, FunctionR2ToR p)
        : base(face, isInGlobal, true)
    { press = p; }

    /// <summary>
    /// Constructor class
    /// </summary>
    /// <param name="face">face of element which be applied load</param>
    /// <param name="isInGlobal">Coresponding on global coordinate (x-y) or local coordinate (n-t)</param>
    /// <param name="p">global coordinate (x-y) : p[0] = px, p[1] = py,  local coordinate (n-t) : p[0] = pn, p[1] = pt</param>
    public PressureSurface(Face face, bool isInGlobal, bool isNaturalCoordinate, FunctionR2ToR p)
        : base(face, isInGlobal, isNaturalCoordinate)
    { press = p; }

    /// <summary>
    /// Compute local load vector
    /// </summary>
    /// <returns></returns>
    public override DoubleVector ComputeLocalLoadVector(double time = 1)
    {
      int d = GetMeshPart().GetNumberOfFields();
      Face face = (Face)GetMeshPart();
      AbstractElement2DPlate elem = (AbstractElement2DPlate)face.GetElement();
      BivariateNURBSBasisFunction basis = (BivariateNURBSBasisFunction)((NURBSSurface)(elem.GetPatch().GetGeometry(0))).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int nen = (p + 1) * (q + 1);
      int n = d * nen;
      DoubleVector re = new DoubleVector(n);

      var paraEndPatchU = elem.GetParameterTwoEndElement(0);
      var paraEndPatchV = elem.GetParameterTwoEndElement(1);

      //double[,] paraEndFace = face.GetParametricEndEdge();
      int numGaussPoint = elem.GetNumberOfGaussPointOnEachDirection();
      for (int k = 0; k < numGaussPoint; k++)
      {
        for (int kk = 0; kk < numGaussPoint; kk++)
        {
          GaussPoints gp = elem.GetGaussPoint(k, kk);
          double xi = gp.location[0];
          double eta = gp.location[1];
          double w1w2 = gp.weight;
          double detJbar = 1.0 / 4.0
                * (paraEndPatchU[1] - paraEndPatchU[0])
                * (paraEndPatchV[1] - paraEndPatchV[0]);
          double[] paraLoad = { 0, 0 };

          if (isNaturalCoordinate)
          {
            paraLoad[0] = xi;
            paraLoad[1] = eta;
          }
          else
          {
            double[] pp = elem.PointAt(xi, eta);
            paraLoad[0] = pp[0];
            paraLoad[1] = pp[1];
          }
          DoubleVector Ni = (DoubleVector)gp.GetValue(DataInGausspoint.Ni);
          double detJ = (double)gp.GetValue(DataInGausspoint.detJ);

          for (int i = 0; i < nen; i++)
          {
            re[i * d + 2] += w1w2 * detJbar * detJ * Ni[i] * press.ValueAt(paraLoad[0], paraLoad[1]);
          }
        }
      }
      return re;
    }
    //public override DoubleVector ComputeLocalLoadVector(double time = 1)
    //{
    //    PatchStructurePlate patch = (PatchStructurePlate)this.GetMeshPart();
    //    int d = patch.GetCountField();
    //    BivariateNURBSBasisFunction basis = (BivariateNURBSBasisFunction)((NURBSSurface)(patch.GetGeometry())).Basis;
    //    int p = basis.GetDegree(0);
    //    int q = basis.GetDegree(1);
    //    int nen = (p + 1) * (q + 1); // number of local basis functions
    //    var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
    //    var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
    //    int idx1 = patch.GetIPN(id, 0);
    //    int idx2 = patch.GetIPN(id, 1);

    //    int d = GetMeshPart().GetNumberOfFields();
    //    Face face = (Face)GetMeshPart();
    //    AbstractElement2DPlate elem = (AbstractElement2DPlate)face.GetElement();
    //    BivariateNURBSBasisFunction basis = (BivariateNURBSBasisFunction)((NURBSSurface)(elem.GetPatch().GetGeometry(0))).Basis;
    //    int p = basis.GetDegree(0);
    //    int q = basis.GetDegree(1);
    //    int nen = (p + 1) * (q + 1);
    //    int n = d * nen;
    //    DoubleVector re = DoubleVector(n);

    //    var paraEndPatchU = elem.GetParameterTwoEndElement(0);
    //    var paraEndPatchV = elem.GetParameterTwoEndElement(1);
    //    double xi = 0;
    //    double eta = 0;
    //    int numGaussPoint = Math.Max(p, q) + 1;
    //    for (int k = 0; k < numGaussPoint; k++)
    //    {
    //        double psi1 = GaussPoints.GetPoint(numGaussPoint, k);
    //        double w1 = GaussPoints.GetWeight(numGaussPoint, k);
    //        for (int kk = 0; kk < numGaussPoint; kk++)
    //        {
    //            double psi2 = GaussPoints.GetPoint(numGaussPoint, kk);
    //            double w2 = GaussPoints.GetWeight(numGaussPoint, kk);
    //            xi = 0.5 * ((paraEndPatchU[1] - paraEndPatchU[0]) * psi1 + paraEndPatchU[1] + paraEndPatchU[0]);
    //            eta = 0.5 * ((paraEndPatchV[1] - paraEndPatchV[0]) * psi2 + paraEndPatchV[1] + paraEndPatchV[0]);
    //            double J2 = 0.25 * (paraEndPatchU[1] - paraEndPatchU[0]) * (paraEndPatchV[1] - paraEndPatchV[0]);

    //            double[] paraLoad = { 0, 0 };

    //            if (isNaturalCoordinate)
    //            {
    //                paraLoad[0] = xi;
    //                paraLoad[1] = eta;
    //            }
    //            else
    //            {
    //                double[] pp = elem.PointAt(xi, eta);
    //                paraLoad[0] = pp[0];
    //                paraLoad[1] = pp[1];
    //            }
    //            DoubleMatrix gradNE = elem.GradBasisFunction(xi, eta);
    //            DoubleMatrix J = elem.JacobianAt(gradNE);
    //            double[,] Ni = face.GetBivariateBasisFunctionOnFace(xi, eta);

    //            for (int i = 0; i < q; i++)
    //            {
    //                for (int j = 0; j < p; j++)
    //                {
    //                    re[(j * (p + 1) + i) * d + 2] += w1 * w2 * J2 * J.Determinant() * press.ValueAt(paraLoad[0], paraLoad[1]);
    //                }
    //            }
    //        }
    //    }
    //    return re;
    //}

    //private double normVectorJacobian(DoubleMatrix J)
    //{
    //    Face face = (Face)GetMeshPart();

    //    double ee = 0;
    //    double gg = 0;
    //    double ff = 0;

    //    switch (face.GetIndexCoordinate())
    //    {
    //        case 0:
    //            ee = J[0, 0] * J[0, 0] + J[1, 0] * J[1, 0] + J[2, 0] * J[2, 0];
    //            gg = J[0, 1] * J[0, 1] + J[1, 1] * J[1, 1] + J[2, 1] * J[2, 1];
    //            ff = J[0, 0] * J[0, 1] + J[1, 0] * J[1, 1] + J[2, 0] * J[2, 1];
    //            break;
    //        case 1:
    //            ee = J[0, 0] * J[0, 0] + J[1, 0] * J[1, 0] + J[2, 0] * J[2, 0];
    //            gg = J[0, 2] * J[0, 2] + J[1, 2] * J[1, 2] + J[2, 2] * J[2, 2];
    //            ff = J[0, 0] * J[0, 2] + J[1, 0] * J[1, 2] + J[2, 0] * J[2, 2];
    //            break;
    //        case 2:
    //            ee = J[0, 1] * J[0, 1] + J[1, 1] * J[1, 1] + J[2, 1] * J[2, 1];
    //            gg = J[0, 2] * J[0, 2] + J[1, 2] * J[1, 2] + J[2, 2] * J[2, 2];
    //            ff = J[0, 1] * J[0, 2] + J[1, 1] * J[1, 2] + J[2, 1] * J[2, 2];
    //            break;
    //    }
    //    return Math.Sqrt(ee * gg - ff * ff);
    //}

    //private DoubleVector GetNormalFace(DoubleMatrix J)
    //{
    //    DoubleVector vec = DoubleVector(3);
    //    Face face = (Face)GetMeshPart();
    //    switch (face.GetIndexCoordinate())
    //    {
    //        case 0:
    //            if (face.GetIndexFrontBack() == 0)
    //            {
    //                vec[0] = J[1, 1] * J[2, 0] - J[2, 1] * J[1, 0];
    //                vec[1] = J[2, 1] * J[0, 0] - J[0, 1] * J[2, 0];
    //                vec[2] = J[0, 1] * J[1, 0] - J[1, 1] * J[0, 0];
    //            }
    //            else
    //            {
    //                vec[0] = -J[1, 1] * J[2, 0] + J[2, 1] * J[1, 0];
    //                vec[1] = -J[2, 1] * J[0, 0] + J[0, 1] * J[2, 0];
    //                vec[2] = -J[0, 1] * J[1, 0] + J[1, 1] * J[0, 0];
    //            }
    //            break;
    //        case 1:
    //            if (face.GetIndexFrontBack() == 0)
    //            {
    //                vec[0] = J[1, 0] * J[2, 2] - J[2, 0] * J[1, 2];
    //                vec[1] = J[2, 0] * J[0, 2] - J[0, 0] * J[2, 2];
    //                vec[2] = J[0, 0] * J[1, 2] - J[1, 0] * J[0, 2];
    //            }
    //            else
    //            {
    //                vec[0] = -J[1, 0] * J[2, 2] + J[2, 0] * J[1, 2];
    //                vec[1] = -J[2, 0] * J[0, 2] + J[0, 0] * J[2, 2];
    //                vec[2] = -J[0, 0] * J[1, 2] + J[1, 0] * J[0, 2];
    //            }
    //            break;
    //        case 2:
    //            if (face.GetIndexFrontBack() == 0)
    //            {
    //                vec[0] = -J[1, 1] * J[2, 2] + J[2, 1] * J[1, 2];
    //                vec[1] = -J[2, 1] * J[0, 2] + J[0, 1] * J[2, 2];
    //                vec[2] = -J[0, 1] * J[1, 2] + J[1, 1] * J[0, 2];
    //            }
    //            else
    //            {
    //                vec[0] = J[1, 1] * J[2, 2] - J[2, 1] * J[1, 2];
    //                vec[1] = J[2, 1] * J[0, 2] - J[0, 1] * J[2, 2];
    //                vec[2] = J[0, 1] * J[1, 2] - J[1, 1] * J[0, 2];
    //            }
    //            break;
    //    }
    //    return vec;
    //}

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
