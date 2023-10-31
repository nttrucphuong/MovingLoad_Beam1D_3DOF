using System;
using DEMSoft.NURBS;
using DEMSoft.Function;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Pressure load on edge of 2D patch
  /// </summary>
  public class PressureEdgeTime : Pressure
  {
    private FunctionRToR[] press;
    private FunctionRToR piecewiseLoad;

    /// <summary>
    /// Constructor class
    /// </summary>
    /// <param name="e">edge of element which be applied load</param>
    /// <param name="isInGlobal">Direction coresponding on global coordinate (x-y) or local coordinate (n-t)</param>
    /// <param name="piecewiseLoad">load function coressponse the time</param>
    /// <param name="p">global coordinate (x-y) : p[0] = px, p[1] = py,  local coordinate (n-t) : p[0] = pn, p[1] = pt</param>
    public PressureEdgeTime(Edge e, bool isInGlobal, FunctionRToR piecewiseLoad, params FunctionRToR[] p)
             : base(e, isInGlobal, true)
    {
      press = p;
      this.piecewiseLoad = piecewiseLoad;
    }

    /// <summary>
    /// Constructor class
    /// </summary>
    /// <param name="e">edge of element which be applied load</param>
    /// <param name="isInGlobal">Coresponding on global coordinate (x-y) or local coordinate (n-t)</param>
    /// <param name="isNaturalCoordinate">Value coresponding on natural coordinate or physical coordinate</param>
    /// <param name="piecewiseLoad">load function coressponse the time</param>
    /// <param name="p">global coordinate (x-y) : p[0] = px, p[1] = py,  local coordinate (n-t) : p[0] = pn, p[1] = pt</param>
    public PressureEdgeTime(Edge e, bool isInGlobal, bool isNaturalCoordinate, FunctionRToR piecewiseLoad, params FunctionRToR[] p)
             : base(e, isInGlobal, isNaturalCoordinate)
    {
      press = p;
      this.piecewiseLoad = piecewiseLoad;
    }

    /// <summary>
    /// Compute local load vector
    /// </summary>
    /// <param name="time">time of load vector to compute</param>
    /// <returns></returns>
    public override DoubleVector ComputeLocalLoadVector(double time)
    {
      Edge e = (Edge)GetMeshPart();
      int d = e.GetFace().GetElement().GetPatch().GetCountDimension();
      int p = e.GetDegree(0);
      int n = d * (p + 1);
      DoubleVector re = new DoubleVector(n);

      var elem = (AbstractElement2D)e.GetFace().GetElement();
      var patch = (AbstractPatch2D)elem.GetPatch();
      var surface = (NURBSSurface)(patch.GetGeometry(0));
      var basis = surface.Basis;
      var paraEndPatchU = e.GetFace().GetElement().GetParameterTwoEndElement(0);
      var paraEndPatchV = e.GetFace().GetElement().GetParameterTwoEndElement(1);
      double xi = 0;
      double eta = 0;
      double loadTime = piecewiseLoad.ValueAt(time);
      //double[] paraEndEdge = e.GetParametricEndEdge(0);
      int numGaussPoint = p + 1;//patch.GetNumberOfGaussPoint();
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
        double paraLoad;
        if (e.GetIndexCoordinate() == 0)
          paraLoad = xi;
        else
          paraLoad = eta;
        double[] Ni = e.GetBivariateBasisFunctionOnEdge(0, xi, eta);
        double J2 = (e.GetIndexCoordinate() == 0) ? (0.5 * (paraEndPatchU[1] - paraEndPatchU[0]))
             : (0.5 * (paraEndPatchV[1] - paraEndPatchV[0]));

        DoubleVector dNi = new DoubleVector(e.GetDerivativeBivariateBasisFunctionOnEdge(0, xi, eta));
        double[] J = JacobianAt(dNi).ToArray();
        double J1 = Math.Sqrt(J[0] * J[0] + J[1] * J[1]);
        double[] unitTangentVector = { J[0] / J1, J[1] / J1 };
        double[] unitNormalVectorSurface = { 0, 0, 1.0 };//Plane surface
        double[] normalVector = { J[1], -J[0] };
        double[] unitNormalVector = { normalVector[0] / J1, normalVector[1] / J1 };
        double a = 1;
        if (((IPatchStructure)patch).StateStress == Structure2DState.Axisymetric)
        {
          double x1 = elem.PointAt(xi, eta)[0];
          a = 2 * Math.PI * x1;
        }
        else if (((IPatchStructure)patch).StateStress == Structure2DState.PlaneStress)
        {
          a = patch.Thickness;
          if (a == 0)
            throw new InvalidArgumentException("Plane stress condition must had thickness");
        }
        for (int i = 0; i <= p; i++)
        {
          if (isInGlobal)//Load distribute on x-y direction
          {
            for (int j = 0; j < d; j++)
            {
              if (press[j] != null)
              {
                re[i * d + j] += a * loadTime * w * J1 * J2 * Ni[i] * press[j].ValueAt(paraLoad);//x y
              }
            }
          }
          else//n-t direction
          {
            for (int j = 0; j < d; j++)
            {
              for (int jj = 0; jj < press.Length; jj++)
              {
                if (press[jj] != null)
                {
                  if (jj == 0)//Normal
                    re[i * d + j] += a * loadTime * w * J1 * J2 * Ni[i] * press[jj].ValueAt(paraLoad) * unitNormalVector[j];//tao, eta
                  else//Tangent
                    re[i * d + j] += a * loadTime * w * J1 * J2 * Ni[i] * press[jj].ValueAt(paraLoad) * unitTangentVector[j];//tao, eta
                }
              }
            }
          }
        }
      }
      return re;
    }

    private DoubleVector JacobianAt(DoubleVector dNi)
    {
      int d = GetMeshPart().CountDimension();
      Edge e = (Edge)GetMeshPart();
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
  }
}
