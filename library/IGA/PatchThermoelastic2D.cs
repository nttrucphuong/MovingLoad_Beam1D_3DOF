using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public class PatchThermoelastic2D : AbstractPatch2D, IPatchStructure
  {
    public PatchThermoelastic2D(NURBSSurface surface, Structure2DState state)
             : base(surface, TypeStructure.Plane)
    {
      StateStress = state;
    }

    public Structure2DState StateStress
    { get; set; }

    public double StrainAt(Result re, params double[] xi)
    {
      throw new NotImplementedException();
    }

    public double StrainOnGlobalAt(Result re, TypeCoordination coord, params double[] x)
    {
      throw new NotImplementedException();
    }

    public double StressOnGlobalAt(Result re, TypeCoordination coord, params double[] x)
    {
      throw new NotImplementedException();
    }

    /////////////////// Doi Luu
    public void ComputeHeatTransferConvectionPatch(ref DoubleVector rGlobal)
    {
      int nel = CountElements();//(n - p) * (m - q);//number of elements
                                //if (!IsParallelProcesing)
                                //{
      for (int i = 0; i < nel; i++)//Loop through elements
      {
        ElementThermoelastic2D elem = (ElementThermoelastic2D)listElement[i];
        if (elem.GetHeatTransferConvection() != null)
        {
          elem.ComputeHeatTransferConvectionLoadVectorElement(ref rGlobal);
        }
      }
    }
    //////////////////////

    public double StressAt(Result re, params double[] xi)
    {
      double disp = 0;
      int numElem = FindIndexOfElementAt(xi[0], xi[1]);
      ElementThermoelastic2D elem = (ElementThermoelastic2D)listElement[numElem];
      DoubleVector stress = elem.StressAt(xi[0], xi[1]);
      switch (re)
      {
        case Result.SIGMAXX:
          disp = stress[0];
          break;
        case Result.SIGMAYY:
          disp = stress[1];
          break;
        case Result.SIGMAZZ:
          disp = stress[2];
          if (StateStress == Structure2DState.PlaneStrain)
          {
            double nu = elem.Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
            disp = nu * (stress[0] + stress[1]);
          }
          break;
        case Result.SIGMAXY:
          disp = stress[3];
          break;
        case Result.SIGMAEQV:
          double sxx = stress[0], syy = stress[1], szz = stress[2], sxy = stress[3];
          if (StateStress == Structure2DState.PlaneStrain)
          {
            double nu = elem.Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
            szz = nu * (sxx + syy);
          }
          disp = Math.Sqrt(1.0 / 2.0 * (Math.Pow(sxx - syy, 2) + Math.Pow(syy - szz, 2) + Math.Pow(sxx - szz, 2) + 6.0 * Math.Pow(sxy, 2)));
          break;
      }
      return disp;
    }


    //public double GetStrainAt(double xi, double eta, Result re)
    //{
    //    double disp = 0;
    //    var surface = GetSurface();
    //    var cps = surface.ControlPoints;
    //    var basis = (BivariateNURBSBasisFunction)surface.Basis;
    //    int p = basis.GetDegree(0);
    //    int uSpan = basis.GetKnotVector(0).FindSpan(xi, p);
    //    int q = basis.GetDegree(1);
    //    int vSpan = basis.GetKnotVector(1).FindSpan(eta, q);
    //    double[,] Nuv = basis.GetValueBivariateBasisFunctions(xi, eta);
    //    double exx = 0, eyy = 0, exy = 0, ezz = 0;
    //    for (int i = 0; i <= p; i++)
    //        for (int j = 0; j <= q; j++)
    //        {
    //            if (re == Result.EPSILONEQV)
    //            {
    //                exx += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.EPSILONXX);
    //                eyy += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.EPSILONYY);
    //                exy += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.EPSILONXY);
    //                ezz += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.EPSILONZZ);
    //            }
    //            else
    //                disp += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(re);
    //        }
    //    if (re == Result.EPSILONEQV)
    //    {
    //        double xx = 2.0 / 3.0 * exx - 1.0 / 3.0 * eyy - 1.0 / 3.0 * ezz;
    //        double yy = -1.0 / 3.0 * exx + 2.0 / 3.0 * eyy - 1.0 / 3.0 * ezz;
    //        double zz = -1.0 / 3.0 * exx - 1.0 / 3.0 * eyy + 2.0 / 3.0 * ezz;
    //        disp = 2.0 / 3.0 * Math.Sqrt(3.0 / 2.0 * (Math.Pow(xx, 2) + Math.Pow(yy, 2) + Math.Pow(zz, 2)) + 3.0 / 4.0 * (Math.Pow(2 * exy, 2)));
    //    }
    //    return disp;
    //}

  }
}
