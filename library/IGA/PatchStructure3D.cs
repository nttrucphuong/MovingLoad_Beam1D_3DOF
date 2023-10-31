using System;
using System.Threading.Tasks;
using DEMSoft.NURBS;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public class PatchStructure3D : AbstractPatch3D, IPatchStructure
  {
    public PatchStructure3D(NURBSVolume volume)
        : base(volume, 3)
    {
      switch (AbstractModel.TypeModel)
      {
        case TypeModelProblem.Piezoelectric:
          //case TypeModelProblem.PhaseField:
          countField[0] = 4;
          break;
      }
    }

    public Structure2DState StateStress { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

    public double StrainAt(Result re, params double[] xi)
    {
      double disp = 0;
      int numElem = FindIndexOfElementAt(xi[0], xi[1], xi[2]);
      AbstractElementStructure3D elem = (AbstractElementStructure3D)listElement[numElem];
      DoubleVector strain = elem.StrainAt(xi[0], xi[1], xi[2]);
      switch (re)
      {
        case Result.EPSILONXX:
          disp = strain[0];
          break;
        case Result.EPSILONYY:
          disp = strain[1];
          break;
        case Result.EPSILONZZ:
          disp = strain[2];
          break;
        case Result.EPSILONXY:
          disp = strain[3];
          break;
        case Result.EPSILONYZ:
          disp = strain[4];
          break;
        case Result.EPSILONXZ:
          disp = strain[5];
          break;
        case Result.EPSILONEQV:
          double sxx = strain[0], syy = strain[1], szz = strain[2], sxy = strain[3], syz = strain[4], sxz = strain[5];

          double xx = 2.0 / 3.0 * sxx - 1.0 / 3.0 * syy - 1.0 / 3.0 * szz;
          double yy = -1.0 / 3.0 * sxx + 2.0 / 3.0 * syy - 1.0 / 3.0 * szz;
          double zz = -1.0 / 3.0 * sxx - 1.0 / 3.0 * syy + 2.0 / 3.0 * szz;
          disp = 2.0 / 3.0 * Math.Sqrt(3.0 / 2.0 * (Math.Pow(xx, 2) + Math.Pow(yy, 2) + Math.Pow(zz, 2)) + 3.0 / 4.0 * (Math.Pow(2 * sxy, 2) + Math.Pow(2 * syz, 2) + Math.Pow(2 * sxz, 2)));
          break;
      }
      return disp;
    }

    public double StrainOnGlobalAt(Result re, TypeCoordination coord, params double[] x)
    {
      throw new NotImplementedException();
    }

    public double StressAt(Result re, params double[] xi)
    {
      throw new NotImplementedException();
    }

    public double StressOnGlobalAt(Result re, TypeCoordination coord, params double[] x)
    {
      throw new NotImplementedException();
    }
  }
}
