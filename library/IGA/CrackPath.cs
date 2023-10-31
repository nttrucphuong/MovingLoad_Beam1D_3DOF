using DEMSoft.NURBS;
using System;

namespace DEMSoft.IGA
{
  public class CrackEdge
  {
    private NURBSCurve pl;
    public int[] IndexPatchHasCrossedCrack
    {
      get;
    }
    private double[,] bound;
    public CrackEdge(NURBSCurve pl, int[] indexPatchHasCrossedCrack)
    {
      this.pl = pl;
      IndexPatchHasCrossedCrack = indexPatchHasCrossedCrack;
    }
    public double DistanceFromPointToCrack(double x, double y, double z, out double lengthFromXiToCrackTip)
    {
      double xiProjection = pl.Projection(x, y, z)[0];
      double[] point = pl.PointAt(xiProjection);
      lengthFromXiToCrackTip = pl.ComputeLength(xiProjection, false);
      return Math.Sqrt(Math.Pow(x - point[0], 2) + Math.Pow(y - point[1], 2) + Math.Pow(z - point[2], 2));
    }
    public double ComputeLength()
    {
      return pl.ComputeLength();
    }
    public double[,] GetBoundBoxCoverCrackEdge(double lengthBigger = 0)
    {

      if (bound == null)
      {
        double[] para = pl.GetParametricOnCurve(5);
        bound = new double[,] { { 100000, -100000 }, { 100000, -100000 } };
        for (int i = 0; i < para.Length; i++)
        {
          double[] point = pl.PointAt(para[i]);
          if (bound[0, 0] > point[0])
            bound[0, 0] = point[0];
          if (bound[0, 1] < point[0])
            bound[0, 1] = point[0];
          if (bound[1, 0] > point[1])
            bound[1, 0] = point[1];
          if (bound[1, 1] < point[1])
            bound[1, 1] = point[1];
        }
      }
      double[,] biggerBound = new double[2, 2];
      biggerBound[0, 0] = bound[0, 0] - lengthBigger;
      biggerBound[0, 1] = bound[0, 1] + lengthBigger;
      biggerBound[1, 0] = bound[1, 0] - lengthBigger;
      biggerBound[1, 1] = bound[1, 1] + lengthBigger;
      return biggerBound;
    }
  }
}
