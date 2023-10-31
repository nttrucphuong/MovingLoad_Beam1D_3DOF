//using DEMSoft.NURBS;
//using System;

//namespace DEMSoft.IGA
//{
//  /// <summary>
//  /// Void of structure using phase field
//  /// </summary>
//  public class Void2D
//  {
//    private NURBSSurface pl;
//    public int[] IndexPatchHasCrossedCrack
//    {
//      get;
//    }
//    public int IndexBoundaryEdge { get; }
//    private double[,] bound;
//    public Void2D(NURBSSurface pl, int[] indexPatchHasCrossedCrack, int indexBoundaryEdge)
//    {
//      this.pl = pl;
//      IndexPatchHasCrossedCrack = indexPatchHasCrossedCrack;
//      IndexBoundaryEdge = indexBoundaryEdge;
//    }
//    public bool IsInside(double x, double y, double z)
//    {
//      return pl.IsInsideGeometry(x, y, z);
//    }
//    public double DistanceFromOutsidePointToVoid(double x, double y, double z)
//    {
//      //double distance = 0;
//      //if (!pl.IsInsideGeometry(x, y, z))
//      //{
//      //  double[] xiProjection = pl.Projection(x, y, z, IndexBoundaryEdge);
//      //  double[] point = pl.PointAt(xiProjection[0], xiProjection[1]);
//      //  return Math.Sqrt(Math.Pow(x - point[0], 2) + Math.Pow(y - point[1], 2) + Math.Pow(z - point[2], 2));
//      //}
//      double[] xiProjection = pl.Projection(x, y, z);
//      double[] point = pl.PointAt(xiProjection[0], xiProjection[1]);
//      //distance 
//      return Math.Sqrt(Math.Pow(x - point[0], 2) + Math.Pow(y - point[1], 2) + Math.Pow(z - point[2], 2));
//    }
//    //public double DistanceWithSignFromPointToVoid(double x, double y, double z)
//    //{
//    //  double xiProjection = pl.Projection(x, y, 0)[0];
//    //  double[] aProjection = pl.TangentUnitVectorCurve(xiProjection);
//    //  double[] pointProjection = pl.PointAt(xiProjection);
//    //  double d = ((pointProjection[0] - p[0]) * aProjection[1] - (pointProjection[1] - p[1]) * aProjection[0]) * Math.Sqrt(Math.Pow(pointProjection[0] - p[0], 2) + Math.Pow(pointProjection[1] - p[1], 2));
//    //  minDistance = Math.Max(minDistance, d);//tính có dấu nên lấy số dương lớn nhất
//    //}
//    public double[,] GetRegionBoxCoverVoid(double lengthBigger = 0)
//    {
//      if (bound == null)
//      {
//        double[][] para = pl.GetParametricOnSurface(5, 5);
//        bound = new double[,] { { 100000, -100000 }, { 100000, -100000 } };
//        for (int j = 0; j < para[0].Length; j++)
//        {
//          double[] point = pl.PointAt(para[0][j], para[1][j]);
//          if (bound[0, 0] > point[0])
//            bound[0, 0] = point[0];
//          if (bound[0, 1] < point[0])
//            bound[0, 1] = point[0];
//          if (bound[1, 0] > point[1])
//            bound[1, 0] = point[1];
//          if (bound[1, 1] < point[1])
//            bound[1, 1] = point[1];
//        }
//      }
//      double[,] biggerBound = new double[2, 2];
//      biggerBound[0, 0] = bound[0, 0] - lengthBigger;
//      biggerBound[0, 1] = bound[0, 1] + lengthBigger;
//      biggerBound[1, 0] = bound[1, 0] - lengthBigger;
//      biggerBound[1, 1] = bound[1, 1] + lengthBigger;
//      return biggerBound;
//    }
//  }
//}

using DEMSoft.NURBS;
using System;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Void of structure using phase field
  /// </summary>
  public class Void2D
  {
    private NURBSCurve curveBoundVoid;
    public int[] IndexPatchHasCrossedCrack
    {
      get;
    }
    public int IndexBoundaryEdge { get; }
    private double[,] bound;
    public Void2D(NURBSCurve pl, int[] indexPatchHasCrossedCrack, int indexBoundaryEdge)
    {
      this.curveBoundVoid = pl;
      IndexPatchHasCrossedCrack = indexPatchHasCrossedCrack;
      IndexBoundaryEdge = indexBoundaryEdge;
    }
    //public bool IsInside(double x, double y, double z)
    //{
    //  return curveBoundVoid.IsInsideGeometry(x, y, z);
    //}
    //public double DistanceFromOutsidePointToVoid(double x, double y, double z)
    //{
    //  //double distance = 0;
    //  //if (!pl.IsInsideGeometry(x, y, z))
    //  //{
    //  //  double[] xiProjection = pl.Projection(x, y, z, IndexBoundaryEdge);
    //  //  double[] point = pl.PointAt(xiProjection[0], xiProjection[1]);
    //  //  return Math.Sqrt(Math.Pow(x - point[0], 2) + Math.Pow(y - point[1], 2) + Math.Pow(z - point[2], 2));
    //  //}
    //  double[] xiProjection = curveBoundVoid.Projection(x, y, z);
    //  double[] point = curveBoundVoid.PointAt(xiProjection[0], xiProjection[1]);
    //  //distance 
    //  return Math.Sqrt(Math.Pow(x - point[0], 2) + Math.Pow(y - point[1], 2) + Math.Pow(z - point[2], 2));
    //}
    public double DistanceWithSignFromPointToVoid(double x, double y, double z)
    {
      double xiProjection = curveBoundVoid.Projection(x, y, 0)[0];
      double[] aProjection = curveBoundVoid.TangentUnitVectorCurve(xiProjection);
      double[] pointProjection = curveBoundVoid.PointAt(xiProjection);
      return ((pointProjection[0] - x) * aProjection[1] - (pointProjection[1] - y) * aProjection[0]) * Math.Sqrt(Math.Pow(pointProjection[0] - x, 2) + Math.Pow(pointProjection[1] - y, 2));
    }
    public double[,] GetRegionBoxCoverVoid(double lengthBigger = 0)
    {
      if (bound == null)
      {
        double[] para = curveBoundVoid.GetParametricOnCurve(5);
        bound = new double[,] { { 100000, -100000 }, { 100000, -100000 } };
        for (int j = 0; j < para.Length; j++)
        {
          double[] point = curveBoundVoid.PointAt(para[j]);
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

