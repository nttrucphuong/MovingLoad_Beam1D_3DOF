using CenterSpace.NMath.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Enumerate of coordination
  /// </summary>
  public enum TypeCoordination
  {
    Cartesian, Cylindrical, Spherical
  }
  //public class Coordination
  //{
  //  public TypeCoordination TypeOfCoordination { get; set; }
  //  public List<ActionCoordination> listOfActionCoordination;
  //  public Coordination(TypeCoordination coord)
  //  {
  //    TypeOfCoordination = coord;
  //    listOfActionCoordination = new List<ActionCoordination>();
  //  }
  //  public void AddActionCoordination(ActionCoordination action)
  //  {
  //    listOfActionCoordination.Add(action);
  //  }
  //  public DoubleMatrix CalculateTransformationMatrix()
  //  {
  //    DoubleMatrix mat = DoubleMatrix.Identity(4);
  //    for (int i = 0; i < listOfActionCoordination.Count; i++)
  //    {
  //      mat = NMathFunctions.Product(listOfActionCoordination[i].CalculateTransformationMatrix(), mat);
  //    }
  //    return mat;
  //  }
  //}

  //public enum TypeActionCoordination
  //{
  //  Translation, Rotation
  //}
  //public enum Direction
  //{
  //  XGlobal, YGlobal, ZGlobal, XLocal, YLocal, ZLocal
  //}
  //public class ActionCoordination
  //{
  //public TypeActionCoordination typeAction { get; }
  //public Direction direction { get; }
  //public double value { get; }
  //public ActionCoordination(TypeActionCoordination typeAction, Direction direction, double value)
  //{
  //  this.typeAction = typeAction;
  //  this.direction = direction;
  //  this.value = value;
  //}
  //public DoubleMatrix CalculateTransformationMatrix()
  //{
  //  DoubleMatrix mat = new DoubleMatrix(4, 4);
  //  mat[3, 3] = 1;
  //  switch (typeAction)
  //  {
  //    case TypeActionCoordination.Translation:
  //      switch (direction)
  //      {
  //        case Direction.XGlobal:
  //          mat[3, 0] = value;
  //          break;
  //        case Direction.YGlobal:
  //          mat[3, 1] = value;
  //          break;
  //        case Direction.ZGlobal:
  //          mat[3, 2] = value;
  //          break;
  //      }
  //      break;
  //    case TypeActionCoordination.Rotation:
  //      switch (direction)
  //      {
  //        case Direction.XGlobal:
  //          mat[0, 0] = 1;
  //          mat[1, 1] = Math.Cos(value);
  //          mat[1, 2] = -Math.Sin(value);
  //          mat[2, 1] = Math.Sin(value);
  //          mat[2, 2] = Math.Cos(value);
  //          break;
  //        case Direction.YGlobal:
  //          mat[1, 1] = 1;
  //          mat[0, 0] = Math.Cos(value);
  //          mat[0, 2] = Math.Sin(value);
  //          mat[2, 0] = -Math.Sin(value);
  //          mat[2, 2] = Math.Cos(value);
  //          break;
  //        case Direction.ZGlobal:
  //          mat[2, 2] = 1;
  //          mat[0, 0] = Math.Cos(value);
  //          mat[0, 1] = -Math.Sin(value);
  //          mat[1, 0] = Math.Sin(value);
  //          mat[1, 1] = Math.Cos(value);
  //          break;
  //      }
  //      break;
  //  }
  //  return mat;
  //}
  //}
}