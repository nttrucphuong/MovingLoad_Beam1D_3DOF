using System;
using System.Collections.Generic;
using DEMSoft.Drawing;
using DEMSoft.Function;
using DEMSoft.NURBS;

namespace DEMSoft.IGA
{
  /// <summary>
  /// A class of Constraint with value of degree of freedom
  /// </summary>
  public abstract class AbstractConstraintValue
  {
    private Object valueConstraint;
    private FunctionRToR piecewiseLoad;
    private int fieldID;
    private Coordination coord;
    /// <summary>
    /// Construct class with value of constraint, if no constraint, so input double.NaN
    /// </summary>
    /// <param name="fieldID">ID of field to be constrainted</param>
    /// <param name="piecewiseLoad">Constraint function coressponse the time</param>
    /// <param name="valueConstraint">Value constraint</param>
    public AbstractConstraintValue(int fieldID, FunctionRToR piecewiseLoad, Object valueConstraint)
    {
      this.valueConstraint = valueConstraint;
      this.piecewiseLoad = piecewiseLoad;
      this.fieldID = fieldID;
      this.coord = new Coordination(Drawing.TypeCoordination.Cartesian);
    }
    public Coordination GetCoordination()
    { return coord; }
    public void SetCoordination(Coordination coord)
    {
      this.coord = coord;
    }
    public int GetIDOfField()
    { return fieldID; }

    /// <summary>
    /// Get constraint at index is value constraint
    /// </summary>
    /// <param name="time">Time of constraint</param>
    /// <returns>value of </returns>
    public double GetValueConstraint(double time)
    {
      if (valueConstraint is double)
        return ((double)valueConstraint) * piecewiseLoad.ValueAt(time);
      else if (valueConstraint is ConstantFunctionRToR)
        return ((ConstantFunctionRToR)valueConstraint).ValueAt(0) * piecewiseLoad.ValueAt(time);
      return 0;
    }

    public Object GetValueConstraint()
    { return valueConstraint; }

    public double GetValueVelocityConstraint(double time)
    {
      if (valueConstraint is double)
        return ((double)valueConstraint) * piecewiseLoad.DerivativeAt(time);
      else if (valueConstraint is ConstantFunctionRToR)
        return ((ConstantFunctionRToR)valueConstraint).ValueAt(0) * piecewiseLoad.DerivativeAt(time);
      return 0;
    }
    public double GetValueAccelerationConstraint(double time)
    {
      if (valueConstraint is double)
        return ((double)valueConstraint) * piecewiseLoad.Derivative().Derivative().ValueAt(time);
      else if (valueConstraint is ConstantFunctionRToR)
        return ((ConstantFunctionRToR)valueConstraint).ValueAt(0) * piecewiseLoad.Derivative().Derivative().ValueAt(time);
      return 0;
    }
    public FunctionRToR GetPiecewiseLoadTime()
    { return piecewiseLoad; }
    public abstract ControlPoint[] GetControlPointsConstrainted();
    public int[] GetTArrayConstrainted()
    {
      ControlPoint[] cpsel = GetControlPointsConstrainted();
      List<int> tArrayConstraint = new List<int>();
      for (int k = 0; k < cpsel.Length; k++)
      {

        int[] tArrayCp = cpsel[k].GetTArrayGlobal();
        int tArray = tArrayCp[GetIDOfField()];
        tArrayConstraint.Add(tArray);
      }
      return tArrayConstraint.ToArray();
    }
    public abstract double[] ComputeValueConstraintOnControlPoints(double time);
  }
}
