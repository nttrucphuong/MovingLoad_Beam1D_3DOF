using DEMSoft.Function;
using System;
using DEMSoft.NURBS;

namespace DEMSoft.IGA
{
  public class ConstraintInitialTemperaturePatch : AbstractConstraintValue
  {
    private AbstractPatch patch;
    public ConstraintInitialTemperaturePatch(AbstractPatch patch, int fieldID, double valueConstraint)
       : base(fieldID, new ConstantFunctionRToR(1.0), valueConstraint)
    {
      this.patch = patch;
    }

    public override double[] ComputeValueConstraintOnControlPoints(double time)
    {
      throw new NotImplementedException();
    }

    public override ControlPoint[] GetControlPointsConstrainted()
    {
      return patch.SelectAllControlPoints(0);
    }
    public AbstractPatch GetPatch()
    { return patch; }
  }
}
