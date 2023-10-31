using DEMSoft.Function;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using DEMSoft.NURBS;

namespace DEMSoft.IGA
{
   public class ConstraintValueArrayOfControlPoints : AbstractConstraintValue
   {
      private ControlPoint[] cps;
      public ConstraintValueArrayOfControlPoints(ControlPoint[] cps, int fieldID, FunctionRToR piecewiseLoad, double valueConstraint)
         : base(fieldID, piecewiseLoad, valueConstraint)
      {
         this.cps = cps;
      }
      public ConstraintValueArrayOfControlPoints(ControlPoint[] cps, int fieldID, FunctionRToR piecewiseLoad, double[] valueConstraint)
         : base(fieldID, piecewiseLoad, valueConstraint)
      {
         if (cps.Length != valueConstraint.Length)
            throw new ArgumentException("Length of array of control points and valueConstraint must be equal");
         this.cps = cps;
      }
      public override ControlPoint[] GetControlPointsConstrainted()
      {
         return cps;
      }
      public override double[] ComputeValueConstraintOnControlPoints(double time)
      {
         double[] value = null;
         if (GetValueConstraint() is double[])
         {
            double[] valueConstraint = (double[])GetValueConstraint();
            for (int i = 0; i < cps.Length; i++)
            {
               value[i] = valueConstraint[i] * GetPiecewiseLoadTime().ValueAt(time);
            }
         }
         return value;
      }
   }
}
