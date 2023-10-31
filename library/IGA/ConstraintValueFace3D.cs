using DEMSoft.Function;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using DEMSoft.NURBS;

namespace DEMSoft.IGA
{
   public class ConstraintValueFace3D : AbstractConstraintValue
   {
      private AbstractPatch3D patch;
      private int indexFace;
      public ConstraintValueFace3D(AbstractPatch3D patch, int indexFace, int fieldID, FunctionRToR piecewiseLoadTime, double valueConstraint)
         : base(fieldID, piecewiseLoadTime, valueConstraint)
      {
         this.patch = patch;
         this.indexFace = indexFace;
      }
      public AbstractPatch3D GetPatch()
      { return patch; }
      public int GetIndexFace()
      { return indexFace; }

      public override ControlPoint[] GetControlPointsConstrainted()
      {
         ControlPoint[,] cpTemp = patch.SelectEndPatchControlPoints(indexFace);
        ControlPoint[] cpsel = new ControlPoint[cpTemp.GetLength(0) * cpTemp.GetLength(1)];
         int count = 0;
         for (int ii = 0; ii < cpTemp.GetLength(0); ii++)
            for (int jj = 0; jj < cpTemp.GetLength(1); jj++)
            {
               cpsel[count++] = cpTemp[ii, jj];
            }
         return cpsel;
      }

      public override double[] ComputeValueConstraintOnControlPoints(double time)
      {
         throw new NotImplementedException();
      }
   }
}
