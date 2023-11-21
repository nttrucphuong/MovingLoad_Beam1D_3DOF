using DEMSoft.Common;
using DEMSoft.Drawing;
using DEMSoft.EngineeringData;
using DEMSoft.Function;
using DEMSoft.IGA;
using DEMSoft.NURBS;
using System;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Reflection;
using System.Text.RegularExpressions;
using System.Windows.Forms;

namespace SampleTesting
{
  internal static class Satic
  {
    /// <summary>
    /// The main entry point for the application.
    /// </summary>
    [STAThread]
    static void Main()
    {
      var targetDir = Path.Combine(@"F:\ppp\CTDT\Nam_tu\LVTN\github\MovingLoad\MovingLoad\bin\Debug\temp");
      if (Directory.Exists(targetDir))
      {
        Directory.Delete(targetDir, true);
      }

      //DOF in each patch: ux,uy, theta

      //double L = 8; //m
      //double h = 1;//4.1602;
      //double b = 1;
      //double E = 100000;

      double L = 34; //m
      double h = 0.2;//4.1602;
      double b = 0.5; //m
      double E = 2.976E+14; //Pa
      double rho = 114000; //kg/m^3
      double I = b * Math.Pow(h, 3) / 12.0;
      double area = b * h;
      //ViewerForm viewer = new ViewerForm(true);
      NURBSCurve curve = GeometryCreator.CreateStraightNURBSCurve(0, 0, 0, L, 0, 0);
      // danh sach cac Curve
      List<NURBSCurve> listCurve = new List<NURBSCurve>();
      curve.pRefinement(3);
      curve.hRefinement(50); //inssert KVector - incerease p, min luoi

      //curve.colorKnot = Color.Red;
      //curve.colorCurve = Color.Blue;
      //curve.colorControlNet = Color.Orange;
      //curve.isDrawOriginal = false;
      //curve.colorControlPoint = Color.Green;
      //curve.Draw(viewer);

      //viewer.UpdateCamera();
      //viewer.Run();


      ////// --------------------Khai bao vat lieu--------------------

      Material steel = new Material("steel");
      steel.AddProperty(new IsotropicElasticity(
          PairOfIsotropicElasticity.YoungModulusAndPoissonRatio, E, 0.2));
      steel.AddProperty(new Density(rho));
      steel.TypeMaterialStructure = TypeMaterialStructure.Elasticity;
      //static
      ModelStructureStatic model = new ModelStructureStatic(Dimension.Beam);

      // dinh nghia trang thai bai toan
      model.AddMaterial(steel);
      model.AddPatch(curve);
      model.AttachMaterialToPatch(0, 0);
      model.SetThicknessBeam(h, 0);
      model.SetTCrossSectionAreaBeam(area, 0);
      model.SetMomentOfInertiaBeam(I, 0);

      // them dai luong can tinh 
      //model.AddComputeResult(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

      ControlPoint[] selAllCPs = ((AbstractPatch1D)model.GetPatch(0)).SelectAllControlPoints(0);
      ControlPoint[] selCPStart = new ControlPoint[] { selAllCPs[0] };
      ControlPoint[] selCPEnd = new ControlPoint[] { selAllCPs[selAllCPs.Length - 1] };

      ///1D-3Dof
      ConstraintValueArrayOfControlPoints cX0 = new ConstraintValueArrayOfControlPoints(selCPStart, 0, new NullFunctionRToR(), 0);
      ConstraintValueArrayOfControlPoints cY0 = new ConstraintValueArrayOfControlPoints(selCPStart, 1, new NullFunctionRToR(), 0);
      ConstraintValueArrayOfControlPoints cThetaZ0 = new ConstraintValueArrayOfControlPoints(selCPStart, 2, new NullFunctionRToR(), 0);

      ConstraintValueArrayOfControlPoints cX1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 0, new NullFunctionRToR(), 0);
      ConstraintValueArrayOfControlPoints cY1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 1, new NullFunctionRToR(), 0);
      //ConstraintValueArrayOfControlPoints cThetaZ1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 2, new NullFunctionRToR(), 0);

      ////apply constraint
      model.AddConstraint(cX0, cY0, cThetaZ0);
      model.AddConstraint(cX1, cY1);

      AbstractModel.IsParallelProcesing = false;

      double centerPos = selAllCPs.Length / 2;
      ControlPoint cpLoad = selAllCPs[(int)Math.Round(centerPos)];
      Force f = new Force(cpLoad, 0, -347000, 0);
      model.AddLoad(f);
      model.IsSaveStepByStep = false;
      model.AddComputeResult(Result.UY);
      model.SetKinematicsFunction(new LinearFunctionRToR(0, 1, 0, 1), 0);
      
      model.InitializePatch();
      model.PreProcessing();
      model.Solve();
      model.PostProcessing();

      //RESULT

      double xi;
      List<double> arrayData = new List<double>();
      NURBSBasisFunction basis = (NURBSBasisFunction)((NURBSCurve)(model.GetPatch(0).GetGeometry(0))).Basis;
      double[] kvNoMulticiply = basis.GetKnotVector().GetKnotVectorNoMultiplicity();
      //double[,] results = ExactSolution(-347000, L, h, b, E, kvNoMulticiply.Length);
      //double[] subData = new double[results.GetLength(0)];
      for (int i = 0; i < kvNoMulticiply.Length; i++)
      {
        xi = kvNoMulticiply[i];
        arrayData.Add(model.GetPatch(0).GetApproximateAt(Result.UY, xi));
        //subData[i] = results[i, 1];
      }

      double centerPointDisplacement = model.GetPatch(0).GetApproximateAt(Result.UY, 0.5);
    }
  }
}