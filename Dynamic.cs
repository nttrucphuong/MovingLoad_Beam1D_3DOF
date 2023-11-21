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
  internal static class Dynamic
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

      //dynamic
      ModelStructureTransient model = new ModelStructureTransient(Dimension.Beam);

      // dinh nghia trang thai bai toan
      model.AddMaterial(steel);
      model.AddPatch(curve);
      model.AttachMaterialToPatch(0, 0);
      model.SetThicknessBeam(h, 0);
      model.SetTCrossSectionAreaBeam(area, 0);
      model.SetMomentOfInertiaBeam(I, 0);

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
      ///

      AbstractModel.IsParallelProcesing = false;

      ///////Dynamic
      //FunctionRToR sinForce = new SinFunctionRToR(1, 500);
      //FunctionRToR concertrateForce = new ConstantFunctionRToR(10);
      FunctionRToR concertrateForce = new ConstantFunctionRToR(347000);

      //ForceTime f = new ForceTime(curve.ControlPoints[curve.ControlPoints.Length - 1], concertrateForce, 0, -1, 0);
      //ForceTime f = new ForceTime(curve.ControlPoints[curve.ControlPoints.Length - 1], sinForce, 0, 1, 0);
      //model.AddLoad(f);

      //double V0 = 0.5;
      double V0 = 136.284/*68.1419*//*13.6*/; //m/s
      double V0_ = V0; //V0*6.97
      FunctionRToR positionLoad = new LinearFunctionRToR(0, 1, 0, V0_);
      ForceMovingOnPatch fMoving = new ForceMovingOnPatch((AbstractPatch1D)model.GetPatch(0), positionLoad, concertrateForce, 0, -1, 0.0);
      model.AddLoad(fMoving);

      //model.IsSaveStepByStep = true;

      double totalTime = 1 / V0_;//L / V0
      int numStep = 100;

      model.SetKinematicsFunction(new LinearFunctionRToR(0, 1, 0, 1), 0);//  new LinearFunctionRToR(0, 1, 0, 1)  new NullFunctionRToR()
      double[] deltaT = { totalTime / numStep };
      int[] numberOfSubStepLoad = { numStep };
      int[] numberOfStepSave = { numStep };
      model.SetTransientSolver(deltaT, numberOfSubStepLoad, numberOfStepSave);
      //=============================================================
      model.InitializePatch();
      model.PreProcessing();
      model.Solve();
      model.PostProcessing();

      //RESULT
      List<double> arrayData = new List<double>();
      for (int i = 1; i <= numStep; i++)
      {
        model.ReadResultByLoadStep(i);
        arrayData.Add(model.GetPatch(0).GetApproximateAt(Result.UY, 0.5));
      }
    }
    private static double[,] ExactSolution(double F, double L, double h, double b, double E, int num)
    {
      double I = b * Math.Pow(h, 3) / 12.0;
      double dx = L / (num - 1);
      double[,] results = new double[num, 3];
      for (int i = 0; i < num; i++)
      {
        results[i, 0] = i * dx;
        double x = L - results[i, 0];
        results[i, 1] = F / (6.0 * E * I) * (x * x * x - 3 * L * L * x + 2 * L * L * L);
      }
      return results;
    }
  }
}