using DEMSoft.Common;
using DEMSoft.Drawing;
using DEMSoft.Drawing.Geometry;
using DEMSoft.EngineeringData;
using DEMSoft.Function;
using DEMSoft.IGA;
using DEMSoft.NURBS;
using DEMSoft.Plot;
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
  internal static class Static_pinned_pinned
  {
    /// <summary>
    /// The main entry point for the application.
    /// </summary>
    [STAThread]
    static void Main()
    {
      //int[] pp = { 2, 3, 4 };
      //int numElem = 6;
      //int[] numElems = new int[numElem];
      //double elem = 1;
      //for (int i = 0; i < numElem; i++)
      //{
      //  elem = Math.Pow(2, i + 1);
      //  numElems[i] = (int)elem;
      //}
      //int[,] arrDof = new int[pp.Length, numElems.Length];
      //double[,] arrError = new double[pp.Length, numElems.Length];
      //double[,] Cdisplacemts = new double[pp.Length, numElems.Length];
      //double Edisplacemts = 0.0;
      //for (int i = 0; i < pp.Length; i++)
      //{
      //  for (int j = 0; j < numElems.Length; j++)
      //  {
      //    double val = Calculate(pp[i], numElems[j], out int dof, out double Cdis, out double Edis);
      //    arrDof[i, j] = dof;
      //    arrError[i, j] = val*100;
      //    Cdisplacemts[i, j] = Math.Round(Cdis,6);
      //    Edisplacemts = Math.Round(Edis,6);
      //  }
      //}

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
          PairOfIsotropicElasticity.YoungModulusAndPoissonRatio, E, 0.2)); //0.01
      steel.AddProperty(new Density(rho));
      steel.TypeMaterialStructure = TypeMaterialStructure.Elasticity;
      //static
      ModelStructureStatic model = new ModelStructureStatic(Dimension.Beam);

      // dinh nghia trang thai bai toan
      model.AddMaterial(steel);
      model.AddPatch(curve);
      model.AttachMaterialToPatch(0, 0);
      model.SetThicknessBeam(h, 0);
      model.SetWidthBeam(b, 0);

      // them dai luong can tinh 
      //model.AddComputeResult(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

      ControlPoint[] selAllCPs = ((AbstractPatch1D)model.GetPatch(0)).SelectAllControlPoints(0);
      ControlPoint[] selCPStart = new ControlPoint[] { selAllCPs[0] };
      ControlPoint[] selCPEnd = new ControlPoint[] { selAllCPs[selAllCPs.Length - 1] };

      ///1D-3Dof
      ConstraintValueArrayOfControlPoints cX0 = new ConstraintValueArrayOfControlPoints(selCPStart, 0, new NullFunctionRToR(), 0);
      ConstraintValueArrayOfControlPoints cY0 = new ConstraintValueArrayOfControlPoints(selCPStart, 1, new NullFunctionRToR(), 0);

      ConstraintValueArrayOfControlPoints cX1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 0, new NullFunctionRToR(), 0);
      ConstraintValueArrayOfControlPoints cY1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 1, new NullFunctionRToR(), 0);

      ////apply constraint
      model.AddConstraint(cX0, cY0, cX1, cY1);
      ///
      AbstractModel.IsParallelProcesing = false;

      //==================================================================
      ////Apply load
      ///
      double[] f = { 0, -347000, 0 };
      double posLoad = 0.5;
      ForcePatch force = new ForcePatch((AbstractPatch1D)model.GetPatch(0), posLoad, f);
      model.AddLoad(force);
      model.IsSaveStepByStep = false;
      model.AddComputeResult(Result.UY);
      model.SetKinematicsFunction(new PolynomialFunctionRToR(0, 1, 0, -4 / (3 * h * h)), 0);
      model.InitializePatch();
      model.PreProcessing();
      model.Solve();
      model.PostProcessing();

      //RESULT

      //===============STATIC===============
      //List<double> arrayData = new List<double>();
      int numPoint = 100;
      double[] arrayData = new double[numPoint];
      double positionLoad = posLoad * L;
      double[,] results = ExactSolution(f[1], L, h, b, E, positionLoad, numPoint);
      double[] subData = new double[results.GetLength(0)];
      double[] xi = new double[subData.Length];
      double dxi = 1.0 / (numPoint - 1);
      double errs = 0.0;
      //ViewerForm viewer = new ViewerForm(true);
      double pointPosMax = Math.Abs(results[(int)numPoint / 2, 1]);
      for (int i = 0; i < xi.Length; i++)
      {
        xi[i] = i * dxi;
        arrayData[i] = model.GetPatch(0).GetApproximateAt(Result.UY, xi[i]);
        subData[i] = results[i, 1];
        if (arrayData[i] != 0)
          errs += Math.Pow((subData[i] - arrayData[i]) / subData[i], 2);
      }
      errs = Math.Sqrt(1.0 / numPoint * errs);
      Console.WriteLine(errs);
      Plotter plotter = new Plotter();
      PlotLine line = new PlotLine();
      line.InputData(xi, arrayData, "IGA");
      line.SetColor(Color.Red);
      line.SetLineWidth(3);
      plotter.AddPlot(line);

      PlotLine line_1 = new PlotLine();
      line_1.InputData(xi, subData, "Exact");
      line_1.SetColor(Color.Blue);
      plotter.AddPlot(line_1);
      plotter.SetShowLegend(true);
      line_1.SetLineType(LineType.DOT1);
      line_1.SetLineWidth(3);
      plotter.SetShowLegend(true);
      plotter.SetXLabel("Normalized distance");
      plotter.SetYLabel("Normalized deflection");
      plotter.Equals(true);
      plotter.Plot();
    }

    public static double Calculate(int p,int num, out int dof, out double Cdis, out double Edis)
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
      curve.pRefinement(p-1);
      curve.hRefinement(num-1); //inssert KVector - incerease p, min luoi

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
          PairOfIsotropicElasticity.YoungModulusAndPoissonRatio, E, 0)); //0.01
      steel.AddProperty(new Density(rho));
      steel.TypeMaterialStructure = TypeMaterialStructure.Elasticity;
      //static
      ModelStructureStatic model = new ModelStructureStatic(Dimension.Beam);

      // dinh nghia trang thai bai toan
      model.AddMaterial(steel);
      model.AddPatch(curve);
      model.AttachMaterialToPatch(0, 0);
      model.SetThicknessBeam(h, 0);
      model.SetWidthBeam(b, 0);

      // them dai luong can tinh 
      //model.AddComputeResult(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

      ControlPoint[] selAllCPs = ((AbstractPatch1D)model.GetPatch(0)).SelectAllControlPoints(0);
      ControlPoint[] selCPStart = new ControlPoint[] { selAllCPs[0] };
      ControlPoint[] selCPEnd = new ControlPoint[] { selAllCPs[selAllCPs.Length - 1] };

      ///1D-3Dof
      ConstraintValueArrayOfControlPoints cX0 = new ConstraintValueArrayOfControlPoints(selCPStart, 0, new NullFunctionRToR(), 0);
      ConstraintValueArrayOfControlPoints cY0 = new ConstraintValueArrayOfControlPoints(selCPStart, 1, new NullFunctionRToR(), 0);

      ConstraintValueArrayOfControlPoints cX1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 0, new NullFunctionRToR(), 0);
      ConstraintValueArrayOfControlPoints cY1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 1, new NullFunctionRToR(), 0);

      ////apply constraint
      model.AddConstraint(cX0, cY0, cX1, cY1);
      ///
      AbstractModel.IsParallelProcesing = false;

      //==================================================================
      ////Apply load
      ///
      double[] f = { 0, -347000, 0 };
      double posLoad = 0.5;
      ForcePatch force = new ForcePatch((AbstractPatch1D)model.GetPatch(0), posLoad, f);
      model.AddLoad(force);
      model.IsSaveStepByStep = false;
      model.AddComputeResult(Result.UY);
      model.SetKinematicsFunction(new PolynomialFunctionRToR(0, 1, 0, -4 / (3 * h * h)), 0);
      model.InitializePatch();
      model.PreProcessing();
      model.Solve();
      model.PostProcessing();

      //RESULT

      //===============STATIC===============
      //List<double> arrayData = new List<double>();
      int numPoint = 100;
      double[] arrayData = new double[numPoint];
      double positionLoad = posLoad * L;
      double[,] results = ExactSolution(f[1], L, h, b, E, positionLoad, numPoint);
      double[] subData = new double[results.GetLength(0)];
      double[] xi = new double[subData.Length];
      double dxi = 1.0 / (numPoint - 1);
      double err = 0.0;

      double pointPosMax = results[(int)numPoint / 2, 1];
      for (int i = 0; i < xi.Length; i++)
      {
        xi[i] = i * dxi;
        arrayData[i] = model.GetPatch(0).GetApproximateAt(Result.UY, xi[i]);
        subData[i] = results[i, 1];
        if (arrayData[i] != 0)
          err += Math.Pow((subData[i] - arrayData[i]), 2);
      }
      err = Math.Sqrt(1.0 / numPoint * err);
      dof = model.CountDOF();
      Cdis = model.GetPatch(0).GetApproximateAt(Result.UY, posLoad);
      Edis = subData[numPoint/2];
      return err;
    }
    private static double[,] ExactSolution(double F, double L, double h, double width, double E, double positionLoad, int num)
    {
      double a = positionLoad;
      double b = L - positionLoad;
      double I = width * Math.Pow(h, 3) / 12.0;
      double dx = L / (num - 1);
      double[,] results = new double[num, 3];
      for (int i = 0; i < num; i++)
      {
        results[i, 0] = i * dx;
        double x = results[i, 0];
        if (x < positionLoad)
        {
          results[i, 1] = F * b * x * (L * L - x * x - b * b) / (6.0 * L * E * I);
        }
        else
        {
          results[i, 1] = F * b * (L / b * Math.Pow((x - a), 3) + (L * L - b * b) * x - Math.Pow(x, 3)) / (6.0 * L * E * I);
        }
      }
      return results;
    }
  }
}