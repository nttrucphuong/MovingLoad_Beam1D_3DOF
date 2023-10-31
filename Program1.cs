//using DEMSoft.Common;
//using DEMSoft.Drawing;
//using DEMSoft.EngineeringData;
//using DEMSoft.Function;
//using DEMSoft.IGA;
//using DEMSoft.NURBS;
//using System;
//using System.Collections.Generic;
//using System.Diagnostics;
//using System.Drawing;

//namespace SampleTesting
//{
//  internal static class Program1
//  {
//    /// <summary>
//    /// The main entry point for the application.
//    /// </summary>
//    [STAThread]
//    static void Main()
//    {
//      double L = 8;
//      double h = 1.0;//4.1602;
//      double b = 1.0;
//      double E = 100000;
//      ViewerForm viewer = new ViewerForm(true);
//      NURBSCurve curve = GeometryCreator.CreateStraightNURBSCurve(0, 0, 0, L, 0, 0);
//      // danh sach cac Curve
//      List<NURBSCurve> listCurve = new List<NURBSCurve>();
//      //curve.SetDegreeOnAllDirections(2);
//      curve.pRefinement(1);
//      curve.hRefinement(4); //inssert KVector - incerease p, min luoi

//      curve.colorKnot = Color.Red;
//      curve.colorCurve = Color.Blue;
//      curve.colorControlNet = Color.Orange;
//      curve.isDrawOriginal = false;
//      curve.colorControlPoint = Color.Green;
//      curve.Draw(viewer);

//      viewer.UpdateCamera();
//      viewer.Run();


//      ////// --------------------Khai bao vat lieu--------------------

//      Material steel = new Material("steel");
//      steel.AddProperty(new IsotropicElasticity(
//          PairOfIsotropicElasticity.YoungModulusAndPoissonRatio, E, 0.3));
//      steel.TypeMaterialStructure = TypeMaterialStructure.Elasticity;


//      ModelStructureStatic model = new ModelStructureStatic(Dimension.Beam);
//      // dinh nghia trang thai bai toan

//      model.AddMaterial(steel);
//      model.AddPatch(curve);
//      model.AttachMaterialToPatch(0, 0);
//      model.SetThicknessBeam(h, 0);

//      // them dai luong can tinh 
//      //model.AddComputeResult(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

//      ControlPoint[] selAllCPs = ((AbstractPatch1D)model.GetPatch(0)).SelectAllControlPoints(0);
//      ControlPoint[] selCP = new ControlPoint[] { selAllCPs[0], selAllCPs[1] };
//      ConstraintValueArrayOfControlPoints cY = new ConstraintValueArrayOfControlPoints(selCP, 0, new NullFunctionRToR(), 0);
//      ////apply constraint
//      model.AddConstraint(cY);

//      double[] fy = { -10 };
//      Force f = new Force(curve.ControlPoints[curve.ControlPoints.Length - 1], fy);
//      model.AddLoad(f);
//      //model.AddLoad();
//      model.IsSaveStepByStep = false;
//      AbstractModel.IsParallelProcesing = false;
//      model.InitializePatch();
//      model.SetKinematicsFunction(new LinearFunctionRToR(0, 1, 0, 1), 0);//  new LinearFunctionRToR(0, 1, 0, 1)  new NullFunctionRToR()
//      model.PreProcessing();
//      model.Solve();
//      model.PostProcessing();
//      Console.WriteLine(model.ToString());
//      Console.ReadLine();


//      //sw.Stop();
//      //double t = sw.Elapsed.TotalMilliseconds;
//      //model.Extrapolation(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

//      viewer = new ViewerForm(true);
//      //model.DrawResult(Result.UY, viewer, 20, false, 20);
//      model.DrawGeometry(viewer, 30, 1);

//      double[,] results = ExactSolution(fy[0], L, h, b, E, 20);
//      PointSet ps = new PointSet(results);
//      ps.SetColor(Color.Purple);
//      ps.SetPointSize(7);
//      viewer.AddObject3D(ps);

//      viewer.UpdateCamera();
//      viewer.Run();
//    }
//    private static double[,] ExactSolution(double F, double L, double h, double b, double E, int num)
//    {
//      double I = b * Math.Pow(h, 3) / 12.0;
//      double dx = L / (num - 1);
//      double[,] results = new double[num, 3];
//      for (int i = 0; i < num; i++)
//      {
//        results[i, 0] = i * dx;
//        double x = L - results[i, 0];
//        //results[i, 1] = F / (E * I) * (3 * L * results[i, 0] - results[i, 0] * results[i, 0] * results[i, 0]);
//        results[i, 1] = F / (6.0 * E * I) * (x * x * x - 3 * L * L * x + 2 * L * L * L);
//        //results[i, 1] = 30 * E * I * results[i, 1] / (F * L * L * L);
//      }
//      return results;
//    }
//  }
//}