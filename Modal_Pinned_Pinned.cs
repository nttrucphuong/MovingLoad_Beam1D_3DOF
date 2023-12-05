//using DEMSoft.Common;
//using DEMSoft.Drawing;
//using DEMSoft.EngineeringData;
//using DEMSoft.Function;
//using DEMSoft.IGA;
//using DEMSoft.NURBS;
//using DEMSoft.Plot;
//using System;
//using System.Collections.Generic;
//using System.Data;
//using System.Diagnostics;
//using System.Drawing;
//using System.IO;
//using System.Reflection;
//using System.Text.RegularExpressions;
//using System.Windows.Forms;

//namespace SampleTesting
//{
//  internal static class Modala_Pinned_Pinned
//  {
//    /// <summary>
//    /// The main entry point for the application.
//    /// </summary>
//    [STAThread]
//    static void Main()
//    {
//      int[] pp = { 2, 3 };
//      int numElem = 20;
//      int[] numElems = new int[numElem];
//      int elem = 5;
//      for (int i = 0; i < numElem; i++)
//      {
//        numElems[i] = elem;
//        elem += 5;
//      }
//      int[,] arrDof = new int[pp.Length, numElems.Length];
//      double[,] arrError = new double[pp.Length, numElems.Length];
//      for (int i = 0; i < pp.Length; i++)
//      {
//        for (int j = 0; j < numElems.Length; j++)
//        {
//          double val = CalculateErr(pp[i], numElems[j], out int dof);
//          arrDof[i, j] = dof;
//          arrError[i, j] = val;
//        }
//      }

//      //var targetDir = Path.Combine(@"F:\ppp\CTDT\Nam_tu\LVTN\github\MovingLoad\MovingLoad\bin\Debug\temp");
//      //if (Directory.Exists(targetDir))
//      //{
//      //  Directory.Delete(targetDir, true);
//      //}

//      ////DOF in each patch: ux,uy, theta

//      ////double L = 8; //m
//      ////double h = 1;//4.1602;
//      ////double b = 1;
//      ////double E = 100000;

//      //double L = 34; //m
//      //double h = 0.2;//4.1602;
//      //double b = 0.5; //m
//      //double E = 2.976E+14; //Pa
//      //double rho = 114000; //kg/m^3
//      //double I = b * Math.Pow(h, 3) / 12.0;
//      //double area = b * h;

//      ////double L = 1; //m
//      ////double h = 0.01;//4.1602;
//      ////double b = 0.01; //m
//      ////double E = 206E9; //Pa
//      ////double rho = 7860; //kg/m^3
//      ////double I = b * Math.Pow(h, 3) / 12.0;
//      ////double area = b * h;

//      ////ViewerForm viewer = new ViewerForm(true);
//      //NURBSCurve curve = GeometryCreator.CreateStraightNURBSCurve(0, 0, 0, L, 0, 0);
//      //// danh sach cac Curve
//      //List<NURBSCurve> listCurve = new List<NURBSCurve>();
//      //curve.pRefinement(3);
//      //curve.hRefinement(10); //inssert KVector - incerease p, min luoi

//      ////curve.colorKnot = Color.Red;
//      ////curve.colorCurve = Color.Blue;
//      ////curve.colorControlNet = Color.Orange;
//      ////curve.isDrawOriginal = false;
//      ////curve.colorControlPoint = Color.Green;
//      ////curve.Draw(viewer);

//      ////viewer.UpdateCamera();
//      ////viewer.Run();


//      //////// --------------------Khai bao vat lieu--------------------

//      //Material steel = new Material("steel");
//      //steel.AddProperty(new IsotropicElasticity(
//      //    PairOfIsotropicElasticity.YoungModulusAndPoissonRatio, E, 0));
//      //steel.AddProperty(new Density(rho));
//      //steel.TypeMaterialStructure = TypeMaterialStructure.Elasticity;
//      ////static
//      //ModelStructureModal model = new ModelStructureModal(Dimension.Beam);

//      //// dinh nghia trang thai bai toan
//      //model.AddMaterial(steel);
//      //model.AddPatch(curve);
//      //model.AttachMaterialToPatch(0, 0);
//      //model.SetThicknessBeam(h, 0);
//      //model.SetWidthBeam(b, 0);

//      //// them dai luong can tinh 
//      ////model.AddComputeResult(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

//      //ControlPoint[] selAllCPs = ((AbstractPatch1D)model.GetPatch(0)).SelectAllControlPoints(0);
//      //ControlPoint[] selCPStart = new ControlPoint[] { selAllCPs[0] };
//      //ControlPoint[] selCPEnd = new ControlPoint[] { selAllCPs[selAllCPs.Length - 1] };

//      /////1D-3Dof
//      //ConstraintValueArrayOfControlPoints cX0 = new ConstraintValueArrayOfControlPoints(selCPStart, 0, new NullFunctionRToR(), 0);
//      //ConstraintValueArrayOfControlPoints cY0 = new ConstraintValueArrayOfControlPoints(selCPStart, 1, new NullFunctionRToR(), 0);

//      //ConstraintValueArrayOfControlPoints cX1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 0, new NullFunctionRToR(), 0);
//      //ConstraintValueArrayOfControlPoints cY1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 1, new NullFunctionRToR(), 0);

//      ////apply constraint
//      //model.AddConstraint(cX0, cY0, cX1, cY1);
//      /////

//      /////1D-1Dof///
//      ////ConstraintValueArrayOfControlPoints cY = new ConstraintValueArrayOfControlPoints(selCP, 0, new NullFunctionRToR(), 0);
//      ////model.AddConstraint(cY);
//      ////model.AddConstraint(cX);
//      //AbstractModel.IsParallelProcesing = false;
//      ////model.SetKinematicsFunction(new LinearFunctionRToR(0, 1, 0, 1), 0);
//      //model.SetKinematicsFunction(new PolynomialFunctionRToR(0, 1, 0, -4 / (3 * h * h)), 0);
//      //int numMode = 4;
//      //model.SetNumberOfMode(numMode);

//      //model.InitializePatch();
//      //model.PreProcessing();
//      //model.Solve();
//      //model.PostProcessing();

//      ////ViewerForm viewer = new ViewerForm();
//      ////model.DrawModeResult(0, viewer, 5, 1);
//      ////viewer.UpdateCamera();
//      ////viewer.Run();

//      ////viewer = new ViewerForm();
//      ////model.DrawModeResult(1, viewer, 5, 1);
//      ////viewer.UpdateCamera();
//      ////viewer.Run();

//      ////viewer = new ViewerForm();
//      ////model.DrawModeResult(2, viewer, 5, 1);
//      ////viewer.UpdateCamera();
//      ////viewer.Run();

//      ////viewer = new ViewerForm();
//      ////model.DrawModeResult(3, viewer, 5, 1);
//      ////viewer.UpdateCamera();
//      ////viewer.Run();

//      ////viewer = new ViewerForm();
//      ////model.DrawModeResult(4, viewer, 5, 1);
//      ////viewer.UpdateCamera();
//      ////viewer.Run();


//      //double[] eigenValue = model.GetEigenValue();
//      //double error = 0;

//      //double[] exactFrequencies = new double[numMode];
//      //double pi = Math.PI; ;
//      //double[] betaNL_pinned_pinned = {pi, 2*pi, 3*pi, 4*pi };

//      //for (int i = 1; i <= numMode; i++)
//      //{
//      //  exactFrequencies[i - 1] = Math.Pow(betaNL_pinned_pinned[i - 1], 2) * Math.Sqrt(E * I / (rho * b * h * Math.Pow(L, 4)));
//      //}
//      //for (int i = 0; i < eigenValue.Length; i++)
//      //{
//      //  error += Math.Pow((exactFrequencies[i] - eigenValue[i]) / exactFrequencies[i], 2);
//      //}
//      //error = 1.0 / eigenValue.Length * Math.Sqrt(error);
//      //Console.WriteLine(error);
//    }

//    public static double CalculateErr(int p, int num, out int dof)
//    {
//      var targetDir = Path.Combine(@"F:\ppp\CTDT\Nam_tu\LVTN\github\MovingLoad\MovingLoad\bin\Debug\temp");
//      if (Directory.Exists(targetDir))
//      {
//        Directory.Delete(targetDir, true);
//      }

//      double L = 34; //m
//      double h = 0.2;//4.1602;
//      double b = 0.5; //m
//      double E = 2.976E+14; //Pa
//      double rho = 114000; //kg/m^3
//      double I = b * Math.Pow(h, 3) / 12.0;
//      NURBSCurve curve = GeometryCreator.CreateStraightNURBSCurve(0, 0, 0, L, 0, 0);
//      // danh sach cac Curve
//      List<NURBSCurve> listCurve = new List<NURBSCurve>();
//      curve.pRefinement(p);
//      curve.hRefinement(num); //inssert KVector - incerease p, min luoi

//      ////// --------------------Khai bao vat lieu--------------------

//      Material steel = new Material("steel");
//      steel.AddProperty(new IsotropicElasticity(
//          PairOfIsotropicElasticity.YoungModulusAndPoissonRatio, E, 0));
//      steel.AddProperty(new Density(rho));
//      steel.TypeMaterialStructure = TypeMaterialStructure.Elasticity;
//      //static
//      ModelStructureModal model = new ModelStructureModal(Dimension.Beam);

//      // dinh nghia trang thai bai toan
//      model.AddMaterial(steel);
//      model.AddPatch(curve);
//      model.AttachMaterialToPatch(0, 0);
//      model.SetThicknessBeam(h, 0);
//      model.SetWidthBeam(b, 0);

//      ControlPoint[] selAllCPs = ((AbstractPatch1D)model.GetPatch(0)).SelectAllControlPoints(0);
//      ControlPoint[] selCPStart = new ControlPoint[] { selAllCPs[0] };
//      ControlPoint[] selCPEnd = new ControlPoint[] { selAllCPs[selAllCPs.Length - 1] };

//      ///1D-3Dof
//      ConstraintValueArrayOfControlPoints cX0 = new ConstraintValueArrayOfControlPoints(selCPStart, 0, new NullFunctionRToR(), 0);
//      ConstraintValueArrayOfControlPoints cY0 = new ConstraintValueArrayOfControlPoints(selCPStart, 1, new NullFunctionRToR(), 0);

//      ConstraintValueArrayOfControlPoints cX1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 0, new NullFunctionRToR(), 0);
//      ConstraintValueArrayOfControlPoints cY1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 1, new NullFunctionRToR(), 0);

//      //apply constraint
//      model.AddConstraint(cX0, cY0, cX1, cY1);
//      ///
//      AbstractModel.IsParallelProcesing = false;
//      //model.SetKinematicsFunction(new LinearFunctionRToR(0, 1, 0, 1), 0);
//      model.SetKinematicsFunction(new PolynomialFunctionRToR(0, 1, 0, -4 / (3 * h * h)), 0);
//      int numMode = 4;
//      model.SetNumberOfMode(numMode);

//      model.InitializePatch();
//      model.PreProcessing();
//      model.Solve();
//      model.PostProcessing();

//      dof = model.CountDOF();

//      double[] eigenValue = model.GetEigenValue();
//      double error = 0;

//      double[] exactFrequencies = new double[numMode];
//      double pi = Math.PI; ;
//      double[] betaNL_pinned_pinned = { pi, 2 * pi, 3 * pi, 4 * pi };

//      for (int i = 1; i <= numMode; i++)
//      {
//        exactFrequencies[i - 1] = Math.Pow(betaNL_pinned_pinned[i - 1], 2) * Math.Sqrt(E * I / (rho * b * h * Math.Pow(L, 4)));
//      }
//      for (int i = 0; i < eigenValue.Length; i++)
//      {
//        error += Math.Pow((exactFrequencies[i] - eigenValue[i]) / exactFrequencies[i], 2);
//      }
//      error = 1.0 / eigenValue.Length * Math.Sqrt(error);
//      return error;
//    }
//  }
//}