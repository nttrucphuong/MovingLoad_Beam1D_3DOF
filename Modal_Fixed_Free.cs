////using DEMSoft.Common;
////using DEMSoft.Drawing;
////using DEMSoft.EngineeringData;
////using DEMSoft.Function;
////using DEMSoft.IGA;
////using DEMSoft.NURBS;
////using DEMSoft.Plot;
////using System;
////using System.Collections.Generic;
////using System.Data;
////using System.Diagnostics;
////using System.Drawing;
////using System.IO;
////using System.Reflection;
////using System.Text.RegularExpressions;
////using System.Windows.Forms;

////namespace SampleTesting
////{
////  internal static class Modal
////  {
////    /// <summary>
////    /// The main entry point for the application.
////    /// </summary>
////    [STAThread]
////    static void Main()
////    {
////      var targetDir = Path.Combine(@"F:\ppp\CTDT\Nam_tu\LVTN\github\MovingLoad\MovingLoad\bin\Debug\temp");
////      if (Directory.Exists(targetDir))
////      {
////        Directory.Delete(targetDir, true);
////      }

////      //DOF in each patch: ux,uy, theta

////      //double L = 8; //m
////      //double h = 1;//4.1602;
////      //double b = 1;
////      //double E = 100000;

////      double L = 34; //m
////      double h = 0.2;//4.1602;
////      double b = 0.5; //m
////      double E = 2.976E+14; //Pa
////      double rho = 114000; //kg/m^3
////      double I = b * Math.Pow(h, 3) / 12.0;
////      double area = b * h;
////      //ViewerForm viewer = new ViewerForm(true);
////      NURBSCurve curve = GeometryCreator.CreateStraightNURBSCurve(0, 0, 0, L, 0, 0);
////      // danh sach cac Curve
////      List<NURBSCurve> listCurve = new List<NURBSCurve>();
////      curve.pRefinement(4);
////      curve.hRefinement(100); //inssert KVector - incerease p, min luoi

////      //curve.colorKnot = Color.Red;
////      //curve.colorCurve = Color.Blue;
////      //curve.colorControlNet = Color.Orange;
////      //curve.isDrawOriginal = false;
////      //curve.colorControlPoint = Color.Green;
////      //curve.Draw(viewer);

////      //viewer.UpdateCamera();
////      //viewer.Run();


////      ////// --------------------Khai bao vat lieu--------------------

////      Material steel = new Material("steel");
////      steel.AddProperty(new IsotropicElasticity(
////          PairOfIsotropicElasticity.YoungModulusAndPoissonRatio, E, 0.2));
////      steel.AddProperty(new Density(rho));
////      steel.TypeMaterialStructure = TypeMaterialStructure.Elasticity;
////      //static
////      ModelStructureModal model = new ModelStructureModal(Dimension.Beam);

////      // dinh nghia trang thai bai toan
////      model.AddMaterial(steel);
////      model.AddPatch(curve);
////      model.AttachMaterialToPatch(0, 0);
////      model.SetThicknessBeam(h, 0);
////      model.SetWidthBeam(b, 0);

////      // them dai luong can tinh 
////      //model.AddComputeResult(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

////      ControlPoint[] selAllCPs = ((AbstractPatch1D)model.GetPatch(0)).SelectAllControlPoints(0);
////      ControlPoint[] selCPStart = new ControlPoint[] { selAllCPs[0] };
////      ControlPoint[] selCPEnd = new ControlPoint[] { selAllCPs[selAllCPs.Length - 1] };

////      ///1D-3Dof
////      ConstraintValueArrayOfControlPoints cX0 = new ConstraintValueArrayOfControlPoints(selCPStart, 0, new NullFunctionRToR(), 0);
////      ConstraintValueArrayOfControlPoints cY0 = new ConstraintValueArrayOfControlPoints(selCPStart, 1, new NullFunctionRToR(), 0);
////      //ConstraintValueArrayOfControlPoints cThetaZ0 = new ConstraintValueArrayOfControlPoints(selCPStart, 2, new NullFunctionRToR(), 0);

////      ConstraintValueArrayOfControlPoints cX1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 0, new NullFunctionRToR(), 0);
////      ConstraintValueArrayOfControlPoints cY1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 1, new NullFunctionRToR(), 0);
////      //ConstraintValueArrayOfControlPoints cThetaZ1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 2, new NullFunctionRToR(), 0);

////      ////apply constraint
////      //model.AddConstraint(cX0, cY0, cThetaZ0);
////      //model.AddConstraint(/*cX1,*/ cY1);

////      model.AddConstraint(cX0, cY0);
////      model.AddConstraint(cX1, cY1);
////      ///

////      ///1D-1Dof///
////      //ConstraintValueArrayOfControlPoints cY = new ConstraintValueArrayOfControlPoints(selCP, 0, new NullFunctionRToR(), 0);
////      //model.AddConstraint(cY);
////      //model.AddConstraint(cX);
////      AbstractModel.IsParallelProcesing = false;
////      model.SetKinematicsFunction(new LinearFunctionRToR(0, 1, 0, 1), 0);
////      int numMode = 5;
////      model.SetNumberOfMode(numMode);

////      model.InitializePatch();
////      model.PreProcessing();
////      model.Solve();
////      model.PostProcessing();

////      double[] eigenValue = model.GetEigenValue();
////      double error = 0;

////      double[] exactFrequencies = new double[numMode];
////      for (int i = 1; i <= numMode; i++)
////      {
////        exactFrequencies[i - 1] = Math.Pow(i * Math.PI, 2) * Math.Sqrt(E * I / (rho * b * h * Math.Pow(L, 4)));
////      }
////      for (int i = 0; i < eigenValue.Length; i++)
////      {
////        error += Math.Pow((exactFrequencies[i] - eigenValue[i]) / eigenValue[i], 2);
////      }
////      error = 1.0 / eigenValue.Length * Math.Sqrt(error);
////      Console.WriteLine(error);
////    }
////  }
////}

/////////////////////////////////////////////////////////////
////using DEMSoft.Drawing;
////using DEMSoft.EngineeringData;
////using DEMSoft.Function;
////using DEMSoft.IGA;
////using DEMSoft.NURBS;
////using System.Collections.Generic;
////using System.IO;
////using System;
///////
////DEMSoft.Common;
////using DEMSoft.Drawing;
////using DEMSoft.EngineeringData;
////using DEMSoft.Function;
////using DEMSoft.IGA;
////using DEMSoft.NURBS;
////using DEMSoft.Plot;
////using System;
////using System.Collections.Generic;
////using System.Data;
////using System.Diagnostics;
////using System.Drawing;
////using System.IO;
////using System.Reflection;
////using System.Text.RegularExpressions;
////using System.Windows.Forms;

////namespace SampleTesting
////{
////  internal static class Modal
////  {
////    /// <summary>
////    /// The main entry point for the application.
////    /// </summary>
////    [STAThread]
////    static void Main()
////    {
////      var targetDir = Path.Combine(@"F:\ppp\CTDT\Nam_tu\LVTN\github\MovingLoad\MovingLoad\bin\Debug\temp");
////      if (Directory.Exists(targetDir))
////      {
////        Directory.Delete(targetDir, true);
////      }

////      //DOF in each patch: ux,uy, theta

////      //double L = 8; //m
////      //double h = 1;//4.1602;
////      //double b = 1;
////      //double E = 100000;

////      double L = 34; //m
////      double h = 0.2;//4.1602;
////      double b = 0.5; //m
////      double E = 2.976E+14; //Pa
////      double rho = 114000; //kg/m^3
////      double I = b * Math.Pow(h, 3) / 12.0;
////      double area = b * h;

////      //double L = 1; //m
////      //double h = 0.01;//4.1602;
////      //double b = 0.01; //m
////      //double E = 206E9; //Pa
////      //double rho = 7860; //kg/m^3
////      //double I = b * Math.Pow(h, 3) / 12.0;
////      //double area = b * h;

////      //ViewerForm viewer = new ViewerForm(true);
////      NURBSCurve curve = GeometryCreator.CreateStraightNURBSCurve(0, 0, 0, L, 0, 0);
////      // danh sach cac Curve
////      List<NURBSCurve> listCurve = new List<NURBSCurve>();
////      curve.pRefinement(2);
////      curve.hRefinement(2); //inssert KVector - incerease p, min luoi

////      //curve.colorKnot = Color.Red;
////      //curve.colorCurve = Color.Blue;
////      //curve.colorControlNet = Color.Orange;
////      //curve.isDrawOriginal = false;
////      //curve.colorControlPoint = Color.Green;
////      //curve.Draw(viewer);

////      //viewer.UpdateCamera();
////      //viewer.Run();


////      ////// --------------------Khai bao vat lieu--------------------

////      Material steel = new Material("steel");
////      steel.AddProperty(new IsotropicElasticity(
////          PairOfIsotropicElasticity.YoungModulusAndPoissonRatio, E, 0.3));
////      steel.AddProperty(new Density(rho));
////      steel.TypeMaterialStructure = TypeMaterialStructure.Elasticity;
////      //static
////      ModelStructureModal model = new ModelStructureModal(Dimension.Beam);

////      // dinh nghia trang thai bai toan
////      model.AddMaterial(steel);
////      model.AddPatch(curve);
////      model.AttachMaterialToPatch(0, 0);
////      model.SetThicknessBeam(h, 0);
////      model.SetWidthBeam(b, 0);

////      // them dai luong can tinh 
////      //model.AddComputeResult(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

////      ControlPoint[] selAllCPs = ((AbstractPatch1D)model.GetPatch(0)).SelectAllControlPoints(0);
////      ControlPoint[] selCPStart = new ControlPoint[] { selAllCPs[0] };
////      ControlPoint[] selCPEnd = new ControlPoint[] { selAllCPs[selAllCPs.Length - 1] };

////      ///1D-3Dof
////      ConstraintValueArrayOfControlPoints cX0 = new ConstraintValueArrayOfControlPoints(selCPStart, 0, new NullFunctionRToR(), 0);
////      ConstraintValueArrayOfControlPoints cY0 = new ConstraintValueArrayOfControlPoints(selCPStart, 1, new NullFunctionRToR(), 0);
////      ConstraintValueArrayOfControlPoints cThetaZ0 = new ConstraintValueArrayOfControlPoints(selCPStart, 2, new NullFunctionRToR(), 0);

////      ConstraintValueArrayOfControlPoints cX1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 0, new NullFunctionRToR(), 0);
////      ConstraintValueArrayOfControlPoints cY1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 1, new NullFunctionRToR(), 0);
////      ConstraintValueArrayOfControlPoints cThetaZ1 = new ConstraintValueArrayOfControlPoints(selCPEnd, 2, new NullFunctionRToR(), 0);

////      //apply constraint
////      //model.AddConstraint(cX0, cY0, cThetaZ0);
////      //model.AddConstraint(cX1, cY1, cThetaZ1);

////      model.AddConstraint(cX0, cY0);
////      model.AddConstraint(cY1);
////      ///

////      ///1D-1Dof///
////      //ConstraintValueArrayOfControlPoints cY = new ConstraintValueArrayOfControlPoints(selCP, 0, new NullFunctionRToR(), 0);
////      //model.AddConstraint(cY);
////      //model.AddConstraint(cX);
////      AbstractModel.IsParallelProcesing = false;
////      model.SetKinematicsFunction(new LinearFunctionRToR(0, 1, 0, 1), 0);
////      int numMode = 4;
////      model.SetNumberOfMode(numMode);

////      model.InitializePatch();
////      model.PreProcessing();
////      model.Solve();
////      model.PostProcessing();
////      ViewerForm viewer = new ViewerForm();

////      model.DrawModeResult(0, viewer, 5, 1);
////      viewer.UpdateCamera();
////      viewer.Run();

////      viewer = new ViewerForm();

////      model.DrawModeResult(1, viewer, 5, 1);
////      viewer.UpdateCamera();
////      viewer.Run();

////      viewer = new ViewerForm();

////      model.DrawModeResult(2, viewer, 5, 1);
////      viewer.UpdateCamera();
////      viewer.Run();

////      double[] eigenValue = model.GetEigenValue();
////      double error = 0;

////      double[] exactFrequencies = new double[numMode];
////      double[] betaNL_fixed_pinned = { 3.926602, 7.068583, 10.210176, 13.351768 };

////      for (int i = 1; i <= numMode; i++)
////      {
////        exactFrequencies[i - 1] = Math.Pow(betaNL_fixed_pinned[i - 1], 2) * Math.Sqrt(E * I / (rho * b * h * Math.Pow(L, 4)));
////      }
////      for (int i = 0; i < eigenValue.Length; i++)
////      {
////        error += Math.Pow((exactFrequencies[i] - eigenValue[i]) / exactFrequencies[i], 2);
////      }
////      error = 1.0 / eigenValue.Length * Math.Sqrt(error);
////      Console.WriteLine(error);
////    }
////  }
////}


//////////////////////////////////////////////////////////////////
//using DEMSoft.NURBS;
//using System;
//using System.IO;///pinned-pinned
//using DEMSoft.Common;
//using DEMSoft.Drawing;
//using DEMSoft.EngineeringData;
//using DEMSoft.Function;
//using DEMSoft.IGA;
//using System.Collections.Generic;
//using System.Data;
//using System.Diagnostics;
//using System.Drawing;
//using System.Reflection;
//using System.Text.RegularExpressions;
//using System.Windows.Forms;

//namespace SampleTesting
//{
//  internal static class Modala_Fixed_Free
//  {
//    /// <summary>
//    /// The main entry point for the application.
//    /// </summary>
//    //[STAThread]
//    static void Main()
//    {

//      int[] pp = { 2, 3, 4 };
//      int numElem = 5;
//      int[] numElems = new int[numElem];
//      int elem = 10;
//      int numMode = 5;
//      for (int i = 0; i < numElem; i++)
//      {
//        numElems[i] = elem;
//        elem += 5;
//      }
//      int[,] arrDof = new int[pp.Length, numElems.Length];
//      double[,] arrError = new double[pp.Length, numElems.Length];
//      double[,] arrAnsysErr = new double[numElems.Length, numMode];
//      List<double> listAnsysErr = new List<double>();
//      double[,][] naturalFrequencies = new double[pp.Length, numElems.Length][];
//      double[] exactFrequencies = new double[numMode];
//      double[][] Ansys = new double[numElem][];
//      double AnsysErr = 0;
//      Ansys[0] = new double[4] { 1.4279, 8.9475, 25.051, 49.101 };
//      Ansys[1] = new double[4] { 1.4279, 8.9473, 25.048, 49.074 };
//      Ansys[2] = new double[4] { 1.4279, 8.9473, 25.047, 49.068 };
//      Ansys[3] = new double[4] { 1.4279, 8.9473, 25.047, 49.065 };
//      Ansys[4] = new double[4] { 1.4279, 8.9473, 25.047, 49.064 };

//      for (int i = 0; i < pp.Length; i++)
//      {
//        for (int j = 0; j < numElems.Length; j++)
//        {
//          double err = CalculateErr(pp[i], numElems[j], out int dof, out double[] eigenValue, out exactFrequencies);
//          arrDof[i, j] = dof;
//          arrError[i, j] = err;
//          naturalFrequencies[i, j] = eigenValue;
//        }
//      }

//      //for (int i = 0; i < numElem; i++)
//      //{
//      //  for (int j = 0; j < numMode; j++)
//      //  {
//      //    AnsysErr += Math.Pow((exactFrequencies[j] - Ansys[i][j] * 2 * Math.PI), 2);
//      //    //arrAnsysErr[i, j] = 1.0 / numMode * AnsysErr; // Mean Square Err
//      //  }
//      //  listAnsysErr.Add(1.0 / numMode * AnsysErr);
//      //}
//      // Mean Square Err

//      ////var targetDir = Path.Combine(@"F:\ppp\CTDT\Nam_tu\LVTN\github\MovingLoad\MovingLoad\bin\Debug\temp");
//      ////if (Directory.Exists(targetDir))
//      ////{
//      ////  Directory.Delete(targetDir, true);
//      ////}

//      ////DOF in each patch: ux,uy, theta

//      ////double L = 8; //m
//      ////double h = 1;//4.1602;
//      ////double b = 1;
//      ////double E = 100000;

//      ////double L = 34; //m
//      ////double h = 0.2;//4.1602;
//      ////double b = 0.5; //m
//      ////double E = 2.976E+14; //Pa
//      ////double nuy = 0.2;
//      ////double rho = 114000; //kg/m^3
//      ////double I = b * Math.Pow(h, 3) / 12.0;
//      ////double area = b * h;

//      ////double L = 1; //m
//      ////double h = 0.01;//4.1602;
//      ////double b = 0.01; //m
//      ////double E = 206E9; //Pa
//      ////double rho = 7860; //kg/m^3
//      ////double I = b * Math.Pow(h, 3) / 12.0;
//      ////double area = b * h;

//      ////ViewerForm viewer = new ViewerForm(true);
//      ////NURBSCurve curve = GeometryCreator.CreateStraightNURBSCurve(0, 0, 0, L, 0, 0);
//      ////danh sach cac Curve
//      ////List<NURBSCurve> listCurve = new List<NURBSCurve>();
//      ////curve.pRefinement(4);
//      ////curve.hRefinement(200); //inssert KVector - incerease p, min luoi

//      ////curve.colorKnot = Color.Red;
//      ////curve.colorCurve = Color.Blue;
//      ////curve.colorControlNet = Color.Orange;
//      ////curve.isDrawOriginal = false;
//      ////curve.colorControlPoint = Color.Green;
//      ////curve.isDrawControlNet = false;
//      ////curve.isDrawControlPoint = false;
//      ////curve.widthCurve = 0.5f;
//      ////curve.opacity = 0.2;
//      ////curve.Draw(viewer);

//      ////viewer.UpdateCamera();
//      ////viewer.Run();


//      //////// --------------------Khai bao vat lieu--------------------

//      ////Material steel = new Material("steel");
//      ////steel.AddProperty(new IsotropicElasticity(
//      ////    PairOfIsotropicElasticity.YoungModulusAndPoissonRatio, E, nuy));
//      ////steel.AddProperty(new Density(rho));
//      ////steel.TypeMaterialStructure = TypeMaterialStructure.Elasticity;
//      ////static
//      ////ModelStructureModal model = new ModelStructureModal(Dimension.Beam);

//      ////dinh nghia trang thai bai toan
//      ////model.AddMaterial(steel);
//      ////model.AddPatch(curve);
//      ////model.AttachMaterialToPatch(0, 0);
//      ////model.SetThicknessBeam(h, 0);
//      ////model.SetWidthBeam(b, 0);

//      ////them dai luong can tinh
//      ////model.AddComputeResult(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

//      ////ControlPoint[] selAllCPs = ((AbstractPatch1D)model.GetPatch(0)).SelectAllControlPoints(0);
//      ////ControlPoint[] selCPStart = new ControlPoint[] { selAllCPs[0], selAllCPs[1] };
//      ////ControlPoint[] selCPNextS = new ControlPoint[] { selAllCPs[1] };

//      ///// 1D - 3Dof
//      ////ConstraintValueArrayOfControlPoints cX0 = new ConstraintValueArrayOfControlPoints(selCPStart, 0, new NullFunctionRToR(), 0);
//      ////ConstraintValueArrayOfControlPoints cY0 = new ConstraintValueArrayOfControlPoints(selCPStart, 1, new NullFunctionRToR(), 0);
//      ////ConstraintValueArrayOfControlPoints cThetaZ0 = new ConstraintValueArrayOfControlPoints(selCPStart, 2, new NullFunctionRToR(), 0);

//      ////ConstraintValueArrayOfControlPoints cY1 = new ConstraintValueArrayOfControlPoints(selCPNextS, 1, new NullFunctionRToR(), 0);

//      ////apply constraint
//      ////model.AddConstraint(cX0, cY0, cThetaZ0);
//      /////
//      ////AbstractModel.IsParallelProcesing = false;
//      ////model.SetKinematicsFunction(new LinearFunctionRToR(0, 1, 0, 1), 0);
//      ////model.SetKinematicsFunction(new PolynomialFunctionRToR(0, 1, 0, -4 / (3 * h * h)), 0);
//      ////int numMode = 10;
//      ////model.SetNumberOfMode(numMode);

//      ////model.InitializePatch();
//      ////model.PreProcessing();
//      ////model.Solve();
//      ////model.PostProcessing();

//      ////ViewerForm viewer = new ViewerForm(true);
//      ////model.DrawGeometry(viewer, 5, 0);
//      ////model.DrawModeResult(0, viewer, 5, 5);
//      ////viewer.UpdateCamera();
//      ////viewer.Run();

//      ////viewer = new ViewerForm(true);
//      ////model.DrawGeometry(viewer, 5, 0);
//      ////model.DrawModeResult(1, viewer, 5, 5);
//      ////viewer.UpdateCamera();
//      ////viewer.Run();

//      ////viewer = new ViewerForm(true);
//      ////model.DrawGeometry(viewer, 5, 0);
//      ////model.DrawModeResult(2, viewer, 5, 5);
//      ////viewer.UpdateCamera();
//      ////viewer.Run();

//      ////viewer = new ViewerForm(true);
//      ////model.DrawGeometry(viewer, 5, 0);
//      ////model.DrawModeResult(3, viewer, 5, 5);
//      ////viewer.UpdateCamera();
//      ////viewer.Run();

//      ////viewer = new ViewerForm(true);
//      ////model.DrawGeometry(viewer, 5, 0);
//      ////model.DrawModeResult(4, viewer, 5, 5);
//      ////viewer.UpdateCamera();
//      ////viewer.Run();


//      ////double[] eigenValue = model.GetEigenValue();
//      ////double[] eigenValueRound = new double[eigenValue.Length];
//      ////for (int i = 0; i < eigenValue.Length; i++)
//      ////{
//      ////  eigenValueRound[i] = Math.Round(eigenValue[i], 4);
//      ////}
//      ////double error = 0;

//      ////double[] exactFrequencies = new double[numMode];
//      ////double[] betaNL_fixed_pinned = { 1.875104, 4.694091, 7.854757, 10.995541 };

//      ////for (int i = 1; i <= numMode; i++)
//      ////{
//      ////  exactFrequencies[i - 1] = Math.Pow(betaNL_fixed_pinned[i - 1], 2) * Math.Sqrt(E * I / (rho * b * h * Math.Pow(L, 4)));
//      ////}
//      ////for (int i = 0; i < eigenValue.Length; i++)
//      ////{
//      ////  error += Math.Pow((exactFrequencies[i] - eigenValue[i]) / exactFrequencies[i], 2);
//      ////}
//      ////error = 1.0 / eigenValue.Length * Math.Sqrt(error);
//      ////Console.WriteLine(error);
//    }

//    public static double CalculateErr(int p, int num, out int dof, out double[] eigenValue, out double[] exactFrequencies)
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
//      double nuy = 0;
//      double I = b * Math.Pow(h, 3) / 12.0;
//      double area = b * h;

//      //ViewerForm viewer = new ViewerForm(true);
//      NURBSCurve curve = GeometryCreator.CreateStraightNURBSCurve(0, 0, 0, L, 0, 0);
//      // danh sach cac Curve
//      List<NURBSCurve> listCurve = new List<NURBSCurve>();
//      curve.pRefinement(p - 1);
//      curve.hRefinement(num - 1); //inssert KVector - incerease p, min luoi


//      //// --------------------Khai bao vat lieu--------------------

//      Material steel = new Material("steel");
//      steel.AddProperty(new IsotropicElasticity(
//          PairOfIsotropicElasticity.YoungModulusAndPoissonRatio, E, nuy));
//      steel.AddProperty(new Density(rho));
//      steel.TypeMaterialStructure = TypeMaterialStructure.Elasticity;
//      //static
//      ModelStructureModal model = new ModelStructureModal(Dimension.Beam);

//      //dinh nghia trang thai bai toan
//      model.AddMaterial(steel);
//      model.AddPatch(curve);
//      model.AttachMaterialToPatch(0, 0);
//      model.SetThicknessBeam(h, 0);
//      model.SetWidthBeam(b, 0);

//      // them dai luong can tinh
//      model.AddComputeResult(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

//      ControlPoint[] selAllCPs = ((AbstractPatch1D)model.GetPatch(0)).SelectAllControlPoints(0);
//      ControlPoint[] selCPStart = new ControlPoint[] { selAllCPs[0] };
//      ControlPoint[] selCPNextS = new ControlPoint[] { selAllCPs[1] };

//      // 1D - 3Dof
//      ConstraintValueArrayOfControlPoints cX0 = new ConstraintValueArrayOfControlPoints(selCPStart, 0, new NullFunctionRToR(), 0);
//      ConstraintValueArrayOfControlPoints cY0 = new ConstraintValueArrayOfControlPoints(selCPStart, 1, new NullFunctionRToR(), 0);
//      ConstraintValueArrayOfControlPoints cThetaZ0 = new ConstraintValueArrayOfControlPoints(selCPStart, 2, new NullFunctionRToR(), 0);

//      ConstraintValueArrayOfControlPoints cY1 = new ConstraintValueArrayOfControlPoints(selCPNextS, 1, new NullFunctionRToR(), 0);

//      // apply constraint
//      model.AddConstraint(cX0, cY0, cThetaZ0, cY1);

//      AbstractModel.IsParallelProcesing = false;
//      model.SetKinematicsFunction(new LinearFunctionRToR(0, 1, 0, 1), 0);
//      model.SetKinematicsFunction(new PolynomialFunctionRToR(0, 1, 0, -4 / (3 * h * h)), 0);
//      int numMode = 5;
//      model.SetNumberOfMode(numMode);

//      model.InitializePatch();
//      model.PreProcessing();
//      model.Solve();
//      model.PostProcessing();

//      dof = model.CountDOF();

//      eigenValue = model.GetEigenValue();
//      for (int i = 0; i < eigenValue.Length; i++)
//      {
//        eigenValue[i] = Math.Round(eigenValue[i], 4);
//      }

//      double error = 0;

//      exactFrequencies = new double[numMode];
//      double[] betaNL_fixed_free_full = new double[10] { 1.875104069, 4.694091133, 7.854757438, 10.99554073, 14.13716839, 17.27875953, 20.42035225, 23.5619449, 26.70353756, 29.84513021 };

//      for (int i = 1; i <= numMode; i++)
//      {

//        exactFrequencies[i - 1] = Math.Round(Math.Pow(betaNL_fixed_free_full[i - 1], 2) * Math.Sqrt(E * I / (rho * b * h * Math.Pow(L, 4))), 4);
//      }
//      for (int i = 0; i < eigenValue.Length; i++)
//      {
//        error += Math.Pow((exactFrequencies[i] - eigenValue[i]), 2);
//      }

//      error = 1.0 / numMode * error; // Mean Square Err
//      //error = Math.Sqrt(1.0 / numMode * error);
//      return error;
//    }
//  }
//}