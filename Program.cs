using DEMSoft.Common;
using DEMSoft.Drawing;
using DEMSoft.EngineeringData;
using DEMSoft.Function;
using DEMSoft.IGA;
using DEMSoft.NURBS;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;

namespace SampleTesting
{
    internal static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            double r = 1;
            double a = 4;

            ViewerForm viewer = new ViewerForm(true);
            NURBSCurve curve = GeometryCreator.CreateStraightNURBSCurve(0, 0, 0, 8, 0, 0);


            // danh sach cac Curve
            List<NURBSCurve> listCurve = new List<NURBSCurve>();

            // Create half Rect surface with Circle Hole
            int p0 = 2;
            double[] w = { 1, 1 };
            double[] kv = { 0, 0, 1, 1 };

            ControlPoint[,] cps0 = new ControlPoint[2, 1];

            cps0[0, 0] = new ControlPoint(0, 0);
            cps0[1, 0] = new ControlPoint(0, 8);


            curve.SetDegreeOnAllDirections(2);
            //curve.hRefinement(1, 0); //inssert KVector - increase p, min luoi
            //curve.isColorfulFace = false;
            //curve.colorControlNet = Color.Blue;
            //curve.colorSurface = Color.Green;
            curve.colorKnot = Color.Red;
            curve.colorCurve = Color.Blue;
            curve.isDrawOriginal = false;
            curve.colorControlPoint = Color.Green;
            //curve.Draw(viewer);

            //viewer.UpdateCamera();
            //viewer.Run();


            ////// --------------------Khai bao vat lieu--------------------

            Material steel = new Material("steel");
            steel.AddProperty(new IsotropicElasticity(
                PairOfIsotropicElasticity.YoungModulusAndPoissonRatio, 1000, 0.3));
            steel.TypeMaterialStructure = TypeMaterialStructure.Elasticity;


            ModelStructureStatic model = new ModelStructureStatic(Dimension.Beam);
            // dinh nghia trang thai bai toan

            model.AddMaterial(steel);
            model.AddPatch(curve);
            model.AttachMaterialToPatch(0, 0);

            // them dai luong can tinh 
            model.AddComputeResult(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

            ControlPoint[] selAllCPs = ((AbstractPatch1D)model.GetPatch(0)).SelectAllControlPoints(0);
            ControlPoint[] selCP = new ControlPoint[] { selAllCPs[0], selAllCPs[1] };
            ConstraintValueArrayOfControlPoints cY = new ConstraintValueArrayOfControlPoints(selCP, 0, new NullFunctionRToR(), 0);
            ////apply constraint
            //ConstraintValueEdge1D cX1 = new ConstraintValueEdge1D(
            //    (AbstractPatch1D)model.GetPatch(0), 0, 0, new NullFunctionRToR(), 0);
            //ConstraintValueEdge1D cY1 = new ConstraintValueEdge1D(
            //  (AbstractPatch1D)model.GetPatch(0), 0, 1, new NullFunctionRToR(), 0);
            //ConstraintValueEdge1D cThetaZ1 = new ConstraintValueEdge1D(
            //  (AbstractPatch1D)model.GetPatch(0), 0, 2, new NullFunctionRToR(), 0);
            //model.AddConstraint(cX1, cY1, cThetaZ1);
            model.AddConstraint(cY);
            Stopwatch sw = new Stopwatch();
            sw.Start();
            model.IsSaveStepByStep = false;
            AbstractModel.IsParallelProcesing = false;
            model.SetKinematicsFunction(new LinearFunctionRToR(0, 1, 0, 1), 0);
            model.InitializePatch();
            model.PreProcessing();
            model.Solve();
            //model.PostProcessing();

            //sw.Stop();
            //double t = sw.Elapsed.TotalMilliseconds;
            //model.Extrapolation(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

            //viewer = new ViewerForm(true);
            //model.DrawResult(Result.SIGMAEQV, viewer, 20, false, 20);
            ////model.DrawResult(Result.USUM, viewer, 5, false, 1);

            //viewer.UpdateCamera();
            //viewer.Run();
        }
    }
}
