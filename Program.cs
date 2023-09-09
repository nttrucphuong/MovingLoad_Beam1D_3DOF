using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using DEMSoft.Drawing;
using DEMSoft.NURBS;
using DEMSoft.Plot;
using DEMSoft.IGA;
using System.Drawing;
using DEMSoft.EngineeringData;
using DEMSoft.Common;
using DEMSoft.Function;
using System.Drawing.Printing;
using System.CodeDom.Compiler;
using System.Diagnostics;

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


            curve.SetDegreeOnAllDirections(3);
            curve.hRefinement(1, 0); //inssert KVector - increase p, min luoi
            //curve.isColorfulFace = false;
            //curve.colorControlNet = Color.Blue;
            //curve.colorSurface = Color.Green;
            curve.colorKnot = Color.Red;
            curve.colorCurve = Color.Blue;
            curve.isDrawOriginal = false;
            curve.colorControlPoint = Color.Green;
            curve.Draw(viewer);

            viewer.UpdateCamera();
            viewer.Run();


            ////// --------------------Khai bao vat lieu--------------------

            Material steel = new Material("steel");
            steel.AddProperty(new IsotropicElasticity(
                PairOfIsotropicElasticity.YoungModulusAndPoissonRatio, 1000, 0.3));
            steel.TypeMaterialStructure = TypeMaterialStructure.Elasticity;


            ModelStructureStatic model = new ModelStructureStatic(Dimension.Beam);
            // dinh nghia trang thai bai toan
            
            model.AddMaterial(steel);
            //model.AddComputeResult(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

            model.AddPatch(curve);

            ////// them dai luong can tinh 
            ////model.AddComputeResult(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAEQV);

            ////apply constraint
            //// clamped edge, patch 6-9
            ConstraintValueEdge2D cX5 = new ConstraintValueEdge2D(
                (PatchStructure2D)model.GetPatch(6), 1, 0, new NullFunctionRToR(), 0);
            ConstraintValueEdge2D cY5 = new ConstraintValueEdge2D(
                (PatchStructure2D)model.GetPatch(6), 1, 1, new NullFunctionRToR(), 0);
            ConstraintValueEdge2D cX6 = new ConstraintValueEdge2D(
                (PatchStructure2D)model.GetPatch(7), 1, 0, new NullFunctionRToR(), 0);
            ConstraintValueEdge2D cY6 = new ConstraintValueEdge2D(
                (PatchStructure2D)model.GetPatch(7), 1, 1, new NullFunctionRToR(), 0);
            ConstraintValueEdge2D cX7 = new ConstraintValueEdge2D(
                (PatchStructure2D)model.GetPatch(8), 1, 0, new NullFunctionRToR(), 0);
            ConstraintValueEdge2D cY7 = new ConstraintValueEdge2D(
                (PatchStructure2D)model.GetPatch(8), 1, 1, new NullFunctionRToR(), 0);
            ConstraintValueEdge2D cX8 = new ConstraintValueEdge2D(
                (PatchStructure2D)model.GetPatch(9), 1, 0, new NullFunctionRToR(), 0);
            ConstraintValueEdge2D cY8 = new ConstraintValueEdge2D(
                (PatchStructure2D)model.GetPatch(9), 1, 1, new NullFunctionRToR(), 0);


            ////// displacement edges, patch 6-9
            //double du = 0.1;
            //ConstraintValueEdge2D cY0 = new ConstraintValueEdge2D(
            //    (PatchStructure2D)model.GetPatch(0), 1, 1, new ConstantFunctionRToR(1.0), du);
            //ConstraintValueEdge2D cY1 = new ConstraintValueEdge2D(
            //    (PatchStructure2D)model.GetPatch(1), 1, 1, new ConstantFunctionRToR(1.0), du);
            //ConstraintValueEdge2D cY2 = new ConstraintValueEdge2D(
            //    (PatchStructure2D)model.GetPatch(2), 1, 1, new ConstantFunctionRToR(1.0), du);
            //ConstraintValueEdge2D cY3 = new ConstraintValueEdge2D(
            //    (PatchStructure2D)model.GetPatch(3), 1, 1, new ConstantFunctionRToR(1.0), du);
            //ConstraintValueEdge2D cX0 = new ConstraintValueEdge2D(
            //    (PatchStructure2D)model.GetPatch(0), 1, 0, new ConstantFunctionRToR(1.0), 0);
            //ConstraintValueEdge2D cX1 = new ConstraintValueEdge2D(
            //    (PatchStructure2D)model.GetPatch(1), 1, 0, new ConstantFunctionRToR(1.0), 0);
            //ConstraintValueEdge2D cX2 = new ConstraintValueEdge2D(
            //    (PatchStructure2D)model.GetPatch(2), 1, 0, new ConstantFunctionRToR(1.0), 0);
            //ConstraintValueEdge2D cX3 = new ConstraintValueEdge2D(
            //    (PatchStructure2D)model.GetPatch(3), 1, 0, new ConstantFunctionRToR(1.0), 0);
            //model.AddConstraint(cX5, cX6, cX7, cX8, cY5, cY6, cY7, cY8,
            //    cY0, cY1, cY2, cY3, cX0, cX1, cX2, cX3);
            //Stopwatch sw = new Stopwatch();
            //sw.Start();
            //model.IsSaveStepByStep = false;
            ////model.NumberOfStepStorageNonlinear = 10;
            //AbstractModel.IsParallelProcesing = true;

            //model.InitializePatch();
            //model.PreProcessing();
            //model.Solve();
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
