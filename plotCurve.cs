//using DEMSoft.Drawing;
//using DEMSoft.EngineeringData;
//using DEMSoft.NURBS;
//using DEMSoft.Plot;
//using System;
//using System.Collections.Generic;
//using System.Drawing;
//using System.Linq;
//using System.Security.Cryptography;
//using System.Text;
//using System.Threading.Tasks;

//namespace PlotBSF
//{
//  internal class plotCurve
//  {
//    static void Main()
//    {
//      //int p = 2;
//      //double[] kv = { 0, 0, 0, 1, 2, 2, 3, 4, 5, 5, 5 };
//      ////double[] kv = { 0, 0, 0, 0.2, 0.4, 0.6, 0.8, 0.8, 1, 1, 1 };
//      ////double[] kv = { 0, 0, 0,1, 1,1 };

//      //ControlPoint cp1 = new ControlPoint(0, 3);
//      //ControlPoint cp2 = new ControlPoint(1.5, 5.5);
//      //ControlPoint cp3 = new ControlPoint(2.5, 4.5);
//      //ControlPoint cp4 = new ControlPoint(4.5, 5.5);
//      //ControlPoint cp5 = new ControlPoint(3, 1.5);
//      //ControlPoint cp6 = new ControlPoint(7.5, 1.5);
//      //ControlPoint cp7 = new ControlPoint(6, 3.5);
//      //ControlPoint cp8 = new ControlPoint(8.5, 4.5);

//      //=============================================
//      //p + n + 1 = numKnot //p-order, n; num cps
//      //int p = 4;
//      //double[] kv_ = { 0, 0, 0, 0, 0, 0.25, 0.5, 1, 1, 1, 1, 1 }; //{ 0, 0, 0, 0, 0, 1 / 3, 2 / 3, 5 / 6, 1, 1, 1, 1, 1 };
//      //KnotVector kv = new KnotVector(kv_);

//      //ControlPoint cp1 = new ControlPoint(1, 0);
//      //ControlPoint cp2 = new ControlPoint(0, 3);
//      //ControlPoint cp3 = new ControlPoint(2, 5);
//      //ControlPoint cp4 = new ControlPoint(4, 4);
//      //ControlPoint cp5 = new ControlPoint(3, 2.5);
//      //ControlPoint cp6 = new ControlPoint(6, 1.5);
//      //ControlPoint cp7 = new ControlPoint(7, 3.5);

//      //ControlPoint cp1 = new ControlPoint(1, 0);
//      //ControlPoint cp2 = new ControlPoint(0, 3);
//      //ControlPoint cp3 = new ControlPoint(2, 5);
//      //ControlPoint cp4 = new ControlPoint(4, 4);
//      //ControlPoint cp5 = new ControlPoint(3, 2.5);
//      //ControlPoint cp6 = new ControlPoint(5.5, 2);
//      //ControlPoint cp7 = new ControlPoint(6.5, 2.5);
//      //ControlPoint cp8 = new ControlPoint(7, 3.5);
//      ///*ControlPoint[] cps = new ControlPoint[]*/ { cp1, cp2, cp3, cp4, cp5, cp6, cp7 }; //ControlPoint[kv.Length - p - 1];


//      //====================================
//      //====================================
//      NURBSCurve curve1 = GeometryCreator.CreateStraightNURBSCurve(0, 0, 0, 1, 0, 0);
//      //NURBSCurve curve1 = new NURBSCurve(p, kv, cps, weights);
//      //BSplineCurve curve1 = new BSplineCurve(p, kv, cps);
//      //====================================================
//      //ViewerForm viewer = new ViewerForm(true);
//      //curve1.hRefinement(1);
//      //curve1.pRefinement(1);
//      //curve1.pRefinement(1);
//      //double[][] data = curve1.GetDataOnCurve();
//      //KnotVector kvNew = curve1.Basis.KnotVector;
//      curve1.colorCurve = Color.Blue;
//      curve1.colorTangent = Color.Black;
//      curve1.colorControlPoint = Color.Red;
//      curve1.resolution = 20;
//      curve1.isDrawKnot = true;
//      curve1.isDrawControlPoint = false;

//      //curve1.Draw(viewer);
//      //viewer.UpdateCamera();
//      //viewer.Run();
//      //===================================
//      //int num = 201;
//      //double[] x = new double[curve1.GetDataOnCurve().Length];
//      //double[] y = new double[curve1.GetDataOnCurve().Length];

//      ////for (int i = 0; i < y.Length; i++)
//      ////{
//      ////  x[i] = new double[curve1.GetDataOnCurve()[0].Length];
//      ////  y[i] = new double[curve1.GetDataOnCurve()[1].Length];
//      ////}
//      //double dx = 4.9999999999 / (num - 1);

//      //for (int i = 0; i < curve1.GetDataOnCurve().Length; i++)
//      //{
//      //  x[i] = curve1.GetDataOnCurve()[i][0];
//      //  y[i] = curve1.GetDataOnCurve()[i][1];
//      //}

//      //Plotter plotter = new Plotter();
//      //for (int i = 0; i < y.Length; i++)
//      //{
//      //  PlotLine l1 = new PlotLine();
//      //  l1.InputData(x, y, "curve");
//      //  l1.SetWidth(3);
//      //  plotter.AddPlot(l1);
//      //}
//      ////plotter.SetShowLegend(true);
//      //plotter.Plot();
//      //==========================
//      curve1.pRefinement(1);
//      curve1.InsertKnot(1.0 / 3.0, 1);
//      curve1.InsertKnot(2.0 / 3.0, 1);
//      //curve1.pRefinement(1);
//      NURBSBasisFunction basicFunction = (NURBSBasisFunction)curve1.Basis;
//      int num = 100;
//      double[] x = new double[num];
//      double[][] y = new double[basicFunction.GetCountBasisFunction()][];
//      for (int i = 0; i < y.Length; i++)
//      {
//        y[i] = new double[num];
//      }
//      double dx = 0.9999999999 / (num - 1);
//      for (int i = 0; i < num; i++)
//      {
//        x[i] = i * dx;
//        for (int j = 0; j < y.Length; j++)
//        {
//          //y[j][i] = basicFunction.GetValueBasisFunction(x[i], j);
//          y[j][i] = basicFunction.GetDerivativeBasisFunction(x[i], j, 0)[0];
//          //y[j][i] = basicFunction.GetValueBasisFunction(x[i], 4);

//        }
//      }
//      ColorsRange colorRange = new ColorsRange(ColorType.Jet, y.Length);
//      Plotter plotter = new Plotter();
//      for (int i = 0; i < y.Length; i++)
//      {
//        PlotLine l1 = new PlotLine();
//        l1.InputData(x, y[i], "N" + i.ToString());
//        double[] c = colorRange.GetColor(i);
//        l1.SetColor(c[0], c[1], c[2]);
//        l1.SetWidth(3);
//        plotter.AddPlot(l1);
//      }
//      //plotter.SetShowLegend(true);
//      plotter.Plot();
//    }
//  }
//}
