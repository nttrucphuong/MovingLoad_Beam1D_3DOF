//using DEMSoft.Drawing;
//using DEMSoft.NURBS;
//using DEMSoft.Plot;
//using System;
//using System.Collections.Generic;
//using System.Drawing;
//using System.Linq;
//using System.Text;
//using System.Threading.Tasks;

//namespace PlotBSF
//{
//  internal class PlotCircle
//  {
//    static void Main()
//    {
//      int p = 2;
//      double a = 1.0 / Math.Sqrt(2.0);
//      double[] weigth = { 1, a, 1, a, 1, a, 1, a, 1 };
//      double[] kv = { 0, 0, 0, 1.0 / 4.0, 1.0 / 4.0, 2.0 / 4.0, 2.0 / 4.0, 3.0 / 4.0, 3.0 / 4.0, 1, 1, 1 };

//      ControlPoint cp1 = new ControlPoint(1, 0);
//      ControlPoint cp2 = new ControlPoint(1, 1);
//      ControlPoint cp3 = new ControlPoint(0, 1);
//      ControlPoint cp4 = new ControlPoint(-1, 1);
//      ControlPoint cp5 = new ControlPoint(-1, 0);
//      ControlPoint cp6 = new ControlPoint(-1, -1);
//      ControlPoint cp7 = new ControlPoint(0, -1);
//      ControlPoint cp8 = new ControlPoint(1, -1);
//      ControlPoint cp9 = new ControlPoint(1, 0);

//      ControlPoint[] cps = new ControlPoint[] { cp1, cp2, cp3, cp4, cp5, cp6, cp7, cp8, cp9 }; //kv.Length - p - 1


//      //BSplineBasisFunction basicFunction = new BSplineBasisFunction(p, kv);
//      NURBSCurve circle = new NURBSCurve(p, kv, cps, weigth);
//      //====================================================
//      ViewerForm viewer = new ViewerForm(true);
//      //curve1.hRefinement(2);
//      //curve1.pRefinement(1);

//      circle.colorCurve = Color.Red;
//      circle.colorControlNet = Color.Blue;
//      circle.colorControlPoint = Color.Green;
//      circle.widthControlNet = (float)2.0;
//      circle.widthCurve = 3;
//      circle.resolution = 20;

//      circle.Draw(viewer);
//      viewer.UpdateCamera();
//      viewer.Run();
//      //===================================
//      //int num = 201;
//      //double[] x = new double[circle.GetDataOnCurve().Length];
//      //double[] y = new double[circle.GetDataOnCurve().Length];

//      ////for (int i = 0; i < y.Length; i++)
//      ////{
//      ////  x[i] = new double[curve1.GetDataOnCurve()[0].Length];
//      ////  y[i] = new double[curve1.GetDataOnCurve()[1].Length];
//      ////}
//      //double dx = 0.9999999999 / (num - 1);

//      //for (int i = 0; i < circle.GetDataOnCurve().Length; i++)
//      //{
//      //  x[i] = circle.GetDataOnCurve()[i][0];
//      //  y[i] = circle.GetDataOnCurve()[i][1];
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
//    }
//  }
//}
