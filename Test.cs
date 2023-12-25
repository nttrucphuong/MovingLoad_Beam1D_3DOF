//using DEMSoft.Plot;
//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using System.Threading.Tasks;

//namespace MovingLoad
//{
//  internal class Test
//  {
//    static void Main()
//    {
//      double L = 34; //m
//      double h = 0.2;//4.1602;
//      double b = 0.5; //m
//      double E = 2.976E+14; //Pa
//      double nuy = 0.2;
//      double rho = 114000; //kg/m^3
//      double I = b * Math.Pow(h, 3) / 12.0;

//      int numMode = 10;
//      double x = 0.0;
//      double pi = Math.PI;
//      double[] beta = new double[10] { 1.875104069, 4.694091133, 7.854757438, 10.99554073, 14.13716839, 17.27875953, 20.42035225, 23.5619449, 26.70353756, 29.84513021 };

//      List<double> listExactFrequenciesCantilever = new List<double>();

//      for (int i = 0; i < numMode; i++)
//      {
//        listExactFrequenciesCantilever.Add(Math.Pow(beta[i], 2) * Math.Sqrt(E * I / (rho * b * h * Math.Pow(L, 4))));
//      }

      
//    }
//  }
//}
