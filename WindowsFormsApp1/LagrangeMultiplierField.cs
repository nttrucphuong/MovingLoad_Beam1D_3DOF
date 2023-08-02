using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOFEM
{
    internal class LagrangeMultiplierField
    {
        public ArrayList lagrangeField { get; set; }
        public LagrangeMultiplierField(ArrayList lagrangeField)
        {
            this.lagrangeField = lagrangeField;
        }
        //private int directionOfField()
        //{
        //    if (nodes[0].GetPosition(0) == nodes[1].GetPosition(0))
        //    {
        //        return 0;
        //    }
        //    else
        //    {
        //        return 1;
        //    }
        //}
        //internal List<double> ComputeJacobian()
        //{
        //    int index = 1;
        //    List<double> J = new List<double>();
        //    for (int i = 0; i < nodes.Count; i++)
        //    {
        //        double x1 = nodes[i].GetPosition(index);
        //        double DxDxi = 0;
        //        if (i + 1 == nodes.Count)
        //        {
        //            DxDxi = (nodes[i - 1].GetPosition(index) - x1) / 2.0;
        //            if (DxDxi > 0)
        //            {
        //                J.Add(DxDxi);
        //            }
        //            else
        //            {
        //                J.Add(0);
        //            }
        //        }
        //        else
        //        {
        //            DxDxi = (nodes[i + 1].GetPosition(index) - x1) / 2.0;
        //            if (DxDxi > 0)
        //            {
        //                J.Add(DxDxi);
        //            }
        //            else
        //            {
        //                J.Add(0);
        //            }
        //        }
        //    }
        //    return J;
        //}
    }
}
