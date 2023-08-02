using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOFEM
{
    internal class Node
    {
        private double[] x;
        private double[] u;

        public int[] TArray { get; internal set; }
        public int[] TArrayGlobal { get; internal set; }


        public Node(double x, double y, double z)
        {
            this.x = new double[] { x, y, z };
            u = new double[3];
        }
        public double GetPosition(int index)
        {
            return this.x[index];
        }
        public double Distance(Node node)
        {
            return Math.Sqrt(Math.Pow(this.x[0] - node.x[0], 2)
                            + Math.Pow(this.x[1] - node.x[1], 2)
                            + Math.Pow(this.x[2] - node.x[2], 2));
        }
        public double GetDisplacement(int index)
        {
            return u[index];
        }
        public void setDisplacement(int index, double v)
        {
            u[index] = v;
        }
    }
}
