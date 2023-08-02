using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOFEM
{
    internal class Force
    {
        public double[] value { get; set; }
        private Node node;
        public Force(Node n, params double[] fValue)
        {
            this.node = n;
            value = fValue;
        }
        public DenseVector ComputeVectorForce()
        {
            return new DenseVector(value);
        }
        public Node GetNode()
        {
            return node;
        }
    }
}
