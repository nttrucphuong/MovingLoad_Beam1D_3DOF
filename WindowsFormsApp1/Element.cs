using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOFEM
{
    internal abstract class Element
    {
        public int[] tArray { get; set; }
        public Node[] nodes { get; set; }
        public double E { get; set; } //Phuong thuc cho user ben ngoai truy cap vao du lieu cua E
        public double nu { get; set; }

        public double A { get; set; }

        internal abstract DenseMatrix ComputeStiffnessMatrix();
    }
}
