using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOFEM
{
    internal class Truss3DElement : Element
    {
        public Truss3DElement(Node n1, Node n2, double E, double A)
        {
            this.E = E;
            this.A = A;
            nodes = new Node[] { n1, n2 };
        }
        public double Length()
        {
            return nodes[0].Distance(nodes[1]);
        }

        internal override DenseMatrix ComputeStiffnessMatrix()
        {
            int countDOF = tArray.Length;
            DenseMatrix Ke = new DenseMatrix(countDOF);
            double x1 = nodes[0].GetPosition(0);
            double y1 = nodes[0].GetPosition(1);
            double z1 = nodes[0].GetPosition(2);
            double x2 = nodes[1].GetPosition(0);
            double y2 = nodes[1].GetPosition(1);
            double z2 = nodes[1].GetPosition(2);
            double L = Length();
            double Cx = (x2 - x1) / L;
            double Cy = (y2 - y1) / L;
            double Cz = (z2 - z1) / L;
            double EAL = E * A / L;
            
            Ke[0, 0] = EAL * Cx * Cx;
            Ke[0, 1] = Ke[1, 0] = EAL * Cx * Cy;
            Ke[0, 2] = Ke[2, 0] = EAL * Cx * Cz;
            Ke[0, 3] = Ke[3, 0] = -EAL * Cx * Cx;
            Ke[0, 4] = Ke[4, 0] = -EAL * Cx * Cy;
            Ke[0, 5] = Ke[5, 0] = -EAL * Cx * Cz;

            Ke[1, 1] = EAL * Cy * Cy;
            Ke[1, 2] = Ke[2, 1] = EAL * Cy * Cz;
            Ke[1, 3] = Ke[3, 1] = -EAL * Cx * Cy;
            Ke[1, 4] = Ke[4, 1] = -EAL * Cy * Cy;
            Ke[1, 5] = Ke[5, 1] = -EAL * Cy * Cz;

            Ke[2, 2] = EAL * Cz * Cz;
            Ke[2, 3] = Ke[3, 2] = -EAL * Cx * Cz;
            Ke[2, 4] = Ke[4, 2] = -EAL * Cy * Cz;
            Ke[2, 5] = Ke[5, 2] = -EAL * Cz * Cz;

            Ke[3, 3] = EAL * Cx * Cx;
            Ke[3, 4] = Ke[4, 3] = EAL * Cx * Cy;
            Ke[3, 5] = Ke[5, 3] = EAL * Cx * Cz;

            Ke[4, 4] = EAL * Cy * Cy;
            Ke[4, 5] = Ke[5, 4] = EAL * Cy * Cz;

            Ke[5, 5] = EAL * Cz * Cz;

            Console.WriteLine("L = " + L);
            Console.WriteLine("Cx = " + Cx);
            Console.WriteLine("Cy = " + Cy);
            Console.WriteLine("Cz = " + Cz);
            Console.WriteLine("EAL = " + EAL);

            Console.WriteLine(Ke);

            Console.WriteLine("\n");

            return Ke;
        }
    }
}
