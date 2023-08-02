using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOFEM
{
    internal class Subdomain
    {
        public int[] tArrayGlobal { get; internal set; }
        public List<Node> nodes { get; set; }
        public List<Element> elements { get; set; }

        public Subdomain(List<Node> nodes, List<Element> element)
        {
            this.nodes = nodes;
            this.elements = element;
        }
        internal DenseMatrix ComputeSubdomainStiffnessMatrix()
        {
            int NOF = nodes[0].TArray.Length;
            int SubdomainDOF = nodes.Count*NOF;
            DenseMatrix K_subdomain = new DenseMatrix(SubdomainDOF);

            foreach (var elem in elements)
            {
                DenseMatrix Ke = elem.ComputeStiffnessMatrix();
                int[] tArrayElem = elem.tArray;
                for (int i = 0; i < Ke.RowCount; i++)
                {
                    for (int j = 0; j < Ke.ColumnCount; j++)
                    {
                        K_subdomain[tArrayElem[i], tArrayElem[j]] += Ke[i, j];
                    }
                }
            }
            Console.WriteLine(K_subdomain);
            return K_subdomain;
        }
        //internal DenseMatrix ComputeCouplingMatrix_Subdomain(List<double> lagrangefield)
        //{
        //    int NOF = nodes[0].TArray.Length;
        //    int numberOfFiels = lagrangefield.Count;
        //    DenseMatrix C_subdomain = new DenseMatrix(nodes.Count * NOF, NOF * numberOfFiels);

        //    return C_subdomain;
        //}
    }
}
