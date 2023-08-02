using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOFEM
{
    internal class Plane2DElement : Element
    {
        public Plane2DElement(Node n1, Node n2, Node n3, Node n4, double E, double nu, double A)
        {
            this.E = E;
            this.nu = nu;
            this.A = A;
            nodes = new Node[] { n1, n2, n3, n4 };
        }
        public double Area()
        {
            return 0;
        }

        internal override DenseMatrix ComputeStiffnessMatrix()
        {
            int countDOF = tArray.Length;
            //int countDOF = 8;

            DenseMatrix Ke = new DenseMatrix(countDOF);

            DenseMatrix C = ComputeMaterialTensor();
            DenseMatrix B = new DenseMatrix(3, countDOF);

            double[] xi = { -1 / Math.Sqrt(3), 1 / Math.Sqrt(3) };


            for (int i = 0; i < 2; i++)
            {
                double xiI = xi[i];
                for (int j = 0; j < 2; j++)
                {
                    double wi = 1;
                    double wj = 1;
                    double etaI = xi[j];

                    double[] shapeDeviationXi = { -1.0 / 4.0 * (1 - etaI),   //N1 
                                                  1.0 / 4.0 * (1 - etaI),   //N2
                                                  1.0 / 4.0 * (1 + etaI),   //N3
                                                 -1.0 / 4.0 * (1 + etaI)};  //N4

                    double[] shapeDeviationEta = { -1.0 / 4.0 * (1 - xiI),   //N1 
                                                   -1.0 / 4.0 * (1 + xiI),   //N2
                                                    1.0 / 4.0 * (1 + xiI),   //N3
                                                    1.0 / 4.0 * (1 - xiI) }; //N4

                    DenseMatrix J = new DenseMatrix(2, 2);

                    for (int k = 0; k < 2; k++)
                    {
                        for (int m = 0; m < 2; m++)
                        {
                            for (int n = 0; n < nodes.Length; n++)
                            {
                                if (k == 0)
                                {
                                    J[k, m] += shapeDeviationXi[n] * nodes[n].GetPosition(m);
                                }
                                else
                                {
                                    J[k, m] += shapeDeviationEta[n] * nodes[n].GetPosition(m);

                                }
                            }
                        }
                    }

                    double detJ = Math.Abs(J.Determinant());

                    /*
                     * B = J.Inverse() * D
                     * 
                     *     | dN1/dXi      0       dN2/dXi      0        ... |
                     * D = |   0       dN1/dEta      0       dN2/dEta   ... |
                     *     | dN1/dEta   dN1/dXi   dN2/dEta   dN2/dXi    ... |
                     * 
                     */

                    DenseMatrix D = new DenseMatrix(3, countDOF);

                    for (int k = 0; k < 3; k++)
                    {
                        for (int m = 0; m < countDOF; m++)
                        {
                            for (int n = 0; n < shapeDeviationEta.Length; n++)
                            {
                                if (m % 2 == 0)
                                {
                                    D[0, m] = D[2, m + 1] = shapeDeviationXi[n];
                                }
                                else
                                {
                                    D[1, m] = D[2, m - 1] = shapeDeviationEta[n];
                                }
                            }
                        }
                    }

                    //var Jinverse = J.Inverse();

                    //B = (DenseMatrix)(Jinverse * D);

                    DenseMatrix dNdxi = new DenseMatrix(2, 4);
                    dNdxi[0, 0] = -1.0 / 4.0 * (1 - xiI);
                    dNdxi[0, 1] = 1.0 / 4.0 * (1 - xiI);
                    dNdxi[0, 2] = 1.0 / 4.0 * (1 + xiI);
                    dNdxi[0, 3] = -1.0 / 4.0 * (1 + xiI);
                    dNdxi[1, 0] = -1.0 / 4.0 * (1 - etaI);
                    dNdxi[1, 1] = -1.0 / 4.0 * (1 + etaI);
                    dNdxi[1, 2] = 1.0 / 4.0 * (1 + etaI);
                    dNdxi[1, 3] = 1.0 / 4.0 * (1 - etaI);

                    for (int kk = 0; kk < nodes.Length; kk++)
                    {
                        B[0, 2 * kk] = J[1, 1] * dNdxi[0, kk] - J[0, 1] * dNdxi[1, kk];
                        B[1, 2 * kk + 1] = -J[1, 0] * dNdxi[0, kk] + J[0, 0] * dNdxi[1, kk];
                        B[2, 2 * kk] = -J[1, 0] * dNdxi[0, kk] + J[0, 0] * dNdxi[1, kk];
                        B[2, 2 * kk + 1] = J[1, 1] * dNdxi[0, kk] - J[0, 1] * dNdxi[1, kk];
                    }

                    B = (DenseMatrix)(B / detJ);


                    Ke += (DenseMatrix)(wi * wj * detJ * B.Transpose() * C * B);
                }
            }

            return Ke;
        }

        //private DenseMatrix ComputeMaterialTensor()
        //{
        //    DenseMatrix C = new DenseMatrix(3);

        //    C[0, 0] = C[1, 1] = E / (1 - nu * nu);
        //    C[0, 1] = C[1, 0] = nu * E / (1 - nu * nu);
        //    C[2, 2] = E / (2 * (1 + nu));
        //    //C[0, 2] = C[1, 2] = C[2, 0] = C[2, 1] = 0;

        //    return C;
        //}

        private DenseMatrix ComputeMaterialTensor()
        {
            DenseMatrix C = new DenseMatrix(3);

            double c = (1 - nu) * E / ((1 + nu) * (1 - 2 * nu));

            C[0, 0] = C[1, 1] = c;
            C[0, 1] = C[1, 0] = c * (nu / (1-nu));
            C[2, 2] = c * ((1 - 2 * nu) / (2 * (1 - nu)));
            //C[0, 2] = C[1, 2] = C[2, 0] = C[2, 1] = 0;

            return C;
        }
    }
}
