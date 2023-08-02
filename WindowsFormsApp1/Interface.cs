using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOFEM
{
    internal class MortarContactElememt
    {
        public Subdomain sub1;
        public Subdomain sub2;
        public List<double[]> lagrangefield;
        public List<Node> master { get; set; }
        public List<Node> slave { get; set; }

        public MortarContactElememt(Subdomain sub1, Subdomain sub2, List<Node> master, List<Node> slave, List<double[]> lagrangefield)
        {
            this.sub1 = sub1;
            this.sub2 = sub2;
            this.master = master;
            this.slave = slave;
            this.lagrangefield = lagrangefield;
        }

        internal DenseMatrix CouplingMatrix()
        {
            int NOF = 2;
            int index = 1;
            double[] gausspoint = { -1.0 / Math.Sqrt(3), 1.0 / Math.Sqrt(3) };
            DenseMatrix C1 = new DenseMatrix(lagrangefield.Count * NOF, (sub2.nodes.Count + sub1.nodes.Count) * NOF);
            for (int i = 0; i < lagrangefield.Count; i++)
            {
                for (int j = 0; j < slave.Count; j++)
                {
                    double c = 0;
                    if (j == 0)
                    {
                        double Jacobi1 = (slave[j + 1].GetPosition(index) - slave[j].GetPosition(index))/2.0;
                        double Jacobi2 = (lagrangefield[i][1] - lagrangefield[i][0]) / (slave[j + 1].GetPosition(index) - slave[j].GetPosition(index));
                        for (int g = 0; g < gausspoint.Length; g++)
                        {
                            double XI = ((lagrangefield[i][1] + lagrangefield[i][0]) - (slave[j + 1].GetPosition(index) + slave[j].GetPosition(index)))/(slave[j + 1].GetPosition(index) - slave[j].GetPosition(index)) +
                                        Jacobi2*gausspoint[g];
                            double x = ((slave[j + 1].GetPosition(index) + slave[j].GetPosition(index)) / 2.0) + Jacobi1 * XI;
                            double[] n1 = { slave[j].GetPosition(index), 1 };
                            double[] n2 = { slave[j+1].GetPosition(index), 0 };
                            double ValueShapeFunction = ComputeShapeFunction(n1, n2, x);
                            c += Jacobi1*Jacobi2*ValueShapeFunction;
                        }
                        C1[NOF * i, (slave[j].TArray)[0]] += -c;
                        C1[NOF * i + 1, (slave[j].TArray)[0] + 1] += -c;
                    }
                    else if (j == slave.Count - 1)
                    {
                        double Jacobi1 = (slave[j].GetPosition(index) - slave[j-1].GetPosition(index)) / 2.0;
                        double Jacobi2 = (lagrangefield[i][1] - lagrangefield[i][0]) / (slave[j].GetPosition(index) - slave[j - 1].GetPosition(index));
                        for (int g = 0; g < gausspoint.Length; g++)
                        {
                            double XI = ((lagrangefield[i][1] + lagrangefield[i][0]) - (slave[j].GetPosition(index) + slave[j - 1].GetPosition(index))) / (slave[j].GetPosition(index) - slave[j - 1].GetPosition(index)) +
                                        Jacobi2 * gausspoint[g];
                            double x = ((slave[j].GetPosition(index) + slave[j - 1].GetPosition(index)) / 2.0) + Jacobi1 * XI;
                            double[] n2 = { slave[j].GetPosition(index), 1 };
                            double[] n1 = { slave[j - 1].GetPosition(index), 0 };
                            double ValueShapeFunction = ComputeShapeFunction(n1, n2, x);
                            c += Jacobi1 * Jacobi2 * ValueShapeFunction;
                        }
                        C1[NOF * i, (slave[j].TArray)[0]] += -c;
                        C1[NOF * i + 1, (slave[j].TArray)[0] + 1] += -c;
                    }
                    else
                    {
                        double Jacobi11 = (slave[j].GetPosition(index) - slave[j - 1].GetPosition(index)) / 2.0;
                        double Jacobi12 = (lagrangefield[i][1] - lagrangefield[i][0]) / (slave[j].GetPosition(index) - slave[j - 1].GetPosition(index));
                        for (int g = 0; g < gausspoint.Length; g++)
                        {
                            double XI = ((lagrangefield[i][1] + lagrangefield[i][0]) - (slave[j].GetPosition(index) + slave[j - 1].GetPosition(index))) / (slave[j].GetPosition(index) - slave[j - 1].GetPosition(index)) +
                                        Jacobi12 * gausspoint[g];
                            double x = ((slave[j].GetPosition(index) + slave[j - 1].GetPosition(index)) / 2.0) + Jacobi11 * XI;
                            double[] n2 = { slave[j].GetPosition(index), 1 };
                            double[] n1 = { slave[j - 1].GetPosition(index), 0 };
                            double ValueShapeFunction = ComputeShapeFunction(n1, n2, x);
                            c += Jacobi11 * Jacobi12 * ValueShapeFunction;
                        }
                        double Jacobi21 = (slave[j + 1].GetPosition(index) - slave[j].GetPosition(index)) / 2.0;
                        double Jacobi22 = (lagrangefield[i][1] - lagrangefield[i][0]) / (slave[j + 1].GetPosition(index) - slave[j].GetPosition(index));
                        for (int g = 0; g < gausspoint.Length; g++)
                        {
                            double XI = ((lagrangefield[i][1] + lagrangefield[i][0]) - (slave[j + 1].GetPosition(index) + slave[j].GetPosition(index))) / (slave[j + 1].GetPosition(index) - slave[j].GetPosition(index)) +
                                        Jacobi22 * gausspoint[g];
                            double x = ((slave[j + 1].GetPosition(index) + slave[j].GetPosition(index)) / 2.0) + Jacobi21 * XI;
                            double[] n1 = { slave[j].GetPosition(index), 1 };
                            double[] n2 = { slave[j + 1].GetPosition(index), 0 };
                            double ValueShapeFunction = ComputeShapeFunction(n1, n2, x);
                            c += Jacobi21 * Jacobi22 * ValueShapeFunction;
                        }
                        C1[NOF * i, (slave[j].TArrayGlobal)[0]] += -c;
                        C1[NOF * i + 1, (slave[j].TArrayGlobal)[0] + 1] += -c;
                    }
                }
            }
            for (int i = 0; i < lagrangefield.Count; i++)
            {
                for (int j = 0; j < master.Count; j++)
                {
                    double c = 0;
                    if (j == 0)
                    {
                        double Jacobi1 = (master[j + 1].GetPosition(index) - master[j].GetPosition(index)) / 2.0;
                        double Jacobi2 = (lagrangefield[i][1] - lagrangefield[i][0]) / (master[j + 1].GetPosition(index) - master[j].GetPosition(index));
                        for (int g = 0; g < gausspoint.Length; g++)
                        {
                            double XI = ((lagrangefield[i][1] + lagrangefield[i][0]) - (master[j + 1].GetPosition(index) + master[j].GetPosition(index))) / (master[j + 1].GetPosition(index) - master[j].GetPosition(index)) +
                                        Jacobi2 * gausspoint[g];
                            double x = ((master[j + 1].GetPosition(index) + master[j].GetPosition(index)) / 2.0) + Jacobi1 * XI;
                            double[] n1 = { master[j].GetPosition(index), 1 };
                            double[] n2 = { master[j + 1].GetPosition(index), 0 };
                            double ValueShapeFunction = ComputeShapeFunction(n1, n2, x);
                            c += Jacobi1 * Jacobi2 * ValueShapeFunction;
                        }
                        C1[NOF * i, (master[j].TArrayGlobal)[0]] += c;
                        C1[NOF * i + 1, (master[j].TArrayGlobal)[0] + 1] += c;
                    }
                    else if (j == master.Count - 1)
                    {
                        double Jacobi1 = (master[j].GetPosition(index) - master[j - 1].GetPosition(index)) / 2.0;
                        double Jacobi2 = (lagrangefield[i][1] - lagrangefield[i][0]) / (master[j].GetPosition(index) - master[j - 1].GetPosition(index));
                        for (int g = 0; g < gausspoint.Length; g++)
                        {
                            double XI = ((lagrangefield[i][1] + lagrangefield[i][0]) - (master[j].GetPosition(index) + master[j - 1].GetPosition(index))) / (master[j].GetPosition(index) - master[j - 1].GetPosition(index)) +
                                        Jacobi2 * gausspoint[g];
                            double x = ((master[j].GetPosition(index) + master[j - 1].GetPosition(index)) / 2.0) + Jacobi1 * XI;
                            double[] n2 = { master[j].GetPosition(index), 1 };
                            double[] n1 = { master[j - 1].GetPosition(index), 0 };
                            double ValueShapeFunction = ComputeShapeFunction(n1, n2, x);
                            c += Jacobi1 * Jacobi2 * ValueShapeFunction;
                        }
                        C1[NOF * i, (master[j].TArrayGlobal)[0]] += c;
                        C1[NOF * i + 1, (master[j].TArrayGlobal)[0] + 1] += c;
                    }
                    else
                    {
                        double Jacobi11 = (master[j].GetPosition(index) - master[j - 1].GetPosition(index)) / 2.0;
                        double Jacobi12 = (lagrangefield[i][1] - lagrangefield[i][0]) / (master[j].GetPosition(index) - master[j - 1].GetPosition(index));
                        for (int g = 0; g < gausspoint.Length; g++)
                        {
                            double XI = ((lagrangefield[i][1] + lagrangefield[i][0]) - (master[j].GetPosition(index) + master[j - 1].GetPosition(index))) / (master[j].GetPosition(index) - master[j - 1].GetPosition(index)) +
                                        Jacobi12 * gausspoint[g];
                            double x = ((master[j].GetPosition(index) + master[j - 1].GetPosition(index)) / 2.0) + Jacobi11 * XI;
                            double[] n2 = { master[j].GetPosition(index), 1 };
                            double[] n1 = { master[j - 1].GetPosition(index), 0 };
                            double ValueShapeFunction = ComputeShapeFunction(n1, n2, x);
                            c += Jacobi11 * Jacobi12 * ValueShapeFunction;
                        }
                        double Jacobi21 = (master[j + 1].GetPosition(index) - master[j].GetPosition(index)) / 2.0;
                        double Jacobi22 = (lagrangefield[i][1] - lagrangefield[i][0]) / (master[j + 1].GetPosition(index) - master[j].GetPosition(index));
                        for (int g = 0; g < gausspoint.Length; g++)
                        {
                            double XI = ((lagrangefield[i][1] + lagrangefield[i][0]) - (master[j + 1].GetPosition(index) + master[j].GetPosition(index))) / (master[j + 1].GetPosition(index) - master[j].GetPosition(index)) +
                                        Jacobi22 * gausspoint[g];
                            double x = ((master[j + 1].GetPosition(index) + master[j].GetPosition(index)) / 2.0) + Jacobi21 * XI;
                            double[] n1 = { master[j].GetPosition(index), 1 };
                            double[] n2 = { master[j + 1].GetPosition(index), 0 };
                            double ValueShapeFunction = ComputeShapeFunction(n1, n2, x);
                            c += Jacobi21 * Jacobi22 * ValueShapeFunction;
                        }
                        C1[NOF * i, (master[j].TArrayGlobal)[0]] += c;
                        C1[NOF * i + 1, (master[j].TArrayGlobal)[0] + 1] += c;
                    }
                }
            }
            Console.WriteLine(C1);
            return C1;
        }


        private double ComputeShapeFunction(double[] n1, double[] n2, double x)
        {
            if (x >= n1[0] && x <= n2[0])
            {
                double a = (n2[1] - n1[1]) / (n2[0] - n1[0]);
                double b = (n2[0] * n1[1] - n1[0] * n2[1]) / (n2[0] - n1[0]);
                return a * x + b;
            }
            else
            {
                return 0;
            }
        }


    }
}
