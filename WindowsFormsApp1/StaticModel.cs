using DEMSoft.Drawing;
using DEMSoft.Drawing.Geometry;
using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOFEM
{
    internal enum Result
    {
        UX,
        UY,
        USUM,
    }
    internal class StaticModel
    {
        private List<Node> nodes;
        private List<Element> elements;
        private List<Subdomain> subdomains;
        private List<Force> forces;
        private List<Constraint> constraints;
        private List<MortarContactElememt> ContactElememt;

        public int NumberOfField { get; set; }
        //public DenseMatrix K_subdomain { get; private set; }

        // NumberOfField: la so bac tu do cua phan tu trong bai toan can giai
        public StaticModel(int numberOfField) //Phuong thuc khoi tao cho moi Class
        {
            nodes = new List<Node>();
            elements = new List<Element>();
            subdomains = new List<Subdomain>();
            forces = new List<Force>();
            constraints = new List<Constraint>();
            ContactElememt = new List<MortarContactElememt>();
            NumberOfField = numberOfField;
        }
        //public StaticModel(int numberOfField) //Phuong thuc khoi tao cho moi Class
        //{
        //    nodes = new List<Node>();
        //    elements = new List<Element>();
        //    subdomains = new List<Subdomain>();
        //    forces = new List<Force>();
        //    constraints = new List<Constraint>();
        //    NumberOfField = numberOfField;
        //}
        public void InsertNode(Node n)
        {
            nodes.Add(n);
        }
        public Node GetNode(int index)
        {
            return nodes[index];
        }
        public void InsertElement(Element elem)
        {
            elements.Add(elem);
        }
        public Element GetElement(int index)
        {
            return elements[index];
        }

        public void InsertSubdomain(Subdomain subdomain)
        {
            subdomains.Add(subdomain);
        }
        public Subdomain GetSubdomain(int index)
        {
            return subdomains[index];
        }

        public void AddMortarContact(MortarContactElememt mortar)
        {
            ContactElememt.Add(mortar);
        }


        internal void DrawNodes(ViewerForm viewer)
        {
            int countNode = nodes.Count;
            double[] x = new double[countNode];
            double[] y = new double[countNode];
            double[] z = new double[countNode];

            for (int i = 0; i < countNode; i++)
            {
                x[i] = nodes[i].GetPosition(0);
                y[i] = nodes[i].GetPosition(1);
                z[i] = nodes[i].GetPosition(2);
            }

            // Nodes la tap hop diem dung PointSet de ve
            PointSet ps = new PointSet(x, y, z);
            ps.SetColor(Color.Black);
            ps.SetPointSize(3);
            viewer.AddObject3D(ps);
        }

        internal void DrawElements(ViewerForm viewer)
        {
            int countElement = elements.Count;
            Random random = new Random();
            foreach (var item in elements)
            {
                if (item is Truss3DElement)
                {
                    Node n0 = item.nodes[0];
                    Node n1 = item.nodes[1];
                    Line line = new Line(n0.GetPosition(0),
                                         n0.GetPosition(1),
                                         n0.GetPosition(2),
                                         n1.GetPosition(0),
                                         n1.GetPosition(1),
                                         n1.GetPosition(2));
                    line.SetColor(Color.Red);
                    line.SetWidth(2);
                    viewer.AddObject3D(line);
                }
                else if (item is Plane2DElement)
                {
                    Quad quad = new Quad(item.nodes[0].GetPosition(0), item.nodes[0].GetPosition(1),
                                         item.nodes[1].GetPosition(0), item.nodes[1].GetPosition(1),
                                         item.nodes[2].GetPosition(0), item.nodes[2].GetPosition(1),
                                         item.nodes[3].GetPosition(0), item.nodes[3].GetPosition(1));
                    quad.SetColor(random.NextDouble(), random.NextDouble(), random.NextDouble());
                    //quad.SetColor(Color.Black);
                    viewer.AddObject3D(quad);
                }
            }
        }
        internal void AddConstraint(Constraint c1)
        {
            constraints.Add(c1);
        }
        internal void DrawConstrain(ViewerForm viewer)
        {
            foreach (var c in constraints)
            {
                double scale = 0.01;
                Node n = c.GetNode();
                for (int i = 0; i < c.IsFixed.Length; i++)
                {
                    if (c.IsFixed[i] == true)
                    {
                        DEMSoft.Drawing.Point p = null;
                        switch (i)
                        {
                            case 0:
                                p = new DEMSoft.Drawing.Point(n.GetPosition(0) - scale, n.GetPosition(1), n.GetPosition(2));
                                break;
                            case 1:
                                p = new DEMSoft.Drawing.Point(n.GetPosition(0), n.GetPosition(1) - scale, n.GetPosition(2));
                                break;
                            case 2:
                                p = new DEMSoft.Drawing.Point(n.GetPosition(0), n.GetPosition(1), n.GetPosition(2) - scale);
                                break;
                        }
                        p.SetColor(Color.Black);
                        p.SetPointSize(10);
                        viewer.AddObject3D(p);
                    }
                }
            }
        }
        internal void AddForce(Force f1)
        {
            forces.Add(f1);
        }

        internal void DrawForce(ViewerForm viewer)
        {
            double scale = 0.0001;
            foreach (var f in forces)
            {
                Node n0 = f.GetNode();
                double fz = 0;
                if (f.value.Length == 3)
                {
                    fz = f.value[2];
                }
                Line line = new Line(n0.GetPosition(0),
                                     n0.GetPosition(1),
                                     n0.GetPosition(2),
                                     n0.GetPosition(0) + scale * f.value[0],
                                     n0.GetPosition(1) + scale * f.value[1],
                                     n0.GetPosition(2) + scale * fz);

                line.SetColor(Color.Red);
                line.SetWidth(2);
                viewer.AddObject3D(line);
            }
        }
        private int[] BCArray;
        private int[] mappingReduceGlobal;
        internal void PreProcessing()
        {
            Enumerate_Mortar();
            BCArray = AppyBC();
            mappingReduceGlobal = new int[countDOF];

            int c = 0;
            //mortar
            for (int i = 0; i < countDOF; i++)
            {
                if (Array.IndexOf(BCArray, i) != -1)
                {
                    mappingReduceGlobal[i] = -1;
                }
                else
                {
                    mappingReduceGlobal[i] = c++;
                }
            }

            //for (int i = 0; i < countDOF; i++)
            //{
            //    if (Array.IndexOf(BCArray, i) != -1)
            //    {
            //        mappingReduceGlobal[i] = -1;
            //    }
            //    else
            //    {
            //        mappingReduceGlobal[i] = c++;
            //    }
            //}
        }

        private int[] AppyBC()
        {
            List<int> bc = new List<int>();
            foreach (var c in constraints)
            {
                Node n = c.GetNode();
                int[] tArrayNode = n.TArray;
                for (int i = 0; i < c.IsFixed.Length; i++)
                {
                    if (c.IsFixed[i] == true)
                    {
                        bc.Add(tArrayNode[i]);
                    }
                }
            }
            return bc.ToArray();
        }

        private int countDOF;
        private void Enumerate()
        {
            int count = 0;
            for (int i = 0; i < nodes.Count; i++)
            {
                Node n = nodes[i];
                nodes[i].TArray = new int[NumberOfField];
                for (int j = 0; j < NumberOfField; j++)
                {
                    n.TArray[j] = count++;
                }
            }

            countDOF = count;
            for (int i = 0; i < elements.Count; i++)
            {
                Element elem = elements[i];
                int numNodeOnElement = elem.nodes.Length;
                elem.tArray = new int[numNodeOnElement * NumberOfField];
                int c = 0;
                for (int j = 0; j < numNodeOnElement; j++)
                {
                    for (int k = 0; k < NumberOfField; k++)
                    {
                        elem.tArray[c] = elem.nodes[j].TArray[k];
                        c++;
                    }
                }
                foreach (var item in elem.tArray)
                {
                    Console.Write(item + "\t");
                }
                Console.WriteLine();
            }
        }

        private int[] countDOF_subdomain;
        private int numberofLagrangeField;
        private void Enumerate_Mortar()
        {
            int count1 = 0;
            for (int i = 0; i < nodes.Count; i++)
            {
                Node n = nodes[i];
                n.TArray = new int[NumberOfField];
                for (int j = 0; j < NumberOfField; j++)
                {
                    n.TArray[j] = count1++;
                }
            }
            countDOF = count1;
            countDOF_subdomain = new int[subdomains.Count];
            int countGlobal = 0;
            int countSubdomainGlobal = 0;
            for (int s = 0; s < subdomains.Count; s++)
            {
                int count = 0;
                Subdomain subdomain = subdomains[s];
                subdomain.tArrayGlobal = new int[NumberOfField * subdomain.nodes.Count];
                for (int i = 0; i < subdomain.nodes.Count; i++)
                {
                    Node ns = subdomain.nodes[i];
                    ns.TArray = new int[NumberOfField];
                    ns.TArrayGlobal = new int[NumberOfField];

                    for (int j = 0; j < NumberOfField; j++)
                    {
                        ns.TArray[j] = count++;
                        ns.TArrayGlobal[j] = countGlobal++;
                    }
                }
                for (int k = 0; k < NumberOfField * subdomain.nodes.Count; k++)
                {
                    subdomain.tArrayGlobal[k] = countSubdomainGlobal++;
                }

                countDOF_subdomain[s] = count;
                for (int i = 0; i < subdomain.elements.Count; i++)
                {
                    Element elem = subdomain.elements[i];
                    int numNodeOnElement = elem.nodes.Length;
                    elem.tArray = new int[numNodeOnElement * NumberOfField];
                    int c = 0;
                    for (int j = 0; j < numNodeOnElement; j++)
                    {
                        for (int k = 0; k < NumberOfField; k++)
                        {
                            elem.tArray[c] = subdomain.elements[i].nodes[j].TArray[k];
                            c++;
                        }
                    }
                }
            }

            for (int i = 0; i < ContactElememt.Count; i++)
            {
                numberofLagrangeField = ContactElememt[i].lagrangefield.Count;
            }
        }


        private double[] uGlobal;
        internal void Solve()
        {
            DenseMatrix K = AssemblyStiffnessMatrix_MortarMethod();
            //DenseMatrix K = AssemblyStiffnessMatrix();
            DenseVector F = AssemblyForceVector_MORTAR();
            DenseMatrix Kreduce = GetSubMatrix(K, BCArray);
            DenseVector Freduce = GetSubVector(F, BCArray);
            double[] u = Kreduce.Solve(Freduce).ToArray();
            for (int i = 0; i < u.Length; i++)
            {
                if (Math.Abs(u[i]) < 1e-6)
                {
                    Console.WriteLine(0);
                }
                else
                {
                    Console.WriteLine(u[i]);
                }
            }


            uGlobal = new double[countDOF];
            for (int i = 0; i < countDOF; i++)
            {
                if (mappingReduceGlobal[i] != -1)
                {
                    uGlobal[i] = u[mappingReduceGlobal[i]];
                }
                else
                {
                    uGlobal[i] = 0;
                }
            }
        }

        private DenseVector GetSubVector(DenseVector f, int[] bCArray)
        {
            int dof = f.Count - BCArray.Length;
            DenseVector RReduce = new DenseVector(dof);
            int ci = 0;
            for (int i = 0; i < f.Count; i++)
            {
                if (Array.IndexOf(bCArray, i) == -1)
                {
                    RReduce[ci] = f[i];
                    ci++;
                }
            }
            return RReduce;
        }

        private DenseMatrix GetSubMatrix(DenseMatrix k, int[] bCArray)
        {
            int dof = k.ColumnCount - BCArray.Length;
            DenseMatrix KReduce = new DenseMatrix(dof);
            int ci = 0;
            for (int i = 0; i < k.RowCount; i++)
            {
                if (Array.IndexOf(bCArray, i) == -1)
                {
                    int cj = 0;
                    for (int j = 0; j < k.ColumnCount; j++)
                    {
                        if (Array.IndexOf(bCArray, j) == -1)
                        {
                            KReduce[ci, cj] = k[i, j];
                            cj++;
                        }
                    }
                    ci++;
                }
            }
            return KReduce;
        }

        private DenseVector AssemblyForceVector_MORTAR()
        {
            int totalDOF = 0;
            for (int i = 0; i < subdomains.Count; i++)
            {
                totalDOF += countDOF_subdomain[i];
            }

            DenseVector R = new DenseVector(totalDOF + (numberofLagrangeField /*+ constraints.Count*/) * NumberOfField);
            foreach (var f in forces)
            {
                Node n = f.GetNode();
                DenseVector re = f.ComputeVectorForce();
                int[] tArray = n.TArrayGlobal;
                for (int i = 0; i < tArray.Length; i++)
                {
                    R[tArray[i]] += re[i];
                }
            }
            Console.WriteLine(R);
            return R;
        }

        private DenseVector AssemblyForceVector()
        {

            DenseVector R = new DenseVector(countDOF);
            foreach (var f in forces)
            {
                Node n = f.GetNode();
                DenseVector re = f.ComputeVectorForce();
                int[] tArray = n.TArray;
                for (int i = 0; i < tArray.Length; i++)
                {
                    R[tArray[i]] += re[i];
                }
            }
            Console.WriteLine(R);
            return R;
        }

        private DenseMatrix AssemblyStiffnessMatrix()
        {
            DenseMatrix K = new DenseMatrix(countDOF);
            foreach (var elem in elements)
            {
                DenseMatrix Ke = elem.ComputeStiffnessMatrix();
                int[] tArrayElem = elem.tArray;
                for (int i = 0; i < Ke.RowCount; i++)
                {
                    for (int j = 0; j < Ke.ColumnCount; j++)
                    {
                        K[tArrayElem[i], tArrayElem[j]] += Ke[i, j];
                    }
                }
            }

            Console.WriteLine(K);

            return K;
        }
        private DenseMatrix AssemblyStiffnessMatrix_MortarMethod()
        {
            int totalDOF = 0;
            for (int i = 0; i < subdomains.Count; i++)
            {
                totalDOF += countDOF_subdomain[i];
            }
            DenseMatrix K = new DenseMatrix(totalDOF + (numberofLagrangeField /*+ constraints.Count*/) * NumberOfField);

            for (int i = 0; i < subdomains.Count; i++)
            {
                DenseMatrix K_subdomain = subdomains[i].ComputeSubdomainStiffnessMatrix();
                int[] tArraySubdomain = subdomains[i].tArrayGlobal;
                for (int j = 0; j < K_subdomain.RowCount; j++)
                {
                    for (int k = 0; k < K_subdomain.ColumnCount; k++)
                    {
                        K[tArraySubdomain[j], tArraySubdomain[k]] += K_subdomain[j, k];
                    }
                }
            }

            //DenseMatrix MM = new DenseMatrix(constraints.Count*2, totalDOF);
            //DenseMatrix MMT = new DenseMatrix(totalDOF, constraints.Count*2);

            //int C_placeholder = 0;
            //foreach (var item in constraints)
            //{
            //    Node n = item.GetNode();
            //    int[] tArrayNode = n.TArray;
            //    for (int i = 0; i < item.IsFixed.Length; i++)
            //    {
            //        if (item.IsFixed[i] == true)
            //        {
            //            MM[i + C_placeholder, tArrayNode[i]] += 1;
            //            MMT[tArrayNode[i], i + C_placeholder] += 1;
            //            K[i + totalDOF + C_placeholder, tArrayNode[i]] += 1;
            //            K[tArrayNode[i], i + totalDOF + C_placeholder] += 1;
            //        }
            //    }
            //    C_placeholder += 2;
            //}

            //Console.WriteLine(MM);


            for (int m = 0; m < ContactElememt.Count; m++)
            {
                DenseMatrix C = ContactElememt[m].CouplingMatrix();

                for (int i = 0; i < C.RowCount; i++)
                {
                    for (int j = 0; j < C.ColumnCount; j++)
                    {
                        if (C[i, j] != 0)
                        {
                            K[i + totalDOF /*+ constraints.Count * 2*/, j] += C[i, j];
                            K[j, i + totalDOF /*+ constraints.Count * 2*/] += C[i, j];
                        }
                    }
                }
            }
            Console.WriteLine(K);
            return K;
        }
        internal void PostProcessing()
        {
            // Dua gia tri u ve node
            foreach (var item in nodes)
            {
                int[] tArray = item.TArray;
                for (int j = 0; j < tArray.Length; j++)
                {
                    item.setDisplacement(j, uGlobal[tArray[j]]);
                }
            }

            //ADDING_mortar
            foreach (var elems in elements)
            {
                Node[] nodelems = elems.nodes;
                foreach (var n in nodelems)
                {
                    int[] indices = n.TArrayGlobal;
                    for (int j = 0; j < indices.Length; j++)
                    {
                        n.setDisplacement(j, uGlobal[indices[j]]);
                    }
                }
            }
            //foreach (var elems in elements)
            //{
            //    Node[] nodelems = elems.nodes;
            //    foreach (var n in nodelems)
            //    {
            //        int[] indices = n.TArray;
            //        for (int j = 0; j < indices.Length; j++)
            //        {
            //            n.setDisplacement(j, uGlobal[indices[j]]);
            //        }
            //    }
            //}
        }
        internal void DrawDefomation(ViewerForm viewer, double scale)
        {
            int countElement = elements.Count;
            foreach (var item in elements)
            {
                if (item is Truss3DElement)
                {
                    Node n0 = item.nodes[0];
                    Node n1 = item.nodes[1];
                    Line line = new Line(n0.GetPosition(0) + scale * n0.GetDisplacement(0),
                                         n0.GetPosition(1) + scale * n0.GetDisplacement(1),
                                         n0.GetPosition(2) + scale * n0.GetDisplacement(2),
                                         n1.GetPosition(0) + scale * n1.GetDisplacement(0),
                                         n1.GetPosition(1) + scale * n1.GetDisplacement(1),
                                         n1.GetPosition(2) + scale * n1.GetDisplacement(2));
                    line.SetColor(Color.Olive);
                    line.SetDottedLine(true);
                    line.SetWidth(2);
                    viewer.AddObject3D(line);
                }
                else if (item is Plane2DElement)
                {
                    Quad quad = new Quad(item.nodes[0].GetPosition(0) + scale * item.nodes[0].GetDisplacement(0), item.nodes[0].GetPosition(1) + scale * item.nodes[0].GetDisplacement(1),
                                         item.nodes[1].GetPosition(0) + scale * item.nodes[1].GetDisplacement(0), item.nodes[1].GetPosition(1) + scale * item.nodes[1].GetDisplacement(1),
                                         item.nodes[2].GetPosition(0) + scale * item.nodes[2].GetDisplacement(0), item.nodes[2].GetPosition(1) + scale * item.nodes[2].GetDisplacement(1),
                                         item.nodes[3].GetPosition(0) + scale * item.nodes[3].GetDisplacement(0), item.nodes[3].GetPosition(1) + scale * item.nodes[3].GetDisplacement(1));
                    quad.SetColor(Color.Red);
                    quad.SetOpacity(0.3);
                    viewer.AddObject3D(quad);
                }
            }
        }

        internal void DrawResult(ViewerForm viewer, double scale, Result result)
        {
            double minVal = 1e100;
            double maxVal = -1e100;
            Contours2D[] contours = new Contours2D[elements.Count];
            for (int lp = 0; lp < elements.Count(); lp++)
            {
                Element elem = elements[lp];
                if (elem is Plane2DElement)
                {
                    double[,] xData = new double[2, 2];
                    double[,] yData = new double[2, 2];
                    double[,] zData = new double[2, 2];
                    double[,] valData = new double[2, 2];

                    for (int i = 0; i < 2; i++)
                        for (int j = 0; j < 2; j++)
                        {
                            int c = 2 * i + j;
                            if (i == 1)
                                c = 2 * i + 1 - j;
                            xData[i, j] = elem.nodes[c].GetPosition(0) + scale * Math.Round(elem.nodes[c].GetDisplacement(0), 3);
                            yData[i, j] = elem.nodes[c].GetPosition(1) + scale * Math.Round(elem.nodes[c].GetDisplacement(1), 3);
                            //zData[i, j] = item.nodes[2 * i + j].GetPosition(2) + scale * item.nodes[2 * i + j].GetDisplacement(2);
                            switch (result)
                            {
                                case Result.UX:
                                    valData[i, j] = elem.nodes[c].GetDisplacement(0);
                                    break;
                                case Result.UY:
                                    valData[i, j] = elem.nodes[c].GetDisplacement(1);
                                    break;
                                case Result.USUM:
                                    valData[i, j] = Math.Sqrt(Math.Pow(elem.nodes[c].GetDisplacement(0), 2) + Math.Pow(elem.nodes[c].GetDisplacement(1), 2));
                                    break;
                            }
                        }
                    contours[lp] = new Contours2D(xData, yData, zData);
                    contours[lp].SetScalarValue(valData);
                    contours[lp].SetOpacity(0.8);
                    contours[lp].ColorType = ColorType.Jet;
                    var val = contours[lp].GetValue();
                    var miVal = val.Min();
                    var maVal = val.Max();
                    if (miVal < minVal)
                        minVal = miVal;
                    if (maVal > maxVal)
                        maxVal = maVal;
                }
            }

            for (int lp = 0; lp < elements.Count; lp++)
            {
                contours[lp].SetMinMaxValue(minVal, maxVal);
                contours[lp].Update(minVal, maxVal);
                viewer.AddObject3D(contours[lp]);
            }
            viewer.SetColormapBarVisible(contours[0], result.ToString());
        }
    }
}
