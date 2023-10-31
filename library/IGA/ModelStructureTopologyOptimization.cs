//using System;
//using System.Threading.Tasks;
//using DEMSoft.NURBS;
//using DEMSoft.EngineeringData;
//using DEMSoft.Drawing;
//using CenterSpace.NMath.Core;

//namespace DEMSoft.IGA
//{
//    /// <summary>
//    /// Structure 2D problem (plane stress, plane strain). Default is plane stress
//    /// </summary>
//    public class ModelStructureTopologyOptimization : AbstractModelStructure
//    {
//        /// <summary>
//        /// Constructor class
//        /// </summary>
//        /// <param name="problem">Define state of stress, plane stress or plane strain</param>
//        public ModelStructureTopologyOptimization(Dimension structureDimension, string pathProject = "temp", string nameProject = "project")
//                     : base(TypeAnalysisModel.TopologyOptimization, structureDimension, pathProject, nameProject)
//        {
//        }

//        /// <summary>
//        /// Solve problem
//        /// </summary>
//        public override void Solve()
//        {
//            if (listPatch.Count == 0)
//                throw new NullReferenceException("Model must be assigned Patch");
//            DoubleMatrix kGlobal = null;
//            DoubleVector uGlobal = null;
//            DoubleVector rGlobal = null;
//            for (int i = 0; i < listPatch.Count; i++)
//            {
//                if (!IsParallelProcesing)
//                    for (int j = 0; j < listPatch[i].GetCountElement(); j++)
//                    {
//                        DoubleMatrix Ke = null;
//                        ElementStructureElastic2D elem = (ElementStructureElastic2D)listPatch[i].GetElement(j);
//                        elem.ComputeStiffnessMatrixElement(out Ke);
//                        elem.SetStiffnessMatrix(Ke);
//                    }
//                else
//                    Parallel.For(0, listPatch[i].GetCountElement(), j =>
//                    {
//                        DoubleMatrix Ke = null;
//                        ElementStructureElastic2D elem = (ElementStructureElastic2D)listPatch[i].GetElement(j);
//                        elem.ComputeStiffnessMatrixElement(out Ke);
//                        elem.SetStiffnessMatrix(Ke);
//                    });
//            }
//            int nelx = ((NURBSSurface)listPatch[0].GetGeometry(0)).Basis.GetKnotVector(0).GetKnotVectorNoMultiplicity().Length - 1;
//            int nely = ((NURBSSurface)listPatch[0].GetGeometry(0)).Basis.GetKnotVector(1).GetKnotVectorNoMultiplicity().Length - 1;
//            //MATERIAL PROPERTIES
//            int ft = 1;
//            double Emin = Math.Pow(10, -9);
//            double rmin = 2.4;
//            double penal = 3;
//            double volfrac = 0.5;
//            //PREPARE FILTER
//            DoubleMatrix iH = new DoubleMatrix((int)(Math.Pow(nelx * nely, 2)), 1);
//            DoubleMatrix jH = new DoubleMatrix(iH.Rows, iH.Cols);
//            DoubleMatrix sH = new DoubleMatrix(iH.Rows, iH.Cols);
//            int k = -1;
//            for (int i1 = 0; i1 < nelx; i1++)
//                for (int j1 = 0; j1 < nely; j1++)
//                {
//                    /////// One patch ///////////////////
//                    int ipn0 = listPatch[0].FindIndexOfElement(i1, j1);
//                    ElementStructureElastic2D curElement = (ElementStructureElastic2D)listPatch[0].GetElement(ipn0);
//                    double[] curCentroidElement = curElement.GetFace().GetCentriodCoordinate();
//                    double e1 = i1 * nely + j1;
//                    for (int i2 = Math.Max((int)(i1 - Math.Round(rmin)), 0); i2 <= Math.Min((int)(i1 + Math.Round(rmin)), nelx - 1); i2++)
//                        for (int j2 = Math.Max((int)(j1 - Math.Round(rmin)), 0); j2 <= Math.Min((int)(j1 + Math.Round(rmin)), nely - 1); j2++)
//                        {
//                            int ipn2 = listPatch[0].FindIndexOfElement(i2, j2);
//                            ElementStructureElastic2D elem = (ElementStructureElastic2D)listPatch[0].GetElement(ipn2);

//                            double[] centroidElement = elem.GetFace().GetCentriodCoordinate();
//                            int e2 = i2 * nely + j2;
//                            k++;
//                            iH[k, 0] = e1;
//                            jH[k, 0] = e2;
//                            sH[k, 0] = Math.Max(0, rmin - Math.Sqrt(Math.Pow(curCentroidElement[0] - centroidElement[0], 2)
//                                         + Math.Pow(curCentroidElement[1] - centroidElement[1], 2)));
//                        }

//                    //for (int i2 = 0; i2 < nelx; i2++)
//                    //    for (int j2 = 0; j2 < nely; j2++)
//                    //    {
//                    //        int ipn2 = listPatch[0].FindIndexOfElement(i2, j2);
//                    //        StructureElastic2DTopologyElement elem = (StructureElastic2DTopologyElement)listPatch[0].GetElement(ipn2);

//                    //        double[] centroidElement = elem.GetFace().GetCentriodCoordinate();
//                    //        double r = Math.Sqrt(Math.Pow(curCentroidElement[0] - centroidElement[0], 2)
//                    //            + Math.Pow(curCentroidElement[1] - centroidElement[1], 2));
//                    //        int e2 = i2 * nely + j2;
//                    //        k++;
//                    //        iH[k, 0] = e1;
//                    //        jH[k, 0] = e2;

//                    //        if (r <= rmin)
//                    //        {
//                    //            sH[k, 0] = Math.Max(0, rmin - r);
//                    //        }
//                    //    }
//                }
//            DoubleMatrix H = MatrixTool.sparse(iH, jH, sH);
//            DoubleMatrix Hs = MatrixTool.sum(H, 1);
//            double Hss = MatrixFunctions.Sum(Hs); //MatrixTool.sum(Hs);
//            //INITIALIZE ITERACTION
//            DoubleMatrix x = new DoubleMatrix(nelx, nely);//(nelx, nely, volfrac);
//                                                      //repmat(volfrac, nely, nelx);
//            DoubleMatrix xPhys = x;
//            DoubleMatrix xnew = new DoubleMatrix(xPhys.Rows, xPhys.Cols);
//            int loop = 0;
//            double change = 1;
//            DoubleMatrix ce = new DoubleMatrix(nelx, nely);
//            DoubleMatrix c = new DoubleMatrix(1, 1);
//            double E0 = new double();

//            //START ITERATION
//            while (change > 0.01)
//            {
//                loop = loop + 1;
//                E0 = ((PatchStructure2D)listPatch[0]).Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();////////////////
//                for (int i = 0; i < listPatch.Count; i++)
//                    for (int j = 0; j < listPatch[i].GetCountElement(); j++)
//                    {
//                        int ii = listPatch[i].GetIPN(j, 0);
//                        int jj = listPatch[i].GetIPN(j, 1);

//                        double Ee = Emin + Math.Pow(xPhys[ii, jj], penal) * (E0 - Emin);
//                        ((ElementStructureElastic2D)GetPatch(i).GetElement(j)).SetCurrentModulus(Ee);
//                    }

//                //ISOGEOMETRIC ANALYSIS
//                kGlobal = null;
//                rGlobal = null;
//                AssemblyStiffnessMatrix(out kGlobal);
//                AssemblyTractionVector(out rGlobal);
//                uGlobal = SolveKuf(kGlobal, rGlobal);
//                SetUGlobal(uGlobal);
//                DoubleVector uu = GetUGlobal(true);
//                //OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
//                int DOF = ((ElementStructureElastic2D)GetPatch(0).GetElement(0)).GetStiffnessMatrix().Cols;
//                DoubleMatrix ce1 = new DoubleMatrix(nelx * nely, DOF);
//                for (int i = 0; i < listPatch.Count; i++)
//                    for (int j = 0; j < listPatch[i].GetCountElement(); j++)
//                    {
//                        DoubleMatrix KE = ((ElementStructureElastic2D)GetPatch(i).GetElement(j)).GetStiffnessMatrix();
//                        DoubleVector uLocalElement = ((ElementStructureElastic2D)GetPatch(i).GetElement(j)).GetDisplacementLocal();
//                        DoubleMatrix a = uLocalElement.ToRowMatrix() * KE;
//                        for (int kk = 0; kk < a.Cols; kk++)
//                            ce1[j, kk] = a[0, kk] * uLocalElement[kk];
//                    }
//                ce = MatrixTool.reshape(MatrixTool.sum(ce1, 1), nelx, nely);
//                //DoubleMatrix c = sum(sum(Mul(Plus(MatPower(xPhys, penal) * (E0 - Emin), Emin), ce), 0), 1);
//                DoubleMatrix dc = (-penal * (E0 - Emin)) * MatrixTool.Mul(MatrixTool.MatPower(xPhys, penal - 1), -ce);//ce);

//                DoubleMatrix dv =new DoubleMatrix(nelx, nely);
//                for (int j = 0; j < listPatch[0].GetCountElement(); j++)
//                {
//                    int ii = listPatch[0].GetIPN(j, 0);
//                    int jj = listPatch[0].GetIPN(j, 1);
//                    dv[ii, jj] = ((ElementStructureElastic2D)GetPatch(0).GetElement(j)).ComputeAreaOfElement();
//                }
//                //for (int i = 0; i < nelx; i++)
//                //    for (int j = 0; j < nely; j++)
//                //    {
//                //        int inp=
//                //        dv[i,j]=
//                //    }
//                //DoubleMatrix dv = ones(nelx, nely);
//                //FILTERING/MODIFICATION OF SENSITIVITIES
//                if (ft == 1)
//                {
//                    DoubleMatrix ma = MatrixTool.Mul(MatrixTool.MatToVec(MatrixTool.max(Math.Pow(10, -3), x.Transpose())), Hs);
//                    DoubleMatrix xx = H * MatrixTool.Mul(MatrixTool.MatToVec(x.Transpose()), MatrixTool.MatToVec(dc.Transpose()));
//                    dc = MatrixTool.VecToMat(MatrixTool.Div(H * MatrixTool.Mul(MatrixTool.MatToVec(x.Transpose()), MatrixTool.MatToVec(dc.Transpose())),
//                                 MatrixTool.Mul(MatrixTool.MatToVec(MatrixTool.max(Math.Pow(10, -3), x.Transpose())), Hs)), dc.Transpose()).Transpose();
//                    //for (int i = 0; i < dc.RowCount; i++)
//                    //    for (int j = 0; j < dc.ColumnCount; j++)
//                    //        dc[i, j] = dc[i, j] / (reshape(Hs, nelx, nely)[i, j] * ma[i, j]);
//                }
//                else if (ft == 2)
//                {
//                    dc = MatrixTool.VecToMat(H * MatrixTool.Div(MatrixTool.MatToVec(dc.Transpose()), Hs), dc.Transpose()).Transpose();
//                    dv = MatrixTool.VecToMat(H * MatrixTool.Div(MatrixTool.MatToVec(dv.Transpose()), Hs), dv.Transpose()).Transpose();
//                }
//                //OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLE AND PHYSICAL DENSITIES
//                double l1 = 0; double l2 = Math.Pow(10, 9); double move = 0.2;
//                while ((l2 - l1) / (l1 + l2) > 0.001)
//                {
//                    double lmid = 0.5 * (l2 + l1);
//                    //DoubleMatrix kk = Mul(x, Sqrt(Div(-dc, MatMul(dv, lmid))));
//                    //DoubleMatrix xn = min(Plus(x, move), Mul(x, Sqrt(Div(-dc, MatMul(dv, lmid)))));
//                    xnew = MatrixTool.max(0, MatrixTool.max(MatrixTool.Plus(x, -move), MatrixTool.min(1, MatrixTool.min(MatrixTool.Plus(x, move), MatrixTool.Mul(x, MatrixTool.Sqrt(MatrixTool.Div(-dc, MatrixTool.MatMul(dv, lmid))))))));
//                    if (ft == 1)
//                        xPhys = xnew;
//                    else if (ft == 2)
//                    {
//                        xPhys = MatrixTool.VecToMat(MatrixTool.Div(H * MatrixTool.MatToVec(xnew.Transpose()), Hs), xPhys.Transpose()).Transpose();
//                    }
//                    if (MatrixTool.sum(MatrixTool.MatToVec(xPhys)) > (volfrac * nelx * nely))
//                        l1 = lmid;
//                    else l2 = lmid;
//                }
//                ///////////////////////////////////////////////
//                //////// Store xPhys into element /////////////
//                ///////////////////////////////////////////////
//                for (int i = 0; i < listPatch.Count; i++)
//                    for (int j = 0; j < listPatch[i].GetCountElement(); j++)
//                    {
//                        int ii = listPatch[i].GetIPN(j, 0);
//                        int jj = listPatch[i].GetIPN(j, 1);
//                        ((ElementStructureElastic2D)GetPatch(i).GetElement(j)).SetDensityFilter(xPhys[ii, jj]);
//                    }
//                //double xphys1 = sum(sum(xPhys, 1));
//                change = MatrixTool.max(MatrixTool.Abs(MatrixTool.MatToVec(xnew) - MatrixTool.MatToVec(x)));
//                x = xnew;
//            }
//            c = MatrixTool.sum(MatrixTool.sum(MatrixTool.Mul(MatrixTool.Plus(MatrixTool.MatPower(xPhys, penal) * (E0 - Emin), Emin), ce), 0), 1);
//        }

//        public void DrawTopologyOptimizationResult(ViewerForm viewer)
//        {
//            if (StructureDimension == Dimension.Plane)
//            {
//                foreach (PatchStructure2D patch in listPatch)
//                {
//                    var surface = patch.GetSurface();
//                    surface.isDrawControlPoint = true;
//                    surface.isDrawControlNet = true;
//                    surface.isColorfulFace = false;
//                    surface.isDrawSurface = false;
//                    surface.resolution1 = 5;
//                    surface.resolution2 = 5;
//                    //surface.Draw(viewer);

//                    double[,][] data = surface.GetDataOnSurface();

//                    double[] knotVectorNoDuplicate1 = surface.Basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
//                    double[] knotVectorNoDuplicate2 = surface.Basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
//                    int numX = surface.resolution1;
//                    int numY = surface.resolution2;

//                    //var colors = new ColorsRange(0, 200);
//                    ColorsRange color = new ColorsRange(ColorType.Gray, 256, 0, 1);
//                    for (int k2 = 0; k2 < knotVectorNoDuplicate2.Length - 1; k2++)
//                    {
//                        for (int k1 = 0; k1 < knotVectorNoDuplicate1.Length - 1; k1++)
//                        {
//                            int idxElement = patch.FindIndexOfElement(k1, k2);
//                            double opacity = ((ElementStructureElastic2D)patch.GetElement(idxElement)).GetDensityFilter();
//                            int sum = 0;

//                            double[] xDataMesh = new double[(numX + 1) * (numY + 1)];
//                            double[] yDataMesh = new double[(numX + 1) * (numY + 1)];
//                            double[] zDataMesh = new double[(numX + 1) * (numY + 1)];
//                            for (int j = 0; j < numY + 1; j++)
//                            {
//                                for (int i = 0; i < numX + 1; i++)
//                                {
//                                    xDataMesh[sum] = data[i + k1 * numX, j + k2 * numY][0];
//                                    yDataMesh[sum] = data[i + k1 * numX, j + k2 * numY][1];
//                                    zDataMesh[sum] = data[i + k1 * numX, j + k2 * numY][2];
//                                    sum++;
//                                }
//                            }

//                            //var mesh = new StructuredGrid(new int[] { numX + 1, numY + 1, 1 }, xDataMesh, yDataMesh, zDataMesh);
//                            //mesh.SetGridOutline(true);
//                            //mesh.SetWidth(widthCurve);
//                            //mesh.SetColor(colorCurve);
//                            //mesh.SetOpacity(opacity);
//                            //form.AddObject3D(mesh);

//                            StructuredGrid2D face = new StructuredGrid2D(new int[] { numY + 1, numX + 1 }, xDataMesh, yDataMesh, zDataMesh);
//                            //face.SetColor(Color.Yellow);
//                            //face.SetRandomColor();
//                            double[] c = color.GetColor(1.0 - opacity);
//                            face.SetColor(c[0], c[1], c[2]);
//                            //face.SetOpacity(opacity);
//                            viewer.AddObject3D(face);
//                        }
//                    }
//                    viewer.SetColormapBarVisible(color, "1.0 - Density Filter");
//                }
//            }
//            else if (StructureDimension == Dimension.Solid)
//            {
//            }
//        }

//        public override void PostProcessing()
//        {
//            throw new NotImplementedException();
//        }
//    }
//}
