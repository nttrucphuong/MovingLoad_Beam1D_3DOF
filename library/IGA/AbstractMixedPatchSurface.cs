//using System.Collections.Generic;
//using DEMSoft.NURBS;

//namespace DEMSoft.IGA
//{
//    /// <summary>
//    /// Abstract class to define NURBS surface to NURBS Patch
//    /// </summary>
//    public abstract class AbstractMixedPatchSurface : AbstractPatchMixedField
//    {
//        //private List<ConnectionInterfaceTwoPatches> connect;

//        /// <summary>
//        /// Constructor class
//        /// </summary>
//        /// <param name="surface">Nurbs surface</param>
//        /// <param name="numberOfFields">number of fields in model</param>
//        public AbstractMixedPatchSurface(NURBSSurface surfaceQ, int numberOfFieldsQ, NURBSSurface surfaceP, int numberOfFieldsP)
//        {
//            //connect = new List<ConnectionInterfaceTwoPatches>();
//            this.geometry = new AbstractParametricGeometry[] { surfaceQ, surfaceP };
//            this.countField = new int[] { numberOfFieldsQ, numberOfFieldsP };
//            this.countDimension = 2;
//            this.INC = new int[2][,];
//            this.IEN = new int[2][,];
//            this.enumeratePatch = new int[2][,];
//            this.enumerateGlobal = new int[2][,];

//            CreateIPN();
//            CreateINC();
//            CreateIEN();
//            //CreateIDArray();
//            //InitialControlPoint();
//        }

//        //private void InitialControlPoint()
//        //{
//        //    int totalField = countField[0] + countField[1];
//        //    for (int k = 0; k < geometry.Length; k++)
//        //    {
//        //        NURBSSurface surface = GetSurface(k);
//        //        ControlPoint[,] cps = surface.ControlPoints;
//        //        for (int i = 0; i < cps.GetLength(0); i++)
//        //            for (int j = 0; j < cps.GetLength(1); j++)
//        //            {
//        //                cps[i, j].SetULocal(new double[4]);
//        //            }
//        //    }
//        //}

//        /// <summary>
//        /// Get NURBS surface
//        /// </summary>
//        /// <returns></returns>
//        public NURBSSurface GetSurface(int idx)
//        {
//            return (NURBSSurface)geometry[idx];
//        }

//        //public void SetConnectionInterface(ConnectionInterfaceTwoPatches con)
//        //{
//        //    connect.Add(con);
//        //}

//        /// <summary>
//        /// Create INC array
//        /// </summary>
//        protected override void CreateINC()
//        {
//            for (int k = 0; k < geometry.Length; k++)
//            {
//                var basis = ((NURBSSurface)geometry[k]).Basis;
//                int n = basis.GetCountBasisFunction(0);
//                int m = basis.GetCountBasisFunction(1);
//                int nen = GetCountLocalBasisFunctions(k);// (p + 1) * (q + 1); // number of local basis functions
//                int nnp = GetCountGlobalBasisFunctions(k);// n * m;//number of global basis functions
//                INC[k] = new int[nnp, 2];
//                int A = 0;

//                for (int j = 0; j < m; j++)
//                {
//                    for (int i = 0; i < n; i++)
//                    {
//                        INC[k][A, 0] = i;
//                        INC[k][A, 1] = j;
//                        A++;
//                    }
//                }
//            }
//        }

//        /// <summary>
//        /// Create IEN array
//        /// </summary>
//        protected override void CreateIEN()
//        {
//            for (int k = 0; k < geometry.Length; k++)
//            {
//                var basis = ((NURBSSurface)geometry[k]).Basis;
//                int n = basis.GetCountBasisFunction(0);
//                int p = basis.GetDegree(0);
//                int m = basis.GetCountBasisFunction(1);
//                int q = basis.GetDegree(1);
//                var knotVectorNoMultiplicity1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
//                var knotVectorNoMultiplicity2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
//                int numelem1 = knotVectorNoMultiplicity1.Length - 1;
//                int numelem2 = knotVectorNoMultiplicity2.Length - 1;
//                int nel = numelem1 * numelem2;// GetNumberOfPatchs();//(n - p) * (m - q);//number of elements
//                int nen = (p + 1) * (q + 1); //GetNumberOfLocalBasisFunctions(); // number of local basis functions

//                IEN[k] = new int[nen, nel];
//                //int e = 0;
//                //int A = 0;
//                for (int ej = 0; ej < numelem2; ej++)
//                {
//                    for (int ei = 0; ei < numelem1; ei++)
//                    {
//                        double mid1 = (knotVectorNoMultiplicity1[ei] + knotVectorNoMultiplicity1[ei + 1]) / 2.0;
//                        double mid2 = (knotVectorNoMultiplicity2[ej] + knotVectorNoMultiplicity2[ej + 1]) / 2.0;
//                        int uspan = basis.FindSpan(mid1, 0);
//                        int vspan = basis.FindSpan(mid2, 1);
//                        int nume = FindIndexOfElement(ei, ej);
//                        int b = 0;
//                        for (int j = 0; j <= q; j++)
//                        {
//                            for (int i = 0; i <= p; i++)
//                            {
//                                int num = FindIndexOfGlobalBasisFunction(k, uspan - p + i, vspan - q + j);
//                                IEN[k][b, nume] = num;
//                                b++;
//                            }
//                        }
//                    }
//                }
//            }
//        }

//        /// <summary>
//        /// Create index of patch on direction
//        /// </summary>
//        protected override void CreateIPN()
//        {
//            var basis = ((NURBSSurface)geometry[0]).Basis;
//            int n = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity().Length - 1;
//            int m = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity().Length - 1;
//            int nnp = n * m;
//            IPN = new int[nnp, 2];
//            int A = 0;
//            for (int j = 0; j < m; j++)
//            {
//                for (int i = 0; i < n; i++)
//                {
//                    IPN[A, 0] = i;
//                    IPN[A, 1] = j;
//                    A++;
//                }
//            }
//        }

//        public int FindIndexOfElementAt(double xi, double eta)
//        {
//            int idxi = -1;
//            int idxj = -1;
//            var kv1 = GetSurface(0).Basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
//            var kv2 = GetSurface(0).Basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
//            for (int i = 0; i < kv1.Length; i++)
//            {
//                if (xi >= kv1[i])
//                    idxi = i;
//                if (xi == kv1[kv1.Length - 1])
//                    idxi = i - 1;
//            }
//            for (int i = 0; i < kv2.Length; i++)
//            {
//                if (eta >= kv2[i])
//                    idxj = i;
//                if (eta == kv2[kv2.Length - 1])
//                    idxj = i - 1;
//            }
//            return FindIndexOfElement(idxi, idxj);
//        }

//        ///// <summary>
//        ///// Get number of patchs
//        ///// </summary>
//        ///// <returns></returns>
//        //public override int GetCountElement()
//        //{
//        //    NURBSSurface surface = GetSurface(0);
//        //    var basis = surface.Basis;
//        //    int n = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity().Length - 1;
//        //    int m = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity().Length - 1;
//        //    return n * m;
//        //}

//        /// <summary>
//        /// Get number of local basis functions
//        /// </summary>
//        /// <returns></returns>
//        public override int GetCountLocalBasisFunctions(int idx)
//        {
//            NURBSSurface surface = GetSurface(idx);
//            int p = surface.Basis.GetDegree(0);
//            int q = surface.Basis.GetDegree(1);
//            return (p + 1) * (q + 1); // number of local basis functions
//        }

//        /// <summary>
//        /// Get number of global basis functions
//        /// </summary>
//        /// <returns></returns>
//        public override int GetCountGlobalBasisFunctions(int idx)
//        {
//            var basis = ((NURBSSurface)geometry[idx]).Basis;
//            int n = basis.GetCountBasisFunction(0);
//            int m = basis.GetCountBasisFunction(1);
//            return n * m; // number of local basis functions
//        }

//        //public override void SetUGlobal(DoubleVector uGlobal)
//        //{
//        //    for (int k = 0; k < geometry.Length; k++)
//        //    {
//        //        var cps = GetSurface(k).ControlPoints;
//        //        for (int i = 0; i < GetCountGlobalBasisFunctions(k); i++)
//        //        {
//        //            var cp = cps[GetINC(k, i, 0), GetINC(k, i, 1)];
//        //            cp.SetUGlobal(uGlobal);
//        //        }
//        //    }
//        //}

//        //public double[] GetULocalAt(double xi, double eta)
//        //{
//        //    double[] disp = new double[countField];
//        //    var surface = GetSurface();
//        //    var cps = surface.ControlPoints;
//        //    var basis = (BivariateNURBSBasisFunction)surface.Basis;
//        //    int p = basis.GetDegree(0);
//        //    int uSpan = basis.GetKnotVector(0).FindSpan(xi, p);
//        //    int q = basis.GetDegree(1);
//        //    int vSpan = basis.GetKnotVector(1).FindSpan(eta, q);
//        //    double[,] Nuv = basis.GetValueBivariateBasisFunctions(xi, eta);

//        //    for (int i = 0; i <= p; i++)
//        //        for (int j = 0; j <= q; j++)
//        //        {
//        //            for (int k = 0; k < countField; k++)
//        //            {
//        //                disp[k] += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetULocal(k);
//        //            }
//        //        }
//        //    return disp;
//        //}

//        public override List<AbstractElement> SelectEndPatchElement(int index)/////////////////////////////////////////////////////////////////
//        {
//            int indexDirection = index / 2;
//            int mod = index % 2;
//            NURBSSurface surface = GetSurface(0);
//            var basis = surface.Basis;
//            int n = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity().Length - 1;
//            int m = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity().Length - 1;
//            List<AbstractElement> selElement = new List<AbstractElement>();
//            for (int i = 0; i < n * m; i++)
//            {
//                if (indexDirection == 0)
//                {
//                    if (mod == 0)
//                    {
//                        if (GetIPN(i, 0) == 0)
//                            selElement.Add(listElement[i]);
//                    }
//                    else
//                    {
//                        if (GetIPN(i, 0) == n - 1)
//                            selElement.Add(listElement[i]);
//                    }
//                }
//                else
//                {
//                    if (mod == 0)
//                    {
//                        if (GetIPN(i, 1) == 0)
//                            selElement.Add(listElement[i]);
//                    }
//                    else
//                    {
//                        if (GetIPN(i, 1) == m - 1)
//                            selElement.Add(listElement[i]);
//                    }
//                }
//            }
//            return selElement;
//        }

//        public ControlPoint[] SelectControlPoints(int idxPatch, int index, int indexCoordinate)
//        {
//            var cps = GetSurface(idxPatch).ControlPoints;
//            ControlPoint[] selCps = null;
//            int n = 0;
//            switch (indexCoordinate)
//            {
//                case 0:
//                    n = GetSurface(idxPatch).Basis.GetCountBasisFunction(1);
//                    selCps = new ControlPoint[n];
//                    for (int i = 0; i < n; i++)
//                        selCps[i] = cps[index, i];
//                    break;
//                case 1:
//                    n = GetSurface(idxPatch).Basis.GetCountBasisFunction(0);
//                    selCps = new ControlPoint[n];
//                    for (int i = 0; i < n; i++)
//                        selCps[i] = cps[i, index];
//                    break;
//            }
//            return selCps;
//        }

//        /// <summary>
//        /// Select serial of patchs from mesh by index on each direction
//        /// </summary>
//        /// <param name="index">index = 0-first column or row of patch on mesh, 1-last column or row of patch on mesh</param>
//        /// <param name="indexDirection">index of normal direction coordinate</param>
//        /// <returns></returns>
//        public ControlPoint[] SelectEndPatchControlPoints(int idxPatch, int index, int indexDirection)
//        {
//            ControlPoint[] selCps = null;
//            int n = 0;
//            var cps = GetSurface(idxPatch).ControlPoints;
//            switch (indexDirection)
//            {
//                case 0://normal xi-coordinate
//                    n = GetSurface(idxPatch).Basis.GetCountBasisFunction(1);
//                    selCps = new ControlPoint[n];

//                    for (int i = 0; i < n; i++)
//                    {
//                        if (index == 0)
//                            selCps[i] = cps[0, i];
//                        else
//                            selCps[i] = cps[GetSurface(idxPatch).Basis.GetCountBasisFunction(0) - 1, i];
//                    }
//                    break;
//                case 1://normal eta-coordinate
//                    n = GetSurface(idxPatch).Basis.GetCountBasisFunction(0);
//                    selCps = new ControlPoint[n];
//                    for (int i = 0; i < n; i++)
//                    {
//                        if (index == 0)
//                            selCps[i] = cps[i, 0];
//                        else
//                            selCps[i] = cps[i, GetSurface(idxPatch).Basis.GetCountBasisFunction(1) - 1];
//                    }
//                    break;
//            }
//            return selCps;
//        }

//        public override int EnumerateInGlobalMultiPatch(int countDof)
//        {
//            for (int k = 0; k < geometry.Length; k++)
//            {
//                var cps = ((NURBSSurface)geometry[k]).ControlPoints;
//                enumerateGlobal[k] = new int[countField[k], cps.GetLength(0) * cps.GetLength(1)];
//                for (int i = 0; i < enumerateGlobal[k].GetLength(1); i++)
//                {
//                    int psi = INC[k][i, 0];
//                    int eta = INC[k][i, 1];
//                    var cp = cps[psi, eta];
//                    if (cp.GetCoupleControlPoint() == null)
//                    {
//                        for (int j = 0; j < countField[k]; j++)
//                        {
//                            enumerateGlobal[k][j, i] = countDof;
//                            countDof++;
//                        }
//                    }
//                    else
//                    {
//                        for (int j = 0; j < countField[k]; j++)
//                        {
//                            enumerateGlobal[k][j, i] = -2;//is coupled
//                        }
//                    }

//                }
//                //countDofU = countDof;
//                for (int j = 0; j < cps.GetLength(1); j++)
//                    for (int i = 0; i < cps.GetLength(0); i++)
//                    {
//                        cps[i, j].SetNumberOfFields(countField[k]);
//                        int[] tArrayGlobal = new int[countField[k]];
//                        for (int kk = 0; kk < countField[k]; kk++)
//                            tArrayGlobal[kk] = enumerateGlobal[k][kk, FindIndexOfGlobalBasisFunction(k, i, j)];
//                        cps[i, j].SetTArrayGlobal(tArrayGlobal);
//                    }
//            }
//            return countDof;
//        }

//        public override int[] EnumerateInPatch()
//        {
//            int id = 0;
//            int[] enumerate = new int[geometry.Length];
//            for (int k = 0; k < geometry.Length; k++)
//            {
//                var cps = ((NURBSSurface)geometry[k]).ControlPoints;
//                enumeratePatch[k] = new int[countField[k], cps.GetLength(0) * cps.GetLength(1)];

//                for (int i = 0; i < enumeratePatch[k].GetLength(1); i++)
//                {
//                    for (int j = 0; j < countField[k]; j++)
//                    {
//                        int psi = INC[k][i, 0];
//                        int eta = INC[k][i, 1];

//                        enumeratePatch[k][j, i] = id;
//                        id++;
//                    }
//                }
//                enumerate[k] = id;


//                /////////////////////////////////////////////////////////////
//                /// Distribute tArray patch into control point //////////////
//                /////////////////////////////////////////////////////////////

//                for (int j = 0; j < cps.GetLength(1); j++)
//                    for (int i = 0; i < cps.GetLength(0); i++)
//                    {
//                        cps[i, j].SetNumberOfFields(countField[k]);
//                        int[] tArray = new int[countField[k]];
//                        for (int kk = 0; kk < countField[k]; kk++)
//                            tArray[kk] = GetEnumerateInPatch(k, kk, FindIndexOfGlobalBasisFunction(k, i, j));
//                        cps[i, j].SetTArray(tArray);
//                    }
//            }
//            countDOF = id;
//            //for (int i = 0; i < GetCountElement(); i++)
//            //{
//            //    //((AbstractStructure2DElement)GetElement(i)).GetFace().EnumerateOnFace();
//            //}
//            return enumerate;
//        }
//    }
//}
