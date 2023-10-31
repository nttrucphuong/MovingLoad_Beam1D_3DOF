using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using DEMSoft.NURBS;
using DEMSoft.Function;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
    /// <summary>
    /// Abstract class to define NURBS surface to NURBS Patch
    /// </summary>
    public abstract class AbstractPatch1D : AbstractPatchOneField
    {
        public double Thickness
        { get; set; }
        public double[][,] ExtractionOperator
        { get; set; }
        public override void ComputeExtractionOperator()
        {
            ExtractionOperator = GetCurve().ComputeExtractionOperator();//GetCurve().ComputeExtractionOperator();
        }

        /// <summary>
        /// Constructor class
        /// </summary>
        /// <param name="curve">Nurbs surface</param>
        /// <param name="numberOfFields">number of fields in model</param>
        public AbstractPatch1D(Abstract1DParametricGeometry curve)
           : base(1)
        {
            geometry = new Abstract1DParametricGeometry[] { curve };
        }

        /// <summary>
        /// Get NURBS curve
        /// </summary>
        /// <returns></returns>
        public NURBSCurve GetCurve()
        {
            return (NURBSCurve)GetGeometry();
        }

        /// <summary>
        /// Create INC array
        /// </summary>
        protected override void CreateINC()
        {
            var basis = GetCurve().Basis;
            int n = basis.GetCountBasisFunction();
            INC = new int[1][,];
            INC[0] = new int[n, 1];
            int A = 0;

            for (int i = 0; i < n; i++)
            {
                INC[0][A, 0] = i;
                A++;
            }
        }

        /// <summary>
        /// Create INC array
        /// </summary>
        protected override void CreateIEN()
        {
            var basis = GetCurve().Basis;
            int n = basis.GetCountBasisFunction();
            int p = basis.GetDegree();
            var knotVectorNoMultiplicity1 = basis.GetKnotVector().GetKnotVectorNoMultiplicity();
            int numelem1 = knotVectorNoMultiplicity1.Length - 1;
            int nel = numelem1;// GetNumberOfPatchs();//(n - p) * (m - q);//number of elements
            int nen = (p + 1); //GetNumberOfLocalBasisFunctions(); // number of local basis functions
            IEN = new int[1][,];
            IEN[0] = new int[nen, nel];
            //int e = 0;
            //int A = 0;
            for (int ei = 0; ei < numelem1; ei++)
            {
                double mid1 = (knotVectorNoMultiplicity1[ei] + knotVectorNoMultiplicity1[ei + 1]) / 2.0;
                int uspan = basis.FindSpan(mid1, 0);
                int nume = FindIndexOfElement(ei);
                int b = 0;
                for (int i = 0; i <= p; i++)
                {
                    int num = FindIndexOfGlobalBasisFunction(uspan - p + i);
                    IEN[0][b, nume] = num;
                    b++;
                }
            }
        }

        /// <summary>
        /// Create index of patch on direction
        /// </summary>
        protected override void CreateIPN()
        {
            var basis = GetCurve().Basis;
            int n = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity().Length - 1;
            int nnp = n;
            IPN = new int[nnp, 1];
            int A = 0;
            for (int i = 0; i < n; i++)
            {
                IPN[A, 0] = i;
                A++;
            }
        }

        public override int FindIndexOfElementAt(params double[] xi)
        {
            int idxi = -1;
            var kv1 = GetCurve().Basis.GetKnotVector().GetKnotVectorNoMultiplicity();
            for (int i = 0; i < kv1.Length; i++)
            {
                if (xi[0] >= kv1[i])
                    idxi = i;
                if (xi[0] == kv1[kv1.Length - 1])
                    idxi = i - 1;
            }
            return FindIndexOfElement(idxi);
        }

        /// <summary>
        /// Get number of patchs
        /// </summary>
        /// <returns></returns>
        public override int CalculateNumberOfElements()
        {
            var basis = GetCurve().Basis;
            int n = basis.GetKnotVector().GetKnotVectorNoMultiplicity().Length - 1;
            return n;
        }

        public override List<AbstractElement> SelectEndPatchElement(int index)/////////////////////////////////////////////////////////////////
        {
            int mod = index % 2;
            NURBSCurve curve = GetCurve();
            var basis = curve.Basis;
            int n = basis.GetKnotVector().GetKnotVectorNoMultiplicity().Length - 1;
            List<AbstractElement> selPatch = new List<AbstractElement>();
            for (int i = 0; i < n; i++)
            {
                if (mod == 0)
                {
                    if (GetIPN(i, 0) == 0)
                        selPatch.Add(listElement[i]);
                }
                else
                {
                    if (GetIPN(i, 0) == n - 1)
                        selPatch.Add(listElement[i]);
                }
            }
            return selPatch;
        }

        public override ControlPoint[] GetAllControlPoints()
        {
            return GetCurve().ControlPoints;
        }

        /// <summary>
        /// Create ID array
        /// </summary>
        public override int[] EnumerateInPatch()//Enumerate DOF
        {
            int d = GetCountField();
            var cps = GetCurve().ControlPoints;
            enumeratePatch = new int[1][,];
            enumeratePatch[0] = new int[d, cps.Length];
            int id = 0;
            for (int i = 0; i < enumeratePatch[0].GetLength(1); i++)
            {
                for (int j = 0; j < d; j++)
                {
                    enumeratePatch[0][j, i] = id;
                    id++;
                }
            }
            countDOF = id;

            /////////////////////////////////////////////////////////////
            /// Distribute tArray patch into control point //////////////
            /////////////////////////////////////////////////////////////
            for (int i = 0; i < cps.Length; i++)
            {
                cps[i].SetDimension(1);
                cps[i].SetNumberOfFields(d);

                int[] tArray = new int[d];
                for (int k = 0; k < d; k++)
                    tArray[k] = GetIDInPatch(k, FindIndexOfGlobalBasisFunction(i));
                cps[i].SetTArray(tArray);

                cps[i].Initialize();
            }
            return new int[] { id };
        }

        public override int EnumerateInGlobalMultiPatch(int countDof)
        {
            int d = GetCountField();
            var cps = GetCurve().ControlPoints;
            enumerateGlobal = new int[1][,];
            enumerateGlobal[0] = new int[d, cps.Length];
            for (int i = 0; i < enumerateGlobal[0].GetLength(1); i++)
            {
                int psi = GetINC(i, 0);
                var cp = cps[psi];

                for (int j = 0; j < d; j++)
                {
                    if (cp.GetCoupleControlPoint() == null)
                    {
                        enumerateGlobal[0][j, i] = countDof;
                        countDof++;
                    }
                    else
                    {
                        int dofUnCoupling = cp.ListDOFUnCoupling.IndexOf(j);
                        if (dofUnCoupling == -1)
                            enumerateGlobal[0][j, i] = -2;//is coupled
                        else
                        {
                            enumerateGlobal[0][j, i] = countDof;
                            countDof++;
                        }
                    }
                }
                ///////////////////////////////////////////////////////////

            }

            /////////////////////////////////////////////////////////////
            /// Distribute tArray patch into control point //////////////
            /////////////////////////////////////////////////////////////
            for (int i = 0; i < cps.Length; i++)
            {
                cps[i].SetDimension(1);
                cps[i].SetNumberOfFields(d);

                int[] tArrayGlobal = new int[d];
                for (int k = 0; k < d; k++)
                    tArrayGlobal[k] = enumerateGlobal[0][k, FindIndexOfGlobalBasisFunction(i)];
                cps[i].SetTArrayGlobal(tArrayGlobal);
            }


            return countDof;
        }

        public override int GetCountLocalBasisFunctions(int idx = 0)
        {
            NURBSCurve curve = GetCurve();
            int p = curve.Basis.GetDegree();
            return (p + 1); // number of local basis functions
        }

        public override int GetCountGlobalBasisFunctions(int idx = 0)
        {
            var basis = GetCurve().Basis;
            int n = basis.GetCountBasisFunction();
            return n; // number of local basis functions
        }

        public double ComputeLength()
        {
            return ((Abstract1DParametricGeometry)geometry[0]).ComputeLength();
        }

        /// <summary>
        ///    NOT FIX
        /// </summary>
        /// <param name="xi"></param>
        /// <returns></returns>
        public double GetMaterialPropertyValueApproximationAt(double xi)
        {
            double disp = 0;
            var curve = GetCurve();
            var cps = curve.ControlPoints;
            var basis = (NURBSBasisFunction)curve.Basis;
            int p = basis.GetDegree(0);
            int uSpan = basis.GetKnotVector(0).FindSpan(xi, p);
            double[] Nuv = basis.GetValueBasisFunctions(xi);

            for (int i = 0; i <= p; i++)
                disp += Nuv[i] * cps[uSpan - p + i].MaterialPropertyValue;
            return disp;
        }

        public double[] ComputeMaterialPropertyValue(int resolution)
        {
            ///////////////////////////////////////////////////////////////////////////////
            ////// Parallel ///////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////
            var curve = GetCurve();
            var data = curve.GetParametricOnCurve(resolution);
            double[] val = new double[data.Length];
            if (!AbstractModel.IsParallelProcesing)
            {
                for (int i = 0; i < data.Length; i++)
                    val[i] = GetMaterialPropertyValueApproximationAt(data[i]);
            }
            else
            {
                //Parallel.For(0, data[1].Length, j =>
                //{
                //  for (int i = 0; i < data[0].Length; i++)
                //  {
                //    val[i, j] = GetMaterialPropertyValueApproximationAt(data[0][i], data[1][j]);
                //  }
                //});
                var degreeOfParallelism = AbstractModel.NumberOfCPUs;
                var tasks = new Task[degreeOfParallelism];
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                         () =>
                         {
                             var max = data.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = data.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 val[i] = GetMaterialPropertyValueApproximationAt(data[i]);
                             }
                         });
                }
                Task.WaitAll(tasks);
            }
            return val;
        }

        public override double GetApproximateAt(Result result, params double[] xi)
        {
            if (this is PatchStructureBeam)
            {
                double gzz = ((PatchStructureBeam)this).Getzz;
                FunctionRToR fz = ((PatchStructureBeam)this).KinematicsFunction;
                var typebeam = ((PatchStructureBeam)this).TypeBeam;
            }
            double disp = 0;
            var curve = GetCurve();
            var cps = curve.ControlPoints;
            var basis = (NURBSBasisFunction)curve.Basis;
            int p = basis.GetDegree(0);
            int uSpan = basis.GetKnotVector(0).FindSpan(xi[0], p);
            //int q = basis.GetDegree(1);
            //int vSpan = basis.GetKnotVector(1).FindSpan(xi[1], q);
            double[] Nuv = basis.GetValueBasisFunctions(xi[0]);
            double sxx = 0, syy = 0, szz = 0, sxy = 0, sxz = 0, syz = 0;
            for (int i = 0; i <= p; i++)
            {
                switch (result)
                {
                    case Result.USUM:
                        sxx += Nuv[i] * cps[uSpan - p + i].GetResult(Result.UX);
                        syy += Nuv[i] * cps[uSpan - p + i].GetResult(Result.UY);
                        break;
                    case Result.SIGMAEQV:
                        sxx += Nuv[i] * cps[uSpan - p + i].GetResult(Result.SIGMAXX);
                        syy += Nuv[i] * cps[uSpan - p + i].GetResult(Result.SIGMAYY);
                        sxy += Nuv[i] * cps[uSpan - p + i].GetResult(Result.SIGMAXY);
                        szz += Nuv[i] * cps[uSpan - p + i].GetResult(Result.SIGMAZZ);
                        break;
                    case Result.EPSILONEQV:
                        sxx += Nuv[i] * cps[uSpan - p + i].GetResult(Result.EPSILONXX);
                        syy += Nuv[i] * cps[uSpan - p + i].GetResult(Result.EPSILONYY);
                        szz += Nuv[i] * cps[uSpan - p + i].GetResult(Result.EPSILONZZ);
                        sxy += Nuv[i] * cps[uSpan - p + i].GetResult(Result.EPSILONXY);
                        break;
                    default:
                        disp += Nuv[i] * cps[uSpan - p + i].GetResult(result);
                        break;
                }
                //disp += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(result);
            }

            switch (result)
            {
                case Result.USUM:
                    disp = Math.Sqrt(sxx * sxx + syy * syy);
                    break;
                //case Result.SIGMAXX:
                //case Result.SIGMAYY:
                //case Result.SIGMAXY:
                //case Result.SIGMAZZ:
                case Result.SIGMAEQV:
                    disp = Math.Sqrt(1.0 / 2.0 * (Math.Pow(sxx - syy, 2) + Math.Pow(syy - szz, 2) + Math.Pow(sxx - szz, 2) + 6.0 * Math.Pow(sxy, 2)));
                    //disp = ((PatchStructure2D)this).StressAt(result, xi);
                    break;
                case Result.EPSILONEQV:
                    double xx = 2.0 / 3.0 * sxx - 1.0 / 3.0 * syy - 1.0 / 3.0 * szz;
                    double yy = -1.0 / 3.0 * sxx + 2.0 / 3.0 * syy - 1.0 / 3.0 * szz;
                    double zz = -1.0 / 3.0 * sxx - 1.0 / 3.0 * syy + 2.0 / 3.0 * szz;
                    disp = 2.0 / 3.0 * Math.Sqrt(3.0 / 2.0 * (Math.Pow(xx, 2) + Math.Pow(yy, 2) + Math.Pow(zz, 2)) + 3.0 / 4.0 * (Math.Pow(2 * sxy, 2)));
                    break;
            }
            return disp;
        }

        public override double GetApproximateOnGlobalAt(Result result, params double[] x)
        {
            double[] param = ((NURBSCurve)geometry[0]).Projection(x[0], x[1], 0);
            return GetApproximateAt(result, param[0], param[1]);
        }

        public override double[] CalculateResult(Result re, int resolution)
        {
            var curve = GetCurve();
            var data = curve.GetParametricOnCurve(resolution);
            double[] val = new double[data.Length * data.Length];
            if (!AbstractModel.IsParallelProcesing)
            {
                int count = 0;
                for (int i = 0; i < data.Length; i++)
                    val[count++] = GetApproximateAt(re, data[i], data[i]);
            }
            else
            {
                //Parallel.For(0, data[0].Length, i =>
                //{
                //  for (int j = 0; j < data[1].Length; j++)
                //    val[i * data[1].Length + j] = GetApproximateAt(re, data[0][i], data[1][j]);
                //});
                var degreeOfParallelism = AbstractModel.NumberOfCPUs;
                var tasks = new Task[degreeOfParallelism];
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                         () =>
                         {
                             var max = data.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = data.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 val[i * data.Length] = GetApproximateAt(re, data[i], data[i]);

                             }
                         });
                }
                Task.WaitAll(tasks);
            }
            return val;
        }

        /// <summary>
        /// Get deformation surface to draw result with scale factor
        /// </summary>
        /// <param name="scale">scale factor</param>
        /// <returns></returns>
        public NURBSCurve GetDeformationCurve(double scale)
        {
            NURBSCurve curve = GetCurve();
            NURBSCurve curveCopy = (NURBSCurve)curve.Clone();
            curveCopy.isDrawControlNet = curve.isDrawControlNet;
            curveCopy.isDrawControlPoint = curve.isDrawControlPoint;
            curveCopy.isDrawCurve = curve.isDrawCurve;
            curveCopy.isDrawKnot = curve.isDrawKnot;
            //curveCopy.isColorfulFace = curve.isColorfulFace;
            //curveCopy.colorSurface = curve.colorSurface;
            curveCopy.colorCurve = curve.colorCurve;
            curveCopy.colorControlPoint = curve.colorControlPoint;
            curveCopy.colorControlNet = curve.colorControlNet;

            curveCopy.resolution = curve.resolution;
            //curveCopy.resolution2 = curve.resolution2;
            curveCopy.opacity = curve.opacity;
            var cps = GetCurve().ControlPoints;
            var cpsCopy = curveCopy.ControlPoints;
            if (!AbstractModel.IsParallelProcesing)
            {
                for (int i = 0; i < cpsCopy.Length; i++)
                {
                    cpsCopy[i].SetCoordinate(1, cpsCopy[i].GetCoordinate(1) + scale * cps[i].GetResult(Result.UY));
                }
            }
            else
            {
                //Parallel.For(0, cpsCopy.GetLength(0), i =>
                //{
                //  for (int j = 0; j < cpsCopy.GetLength(1); j++)
                //  {
                //      cpsCopy[i, j].SetCoordinate(0, cpsCopy[i, j].GetCoordinate(0) + scale * cps[i, j].GetResult(Result.UX));
                //    cpsCopy[i, j].SetCoordinate(1, cpsCopy[i, j].GetCoordinate(1) + scale * cps[i, j].GetResult(Result.UY));
                //    cpsCopy[i, j].SetCoordinate(2, cpsCopy[i, j].GetCoordinate(2) + scale * cps[i, j].GetResult(Result.UZ));
                //  }
                //});
                var degreeOfParallelism = AbstractModel.NumberOfCPUs;
                var tasks = new Task[degreeOfParallelism];
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                         () =>
                         {
                             var max = cpsCopy.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = cpsCopy.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 cpsCopy[i].SetCoordinate(1, cpsCopy[i].GetCoordinate(1) + scale * cps[i].GetResult(Result.UY));
                             }
                         });
                }
                Task.WaitAll(tasks);
            }
            return curveCopy;
        }
    }
}
