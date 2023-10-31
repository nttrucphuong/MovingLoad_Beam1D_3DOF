using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
    /// <summary>
    /// Abstract of NURBS patch mesh
    /// </summary>
    public abstract class AbstractPatch : AbstractMeshPart
    {
        /// <summary>
        /// Nurbs surface 
        /// [0] ux - velocity x, uy - velocity y
        /// [1] theta - volume change, p - pressure
        /// </summary>
        protected AbstractParametricGeometry[] geometry;
        private int id = -1;

        public Material Material
        { get; set; }

        /// <summary>
        /// List of patch
        /// </summary>
        protected List<AbstractElement> listElement;

        //protected List<ConnectionInterfaceTwoPatches> connect;

        /// <summary>
        /// Number of fields of problem 
        /// [0] ux - velocity x, uy - velocity y
        /// [1] theta - volume change, p - pressure
        /// </summary>
        protected int[] countField;

        /// <summary>
        /// A corresponding NURBS coordinate 
        /// [0] ux - velocity x, uy - velocity y
        /// [1] theta - volume change, p - pressure
        /// [][number of global basis functions, index coordinate]
        /// </summary>
        protected int[][,] INC;

        /// <summary>
        /// A corresponding global basis function number
        /// [0] ux - velocity x, uy - velocity y
        /// [1] theta - volume change, p - pressure
        /// [][number of local basis functions, number of elements]
        /// </summary>
        protected int[][,] IEN;

        /// <summary>
        /// A corresponding equation number
        /// [0] ux - velocity x, uy - velocity y
        /// [1] theta - volume change, p - pressure
        /// [][number Of Fields, number of global basis functions]
        /// </summary>
        protected int[][,] enumeratePatch;

        protected int[][,] enumerateGlobal;

        /// <summary>
        /// A corresponding patch number [index of patch, index coordinate]
        /// </summary>
        protected int[,] IPN;

        /// <summary>
        /// number of DOFs of problem
        /// </summary>
        protected int countDOF;

        /// <summary>
        /// number of dimension (1-curve, 2-surface, 3-solid)
        /// </summary>
        protected int countDimension;
        public AbstractPatch(int countDimension)
        {
            this.countDimension = countDimension;
        }
        public List<TypeFields> DisableField { get; set; }

        /// <summary>
        /// Get surface 
        /// </summary>
        /// <param name="idx">
        /// [0] ux - velocity x, uy - velocity y -----
        /// [1] theta - volume change, p - pressure</param>
        /// <returns></returns>
        public AbstractParametricGeometry GetGeometry(int idx)
        {
            return geometry[idx];
        }

        /// <summary>
        /// Get number of DOFs in patch
        /// </summary>
        /// <returns></returns>
        public int GetCountDOF()
        {
            return countDOF;
        }

        public List<TypeFields> TypeOfFieldsInMultifield { get; set; }

        /// <summary>
        /// Get number of fields
        /// </summary>
        /// <returns></returns>
        public int GetCountField(int idx)
        {
            return countField[idx];
        }

        /// <summary>
        /// Get number of dimension
        /// </summary>
        /// <returns></returns>
        public int GetCountDimension()
        {
            return countDimension;
        }

        /// <summary>
        /// Create INC array
        /// </summary>
        protected abstract void CreateINC();

        /// <summary>
        /// Create INC array
        /// </summary>
        protected abstract void CreateIEN();

        /// <summary>
        /// Create IPN array (index of element)
        /// </summary>
        protected abstract void CreateIPN();

        /// <summary>
        /// Create ID array
        /// </summary>
        public abstract int[] EnumerateInPatch();//Enumerate DOF

        public abstract int EnumerateInGlobalMultiPatch(int countDof);

        /// <summary>
        /// Get corresponding NURBS coordinate with a global basis function number and a parametric direction number
        /// </summary>
        /// <param name="idx">
        /// [0] ux - velocity x, uy - velocity y -----
        /// [1] theta - volume change, p - pressure</param>
        /// <param name="indexGlobalBasisFunction">Index of global basis function in 1D array</param>
        /// <param name="indexCoordinate">a parametric direction number</param>
        /// <returns></returns>
        public int GetINC(int idx, int indexGlobalBasisFunction, int indexCoordinate)
        {
            return INC[idx][indexGlobalBasisFunction, indexCoordinate];
        }

        /// <summary>
        /// Find index of global basis function corresponding global coordinate
        /// </summary>
        /// <param name="index">
        /// [0] ux - velocity x, uy - velocity y -----
        /// [1] theta - volume change, p - pressure</param>
        /// <param name="idx">index of global basis functions in global coordinate</param>
        /// <returns></returns>
        public int FindIndexOfGlobalBasisFunction(int index, params int[] idx)
        {
            int indexout = -1;
            ////////////////////////////////////////////////////////
            ///////// serial ////////////////////////////////////
            ///////////////////////////////////////////////////////
            if (!AbstractModel.IsParallelProcesing)
            {
                for (int i = 0; i < GetCountGlobalBasisFunctions(index); i++)
                {
                    bool flag = true;
                    for (int k = 0; k < idx.Length; k++)
                        if (INC[index][i, k] != idx[k])
                        {
                            flag = false;
                            break;
                        }
                    if (flag)
                    {
                        indexout = i;
                        break;
                    }
                }
            }
            else
            {
                ////////////////////////////////////////////////////////
                ///////// parallel ////////////////////////////////////
                ///////////////////////////////////////////////////////
                //Parallel.For(0, GetCountGlobalBasisFunctions(index), (i, loopState) =>
                //{
                //  bool flag = true;
                //  for (int k = 0; k < idx.Length; k++)
                //    if (INC[index][i, k] != idx[k])
                //    {
                //      flag = false;
                //      break;
                //    }
                //  if (flag)
                //  {
                //    indexout = i;
                //    loopState.Stop();
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
                             var max = GetCountGlobalBasisFunctions(index) * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = GetCountGlobalBasisFunctions(index) * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 bool flag = true;
                                 for (int k = 0; k < idx.Length; k++)
                                     if (INC[index][i, k] != idx[k])
                                     {
                                         flag = false;
                                         break;
                                     }
                                 if (flag)
                                 {
                                     indexout = i;
                                     break;
                                 }
                             }
                         });
                }
                Task.WaitAll(tasks);
            }
            return indexout;
        }

        /// <summary>
        /// Get corresponding global basis function number with a local basis function number and an element number
        /// </summary>
        /// <param name="idx">
        /// [0] ux - velocity x, uy - velocity y -----
        /// [1] theta - volume change, p - pressure</param>
        /// <param name="indexElement">index of element</param>
        /// <param name="indexLocalBasisFunction">a local basis function number in element</param>
        /// <returns></returns>
        public int GetIEN(int idx, int indexElement, int indexLocalBasisFunction)
        {
            return IEN[idx][indexLocalBasisFunction, indexElement];
        }

        /// <summary>
        /// Get corresponding equation number with a global basis function number and a dof number.
        /// Return -1 then this DOF was constrainted.
        /// </summary>
        /// <param name="idx">
        /// [0] ux - velocity x, uy - velocity y -----
        /// [1] theta - volume change, p - pressure</param>
        /// <param name="indexFields">index of fields</param>
        /// <param name="indexGlobalBasisFunction">index of a global basis function</param>
        /// <returns></returns>
        public int GetEnumerateInPatch(int idx, int indexFields, int indexGlobalBasisFunction)
        {
            return enumeratePatch[idx][indexFields, indexGlobalBasisFunction];
        }

        public int GetEnumerateInGlobal(int idx, int indexFields, int indexGlobalBasisFunction)
        {
            return enumerateGlobal[idx][indexFields, indexGlobalBasisFunction];
        }

        public void SetEnumerateGlobal(int idx, int idDof, int indexFields, int indexGlobalBasisFunction)
        {
            enumerateGlobal[idx][indexFields, indexGlobalBasisFunction] = idDof;
        }

        /// <summary>
        /// Get corresponding patch(element) number with an patch(element) number and direction number
        /// </summary>
        /// <param name="indexElement">id of element</param>
        /// <param name="indexCoordinate">a parametric direction number</param>
        /// <returns></returns>
        public int GetIPN(int indexElement, int indexCoordinate)
        {
            return IPN[indexElement, indexCoordinate];
        }

        /// <summary>
        /// Find index of global basis function corresponding global coordinate
        /// </summary>
        /// <param name="idx">index of global basis functions in global coordinate</param>
        /// <returns></returns>
        public int FindIndexOfElement(params int[] idx)
        {
            int indexOfElement = -1;
            if (!AbstractModel.IsParallelProcesing)
            {
                for (int i = 0; i < CalculateNumberOfElements(); i++)
                {
                    bool flag = true;
                    for (int k = 0; k < idx.Length; k++)
                        if (IPN[i, k] != idx[k])
                        {
                            flag = false;
                            break;
                        }
                    if (flag)
                    {
                        indexOfElement = i;
                        break;
                    }
                }
            }
            else
            {
                var degreeOfParallelism = AbstractModel.NumberOfCPUs;
                var tasks = new Task[degreeOfParallelism];
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                         () =>
                         {
                             var max = CalculateNumberOfElements() * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = CalculateNumberOfElements() * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 bool flag = true;
                                 for (int k = 0; k < idx.Length; k++)
                                     if (IPN[i, k] != idx[k])
                                     {
                                         flag = false;
                                         break;
                                     }
                                 if (flag)
                                 {
                                     indexOfElement = i;
                                     break;
                                 }
                             }
                         });
                }
                Task.WaitAll(tasks);
            }
            return indexOfElement;
        }

        /// <summary>
        /// Get corresponding equation number with a dof number, a local basis function number and  an element number
        /// </summary>
        /// <param name="idx">
        /// [0] ux - velocity x, uy - velocity y -----
        /// [1] theta - volume change, p - pressure</param>
        /// <param name="indexLocalBasisFunction">index of local basis function</param>
        /// <param name="indexElement">index of element</param>
        /// <param name="indexDOF">index of DOF</param>
        /// <returns></returns>
        public int GetLM(int idx, int indexLocalBasisFunction, int indexElement, int indexDOF)
        {
            return enumeratePatch[idx][indexDOF, GetIEN(idx, indexElement, indexLocalBasisFunction)];
        }

        /// <summary>
        /// Select control points in region
        /// </summary>
        /// <param name="idx">
        /// [0] ux - velocity x, uy - velocity y -----
        /// [1] theta - volume change, p - pressure</param>
        ///<param name="loc">location is a region object</param>
        /// <returns></returns>
        public ControlPoint[] SelectControlPoints(int idx, IRegion loc)
        {
            return geometry[idx].SelectControlPoints(loc);
        }

        /// <summary>
        /// Select all of control points on patch
        /// </summary>
        /// <param name="idx"></param>
        /// <returns></returns>
        public ControlPoint[] SelectAllControlPoints(int idx)
        {
            return geometry[idx].SelectAllControlPoints();
        }

        /// <summary>
        /// Get number of patchs
        /// </summary>
        /// <returns></returns>
        public abstract int CalculateNumberOfElements();
        public int CountElements()
        {
            return listElement.Count;
        }

        /// <summary>
        /// Get number of local basis functions
        /// </summary>
        /// <returns></returns>
        public abstract int GetCountLocalBasisFunctions(int idx);

        /// <summary>
        /// Get number of global basis functions
        /// </summary>
        /// <returns></returns>
        public abstract int GetCountGlobalBasisFunctions(int idx);

        //public abstract void SetUGlobal(DoubleVector uGlobal);

        /// <summary>
        /// Select serial of patchs from mesh by index on each direction
        /// </summary>
        /// <param name="index">[0], [1]: u-direction, [2], [3]: v-direction, [4], [5]: w-direction</param>
        /// <returns></returns>
        public abstract List<AbstractElement> SelectEndPatchElement(int index);

        /// <summary>
        /// Select serial of patchs from mesh by location
        /// </summary>
        /// <param name="idx">
        /// [0] ux - velocity x, uy - velocity y -----
        /// [1] theta - volume change, p - pressure</param>
        /// <param name="xLower">x-lower</param>
        /// <param name="xUpper">x-upper</param>
        /// <param name="yLower">y-lower</param>
        /// <param name="yUpper">y-upper</param>
        /// <param name="zLower">z-lower</param>
        /// <param name="zUpper">z-upper</param>
        /// <returns></returns>
        public List<AbstractElement> SelectElement(double xLower, double xUpper, double yLower, double yUpper, double zLower, double zUpper)
        {
            List<AbstractElement> selElement = new List<AbstractElement>();
            int nel = listElement.Count;
            //for (int i = 0; i < listElement.Count; i++)
            //{
            //  if (listElement[i].IsInsideRegion(new RegionBox(xLower, xUpper, yLower, yUpper, zLower, zUpper)))
            //    selElement.Add(listElement[i]);
            //}

            if (!AbstractModel.IsParallelProcesing)
            {
                for (int i = 0; i < nel; i++)//Loop through elements
                {
                    AbstractElement elem = listElement[i];
                    if (elem.IsInsideRegion(new RegionBox(xLower, xUpper, yLower, yUpper, zLower, zUpper)))
                        selElement.Add(elem);
                }
            }
            else
            {
                var degreeOfParallelism = Environment.ProcessorCount;
                object monitor = new object();
                var tasks = new Task[degreeOfParallelism];

                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                        () =>
                        {
                            var max = nel * (taskNumberCopy + 1) / degreeOfParallelism;
                            for (int i = nel * taskNumberCopy / degreeOfParallelism; i < max; i++)
                            {
                                AbstractElement elem = listElement[i];
                                if (elem.IsInsideRegion(new RegionBox(xLower, xUpper, yLower, yUpper, zLower, zUpper)))
                                {
                                    lock (monitor)
                                    {
                                        selElement.Add(elem);
                                    }
                                }
                            }
                        });
                }
                Task.WaitAll(tasks);
            }
            return selElement;
        }


        public int GetCountGeometry()
        {
            return geometry.Length;
        }

        public AbstractElement GetElement(int idx)
        {
            return listElement[idx];
        }

        //public void SetConnectionInterface(ConnectionInterfaceTwoPatches con)
        //{
        //    connect.Add(con);
        //}

        /// <summary>
        /// ID of patch
        /// </summary>
        public int ID
        {
            get { return id; }
            set { id = value; }
        }

        public void SetMaterial(Material ma)
        {
            this.Material = ma;
            //for (int i = 0; i < listElement.Count; i++)
            //    listElement[i].Material = ma;
        }

        public void AssignMaterialToElements()
        {
            int countElement = listElement.Count;
            if (!AbstractModel.IsParallelProcesing)
            {
                for (int i = 0; i < countElement; i++)
                    listElement[i].Material = Material;
            }
            else
            {
                var degreeOfParallelism = AbstractModel.NumberOfCPUs;
                var tasks = new Task[degreeOfParallelism];
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                         () =>
                         {
                             var max = countElement * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = countElement * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 listElement[i].Material = Material;
                             }
                         });
                }
                Task.WaitAll(tasks);
            }
        }
        public GaussPoints[] GaussPointOnBezierElement
        { get; set; }
        public void ComputeDataDrawMaterialDistribution(ref List<double> x, ref List<double> y, ref List<double> z, ref List<double> val)
        {
            for (int i = 0; i < listElement.Count; i++)
            {
                if (this is AbstractPatch2D)
                {
                    for (int ii = 0; ii < ((AbstractElement2D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); ii++)
                        for (int jj = 0; jj < ((AbstractElement2D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); jj++)
                        {
                            GaussPoints gps = ((AbstractElement2D)listElement[i]).GetGaussPoint(ii, jj);
                            double[] pointAt = ((AbstractPatch2D)this).GetSurface().PointAt(gps.location[0], gps.location[1]);
                            x.Add(pointAt[0]);
                            y.Add(pointAt[1]);
                            z.Add(0);
                            val.Add((double)gps.GetValue(DataInGausspoint.MaterialPropertyValue));//gps.materialPropertyValue);
                        }
                }
                else if (this is AbstractPatch3D)
                {
                    for (int ii = 0; ii < ((AbstractElement3D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); ii++)
                        for (int jj = 0; jj < ((AbstractElement3D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); jj++)
                            for (int kk = 0; kk < ((AbstractElement3D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); kk++)
                            {
                                GaussPoints gps = ((AbstractElement3D)listElement[i]).GetGaussPoint(ii, jj, kk);
                                double[] pointAt = ((AbstractPatch3D)this).GetVolume().PointAt(gps.location[0], gps.location[1], gps.location[2]);
                                x.Add(pointAt[0]);
                                y.Add(pointAt[1]);
                                z.Add(pointAt[2]);
                                val.Add((double)gps.GetValue(DataInGausspoint.MaterialPropertyValue));
                            }
                }
            }
        }

        public void ComputeDataDrawValueDistribution(ref List<double> x, ref List<double> y, ref List<double> z, ref List<double> val)
        {
            for (int i = 0; i < listElement.Count; i++)
            {
                if (this is AbstractPatch2D)
                {
                    for (int ii = 0; ii < ((AbstractElement2D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); ii++)
                        for (int jj = 0; jj < ((AbstractElement2D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); jj++)
                        {
                            GaussPoints gps = ((AbstractElement2D)listElement[i]).GetGaussPoint(ii, jj);
                            double[] pointAt = ((AbstractPatch2D)this).GetSurface().PointAt(gps.location[0], gps.location[1]);
                            x.Add(pointAt[0]);
                            y.Add(pointAt[1]);
                            z.Add(0);
                            val.Add((double)gps.GetValue(DataInGausspoint.DrawValue));
                        }
                }
                else if (this is AbstractPatch3D)
                {
                    for (int ii = 0; ii < ((AbstractElement3D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); ii++)
                        for (int jj = 0; jj < ((AbstractElement3D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); jj++)
                            for (int kk = 0; kk < ((AbstractElement3D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); kk++)
                            {
                                GaussPoints gps = ((AbstractElement3D)listElement[i]).GetGaussPoint(ii, jj, kk);
                                double[] pointAt = ((AbstractPatch3D)this).GetVolume().PointAt(gps.location[0], gps.location[1], gps.location[2]);
                                x.Add(pointAt[0]);
                                y.Add(pointAt[1]);
                                z.Add(pointAt[2]);
                                val.Add((double)gps.GetValue(DataInGausspoint.DrawValue));
                            }
                }
            }
        }
        public void ComputeResultToGaussPoint(Result result, ref List<double> x, ref List<double> y, ref List<double> z, ref List<double> val)
        {
            for (int i = 0; i < listElement.Count; i++)
            {
                if (this is AbstractPatch2D)
                {
                    for (int ii = 0; ii < ((AbstractElement2D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); ii++)
                        for (int jj = 0; jj < ((AbstractElement2D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); jj++)
                        {
                            GaussPoints gps = ((AbstractElement2D)listElement[i]).GetGaussPoint(ii, jj);
                            double[] pointAt = ((AbstractPatch2D)this).GetSurface().PointAt(gps.location[0], gps.location[1]);
                            x.Add(pointAt[0]);
                            y.Add(pointAt[1]);
                            z.Add(0);
                            val.Add(GetApproximateAt(result, gps.location[0], gps.location[1]));
                        }
                }
                else if (this is AbstractPatch3D)
                {
                    for (int ii = 0; ii < ((AbstractElement3D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); ii++)
                        for (int jj = 0; jj < ((AbstractElement3D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); jj++)
                            for (int kk = 0; kk < ((AbstractElement3D)listElement[i]).GetNumberOfGaussPointOnEachDirection(); kk++)
                            {
                                GaussPoints gps = ((AbstractElement3D)listElement[i]).GetGaussPoint(ii, jj, kk);
                                double[] pointAt = ((AbstractPatch3D)this).GetVolume().PointAt(gps.location[0], gps.location[1], gps.location[2]);
                                x.Add(pointAt[0]);
                                y.Add(pointAt[1]);
                                z.Add(pointAt[2]);
                                val.Add(GetApproximateAt(result, gps.location[0], gps.location[1], gps.location[2]));
                            }
                }
            }
        }
        public abstract double[] CalculateResult(Result re, int resolution);
        public abstract void ComputeExtractionOperator();
        //public abstract void ComputeStiffnessMatrixPatch(ref SparseMatrixBuilder<double> kGlobal);
        //public abstract void ComputeStiffnessMatrixPatch(ref ConcurrentDictionary<IntPair, double> kGlobal);
        public abstract void ComputeStiffnessMatrixPatch(ref DoubleCsrSparseMatrix kGlobal);
        public abstract void ComputeStiffnessMatrixPatch(ref DoubleMatrix kGlobal);
        public abstract void ComputeGeometricStiffnessMatrixPatch(double[,] stressTensor, ref DoubleMatrix kGGlobal);
        public abstract void ComputeGeometricStiffnessMatrixPatch(double[,] stressTensor, ref DoubleCsrSparseMatrix kGGlobal);
        public abstract void ComputeInternalForcePatch(ref DoubleVector residualGlobal);
        public virtual double ComputeSharpCrackSurface()
        { return 0; }
        public abstract void ComputeMassMatrixPatch(ref DoubleCsrSparseMatrix mGlobal);
        public abstract void ComputeMassMatrixPatch(ref DoubleMatrix mGlobal);
        public void Initialize()
        {
            int nFields = 0;
            for (int i = 0; i < TypeOfFieldsInMultifield.Count; i++)
            {
                switch (TypeOfFieldsInMultifield[i])
                {
                    case TypeFields.Structural:
                        nFields += countDimension;
                        if (this is PatchStructureBeam)
                        {
                            nFields = 3;
                        }
                        break;
                    case TypeFields.Thermal:
                        nFields += 1;
                        break;
                    case TypeFields.PhaseField:
                        nFields += 1;
                        break;
                    case TypeFields.Electric:
                        nFields += 1;
                        break;
                }
            }
            countField = new int[] { nFields };

            CreateIPN();
            CreateINC();
            CreateIEN();

            int nel = CalculateNumberOfElements();// (n - p) * (m - q);//number of elements
            listElement = new List<AbstractElement>();
            for (int i = 0; i < nel; i++)
            {
                AbstractElement elem = null;
                switch (AbstractModel.TypeModel)
                {
                    case TypeModelProblem.Structural:
                        if (this is PatchStructurePlate)
                        {
                            elem = new ElementStructureElasticPlate((AbstractPatch2D)this, i);
                        }
                        else if (this is PatchStructureBeam)
                        {
                            elem = new ElementStructureElasticBeam((AbstractPatch1D)this, i);
                        }
                        else
                        {
                            if (countDimension == 2)
                            {
                                switch (Material.TypeMaterialStructure)
                                {
                                    case TypeMaterialStructure.Elasticity:
                                        elem = new ElementStructureElastic2D((AbstractPatch2D)this, i);
                                        break;
                                    case TypeMaterialStructure.Plasticity:
                                        elem = new ElementStructurePlastic2D((AbstractPatch2D)this, i);
                                        break;
                                    case TypeMaterialStructure.Hyperelasticity:
                                        elem = new ElementStructureHyperElastic2D((AbstractPatch2D)this, i);
                                        break;
                                }
                            }
                            else if (countDimension == 3)
                            {
                                switch (Material.TypeMaterialStructure)
                                {
                                    case TypeMaterialStructure.Elasticity:
                                        elem = new ElementStructureElastic3D((AbstractPatch3D)this, i);
                                        break;
                                    case TypeMaterialStructure.Plasticity:
                                        elem = new ElementStructurePlastic3D((AbstractPatch3D)this, i);
                                        break;
                                    case TypeMaterialStructure.Hyperelasticity:
                                        elem = new ElementStructureHyperElastic3D((AbstractPatch3D)this, i);
                                        break;
                                }
                            }
                        }
                        break;
                    case TypeModelProblem.Piezoelectric:
                        if (countDimension == 2)
                            elem = new ElementPiezoelectric2D((AbstractPatch2D)this, i);
                        else if (countDimension == 3)
                            elem = new ElementPiezoelectric3D((AbstractPatch3D)this, i);
                        break;
                    case TypeModelProblem.Thermal:
                        ///////////// Tien
                        if (countDimension == 2)
                            elem = new ElementThermal2D((AbstractPatch2D)this, i);
                        else if (countDimension == 3)
                            elem = new ElementThermal3D((AbstractPatch3D)this, i);
                        break;
                    ////////////
                    case TypeModelProblem.StructuralThermal:
                        if (countDimension == 2)
                            elem = new ElementThermoelastic2D((AbstractPatch2D)this, i);
                        else if (countDimension == 3)
                            elem = new ElementThermoelastic3D((AbstractPatch3D)this, i);
                        break;
                    case TypeModelProblem.PhaseField:
                        if (countDimension == 2)
                        {
                            if (TypeOfFieldsInMultifield.IndexOf(TypeFields.PhaseField) != -1)
                                elem = new ElementStructurePhaseField2D((AbstractPatch2D)this, i);
                            else
                                elem = new ElementStructureElastic2D((AbstractPatch2D)this, i);
                        }
                        else if (countDimension == 3)
                        {
                            if (TypeOfFieldsInMultifield.IndexOf(TypeFields.PhaseField) != -1)
                                elem = new ElementStructurePhaseField3D((AbstractPatch3D)this, i);
                            else
                                elem = new ElementStructureElastic3D((AbstractPatch3D)this, i);
                        }
                        break;
                    default:
                        break;
                }
                listElement.Add(elem);
            }
        }
        public abstract double GetApproximateAt(Result result, params double[] xi);
        public abstract double GetApproximateOnGlobalAt(Result result, params double[] x);
    }
}
