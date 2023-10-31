using System;
using System.Collections.Generic;
using System.Linq;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using DEMSoft.Drawing;
using DEMSoft.Function;
using System.IO;

using System.Drawing;
using System.Threading.Tasks;
using CenterSpace.NMath.Core;
using System.Diagnostics;
using static DEMSoft.IGA.MonitorDataHistorism;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
    #region enumTools
    /// <summary>
    /// Enumerate of dimension of model
    /// </summary>
    public enum Dimension
    {
        /// <summary>
        /// In two-dimension
        /// </summary>
        Plane = 2,
        /// <summary>
        /// In three-dimension
        /// </summary>
        Solid = 3,
        Plate = 4,// Bo sung/
        Shell,
        Beam = 1//
    };
    /// <summary>
    /// Type of model
    /// </summary>
    public enum TypeModelProblem
    {
        /// <summary>
        /// Structural analysis.
        /// DOF: UX, UY, UZ
        /// Force label: FX, FY, FZ
        /// Reaction solution: Force
        /// </summary>
        Structural,
        /// <summary>
        /// Thermal analysis
        /// DOF: TEMP
        /// Force label: HEAT
        /// Reaction solution: Heat flow
        /// </summary>
        Thermal,
        /// <summary>
        /// Couple-field analysis piezoelectric
        /// DOF: UX, UY, UZ, VOLT
        /// Force label: FX, FY, FZ, CHRG
        /// Reaction solution: Force, electric charge (negative)
        /// </summary>
        Piezoelectric,
        /// <summary>
        /// Couple-field analysis structural-thermal
        /// DOF: UX, UY, UZ, TEMP
        /// Force label: FX, FY, FZ, HEAT
        /// Reaction solution: Force, Heat flow (negative)
        /// </summary>
        StructuralThermal,
        /// <summary>
        /// Couple-field analysis phase-field
        /// DOF: UX, UY, UZ, PHASEFIELD
        /// Force label: FX, FY, FZ
        /// Reaction solution: Force
        /// </summary>
        PhaseField,
    };
    public enum TypeFields
    {
        /// <summary>
        /// Structural analysis.
        /// DOF: UX, UY, UZ
        /// Force label: FX, FY, FZ
        /// Reaction solution: Force
        /// </summary>
        Structural,
        /// <summary>
        /// Thermal analysis
        /// DOF: TEMP
        /// Force label: HEAT
        /// Reaction solution: Heat flow (negative)
        /// </summary>
        Thermal,
        /// <summary>
        /// Electric analysis
        /// DOF: VOLT
        /// Force label: CHRG
        /// Reaction solution: electric charge (negative)
        /// </summary>
        Electric,
        /// <summary>
        /// Phase-field analysis 
        /// DOF:PHASEFIELD
        /// Force label: none
        /// Reaction solution: none
        /// </summary>
        PhaseField,
    };
    /// <summary>
    /// Enumerate type of analysis model
    /// </summary>
    public enum TypeAnalysisModel
    {
        Static, Transient, Modal, TopologyOptimization, Buckling
    };
    #endregion

    /// <summary>
    /// Abstract model class
    /// </summary>
    public abstract class AbstractModel
    {
        protected double[,] bcMatrix;
        /// <summary>
        /// List of patches in model
        /// </summary>
        protected List<AbstractPatch> listPatch;
        /// <summary>
        /// List of loads in model
        /// </summary>
        protected List<AbstractLoad> loads;
        /// <summary>
        /// List of materials in model
        /// </summary>
        protected List<Material> materials;
        /// <summary>
        /// List of connection between two patches, connection 1:1
        /// </summary>
        protected List<ICoupling> listConnection;
        private List<int[]> listIndexConnectionPatch;
        private List<List<ControlPoint>> listControlPointCoupling;
        private List<List<int>> listIndexPatchCoupling;
        internal SparseMatrixBuilder<double> TSparseBuilder;//matrix transform
        internal DoubleMatrix TDense;//matrix transform
        internal DoubleMatrix TDenseTranspose;//matrix transform
        internal DoubleCsrSparseMatrix TSparse;
        internal DoubleCsrSparseMatrix TSparseTranspose;
        internal DoubleVector g;//inhomogenous
                                //internal DoubleMatrix TReduce;//matrix transform
                                //internal DoubleVector gReduce;//inhomogenous
        internal List<int> tArrayUncouple;
        internal List<int> tArrayIndependent;
        internal List<int> tArrayDependent;
        //internal bool isUseTransformMatrix;
        /// <summary>
        /// Store results of field variable of analysis and results of component (stress, strain, total displacement,...)
        /// </summary>
        public List<Result> listComputeResult { get; set; }
        /// <summary>
        /// Store constrain
        /// </summary>
        internal List<AbstractConstraintValue> listConstraint;
        /// <summary>
        /// Store initial constrain
        /// </summary>
        internal List<AbstractConstraintValue> listInitialConstraint;
        /// <summary>
        /// Number of degree of freedoms in model
        /// </summary>
        internal static int countDOF = -1;
        /// <summary>
        /// Store enum Dimension of structure in model, 2D: Dimension.Plane, 3D: Dimension.Solid
        /// </summary>
        internal Dimension StructureDimension
        { get; set; }
        public static int NumberOfCPUs
        { get; set; }
        ///// <summary>
        ///// Determine model processing is parallel or serial
        ///// </summary>
        public static bool IsParallelProcesing
        { get; set; }
        /// <summary>
        /// Compute Ke in parallel or not
        /// </summary>
        public static bool IsComputeKeParallel
        { get; set; }
        /// <summary>
        /// Use bezier extraction or not
        /// </summary>
        public bool IsBezierExtraction
        { get; set; }
        /// <summary>
        /// Store list of monitor data at points
        /// </summary>
        protected MonitorDataHistorism monitorDataHistorism;

        public void TurnGraphRealTimeOn()
        {
            monitorDataHistorism.IsGraphRealTime = true;
        }
        /// <summary>
        /// Matrix or Vector in Sparse or Dense
        /// </summary>
        public bool IsSparseData
        { get; set; }

        /// <summary>
        /// Save result in step by step or in the end of solving
        /// </summary>
        public bool IsSaveStepByStep
        { get; set; }

        public bool IsMatlabSolver
        { get; set; }

        public MLApp.MLApp matlab
        { get; set; }

        public Convergence convergenceOption;
        /// <summary>
        /// Draw convergence error results in iteration problem
        /// </summary>
        public bool IsDrawConvergencePlot { get; set; }
        public int[] IndexRowColPlot { get; set; }
        /// <summary>
        /// Draw converged result in monitor data store
        /// </summary>
        public bool IsDrawMonitorPlot { get; set; }

        public int NumberOfStepStorageNonlinear { get; set; }
        /// <summary>
        /// Store type of analysis of model, such as Static, Transient, Modal, TopologyOptimization
        /// </summary>
        internal static TypeAnalysisModel TypeAnalysisModel
        { get; set; }
        /// <summary>
        /// Store type of model Structural, Thermal, Piezoelectric, StructuralThermal, PhaseField
        /// </summary>
        internal static TypeModelProblem TypeModel
        { get; set; }
        /// <summary>
        /// Store variables in the model in time step
        /// </summary>
        //internal List<DoubleVector> DisplacementTime;
        internal string FileNameDataTime
        { get; set; }
        internal string FileNameTime
        { get; set; }
        /// <summary>
        /// Store time step
        /// </summary>
        internal List<double> Time;
        /// <summary>
        /// Time step
        /// </summary>
        internal double[] deltaT;
        /// <summary>
        /// Check process is finished or not
        /// </summary>
        //internal bool IsFinished
        //{ get; set; }
        public bool IsRunFromInitial
        { get; set; }
        //internal string[] filenameCheckpoint;
        /// <summary>
        /// If this is true, program will apply constraint>0, else apply zero
        /// </summary>
        internal bool IsApplyBoundaryValue;
        /// <summary>
        /// One field displacement:1, coupled-fields: 2...
        /// </summary>
        internal List<TypeFields> typeOfFieldsInMultifield;
        /// <summary>
        /// Write log file storing all parameters of problem
        /// </summary>
        public bool IsWriteLogFile = true;
        /// <summary>
        /// Stop watch to monitor time computing whole process
        /// </summary>
        internal Stopwatch stopWatchWholeModel;
        /// <summary>
        /// Contain log file storing all parameters of problem
        /// </summary>
        internal List<string> logFile;
        ///// <summary>
        ///// Use for multi-load function
        ///// </summary>
        //protected int numberOfStepLoad;
        /// <summary>
        /// Use for iteration problems in 1 Load step 
        /// </summary>
        internal int[] numberOfSubstepLoad;
        /// <summary>
        /// Save number of Step load (&lt;numberOfStepLoad)
        /// </summary>
        internal int[] numberOfStepSave;
        /// <summary>
        /// Store all t-array, including constraint and unconstraint
        /// [number of field][count DOF]
        /// Number of field of One field problem is [0], coupled-field problem: [0]: mechanical, [1]: thermal or [1] electric
        /// </summary>
        internal int[][] tArray;

        /// <summary>
        /// Store unconstraint t-array
        /// [number of field][count DOF unconstraint]
        /// Number of field of One field problem is [0], coupled-field problem: [0]: mechanical, [1]: thermal or [1] electric
        /// </summary>
        internal int[][] tArrayUnconstraint;

        /// <summary>
        /// Store constraint t-array
        /// [number of field][count DOF constraint]
        /// Number of field of One field problem is [0], coupled-field problem: [0]: mechanical, [1]: thermal or [1] electric
        /// </summary>
        internal int[][] tArrayConstraint;
        public static TypeSchemeNonlinearSolver TypeSchemeNonlinearSolver { get; set; }
        public TypeNonlinearSolver type;
        internal DoubleCsrSparseMatrix kGlobalSparse, mGlobalSparse;
        public DoubleVector residualForce;
        //public DataInGausspoint Variable;
        public ColorType colorType = ColorType.Jet;
        public int numberOfDecimal = 12;
        public bool IsMergeControlpoints1_1Constraint { get; set; }
        public void DrawResultAtGaussPoints(ViewerForm viewer, Result re, Orientation orient = Orientation.Horizontal, ColorType colorType = ColorType.Jet)
        {
            List<double> x = new List<double>();
            List<double> y = new List<double>();
            List<double> z = new List<double>();
            List<double> val = new List<double>();
            PointSetColor ps;
            if (StructureDimension == Dimension.Beam)
            {
                foreach (AbstractPatch1D patch in listPatch)
                {
                    var curve = patch.GetCurve();
                    curve.isDrawControlPoint = false;
                    curve.isDrawControlNet = false;
                    curve.isDrawKnot = false;
                    curve.colorCurve = Color.Black;
                    curve.widthCurve = 1;
                    curve.resolution = 5;
                    //curve.opacity = 0.2;
                    curve.Draw(viewer);
                    patch.ComputeResultToGaussPoint(re, ref x, ref y, ref z, ref val);
                }
            }
            else if (StructureDimension == Dimension.Plane)
            {
                foreach (AbstractPatch2D patch in listPatch)
                {
                    var surface = patch.GetSurface();
                    surface.isDrawControlPoint = false;
                    surface.isDrawControlNet = false;
                    surface.isColorfulFace = false;
                    surface.isDrawSurface = false;
                    surface.isDrawCurve = false;
                    surface.isDrawKnot = false;
                    surface.colorCurve = Color.Black;
                    surface.widthCurve = 1;
                    surface.resolution1 = 5;
                    surface.resolution2 = 5;
                    //surface.opacity = 0.2;
                    surface.Draw(viewer);
                    patch.ComputeResultToGaussPoint(re, ref x, ref y, ref z, ref val);
                }
            }
            else if (StructureDimension == Dimension.Solid)
            {
                foreach (AbstractPatch3D patch in listPatch)
                {
                    var vol = patch.GetVolume();
                    vol.isDrawGrid = true;
                    vol.isDrawVolume = false;
                    vol.isDrawControlPoint = false;
                    vol.isDrawControlNet = false;
                    vol.isColorfulFace = false;
                    vol.isDrawCurve = false;
                    vol.isDrawKnot = false;
                    vol.colorGrid = Color.Black;
                    vol.widthGrid = 1;
                    vol.resolution1 = 5;
                    vol.resolution2 = 5;
                    vol.resolution3 = 5;
                    //surface.opacity = 0.2;
                    vol.Draw(viewer);
                    patch.ComputeResultToGaussPoint(re, ref x, ref y, ref z, ref val);
                }
            }
            ps = new PointSetColor(x.ToArray(), y.ToArray(), z.ToArray(), val.ToArray());
            ps.SetPointSize(10);
            ps.SetOpacity(0.6);
            viewer.AddObject3D(ps);
            viewer.SetColormapBarVisible(ps, re.ToString());
            viewer.Run();
        }


        public static bool IsCreateSparseMatrixForm;
        /// <summary>
        /// Constructor of abstract model
        /// </summary>
        /// <param name="typeProblem">Type of model</param>
        /// <param name="typeAnalysis">Type of analysis model</param>
        /// <param name="structureDimension">Dimension of structure</param>
        /// <param name="pathProject">Path of project</param>
        /// <param name="nameProject">Name of project</param>
        public AbstractModel(TypeModelProblem typeProblem, TypeAnalysisModel typeAnalysis, Dimension structureDimension, string pathProject, string nameProject)
        {
            listPatch = new List<AbstractPatch>();
            loads = new List<AbstractLoad>();
            materials = new List<Material>();
            listConnection = new List<ICoupling>();
            listComputeResult = new List<Result>();
            listConstraint = new List<AbstractConstraintValue>();
            if (pathProject.Substring(pathProject.Length - 2) == "\\")
                PathProject = pathProject;
            else
                PathProject = pathProject + "\\";
            NameProject = nameProject;
            FileNameDataTime = PathProject + "data_time.msgpack";
            FileNameTime = PathProject + "time.msgpack";
            //filenameCheckpoint = new string[] { PathProject + "restore.xml" };

            TypeModel = typeProblem;
            TypeAnalysisModel = typeAnalysis;
            StructureDimension = structureDimension;
            monitorDataHistorism = new MonitorDataHistorism();
            IsParallelProcesing = true;
            IsComputeKeParallel = false;
            //IsObserveForm = false;
            IsSparseData = false;
            IsSaveStepByStep = true;
            IsMatlabSolver = false;
            //IsFinished = false;
            IsRunFromInitial = true;
            IsCreateSparseMatrixForm = true;
            NumberOfStepStorageNonlinear = 10;
            IsDrawConvergencePlot = false;
            IsDrawMonitorPlot = false;

            IndexRowColPlot = new int[] { 0, 1 };
            typeOfFieldsInMultifield = new List<TypeFields>();
            //IsUseOnly1_1Connection = true;
            NumberOfCPUs = Environment.ProcessorCount;
            IsMergeControlpoints1_1Constraint = false;

            if (!Directory.Exists(".\\NMathPackage"))
            {
                MatrixTool.DirectoryCopy("C:\\NMathPackage", ".\\NMathPackage", true);
            }

            NMathConfiguration.LogLocation = ".\\NMathPackage";
            NMathConfiguration.NativeLocation = ".\\NMathPackage";
            NMathConfiguration.LicenseKey = "2DE2E058791B92A";


            //"2DB1884705605CA" "2DAB8554DF60F9F";

            //  NMathConfiguration.LicenseKey = "2DB1884705605CA";//"2DAB8554DF60F9F";
            //  NMathConfiguration.NativeLocation = @"C:\NMathPackage";
            //  NMathConfiguration.LogLocation = @"C:\NMathPackage";
        }
        /// <summary>
        /// Initialize of gauss point of elements
        /// </summary>
        protected void Initialize()
        {
            #region Bezier extraction
            if (IsBezierExtraction)
            {
                int[] p = new int[(int)StructureDimension];
                DoubleVector Bb = null;

                for (int i = 0; i < listPatch.Count; i++)
                {
                    AbstractPatch patchi = listPatch[i];
                    for (int ii = 0; ii < p.Length; ii++)
                        p[ii] = patchi.GetGeometry(0).Basis.GetDegree(ii);
                    int numberOfGaussPointOnEachDirection = (int)Math.Ceiling((2.0 * p.Max() + 1.0) / 2.0);
                    int totalNumberOfGaussPoint = (int)Math.Pow(numberOfGaussPointOnEachDirection, (int)StructureDimension);
                    GaussPoints[] gaussPointOnBezierElement = new GaussPoints[totalNumberOfGaussPoint];
                    //int nen = patchi.GetCountLocalBasisFunctions(0);//p * q
                    AbstractParametricBasisFunction bern = null;
                    int[] indexN = new int[(int)StructureDimension];
                    double[] xi = new double[(int)StructureDimension];
                    for (int k = 0; k < totalNumberOfGaussPoint; k++)
                    {
                        DoubleMatrix dB = null;
                        switch (StructureDimension)
                        {
                            case Dimension.Beam:
                                bern = new BernsteinBasisFunction(p[0]);
                                ///////////////////////

                                indexN[0] = k % numberOfGaussPointOnEachDirection;
                                // xi, eta=[-1, 1]
                                for (int kk = 0; kk < (int)StructureDimension; kk++)
                                    xi[kk] = (GaussPoints.GetPoint(numberOfGaussPointOnEachDirection, indexN[kk]) + 1) / 2; //change scale
                                                                                                                            // xiNew, etaNew=[0, 1]
                                                                                                                            //double xiNew = (xi[0] - (-1)) / 2;
                                                                                                                            //double etaNew = (xi[1] - (-1)) / 2;
                                gaussPointOnBezierElement[k] = new GaussPoints(null, 0);
                                double[] matBb1 = ((BernsteinBasisFunction)bern).GetValueBasisFunctions(xi[0]);
                                double[,] matdBn1 = ((BernsteinBasisFunction)bern).GetDerivativeBasisFunctions(xi[0], 1);
                                Bb = new DoubleVector(matBb1.Length);
                                int c1 = 0;
                                for (int ii = 0; ii < matBb1.Length; ii++)
                                {
                                    Bb[c1] = matBb1[ii];
                                    c1++;
                                }

                                dB = new DoubleMatrix(matdBn1.Length, 2);
                                c1 = 0;
                                for (int ii = 0; ii < matdBn1.Length; ii++)
                                {
                                    dB[c1, 0] = matdBn1[ii, 1];
                                    dB[c1, 1] = matdBn1[ii, 0];
                                    c1++;
                                }
                                break;

                            case Dimension.Plane:
                                bern = new BivariateBernsteinBasisFunction(p[0], p[1]);
                                ///////////////////////

                                indexN[0] = k % numberOfGaussPointOnEachDirection;
                                indexN[1] = k / numberOfGaussPointOnEachDirection;
                                // xi, eta=[-1, 1]
                                for (int kk = 0; kk < (int)StructureDimension; kk++)
                                    xi[kk] = (GaussPoints.GetPoint(numberOfGaussPointOnEachDirection, indexN[kk]) + 1) / 2; //change scale
                                                                                                                            // xiNew, etaNew=[0, 1]
                                                                                                                            //double xiNew = (xi[0] - (-1)) / 2;
                                                                                                                            //double etaNew = (xi[1] - (-1)) / 2;
                                gaussPointOnBezierElement[k] = new GaussPoints(null, 0);
                                double[,] matBb2 = ((BivariateBernsteinBasisFunction)bern).GetValueBivariateBasisFunctions(xi[0], xi[1]);
                                double[,][,] matdBn2 = ((BivariateBernsteinBasisFunction)bern).GetDerivativeBivariateBasisFunctions(xi[0], xi[1], 1);
                                Bb = new DoubleVector(matBb2.GetLength(0) * matBb2.GetLength(1));
                                int c = 0;
                                for (int jj = 0; jj < matBb2.GetLength(1); jj++)
                                    for (int ii = 0; ii < matBb2.GetLength(0); ii++)
                                    {
                                        Bb[c] = matBb2[ii, jj];
                                        c++;
                                    }

                                dB = new DoubleMatrix(matdBn2.GetLength(0) * matdBn2.GetLength(1), 2);
                                c = 0;
                                for (int jj = 0; jj < matdBn2.GetLength(1); jj++)
                                    for (int ii = 0; ii < matdBn2.GetLength(0); ii++)
                                    {
                                        dB[c, 0] = matdBn2[ii, jj][1, 0];
                                        dB[c, 1] = matdBn2[ii, jj][0, 1];
                                        c++;
                                    }
                                break;
                            case Dimension.Solid:
                                bern = new TrivariateBernsteinBasisFunction(p[0], p[1], p[2]);
                                ///////////////////////

                                indexN[0] = (k % (int)Math.Pow(numberOfGaussPointOnEachDirection, (int)StructureDimension - 1)) % numberOfGaussPointOnEachDirection;
                                indexN[1] = (k % (int)Math.Pow(numberOfGaussPointOnEachDirection, (int)StructureDimension - 1)) / numberOfGaussPointOnEachDirection;
                                indexN[2] = k / (int)Math.Pow(numberOfGaussPointOnEachDirection, (int)StructureDimension - 1);
                                // xi, eta=[0, 1]
                                for (int kk = 0; kk < (int)StructureDimension; kk++)
                                    xi[kk] = (GaussPoints.GetPoint(numberOfGaussPointOnEachDirection, indexN[kk]) + 1) / 2; //change scale
                                gaussPointOnBezierElement[k] = new GaussPoints(null, 0);
                                double[,,] matBb3 = ((TrivariateBernsteinBasisFunction)bern).GetValueTrivariateBasisFunctions(xi[0], xi[1], xi[2]);
                                double[,,][,,] matdBn3 = ((TrivariateBernsteinBasisFunction)bern).GetDerivativeTrivariateBasisFunctions(xi[0], xi[1], xi[2], 1);
                                Bb = new DoubleVector(matBb3.GetLength(0) * matBb3.GetLength(1) * matBb3.GetLength(2));
                                int cc = 0;
                                for (int kk = 0; kk < matBb3.GetLength(2); kk++)
                                    for (int jj = 0; jj < matBb3.GetLength(1); jj++)
                                        for (int ii = 0; ii < matBb3.GetLength(0); ii++)
                                        {
                                            Bb[cc] = matBb3[ii, jj, kk];
                                            cc++;
                                        }

                                dB = new DoubleMatrix(matdBn3.GetLength(0) * matdBn3.GetLength(1) * matdBn3.GetLength(2), 3);
                                cc = 0;
                                for (int kk = 0; kk < matdBn3.GetLength(2); kk++)
                                    for (int jj = 0; jj < matdBn3.GetLength(1); jj++)
                                        for (int ii = 0; ii < matdBn3.GetLength(0); ii++)
                                        {
                                            dB[cc, 0] = matdBn3[ii, jj, kk][1, 0, 0];
                                            dB[cc, 1] = matdBn3[ii, jj, kk][0, 1, 0];
                                            dB[cc, 2] = matdBn3[ii, jj, kk][0, 0, 1];
                                            cc++;
                                        }
                                break;
                        }
                        gaussPointOnBezierElement[k].SetValue(DataInGausspoint.NBernsteini, Bb);
                        gaussPointOnBezierElement[k].SetValue(DataInGausspoint.dNBernsteindxi, dB);
                        //gaussPointOnBezierElement[k].NBernsteini = Bb;
                        //gaussPointOnBezierElement[k].dNBernsteindxi = dB;
                    }
                    patchi.ComputeExtractionOperator();
                    patchi.GaussPointOnBezierElement = gaussPointOnBezierElement;
                }
            }
            #endregion
            foreach (AbstractPatch patch in listPatch)
            {
                int nen = patch.CountElements();
                if (!IsParallelProcesing)
                {
                    for (int i = 0; i < nen; i++)
                    {
                        patch.GetElement(i).InitializeGaussPoint();
                    }
                }
                else
                {
                    var degreeOfParallelism = NumberOfCPUs;
                    var tasks = new Task[degreeOfParallelism];
                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;
                        tasks[taskNumber] = Task.Factory.StartNew(
                           () =>
                           {
                               var max = nen * (taskNumberCopy + 1) / degreeOfParallelism;
                               for (int i = nen * taskNumberCopy / degreeOfParallelism; i < max; i++)
                               {
                                   patch.GetElement(i).InitializeGaussPoint();
                               }
                           });
                    }
                    Task.WaitAll(tasks);


                    //Parallel.For(0, nen, i =>
                    //{
                    //  patch.GetElement(i).InitializeGaussPoint();
                    //});
                }
            }
        }

        /// <summary>
        /// Get monitor history data
        /// </summary>
        /// <returns></returns>
        public MonitorDataHistorism GetMonitorDataHistorism()
        { return monitorDataHistorism; }
        /// <summary>
        /// Add monitor data of point
        /// </summary>
        /// <param name="resultType">type of result to store</param>
        /// <param name="indexPatch">index of patch containt point to get result</param>
        /// <param name="xi">parametric coordinations</param>
        public void AddMonitorData(Result resultType, int indexPatch, params double[] xi)
        {
            if (resultType != Result.REACFORCEX && resultType != Result.REACFORCEY && resultType != Result.REACFORCEZ)
                monitorDataHistorism.AddMonitorDataNode(resultType, listPatch[indexPatch], xi[0], xi[1], (StructureDimension == Dimension.Plane) ? 0 : xi[2]);
            else
                throw new ArgumentException("Only use this method excepting reaction forces");
        }

        public void AddMonitorData(Result resultType, int indexConstraint)
        {
            if (resultType == Result.REACFORCEX || resultType == Result.REACFORCEY || resultType == Result.REACFORCEZ)
                monitorDataHistorism.AddMonitorDataNode(resultType, indexConstraint);
            else
                throw new ArgumentException("Only use this method in case of reaction forces");
        }

        public void AddMonitorData(Result resultType, DoFormulationOnMonitorData doFormulation, params int[] indexCol)
        {
            if (resultType == Result.USERDEFINED)
                monitorDataHistorism.AddMonitorDataNode(resultType, doFormulation, indexCol);
            else
                throw new ArgumentException("Only use this method for user defined");
        }
        public void AddIndexColOutput(params int[] indexColMonitorOutput)
        {
            monitorDataHistorism.IndexColOutputMonitor = indexColMonitorOutput;
        }
        ///// <summary>
        ///// Get observe windows form
        ///// </summary>
        ///// <returns></returns>
        //public ObserveForm GetObserveForm()
        //{
        //    return observe;
        //}

        /// <summary>
        /// Get patch by index of patch in list patch
        /// </summary>
        /// <param name="index">index</param>
        /// <returns></returns>
        public AbstractPatch GetPatch(int index)
        {
            return listPatch[index];
        }

        /// <summary>
        /// Add patch into model
        /// </summary>
        /// <param name="geo">Geometry of patch</param>
        public void AddPatch(AbstractParametricGeometry geo)
        {
            AbstractPatch m = null;
            switch (StructureDimension)
            {
                case Dimension.Plane:
                    switch (TypeModel)
                    {
                        case TypeModelProblem.Structural:
                        case TypeModelProblem.Piezoelectric:
                        case TypeModelProblem.PhaseField:
                            m = new PatchStructure2D((NURBSSurface)geo, Structure2DState.PlaneStrain);
                            break;
                        case TypeModelProblem.Thermal:
                            m = new PatchThermal2D((Abstract2DParametricGeometry)geo);
                            break;
                        case TypeModelProblem.StructuralThermal:
                            m = new PatchThermoelastic2D((NURBSSurface)geo, Structure2DState.PlaneStrain);
                            break;
                    }
                    break;
                case Dimension.Solid:
                    switch (TypeModel)
                    {
                        case TypeModelProblem.Structural:
                        case TypeModelProblem.Piezoelectric:
                        case TypeModelProblem.PhaseField:
                            m = new PatchStructure3D((NURBSVolume)geo);
                            break;
                        case TypeModelProblem.Thermal:
                            m = new PatchThermal3D((NURBSVolume)geo);
                            break;
                        case TypeModelProblem.StructuralThermal:
                            m = new PatchThermoelastic3D((NURBSVolume)geo);
                            break;
                    }
                    break;
                case Dimension.Plate:
                    m = new PatchStructurePlate((NURBSSurface)geo, 1, ((AbstractModelStructure)this).TypePlate);
                    break;
                case Dimension.Beam:
                    m = new PatchStructureBeam((NURBSCurve)geo, 1, ((AbstractModelStructure)this).TypeBeam);
                    break;
            }
            m.ID = listPatch.Count;
            listPatch.Add(m);
        }

        /// <summary>
        /// Count number of patch in model
        /// </summary>
        /// <returns></returns>
        public int CountPatch()
        {
            return listPatch.Count;
        }

        /// <summary>
        /// Add material into model
        /// </summary>
        /// <param name="ma">material</param>
        public void AddMaterial(Material ma)
        {
            materials.Add(ma);
        }

        /// <summary>
        /// Get material by index
        /// </summary>
        /// <param name="index">index material</param>
        /// <returns></returns>
        public Material GetMaterial(int index)
        {
            return materials[index];
        }

        /// <summary>
        /// Count number of material
        /// </summary>
        /// <returns></returns>
        public int CountMaterial()
        {
            return materials.Count;
        }

        /// <summary>
        /// Attach material to patch
        /// </summary>
        /// <param name="indexMaterial">index of material in model</param>
        /// <param name="indexPatch">index of patch in model</param>
        public void AttachMaterialToPatch(int indexMaterial, int indexPatch)
        {
            listPatch[indexPatch].SetMaterial(materials[indexMaterial]);
            //switch (StructureDimension)
            //{
            //    case Dimension.Plane:
            //        ((AbstractPatch2D)listPatch[indexPatch]).SetMaterial(materials[indexMaterial]);
            //        break;
            //    case Dimension.Solid:
            //        ((AbstractPatch3D)listPatch[indexPatch]).SetMaterial(materials[indexMaterial]);
            //        break;
            //}
            //listPatch[indexPatch].Initialize();
        }
        /// <summary>
        /// Initializa patch to create elements
        /// </summary>
        public void InitializePatch()
        {
            if (!IsParallelProcesing)
            {
                for (int i = 0; i < listPatch.Count; i++)
                {
                    listPatch[i].TypeOfFieldsInMultifield = typeOfFieldsInMultifield.ToList();
                    if (listPatch[i].DisableField != null)
                    {
                        foreach (var field in listPatch[i].DisableField)
                        {
                            listPatch[i].TypeOfFieldsInMultifield.Remove(field);
                        }
                    }

                    listPatch[i].Initialize();
                }
            }
            else
            {
                var degreeOfParallelism = NumberOfCPUs;
                var tasks = new Task[degreeOfParallelism];
                object monitor = new object();
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                         () =>
                         {
                             var max = listPatch.Count * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = listPatch.Count * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 listPatch[i].TypeOfFieldsInMultifield = typeOfFieldsInMultifield.ToList();
                                 if (listPatch[i].DisableField != null)
                                 {
                                     foreach (var field in listPatch[i].DisableField)
                                     {
                                         listPatch[i].TypeOfFieldsInMultifield.Remove(field);
                                     }
                                 }
                                 listPatch[i].Initialize();
                             }
                         });
                }
                Task.WaitAll(tasks);
            }
        }
        /// <summary>
        /// Add load of model
        /// </summary>
        /// <param name="load">load</param>
        public void AddLoad(AbstractLoad load)
        {
            loads.Add(load);
        }

        /// <summary>
        /// Get load by index
        /// </summary>
        /// <param name="index">index</param>
        /// <returns></returns>
        public AbstractLoad GetLoad(int index)
        {
            return loads[index];
        }

        /// <summary>
        /// Count number of load
        /// </summary>
        /// <returns></returns>
        public int CountLoad()
        {
            return loads.Count;
        }
        /// <summary>
        /// Get number of degree of freedom in model
        /// </summary>
        /// <returns></returns>
        public int CountDOF()
        {
            return countDOF;
        }
        public int CountMaxNumberOfField()
        {
            int maxCountField = int.MinValue;
            foreach (var item in listPatch)
            {
                maxCountField = Math.Max(maxCountField, item.GetCountField(0));
            }
            return maxCountField;
        }
        /// <summary>
        /// Solve problem
        /// </summary>
        public abstract void Solve();// int numberOfLoadStep = 1, int numberOfSubstepLoad = 1, double TOL = 0.000001, int maxIter = 1000);

        /// <summary>
        /// PreProcessing to initialize and enumerate model
        /// </summary>
        public virtual void PreProcessing()
        {
            if (IsWriteLogFile)
            {
                stopWatchWholeModel = new Stopwatch();
                stopWatchWholeModel.Start();
                logFile = new List<string>();
            }
            foreach (AbstractPatch patch in listPatch)
            {
                patch.AssignMaterialToElements();
                if (this is IModelStructure)
                    if (StructureDimension == Dimension.Plane)
                    {
                        ((IPatchStructure)patch).StateStress = ((IModelStructure)this).StressState;
                    }
            }
            Initialize();
            Enumerate();//enumerate in global
            if (IsSparseData)
            {
                kGlobalSparse = CreateFormOfStiffnessMatrix();
                if (this is ModelStructureModal || this is ModelPiezoelectricModal || this is ModelStructureBuckling)
                    mGlobalSparse = (DoubleCsrSparseMatrix)kGlobalSparse.Clone();
            }
            bcMatrix = StoreBoundaryCondition();
            DetachedCoincidentUncoupling();//Detached coupling at point meeting but they are separated
            ApplyCouplingInterfaces();//Create transform matrix
            ApplyUnchangeBoundaryCondition();//Create T and g with unchange boundary conditions
            CalculateMasterSlaveTMatrixAndgVector();
            CreateDependentIndependentUncouple();

            //if (!IsParallelProcesing)
            //{
            RemoveInIndependentByDependent();
            TMatrixPlus1OnDiagonal();
            //}
            //else
            //{
            //  Task[] tasks = new Task[2]{
            //    Task.Factory.StartNew(() => RemoveInIndependentByDependent()),
            //    Task.Factory.StartNew(() => TMatrixPlus1OnDiagonal())
            //  };
            //  //Block until all tasks complete.
            //  Task.WaitAll(tasks);
            //}
            UpdateTMatrixAndTMatrixTranspose();
        }
        private void CalculateMasterSlaveTMatrixAndgVector()
        {

            if (!IsSparseData)
            {
                #region run_good
                ////////run good//////
                ////Substitute to masterDOF of slaveDOF to have one slaveDOF has many last masterDOFs with multi-level masters
                // On diagonal must be equal -1
                for (int i = 0; i < countDOF; i++)////////////////////////////////////
                {
                    double Tii = TDense[i, i];
                    if (Tii != 0)
                    {
                        if (Tii != -1)
                        {
                            for (int j = 0; j < countDOF; j++)
                            {
                                TDense[i, j] = Math.Round(TDense[i, j] / (-Tii), numberOfDecimal);
                            }
                            g[i] = Math.Round(g[i] / (-Tii), numberOfDecimal);
                        }
                        List<int> indexRow = MatrixTool.FindIndexItemNonZeroOnCol(TDense, i, AbstractModel.NumberOfCPUs);
                        List<int> indexCol = MatrixTool.FindIndexItemNonZeroOnRow(TDense, i, AbstractModel.NumberOfCPUs);
                        for (int j = 0; j < indexRow.Count; j++)
                        {
                            if (i != indexRow[j])
                            {
                                double f = TDense[indexRow[j], i];

                                for (int k = 0; k < indexCol.Count; k++)
                                {
                                    if (i != indexCol[k])
                                    {
                                        TDense[indexRow[j], indexCol[k]] = Math.Round(TDense[indexRow[j], indexCol[k]] + f * TDense[i, indexCol[k]], numberOfDecimal);
                                    }
                                }
                                g[indexRow[j]] = Math.Round(g[indexRow[j]] + f * g[i], numberOfDecimal);
                                TDense[indexRow[j], i] = 0;
                            }
                        }
                    }
                }
                #endregion
            }
            else
            {
                ////////////////////////////////////////////////////////////
                #region run_good
                ////////run good//////
                ////Substitute to masterDOF of slaveDOF to have one slaveDOF has many last masterDOFs with multi-level masters
                // On diagonal must be equal -1
                for (int i = 0; i < countDOF; i++)////////////////////////////////////
                {
                    if (TSparseBuilder[i, i] != 0)
                    //if (Math.Abs(TSparseBuilder[i, i]) > 1e-9)
                    {
                        double Tii = TSparseBuilder[i, i];
                        if (Tii != -1)
                        {
                            for (int j = 0; j < countDOF; j++)
                            {
                                TSparseBuilder[i, j] = Math.Round(TSparseBuilder[i, j] / (-Tii), numberOfDecimal);
                            }
                            g[i] = Math.Round(g[i] / (-Tii), numberOfDecimal);
                        }
                        List<int> indexRow = MatrixTool.FindIndexItemNonZeroOnCol(TSparseBuilder, i, AbstractModel.countDOF, AbstractModel.NumberOfCPUs);
                        List<int> indexCol = MatrixTool.FindIndexItemNonZeroOnRow(TSparseBuilder, i, AbstractModel.countDOF, AbstractModel.NumberOfCPUs);
                        for (int j = 0; j < indexRow.Count; j++)
                        {
                            if (i != indexRow[j])
                            {
                                double f = TSparseBuilder[indexRow[j], i];

                                for (int k = 0; k < indexCol.Count; k++)
                                {
                                    if (i != indexCol[k])
                                    {
                                        TSparseBuilder[indexRow[j], indexCol[k]] = Math.Round(TSparseBuilder[indexRow[j], indexCol[k]] + f * TSparseBuilder[i, indexCol[k]], numberOfDecimal);
                                    }
                                }
                                //g[indexRow[j]] = g[indexRow[j]] + f * g[i];
                                g[indexRow[j]] = Math.Round(g[indexRow[j]] + f * g[i], numberOfDecimal);
                                TSparseBuilder[indexRow[j], i] = 0;
                            }
                        }
                    }
                }
                #endregion
            }
        }
        protected void ApplyChangeBoundaryConditionAndUpdateIMatrix(double time)
        {
            ApplyChangeBoundaryCondition(time);
            TMatrixPlus1OnDiagonal();
            UpdateTMatrixAndTMatrixTranspose();
        }
        protected void UpdateTMatrixAndTMatrixTranspose()
        {
            if (!IsSparseData)
            {
                TDenseTranspose = MatrixFunctions.Transpose(TDense);
            }
            else
            {
                MatrixTool.CorrectSparseMatrixBuilderZeroRows(TSparseBuilder, countDOF, countDOF);
                TSparse = new DoubleCsrSparseMatrix(TSparseBuilder, countDOF);
                //Stopwatch sw = new Stopwatch();
                //sw.Start();
                //bool a = true;
                //if (a)
                TSparseTranspose = MatrixTool.Transpose(TSparse);//2440
                                                                 //else
                                                                 //{
                                                                 //  SparseMatrixBuilder<double> sparseMatrixBuilder = MatrixTool.Transpose(TSparseBuilder);
                                                                 //  TSparseTranspose = new DoubleCsrSparseMatrix(sparseMatrixBuilder, countDOF);
                                                                 //  ////2737.8567000000003
                                                                 //}
                                                                 //sw.Stop();
                                                                 //double t = sw.Elapsed.TotalMilliseconds;
            }
        }

        protected void TMatrixPlus1OnDiagonal()
        {
            if (!IsSparseData)
            {
                if (!IsParallelProcesing)
                {
                    for (int i = 0; i < countDOF; i++)
                    {
                        TDense[i, i] += 1;
                    }
                }
                else
                {
                    var degreeOfParallelism = NumberOfCPUs;
                    var tasks = new Task[degreeOfParallelism];
                    object monitor = new object();
                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = countDOF * (taskNumberCopy + 1) / degreeOfParallelism;
                                 //lock (monitor)
                                 //{
                                 for (int i = countDOF * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     TDense[i, i] += 1;
                                 }
                                 //}
                             });
                    }
                    Task.WaitAll(tasks);
                }
            }
            else
            {
                if (!IsParallelProcesing)
                {
                    for (int i = 0; i < countDOF; i++)
                    {
                        TSparseBuilder[i, i] += 1;
                    }
                }
                else
                {
                    var degreeOfParallelism = NumberOfCPUs;
                    var tasks = new Task[degreeOfParallelism];
                    object monitor = new object();
                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = countDOF * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = countDOF * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     lock (monitor)
                                     {
                                         TSparseBuilder[i, i] += 1;
                                     }
                                 }
                             });
                    }
                    Task.WaitAll(tasks);
                }
            }

            //if (!IsParallelProcesing)
            //{
            //  if (!IsSparseData)
            //  {
            //    for (int i = 0; i < countDOF; i++)
            //    {
            //      TDense[i, i] += 1;
            //    }
            //  }
            //  else
            //  {
            //    for (int i = 0; i < countDOF; i++)
            //    {
            //      TSparseBuilder[i, i] += 1;
            //    }
            //  }
            //}
            //else
            //{
            //  var degreeOfParallelism = NumberOfCPUs;
            //  var tasks = new Task[degreeOfParallelism];
            //  object monitor = new object();
            //  for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
            //  {
            //    // capturing taskNumber in lambda wouldn't work correctly
            //    int taskNumberCopy = taskNumber;

            //    tasks[taskNumber] = Task.Factory.StartNew(
            //         () =>
            //         {
            //           var max = countDOF * (taskNumberCopy + 1) / degreeOfParallelism;
            //           lock (monitor)
            //           {
            //             for (int i = countDOF * taskNumberCopy / degreeOfParallelism; i < max; i++)
            //             {
            //               if (!IsSparseData)
            //               {
            //                 TDense[i, i] += 1;
            //               }
            //               else
            //               {
            //                 TSparseBuilder[i, i] += 1;
            //               }
            //             }
            //           }
            //         });
            //  }
            //  Task.WaitAll(tasks);
            //}
        }

        private void CreateDependentIndependentUncouple()
        {
            //string filenameSaveTArray = PathProject + "restore_tArray.soap";
            string filenameSaveTArray = PathProject + "restore_tArray.msgpack";
            string directory = Path.GetDirectoryName(filenameSaveTArray);
            if (!Directory.Exists(directory) || !File.Exists(filenameSaveTArray))
            {
                //if (!IsParallelProcesing)
                //{
                //Stopwatch sw = new Stopwatch();
                //sw.Start();
                CreateDependent();
                CreateIndependent();
                CreateUnCouple();
                //}
                //else
                //{
                //  Task[] tasks = new Task[3]{
                //    Task.Factory.StartNew(() => CreateDependent()),
                //    Task.Factory.StartNew(() => CreateIndependent()),
                //    Task.Factory.StartNew(() => CreateUnCouple())
                //  };
                //  //Block until all tasks complete.
                //  Task.WaitAll(tasks);
                //}
                //SaveSoapData(filenameSaveTArray, tArrayDependent.ToArray(), tArrayIndependent.ToArray(), tArrayUncouple.ToArray());
                //sw.Stop();
                //double aa = sw.Elapsed.TotalMilliseconds;//12080.7625
                //sw.Reset();
                //sw.Start();
                //var tarray = LoadSoapData<int[]>(filenameSaveTArray);
                //List<int> tArrayDependent1 = tarray[0].ToList<int>();
                //List<int> tArrayIndependent1 = tarray[1].ToList<int>();
                //List<int> tArrayUncouple1 = tarray[2].ToList<int>();
                //sw.Stop();
                //double bb = sw.Elapsed.TotalMilliseconds;//28.250700000000002
                IOIGA.WriteMessagePackData(filenameSaveTArray, tArrayDependent, tArrayIndependent, tArrayUncouple);
                //sw.Stop();
                //double aa = sw.Elapsed.TotalMilliseconds;//14377.7263
                //sw.Reset();
                //sw.Start();
                //List<int>[] temp = LoadMessagePackData<List<int>>(PathProject + "restore_tArray.msgpak");
                //List<int> tArrayDependent1 = temp[0];
                //List<int> tArrayIndependent1 = temp[1];
                //List<int> tArrayUncouple1 = temp[2];
                //sw.Stop();
                //double bb = sw.Elapsed.TotalMilliseconds;//4.7027
            }
            else
            {
                //var tarray = LoadSoapData<int[]>(filenameSaveTArray);
                //tArrayDependent = tarray[0].ToList<int>();
                //tArrayIndependent = tarray[1].ToList<int>();
                //tArrayUncouple = tarray[2].ToList<int>();
                List<int>[] temp = IOIGA.ReadMessagePackData<List<int>>(filenameSaveTArray);
                tArrayDependent = temp[0];
                tArrayIndependent = temp[1];
                tArrayUncouple = temp[2];
            }
        }

        protected void RemoveInIndependentByDependent()
        {
            for (int i = 0; i < tArrayDependent.Count; i++)
            {
                bool a = tArrayIndependent.Remove(tArrayDependent[i]);
                if (!a)
                    throw new Exception("Error");
            }
        }

        private void CreateIndependent()
        {
            if (!IsSparseData)
            {
                tArrayIndependent = MatrixTool.FindIndexColAllNotEqualZero(TDense, AbstractModel.NumberOfCPUs);
            }
            else
            {
                tArrayIndependent = MatrixTool.FindIndexColAllNotEqualZero(TSparseBuilder, AbstractModel.countDOF, AbstractModel.NumberOfCPUs);
            }
        }

        private void CreateUnCouple()
        {
            if (!IsSparseData)
            {
                tArrayUncouple = MatrixTool.FindIndexColAllEqualZero(TDense, AbstractModel.NumberOfCPUs);
            }
            else
            {
                tArrayUncouple = MatrixTool.FindIndexColAllEqualZero(TSparseBuilder, AbstractModel.countDOF, AbstractModel.NumberOfCPUs);
            }
        }

        private void CreateDependent()
        {
            tArrayDependent = new List<int>();
            if (!IsParallelProcesing)
            {
                if (!IsSparseData)
                {
                    for (int i = 0; i < countDOF; i++)
                    {
                        if (TDense[i, i] == -1)
                        {
                            tArrayDependent.Add(i);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < countDOF; i++)
                    {
                        if (TSparseBuilder[i, i] == -1)
                        {
                            tArrayDependent.Add(i);
                        }
                    }
                }
            }
            else
            {
                var degreeOfParallelism = NumberOfCPUs;
                object monitor = new object();
                var tasks = new Task[degreeOfParallelism];
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                        () =>
                        {
                            var max = countDOF * (taskNumberCopy + 1) / degreeOfParallelism;
                            for (int i = countDOF * taskNumberCopy / degreeOfParallelism; i < max; i++)
                            {
                                lock (monitor)
                                {
                                    if (!IsSparseData)
                                    {
                                        if (TDense[i, i] == -1)
                                        {
                                            tArrayDependent.Add(i);
                                        }
                                    }
                                    else
                                    {
                                        if (TSparseBuilder[i, i] == -1)
                                        {
                                            tArrayDependent.Add(i);
                                        }
                                    }
                                }
                            }
                        });
                }
                Task.WaitAll(tasks);
            }
        }

        /// <summary>
        /// Post processing to get last result uLocal to control points
        /// </summary>
        public virtual void PostProcessing()
        {
            //if (DisplacementTime == null)
            //  DisplacementTime = IOIGA.ReadListDoubleVectorMessagePackFormatter(FileNameDataTime);
            string fileNameDataTime = FileNameDataTime;
            try
            {
                Time = IOIGA.ReadMessagePackData<List<double>>(FileNameTime)[0];
            }
            catch
            {
                Time = null;
            }

            if (IsSaveStepByStep)
            {
                string[] nameFiles = Directory.GetFiles(PathProject, "data_*.msgpack");
                fileNameDataTime = PathProject + "data_" + nameFiles.Length + ".msgpack";
            }
            DoubleVector uLast = IOIGA.ReadMessagePackDoubleVector(fileNameDataTime);// DisplacementTime[DisplacementTime.Count - 1];
            SetUGlobal(uLast);
            if (IsWriteLogFile)
            {
                stopWatchWholeModel.Stop();
                IOIGA.WriteLogFileAndConsole(logFile, true, "-----------------------------------------------------------");
                IOIGA.WriteLogFileAndConsole(logFile, true, "Finished!");
                IOIGA.WriteLogFileAndConsole(logFile, true, "Elapsed time = " + stopWatchWholeModel.Elapsed.TotalSeconds + "[s]");
                IOIGA.WriteLogFile(logFile, PathProject + NameProject + ".log");
                IOIGA.WriteLogFileAndConsole(null, true, "Saved log file at " + Directory.GetCurrentDirectory() + "\\" + PathProject + NameProject + ".log");
            }
        }
        public void Extrapolation(params Result[] results)
        {
            object monitor = new object();

            foreach (Result re in results)
            {
                for (int iji = 0; iji < listPatch.Count; iji++)
                {
                    AbstractPatch patch = listPatch[iji];
                    int nen = patch.GetCountLocalBasisFunctions(0);//(p + 1) * (q + 1);
                    int nnp1 = patch.GetCountGlobalBasisFunctions(0);
                    SparseMatrixBuilder<double> AsBuilder = new SparseMatrixBuilder<double>();
                    DoubleCsrSparseMatrix AsSparse = null;// ComputeAsCsrSparseMatrix(patch);
                    DoubleMatrix As = null;// new DoubleMatrix(nnp1, nnp1);
                    if (!IsSparseData)
                        As = new DoubleMatrix(nnp1, nnp1);
                    else
                        AsSparse = ComputeAsCsrSparseMatrix(patch);
                    DoubleVector Bs = new DoubleVector(nnp1);
                    if (!IsParallelProcesing)
                    {
                        for (int iel = 0; iel < patch.CountElements(); iel++)
                        {
                            CreateAsBs(re, patch, ref As, ref AsSparse, ref Bs, iel, monitor);
                        }
                    }
                    else
                    {
                        var degreeOfParallelism = NumberOfCPUs;

                        var tasks = new Task[degreeOfParallelism];
                        for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                        {
                            // capturing taskNumber in lambda wouldn't work correctly
                            int taskNumberCopy = taskNumber;

                            tasks[taskNumber] = Task.Factory.StartNew(
                                () =>
                                {
                                    var max = patch.CountElements() * (taskNumberCopy + 1) / degreeOfParallelism;
                                    for (int iel = patch.CountElements() * taskNumberCopy / degreeOfParallelism; iel < max; iel++)
                                    {
                                        CreateAsBs(re, patch, ref As, ref AsSparse, ref Bs, iel, monitor);
                                    }
                                });
                        }
                        Task.WaitAll(tasks);
                    }
                    DoubleVector Ps = (!IsSparseData) ? MatrixFunctions.Solve(As, Bs) : MatrixFunctions.Solve(AsSparse, Bs);
                    //if (!IsParallelProcesing)
                    //{
                    for (int ii = 0; ii < nnp1; ii++)
                    {
                        if (patch is AbstractPatch2D)
                        {
                            int ik = patch.GetINC(0, ii, 0);
                            int jk = patch.GetINC(0, ii, 1);
                            ((NURBSSurface)patch.GetGeometry(0)).ControlPoints[ik, jk].AddResult(re, Ps[ii]);
                        }
                        else if (patch is AbstractPatch3D)
                        {
                            int ik = patch.GetINC(0, ii, 0);
                            int jk = patch.GetINC(0, ii, 1);
                            int kk = patch.GetINC(0, ii, 2);
                            ((NURBSVolume)patch.GetGeometry(0)).ControlPoints[ik, jk, kk].AddResult(re, Ps[ii]);
                        }
                    }
                    //}
                    //else
                    //{
                    //  var degreeOfParallelism = NumberOfCPUs;

                    //  var tasks = new Task[degreeOfParallelism];
                    //  for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    //  {
                    //    // capturing taskNumber in lambda wouldn't work correctly
                    //    int taskNumberCopy = taskNumber;

                    //    tasks[taskNumber] = Task.Factory.StartNew(
                    //        () =>
                    //        {
                    //          var max = nnp1 * (taskNumberCopy + 1) / degreeOfParallelism;
                    //          for (int ii = nnp1 * taskNumberCopy / degreeOfParallelism; ii < max; ii++)
                    //          {
                    //            if (patch is AbstractPatch2D)
                    //            {
                    //              int ik = patch.GetINC(0, ii, 0);
                    //              int jk = patch.GetINC(0, ii, 1);
                    //              ((NURBSSurface)patch.GetGeometry(0)).ControlPoints[ik, jk].AddResult(re, Ps[ii]);
                    //            }
                    //            else if (patch is AbstractPatch3D)
                    //            {
                    //              int ik = patch.GetINC(0, ii, 0);
                    //              int jk = patch.GetINC(0, ii, 1);
                    //              int kk = patch.GetINC(0, ii, 2);
                    //              ((NURBSVolume)patch.GetGeometry(0)).ControlPoints[ik, jk, kk].AddResult(re, Ps[ii]);
                    //            }
                    //          }
                    //        });
                    //  }
                    //  Task.WaitAll(tasks);
                    //}
                }
            }
        }

        //private void CreateAsBsPlate(Result re, int iji, int nen, DoubleMatrix As, DoubleVector Bs, int iel, object monitor)
        //{
        //  AbstractElementStructurePlate elem = (AbstractElementStructurePlate)listPatch[iji].GetElement(iel);
        //  int[] econn = new int[nen];

        //  for (int k = 0; k < nen; k++)
        //  {
        //    econn[k] = ((AbstractPatch2D)listPatch[iji]).GetIEN(iel, k);
        //  }

        //  DoubleMatrix Ae = new DoubleMatrix(nen, nen);
        //  DoubleVector Be = new DoubleVector(nen);

        //  ////change order///////
        //  for (int j = 0; j < elem.GetNumberOfGaussPointOnEachDirection(); j++)
        //  {
        //    for (int i = 0; i < elem.GetNumberOfGaussPointOnEachDirection(); i++)
        //    {
        //      GaussPoints gpsij = elem.GetGaussPoint(i, j);
        //      DoubleVector Nii = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);//new DoubleVector(nen);
        //      //for (int k = 0; k < nen; k++)
        //      //  Nii[k] = gpsij.Ni[k];
        //      Ae += MatrixFunctions.OuterProduct(Nii, Nii);
        //      double reExtra = 0;
        //      switch (re)
        //      {
        //        case Result.SIGMAXX:
        //          if ((listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity))
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[0];
        //          }
        //          else if (listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[0];//gpsij.lastStress[0];
        //          }
        //          break;
        //        case Result.SIGMAYY:
        //          if ((listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity))
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[1];
        //          }
        //          else if (listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[1];
        //          }
        //          break;
        //        case Result.SIGMAXY:
        //          if ((listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity))
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[2];
        //          }
        //          else if (listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[2];
        //          }
        //          break;
        //        case Result.SIGMAXZ:
        //          if ((listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity))
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[3];
        //          }
        //          else if (listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[3];
        //          }
        //          break;
        //        case Result.SIGMAYZ:
        //          if ((listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity))
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[4];
        //          }
        //          else if (listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[4];
        //          }
        //          break;
        //        case Result.EPSILONXX:
        //          if ((listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity))
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[0];
        //          }
        //          else if (listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[0];//gpsij.lastStrain[0];
        //          }
        //          break;
        //        case Result.EPSILONYY:
        //          if ((listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity))
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[1];
        //          }
        //          else if (listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[1];
        //          }
        //          break;
        //        case Result.EPSILONXY:
        //          if ((listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity))
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[1];
        //          }
        //          else if (listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[2];
        //          }
        //          break;
        //        case Result.EPSILONXZ:
        //          if ((listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity))
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[3];
        //          }
        //          else if (listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[3];
        //          }
        //          break;
        //        case Result.EPSILONYZ:
        //          if ((listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity))
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[4];
        //          }
        //          else if (listPatch[iji].Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[4];
        //          }
        //          break;
        //        case Result.EPSILONPEQV:
        //          reExtra = (double)gpsij.GetValue(DataInGausspoint.lastAlpha);//gpsij.lastAlpha;
        //          break;
        //      }
        //      Be += reExtra * Nii;
        //    }
        //  }
        //  lock (monitor)
        //  {
        //    for (int jj = 0; jj < nen; jj++)
        //    {
        //      for (int ii = 0; ii < nen; ii++)
        //      {
        //        As[econn[ii], econn[jj]] += Ae[ii, jj];
        //      }
        //      Bs[econn[jj]] += Be[jj];
        //    }
        //  }
        //}
        private void CreateAsBs(Result re, AbstractPatch patch, ref DoubleMatrix As, ref DoubleCsrSparseMatrix AsSparse, ref DoubleVector Bs, int iel, object monitor)
        {
            int nen = patch.GetCountLocalBasisFunctions(0);
            int[] econn = new int[nen];

            for (int k = 0; k < nen; k++)
            {
                econn[k] = patch.GetIEN(0, iel, k);
            }
            DoubleMatrix Ae = new DoubleMatrix(nen, nen);
            DoubleVector Be = new DoubleVector(nen);
            if (patch is AbstractPatch2D)
            {
                AbstractElementStructure2D elem = (AbstractElementStructure2D)patch.GetElement(iel);
                ////change order///////
                for (int j = 0; j < elem.GetNumberOfGaussPointOnEachDirection(); j++)
                {
                    for (int i = 0; i < elem.GetNumberOfGaussPointOnEachDirection(); i++)
                    {
                        GaussPoints gpsij = elem.GetGaussPoint(i, j);
                        DoubleVector Nii = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni); //new DoubleVector(nen);
                        Ae += MatrixFunctions.OuterProduct(Nii, Nii);
                        double reExtra = 0;
                        switch (re)
                        {
                            case Result.SIGMAXX:

                                if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                                {
                                    reExtra = elem.StressAt(gpsij.location)[0];
                                }
                                else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                                {
                                    reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[0];
                                }
                                break;
                            case Result.SIGMAYY:

                                if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                                {
                                    reExtra = elem.StressAt(gpsij.location)[1];
                                }
                                else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                                {
                                    reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[1];
                                }
                                break;
                            case Result.SIGMAZZ:
                                if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                                {
                                    reExtra = elem.StressAt(gpsij.location)[2];
                                }
                                else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                                {
                                    reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[2];
                                }
                                break;
                            case Result.SIGMAXY:

                                if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                                {
                                    reExtra = elem.StressAt(gpsij.location)[3];
                                }
                                else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                                {
                                    reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[3];
                                }
                                break;
                            case Result.EPSILONXX:
                                if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                                {
                                    reExtra = elem.StrainAt(gpsij.location)[0];
                                }
                                else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                                {
                                    reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[0];
                                }
                                break;
                            case Result.EPSILONYY:
                                if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                                {
                                    reExtra = elem.StrainAt(gpsij.location)[1];
                                }
                                else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                                {
                                    reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[1];
                                }
                                break;
                            case Result.EPSILONZZ:
                                if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                                {
                                    reExtra = elem.StrainAt(gpsij.location)[2];
                                }
                                else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                                {
                                    reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[2];
                                }
                                break;
                            case Result.EPSILONXY:
                                if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                                {
                                    reExtra = elem.StrainAt(gpsij.location)[3];
                                }
                                else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                                {
                                    reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[3];
                                }
                                break;
                            case Result.EPSILONTHXX:
                                reExtra = ((ElementThermoelastic2D)elem).StrainThermoAt(gpsij.location)[0];
                                break;
                            case Result.EPSILONTHYY:
                                reExtra = ((ElementThermoelastic2D)elem).StrainThermoAt(gpsij.location)[1];
                                break;
                            case Result.EPSILONTHZZ:
                                reExtra = ((ElementThermoelastic2D)elem).StrainThermoAt(gpsij.location)[2];
                                break;
                            case Result.EPSILONEXX:
                                reExtra = ((ElementThermoelastic2D)elem).StrainElasticAt(gpsij.location)[0];
                                break;
                            case Result.EPSILONEYY:
                                reExtra = ((ElementThermoelastic2D)elem).StrainElasticAt(gpsij.location)[1];
                                break;
                            case Result.EPSILONEZZ:
                                reExtra = ((ElementThermoelastic2D)elem).StrainElasticAt(gpsij.location)[2];
                                break;
                            case Result.EPSILONEXY:
                                reExtra = ((ElementThermoelastic2D)elem).StrainElasticAt(gpsij.location)[3];
                                break;
                            case Result.EPSILONPEQV:
                                reExtra = (double)gpsij.GetValue(DataInGausspoint.lastAlpha);//gpsij.lastAlpha;
                                break;
                            default:
                                reExtra = 0;
                                break;
                        }
                        Be += reExtra * Nii;
                    }
                }
            }
            else if (patch is AbstractPatch3D)
            {
                AbstractElementStructure3D elem = (AbstractElementStructure3D)patch.GetElement(iel);
                for (int i = 0; i < elem.GetNumberOfGaussPointOnEachDirection(); i++)
                {
                    for (int j = 0; j < elem.GetNumberOfGaussPointOnEachDirection(); j++)
                    {
                        for (int k = 0; k < elem.GetNumberOfGaussPointOnEachDirection(); k++)
                        {
                            GaussPoints gpsijk = elem.GetGaussPoint(i, j, k);
                            DoubleVector Nii = (DoubleVector)gpsijk.GetValue(DataInGausspoint.Ni); //new DoubleVector(nen);
                            Ae += MatrixFunctions.OuterProduct(Nii, Nii);
                            double reExtra = ComputeResultValueExtrapolation(re, patch, elem, gpsijk);
                            Be += reExtra * Nii;
                        }
                    }
                }
            }
            lock (monitor)
            {
                for (int jj = 0; jj < nen; jj++)
                {
                    for (int ii = 0; ii < nen; ii++)
                    {
                        if (!IsSparseData)
                            As[econn[ii], econn[jj]] += Ae[ii, jj];
                        else
                            AsSparse[econn[ii], econn[jj]] += Ae[ii, jj];
                    }
                    Bs[econn[jj]] += Be[jj];
                }
            }
        }

        private static double ComputeResultValueExtrapolation(Result re, AbstractPatch patch, AbstractElementStructure3D elem, GaussPoints gpsijk)
        {
            double reExtra = 0;
            switch (re)
            {
                case Result.SIGMAXX:
                    if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                    {
                        reExtra = elem.StressAt(gpsijk.location)[0];
                    }
                    else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                    {
                        reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStress))[0];
                    }
                    break;

                case Result.SIGMAYY:
                    if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                    {
                        reExtra = elem.StressAt(gpsijk.location)[1];
                    }
                    else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                    {
                        reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStress))[1];
                    }

                    break;
                case Result.SIGMAZZ:
                    if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                    {
                        reExtra = elem.StressAt(gpsijk.location)[2];
                    }
                    else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                    {
                        reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStress))[2];
                    }
                    break;
                case Result.SIGMAXY:
                    if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                    {
                        reExtra = elem.StressAt(gpsijk.location)[3];
                    }
                    else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                    {
                        reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStress))[3];
                    }

                    break;
                case Result.SIGMAYZ:
                    if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                    {
                        reExtra = elem.StressAt(gpsijk.location)[4];
                    }
                    else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                    {
                        reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStress))[4];
                    }

                    break;
                case Result.SIGMAXZ:
                    if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                    {
                        reExtra = elem.StressAt(gpsijk.location)[5];
                    }
                    else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                    {
                        reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStress))[5];
                    }

                    break;
                case Result.EPSILONXX:
                    if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                    {
                        reExtra = elem.StrainAt(gpsijk.location)[0];
                    }
                    else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                    {
                        reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStrain))[0];
                    }
                    break;
                case Result.EPSILONYY:
                    if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                    {
                        reExtra = elem.StrainAt(gpsijk.location)[1];
                    }
                    else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                    {
                        reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStrain))[1];
                    }
                    break;
                case Result.EPSILONZZ:
                    if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                    {
                        reExtra = elem.StrainAt(gpsijk.location)[2];
                    }
                    else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                    {
                        reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStrain))[2];
                    }
                    break;
                case Result.EPSILONXY:
                    if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                    {
                        reExtra = elem.StrainAt(gpsijk.location)[3];
                    }
                    else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                    {
                        reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStrain))[3];
                    }
                    break;
                case Result.EPSILONYZ:
                    if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                    {
                        reExtra = elem.StrainAt(gpsijk.location)[4];
                    }
                    else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                    {
                        reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStrain))[4];
                    }
                    break;
                case Result.EPSILONXZ:
                    if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
                    {
                        reExtra = elem.StrainAt(gpsijk.location)[5];
                    }
                    else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
                    {
                        reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStrain))[5];
                    }
                    break;
                case Result.EPSILONTHXX:
                    reExtra = ((ElementThermoelastic3D)elem).StrainThermoAt(gpsijk.location)[0];
                    break;
                case Result.EPSILONTHYY:
                    reExtra = ((ElementThermoelastic3D)elem).StrainThermoAt(gpsijk.location)[1];
                    break;
                case Result.EPSILONTHZZ:
                    reExtra = ((ElementThermoelastic3D)elem).StrainThermoAt(gpsijk.location)[2];
                    break;
                case Result.EPSILONEXX:
                    reExtra = ((ElementThermoelastic3D)elem).StrainElasticAt(gpsijk.location)[0];
                    break;
                case Result.EPSILONEYY:
                    reExtra = ((ElementThermoelastic3D)elem).StrainElasticAt(gpsijk.location)[1];
                    break;
                case Result.EPSILONEZZ:
                    reExtra = ((ElementThermoelastic3D)elem).StrainElasticAt(gpsijk.location)[2];
                    break;
                case Result.EPSILONEXY:
                    reExtra = ((ElementThermoelastic3D)elem).StrainElasticAt(gpsijk.location)[3];
                    break;
                case Result.EPSILONEYZ:
                    reExtra = ((ElementThermoelastic3D)elem).StrainElasticAt(gpsijk.location)[4];
                    break;
                case Result.EPSILONEXZ:
                    reExtra = ((ElementThermoelastic3D)elem).StrainElasticAt(gpsijk.location)[5];
                    break;
                case Result.EPSILONPEQV:
                    reExtra = (double)gpsijk.GetValue(DataInGausspoint.lastAlpha);//gpsijk.lastAlpha;
                    break;
                default:
                    reExtra = 0;
                    break;
            }
            return reExtra;
        }

        //private void CreateAsBsPlane(Result re, AbstractPatch patch, ref DoubleMatrix As, ref DoubleVector Bs, int iel, object monitor)
        //{
        //  AbstractElementStructure2D elem = (AbstractElementStructure2D)patch.GetElement(iel);
        //  int nen = patch.GetCountLocalBasisFunctions(0);
        //  int[] econn = new int[nen];

        //  for (int k = 0; k < nen; k++)
        //  {
        //    econn[k] = ((AbstractPatch2D)patch).GetIEN(iel, k);
        //  }

        //  DoubleMatrix Ae = new DoubleMatrix(nen, nen);
        //  DoubleVector Be = new DoubleVector(nen);

        //  ////change order///////
        //  for (int j = 0; j < elem.GetNumberOfGaussPointOnEachDirection(); j++)
        //  {
        //    for (int i = 0; i < elem.GetNumberOfGaussPointOnEachDirection(); i++)
        //    {
        //      GaussPoints gpsij = elem.GetGaussPoint(i, j);
        //      DoubleVector Nii = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni); //new DoubleVector(nen);
        //                                                                            //for (int k = 0; k < nen; k++)
        //                                                                            //  Nii[k] = gpsij.Ni[k];
        //      Ae += MatrixFunctions.OuterProduct(Nii, Nii);
        //      double reExtra = 0;
        //      switch (re)
        //      {
        //        case Result.SIGMAXX:

        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[0];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[0];
        //          }
        //          break;
        //        case Result.SIGMAYY:

        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[1];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[1];
        //          }
        //          break;
        //        case Result.SIGMAZZ:
        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[2];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[2];
        //          }
        //          break;
        //        case Result.SIGMAXY:

        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[3];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[3];
        //          }
        //          break;
        //        case Result.EPSILONXX:
        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[0];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[0];
        //          }
        //          break;
        //        case Result.EPSILONYY:
        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[1];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[1];
        //          }
        //          break;
        //        case Result.EPSILONZZ:
        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[2];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[2];
        //          }
        //          break;
        //        case Result.EPSILONXY:
        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[3];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[3];
        //          }
        //          break;
        //        case Result.EPSILONTHXX:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainThermoAt(gpsij.location)[0];
        //          break;
        //        case Result.EPSILONTHYY:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainThermoAt(gpsij.location)[1];
        //          break;
        //        case Result.EPSILONTHZZ:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainThermoAt(gpsij.location)[2];
        //          break;
        //        case Result.EPSILONEXX:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainElasticAt(gpsij.location)[0];
        //          break;
        //        case Result.EPSILONEYY:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainElasticAt(gpsij.location)[1];
        //          break;
        //        case Result.EPSILONEZZ:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainElasticAt(gpsij.location)[2];
        //          break;
        //        case Result.EPSILONEXY:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainElasticAt(gpsij.location)[3];
        //          break;
        //        case Result.EPSILONPEQV:
        //          reExtra = (double)gpsij.GetValue(DataInGausspoint.lastAlpha);//gpsij.lastAlpha;
        //          break;
        //      }
        //      Be += reExtra * Nii;
        //    }
        //  }
        //  lock (monitor)
        //  {
        //    for (int jj = 0; jj < nen; jj++)
        //    {
        //      for (int ii = 0; ii < nen; ii++)
        //      {
        //        As[econn[ii], econn[jj]] += Ae[ii, jj];
        //      }
        //      Bs[econn[jj]] += Be[jj];
        //    }
        //  }
        //}

        //private void CreateAsBsPlane(Result re, AbstractPatch patch, ref DoubleCsrSparseMatrix As, ref DoubleVector Bs, int iel, object monitor)
        //{
        //  AbstractElementStructure2D elem = (AbstractElementStructure2D)patch.GetElement(iel);
        //  int nen = patch.GetCountLocalBasisFunctions(0);
        //  int[] econn = new int[nen];

        //  for (int k = 0; k < nen; k++)
        //  {
        //    econn[k] = ((AbstractPatch2D)patch).GetIEN(iel, k);
        //  }

        //  DoubleMatrix Ae = new DoubleMatrix(nen, nen);
        //  DoubleVector Be = new DoubleVector(nen);

        //  ////change order///////
        //  for (int j = 0; j < elem.GetNumberOfGaussPointOnEachDirection(); j++)
        //  {
        //    for (int i = 0; i < elem.GetNumberOfGaussPointOnEachDirection(); i++)
        //    {
        //      GaussPoints gpsij = elem.GetGaussPoint(i, j);
        //      DoubleVector Nii = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni); //new DoubleVector(nen);
        //                                                                            //for (int k = 0; k < nen; k++)
        //                                                                            //  Nii[k] = gpsij.Ni[k];
        //      Ae += MatrixFunctions.OuterProduct(Nii, Nii);
        //      double reExtra = 0;
        //      switch (re)
        //      {
        //        case Result.SIGMAXX:

        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[0];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[0];
        //          }
        //          break;
        //        case Result.SIGMAYY:

        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[1];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[1];
        //          }
        //          break;
        //        case Result.SIGMAZZ:
        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[2];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[2];
        //          }
        //          break;
        //        case Result.SIGMAXY:

        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StressAt(gpsij.location)[3];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStress))[3];
        //          }
        //          break;
        //        case Result.EPSILONXX:
        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[0];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[0];
        //          }
        //          break;
        //        case Result.EPSILONYY:
        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[1];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[1];
        //          }
        //          break;
        //        case Result.EPSILONZZ:
        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[2];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[2];
        //          }
        //          break;
        //        case Result.EPSILONXY:
        //          if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //          {
        //            reExtra = elem.StrainAt(gpsij.location)[3];
        //          }
        //          else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //          {
        //            reExtra = ((DoubleVector)gpsij.GetValue(DataInGausspoint.lastStrain))[3];
        //          }
        //          break;
        //        case Result.EPSILONTHXX:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainThermoAt(gpsij.location)[0];
        //          break;
        //        case Result.EPSILONTHYY:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainThermoAt(gpsij.location)[1];
        //          break;
        //        case Result.EPSILONTHZZ:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainThermoAt(gpsij.location)[2];
        //          break;
        //        case Result.EPSILONEXX:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainElasticAt(gpsij.location)[0];
        //          break;
        //        case Result.EPSILONEYY:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainElasticAt(gpsij.location)[1];
        //          break;
        //        case Result.EPSILONEZZ:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainElasticAt(gpsij.location)[2];
        //          break;
        //        case Result.EPSILONEXY:
        //          reExtra = ((ElementThermoelastic2D)elem).StrainElasticAt(gpsij.location)[3];
        //          break;
        //        case Result.EPSILONPEQV:
        //          reExtra = (double)gpsij.GetValue(DataInGausspoint.lastAlpha);//gpsij.lastAlpha;
        //          break;
        //      }
        //      Be += reExtra * Nii;
        //    }
        //  }
        //  lock (monitor)
        //  {
        //    for (int jj = 0; jj < nen; jj++)
        //    {
        //      for (int ii = 0; ii < nen; ii++)
        //      {
        //        As[econn[ii], econn[jj]] += Ae[ii, jj];
        //      }
        //      Bs[econn[jj]] += Be[jj];
        //    }
        //  }
        //}
        protected DoubleCsrSparseMatrix ComputeAsCsrSparseMatrix(AbstractPatch patch)
        {
            DoubleCsrSparseMatrix ma = null;
            SparseMatrixBuilder<double> AsBuilder = new SparseMatrixBuilder<double>();
            int nel = patch.CountElements();
            int nen = patch.GetCountLocalBasisFunctions(0);
            //if (!IsParallelProcesing)
            //{
            for (int j = 0; j < nel; j++)
            {
                int[] econn = new int[nen];
                for (int k = 0; k < nen; k++)
                {
                    econn[k] = ((AbstractPatchOneField)patch).GetIEN(j, k);
                }
                for (int jj = 0; jj < nen; jj++)
                {
                    for (int ii = 0; ii < nen; ii++)
                    {
                        IntPair key = new IntPair(econn[ii], econn[jj]);
                        //lock (monitor)
                        //{
                        if (!AsBuilder.ContainsKey(key))
                            AsBuilder.Add(key, 0);
                        //}
                    }
                }
            }
            //}
            //else
            //{
            //  var degreeOfParallelism = NumberOfCPUs;
            //  object monitor = new object();
            //  var tasks = new Task[degreeOfParallelism];
            //  for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
            //  {
            //    // capturing taskNumber in lambda wouldn't work correctly
            //    int taskNumberCopy = taskNumber;

            //    tasks[taskNumber] = Task.Factory.StartNew(
            //        () =>
            //        {
            //          var max = nel * (taskNumberCopy + 1) / degreeOfParallelism;
            //          for (int j = nel * taskNumberCopy / degreeOfParallelism; j < max; j++)
            //          {
            //            int[] econn = new int[nen];
            //            for (int k = 0; k < nen; k++)
            //            {
            //              econn[k] = ((AbstractPatch2D)patch).GetIEN(j, k);
            //            }
            //            for (int jj = 0; jj < nen; jj++)
            //            {
            //              for (int ii = 0; ii < nen; ii++)
            //              {
            //                IntPair key = new IntPair(econn[ii], econn[jj]);
            //                lock (monitor)
            //                {
            //                  if (!AsBuilder.ContainsKey(key))
            //                    AsBuilder.Add(key, 0);
            //                }
            //              }
            //            }
            //          }
            //        });
            //  }
            //  Task.WaitAll(tasks);
            //}
            int nnp1 = patch.GetCountGlobalBasisFunctions(0);
            ma = new DoubleCsrSparseMatrix(AsBuilder, nnp1);
            return ma;
        }
        //private void CreateAsBsSolid(Result re, int iji, int nen, DoubleMatrix Asxx, DoubleVector Bsxx, int iel, object monitor)
        //{
        //  int[] econn = new int[nen];
        //  AbstractPatch patch = listPatch[iji];
        //  for (int k = 0; k < nen; k++)
        //  {
        //    econn[k] = ((AbstractPatch3D)patch).GetIEN(iel, k);
        //  }
        //  DoubleMatrix Ae = new DoubleMatrix(nen, nen);
        //  DoubleVector Bexx = new DoubleVector(nen);

        //  AbstractElementStructure3D elem = (AbstractElementStructure3D)patch.GetElement(iel);
        //  for (int i = 0; i < elem.GetNumberOfGaussPointOnEachDirection(); i++)
        //  {
        //    for (int j = 0; j < elem.GetNumberOfGaussPointOnEachDirection(); j++)
        //    {
        //      for (int k = 0; k < elem.GetNumberOfGaussPointOnEachDirection(); k++)
        //      {
        //        GaussPoints gpsijk = elem.GetGaussPoint(i, j, k);
        //        DoubleVector Nii = (DoubleVector)gpsijk.GetValue(DataInGausspoint.Ni); //new DoubleVector(nen);

        //        //for (int kk = 0; kk < nen; kk++)
        //        //  Nii[kk] = gpsijk.Ni[kk];
        //        Ae += MatrixFunctions.OuterProduct(Nii, Nii);
        //        double reExtra = 0;
        //        switch (re)
        //        {
        //          case Result.SIGMAXX:
        //            if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //            {
        //              reExtra = elem.StressAt(gpsijk.location)[0];
        //            }
        //            else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //            {
        //              reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStress))[0];
        //            }
        //            break;

        //          case Result.SIGMAYY:
        //            if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //            {
        //              reExtra = elem.StressAt(gpsijk.location)[1];
        //            }
        //            else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //            {
        //              reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStress))[1];
        //            }

        //            break;
        //          case Result.SIGMAZZ:
        //            if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //            {
        //              reExtra = elem.StressAt(gpsijk.location)[2];
        //            }
        //            else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //            {
        //              reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStress))[2];
        //            }
        //            break;
        //          case Result.SIGMAXY:
        //            if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //            {
        //              reExtra = elem.StressAt(gpsijk.location)[3];
        //            }
        //            else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //            {
        //              reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStress))[3];
        //            }

        //            break;
        //          case Result.SIGMAYZ:
        //            if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //            {
        //              reExtra = elem.StressAt(gpsijk.location)[4];
        //            }
        //            else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //            {
        //              reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStress))[4];
        //            }

        //            break;
        //          case Result.SIGMAXZ:
        //            if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //            {
        //              reExtra = elem.StressAt(gpsijk.location)[5];
        //            }
        //            else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //            {
        //              reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStress))[5];
        //            }

        //            break;
        //          case Result.EPSILONXX:
        //            if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //            {
        //              reExtra = elem.StrainAt(gpsijk.location)[0];
        //            }
        //            else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //            {
        //              reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStrain))[0];
        //            }
        //            break;
        //          case Result.EPSILONYY:
        //            if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //            {
        //              reExtra = elem.StrainAt(gpsijk.location)[1];
        //            }
        //            else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //            {
        //              reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStrain))[1];
        //            }
        //            break;
        //          case Result.EPSILONZZ:
        //            if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //            {
        //              reExtra = elem.StrainAt(gpsijk.location)[2];
        //            }
        //            else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //            {
        //              reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStrain))[2];
        //            }
        //            break;
        //          case Result.EPSILONXY:
        //            if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //            {
        //              reExtra = elem.StrainAt(gpsijk.location)[3];
        //            }
        //            else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //            {
        //              reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStrain))[3];
        //            }
        //            break;
        //          case Result.EPSILONYZ:
        //            if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //            {
        //              reExtra = elem.StrainAt(gpsijk.location)[4];
        //            }
        //            else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //            {
        //              reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStrain))[4];
        //            }
        //            break;
        //          case Result.EPSILONXZ:
        //            if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Elasticity)
        //            {
        //              reExtra = elem.StrainAt(gpsijk.location)[5];
        //            }
        //            else if (patch.Material.TypeMaterialStructure == TypeMaterialStructure.Plasticity)
        //            {
        //              reExtra = ((DoubleVector)gpsijk.GetValue(DataInGausspoint.lastStrain))[5];
        //            }
        //            break;
        //          case Result.EPSILONTHXX:
        //            reExtra = ((ElementThermoelastic3D)elem).StrainThermoAt(gpsijk.location)[0];
        //            break;
        //          case Result.EPSILONTHYY:
        //            reExtra = ((ElementThermoelastic3D)elem).StrainThermoAt(gpsijk.location)[1];
        //            break;
        //          case Result.EPSILONTHZZ:
        //            reExtra = ((ElementThermoelastic3D)elem).StrainThermoAt(gpsijk.location)[2];
        //            break;
        //          case Result.EPSILONEXX:
        //            reExtra = ((ElementThermoelastic3D)elem).StrainElasticAt(gpsijk.location)[0];
        //            break;
        //          case Result.EPSILONEYY:
        //            reExtra = ((ElementThermoelastic3D)elem).StrainElasticAt(gpsijk.location)[1];
        //            break;
        //          case Result.EPSILONEZZ:
        //            reExtra = ((ElementThermoelastic3D)elem).StrainElasticAt(gpsijk.location)[2];
        //            break;
        //          case Result.EPSILONEXY:
        //            reExtra = ((ElementThermoelastic3D)elem).StrainElasticAt(gpsijk.location)[3];
        //            break;
        //          case Result.EPSILONEYZ:
        //            reExtra = ((ElementThermoelastic3D)elem).StrainElasticAt(gpsijk.location)[4];
        //            break;
        //          case Result.EPSILONEXZ:
        //            reExtra = ((ElementThermoelastic3D)elem).StrainElasticAt(gpsijk.location)[5];
        //            break;
        //          case Result.EPSILONPEQV:
        //            reExtra = (double)gpsijk.GetValue(DataInGausspoint.lastAlpha);//gpsijk.lastAlpha;
        //            break;
        //        }
        //        Bexx += reExtra * Nii;
        //      }
        //    }
        //  }

        //  lock (monitor)
        //  {
        //    for (int jj = 0; jj < nen; jj++)
        //    {
        //      for (int ii = 0; ii < nen; ii++)
        //      {
        //        Asxx[econn[ii], econn[jj]] += Ae[ii, jj];
        //      }
        //      Bsxx[econn[jj]] += Bexx[jj];
        //    }
        //  }
        //}
        /// <summary>
        /// 
        /// </summary>
        /// <param name="variable">only for currentStress and currentStrain</param>
        protected void ComputeValueAtGausspoints(DataInGausspoint variable)
        {
            foreach (AbstractPatchOneField patch in listPatch)
            {
                int countElement = patch.CountElements();
                for (int i = 0; i < countElement; i++)
                {
                    patch.GetElement(i).ComputeValueAtGaussPoint(variable);
                }
            }
        }
        public double[,] ComputeResultAtGausspoints(Result re, params int[] indexPatch)
        {
            List<double> x = new List<double>();
            List<double> y = new List<double>();
            List<double> z = new List<double>();
            List<double> val = new List<double>();
            if (indexPatch != null)
            {
                for (int i = 0; i < indexPatch.Length; i++)
                {
                    int index = indexPatch[i];
                    AbstractPatchOneField patch = (AbstractPatchOneField)listPatch[index];
                    patch.ComputeResultToGaussPoint(re, ref x, ref y, ref z, ref val);
                }
            }
            else
            {
                foreach (AbstractPatchOneField patch in listPatch)
                {
                    patch.ComputeResultToGaussPoint(re, ref x, ref y, ref z, ref val);
                }
            }
            double[,] matrix = new double[val.Count, 4];
            for (int i = 0; i < val.Count; i++)
            {
                matrix[i, 0] = x[i];
                matrix[i, 1] = y[i];
                matrix[i, 2] = z[i];
                matrix[i, 3] = val[i];
            }
            return matrix;
        }
        protected virtual void Enumerate()
        {
            int count = 0;
            if (IsMergeControlpoints1_1Constraint)
            {
                //Enumerate Merged Controlpoints
                int numberOfField = CountMaxNumberOfField();
                int dimension = -1;
                switch (StructureDimension)
                {
                    case Dimension.Beam:
                        dimension = 1;
                        break;
                    case Dimension.Plane:
                        dimension = 2;
                        break;
                    case Dimension.Solid:
                        dimension = 3;
                        break;
                }
                for (int i = 0; i < listControlPointCoincidented.Count; i++)
                {
                    listControlPointCoincidented[i].SetDimension(dimension);
                    listControlPointCoincidented[i].SetNumberOfFields(numberOfField);

                    int[] tArrayGlobal = new int[numberOfField];
                    for (int k = 0; k < numberOfField; k++)
                        tArrayGlobal[k] = count++;
                    listControlPointCoincidented[i].SetTArrayGlobal(tArrayGlobal);
                }
            }

            //Enumerate all patch
            foreach (AbstractPatchOneField patch in listPatch)
            {
                patch.EnumerateInPatch();
                count = patch.EnumerateInGlobalMultiPatch(count);
            }
            countDOF = count;
            //bool isMergeControlpoints1_1Constraint = true;
            //Modified coupled control point -2 into tArrayGlobal of coupleControlPoint
            foreach (CouplingTwoPatches con in listConnection)
            {
                if (con.Is1_1Constraint || IsMergeControlpoints1_1Constraint)
                {
                    var cpsMaster = con.GetControlPointMasterObjectInterface();
                    foreach (ControlPoint cp in cpsMaster)
                    {
                        ControlPoint couplingCp = cp.GetCoupleControlPoint();
                        if (couplingCp != null)
                        {
                            int[] tArrayMaster = couplingCp.GetTArrayGlobal();
                            int[] tArraySlave = cp.GetTArrayGlobal();
                            for (int i = 0; i < tArrayMaster.Length; i++)
                            {
                                if (i < tArraySlave.Length)//Master [ux uy T] Slave [ux uy]
                                {
                                    if (tArraySlave[i] == -2)
                                    {
                                        tArraySlave[i] = tArrayMaster[i];
                                    }
                                }
                            }
                            cp.SetTArrayGlobal(tArraySlave);
                        }
                    }



                    var cpsSlave = con.GetControlPointSlaveObjectInterface();
                    foreach (ControlPoint cp in cpsSlave)
                    {
                        ControlPoint couplingCp = cp.GetCoupleControlPoint();
                        if (couplingCp != null)
                        {
                            int[] tArrayMaster = couplingCp.GetTArrayGlobal();
                            int[] tArraySlave = cp.GetTArrayGlobal();
                            for (int i = 0; i < tArrayMaster.Length; i++)
                            {
                                if (i < tArraySlave.Length)//Master [ux uy T] Slave [ux uy]
                                {
                                    if (tArraySlave[i] == -2)
                                    {
                                        tArraySlave[i] = tArrayMaster[i];
                                    }
                                }
                            }
                            cp.SetTArrayGlobal(tArraySlave);
                        }
                    }
                }
            }

            for (int i = 0; i < listPatch.Count; i++)
            {
                AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
                int d = pa.GetCountField();
                int nen = pa.GetCountGlobalBasisFunctions();
                //if ((this is ModelStructure2D) || (this is ModelPiezoelectric2D))
                if (StructureDimension == Dimension.Beam)
                {
                    var cpsCurve = ((AbstractPatch1D)pa).GetCurve().ControlPoints;

                    for (int j = 0; j < d; j++)
                    {
                        if (!IsParallelProcesing)
                        {
                            for (int k = 0; k < nen; k++)
                            {
                                var idx = pa.GetINC(k, 0);
                                cpsCurve[idx].AddResult(listComputeResult[j], 0);
                                if (pa.GetIDInGlobal(j, k) == -2)
                                {
                                    pa.SetIDGlobal(cpsCurve[idx].GetTArrayGlobal()[j], j, k);
                                }
                            }
                        }
                        else
                        {
                            var degreeOfParallelism = NumberOfCPUs;
                            var tasks = new Task[degreeOfParallelism];
                            for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                            {
                                // capturing taskNumber in lambda wouldn't work correctly
                                int taskNumberCopy = taskNumber;

                                tasks[taskNumber] = Task.Factory.StartNew(
                                     () =>
                                     {
                                         var max = nen * (taskNumberCopy + 1) / degreeOfParallelism;
                                         for (int k = nen * taskNumberCopy / degreeOfParallelism; k < max; k++)
                                         {
                                             var idx1 = pa.GetINC(k, 0);
                                             cpsCurve[idx1].AddResult(listComputeResult[j], 0);
                                             if (pa.GetIDInGlobal(j, k) == -2)
                                             {
                                                 pa.SetIDGlobal(cpsCurve[idx1].GetTArrayGlobal()[j], j, k);
                                             }
                                         }
                                     });
                            }
                            Task.WaitAll(tasks);
                            //Parallel.For(0, nen, k =>
                            //{
                            //  var idx1 = pa.GetINC(k, 0);
                            //  var idx2 = pa.GetINC(k, 1);
                            //  cpsSurface[idx1, idx2].AddResult(listComputeResult[j], 0);
                            //  if (pa.GetIDInGlobal(j, k) == -2)
                            //  {
                            //    pa.SetIDGlobal(cpsSurface[idx1, idx2].GetTArrayGlobal()[j], j, k);
                            //  }
                            //});
                        }
                    }
                }
                else if (StructureDimension == Dimension.Plane)
                {
                    var cpsSurface = ((AbstractPatch2D)pa).GetSurface().ControlPoints;

                    for (int j = 0; j < d; j++)
                    {
                        if (!IsParallelProcesing)
                        {
                            for (int k = 0; k < nen; k++)
                            {
                                var idx1 = pa.GetINC(k, 0);
                                var idx2 = pa.GetINC(k, 1);
                                cpsSurface[idx1, idx2].AddResult(listComputeResult[j], 0);
                                if (pa.GetIDInGlobal(j, k) == -2)
                                {
                                    pa.SetIDGlobal(cpsSurface[idx1, idx2].GetTArrayGlobal()[j], j, k);
                                }
                            }
                        }
                        else
                        {
                            var degreeOfParallelism = NumberOfCPUs;
                            var tasks = new Task[degreeOfParallelism];
                            for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                            {
                                // capturing taskNumber in lambda wouldn't work correctly
                                int taskNumberCopy = taskNumber;

                                tasks[taskNumber] = Task.Factory.StartNew(
                                     () =>
                                     {
                                         var max = nen * (taskNumberCopy + 1) / degreeOfParallelism;
                                         for (int k = nen * taskNumberCopy / degreeOfParallelism; k < max; k++)
                                         {
                                             var idx1 = pa.GetINC(k, 0);
                                             var idx2 = pa.GetINC(k, 1);
                                             cpsSurface[idx1, idx2].AddResult(listComputeResult[j], 0);
                                             if (pa.GetIDInGlobal(j, k) == -2)
                                             {
                                                 pa.SetIDGlobal(cpsSurface[idx1, idx2].GetTArrayGlobal()[j], j, k);
                                             }
                                         }
                                     });
                            }
                            Task.WaitAll(tasks);
                            //Parallel.For(0, nen, k =>
                            //{
                            //  var idx1 = pa.GetINC(k, 0);
                            //  var idx2 = pa.GetINC(k, 1);
                            //  cpsSurface[idx1, idx2].AddResult(listComputeResult[j], 0);
                            //  if (pa.GetIDInGlobal(j, k) == -2)
                            //  {
                            //    pa.SetIDGlobal(cpsSurface[idx1, idx2].GetTArrayGlobal()[j], j, k);
                            //  }
                            //});
                        }
                    }
                }
                else if (StructureDimension == Dimension.Solid)
                {
                    var cpsVolume = ((AbstractPatch3D)pa).GetVolume().ControlPoints;

                    for (int j = 0; j < d; j++)
                    {
                        if (!IsParallelProcesing)
                        {
                            for (int k = 0; k < nen; k++)
                            {
                                var idx1 = pa.GetINC(k, 0);
                                var idx2 = pa.GetINC(k, 1);
                                var idx3 = pa.GetINC(k, 2);
                                cpsVolume[idx1, idx2, idx3].AddResult(listComputeResult[j], 0);
                                if (pa.GetIDInGlobal(j, k) == -2)
                                {
                                    pa.SetIDGlobal(cpsVolume[idx1, idx2, idx3].GetTArrayGlobal()[j], j, k);
                                }
                            }
                        }
                        else
                        {
                            var degreeOfParallelism = NumberOfCPUs;
                            var tasks = new Task[degreeOfParallelism];
                            for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                            {
                                // capturing taskNumber in lambda wouldn't work correctly
                                int taskNumberCopy = taskNumber;

                                tasks[taskNumber] = Task.Factory.StartNew(
                                     () =>
                                     {
                                         var max = nen * (taskNumberCopy + 1) / degreeOfParallelism;
                                         for (int k = nen * taskNumberCopy / degreeOfParallelism; k < max; k++)
                                         {
                                             var idx1 = pa.GetINC(k, 0);
                                             var idx2 = pa.GetINC(k, 1);
                                             var idx3 = pa.GetINC(k, 2);
                                             cpsVolume[idx1, idx2, idx3].AddResult(listComputeResult[j], 0);
                                             if (pa.GetIDInGlobal(j, k) == -2)
                                             {
                                                 pa.SetIDGlobal(cpsVolume[idx1, idx2, idx3].GetTArrayGlobal()[j], j, k);
                                             }
                                         }
                                     });
                            }
                            Task.WaitAll(tasks);
                            //Parallel.For(0, nen, k =>
                            //{
                            //  var idx1 = pa.GetINC(k, 0);
                            //  var idx2 = pa.GetINC(k, 1);
                            //  var idx3 = pa.GetINC(k, 2);
                            //  cpsVolume[idx1, idx2, idx3].AddResult(listComputeResult[j], 0);
                            //  if (pa.GetIDInGlobal(j, k) == -2)
                            //  {
                            //    pa.SetIDGlobal(cpsVolume[idx1, idx2, idx3].GetTArrayGlobal()[j], j, k);
                            //  }
                            //});
                        }
                    }
                }
            }
        }

        public void AssemblyStiffnessMatrix(ref DoubleCsrSparseMatrix kGlobal)
        {
            foreach (AbstractPatch patch in listPatch)
                patch.ComputeStiffnessMatrixPatch(ref kGlobal);//serial  7.4 ---32,   parallel 16,
        }
        public void AssemblyStiffnessMatrix(out DoubleMatrix kGlobal)
        {
            kGlobal = new DoubleMatrix(countDOF, countDOF);
            foreach (AbstractPatch patch in listPatch)
                patch.ComputeStiffnessMatrixPatch(ref kGlobal);
        }

        public void AssemblyMassMatrix(ref DoubleCsrSparseMatrix mGlobal)
        {
            foreach (AbstractPatch patch in listPatch)
                patch.ComputeMassMatrixPatch(ref mGlobal);//serial  7.4 ---32,   parallel 16,
        }
        public void AssemblyMassMatrix(out DoubleMatrix mGlobal)
        {
            mGlobal = new DoubleMatrix(countDOF, countDOF);
            for (int i = 0; i < listPatch.Count; i++)
            {
                listPatch[i].ComputeMassMatrixPatch(ref mGlobal);
            }
        }

        public void AssemblyTractionVector(out DoubleVector rGlobal, double time = 1)
        {
            //if (IsSparseData)
            //    rGlobal = DoubleVector.Build.Sparse(countDOF);
            //else
            rGlobal = new DoubleVector(countDOF);
            for (int i = 0; i < CountLoad(); i++)
            {
                AbstractLoad load = loads[i];
                var mp = load.GetMeshPart();
                int[] tArrayGlobal = mp.GetTArrayGlobal();
                DoubleVector rLocal = load.ComputeLocalLoadVector(time);
                for (int j = 0; j < rLocal.Length; j++)
                {
                    int n = -1;
                    if (this is ModelThermoelasticStatic && load is HeatFlux)
                    {
                        switch (StructureDimension)
                        {
                            case Dimension.Beam:
                                n = tArrayGlobal[2 * j + 1];
                                break;
                            case Dimension.Plane:
                                n = tArrayGlobal[3 * j + 2];
                                break;
                            case Dimension.Solid:
                                n = tArrayGlobal[4 * j + 3];
                                break;
                        }
                    }
                    else if (this is ModelPiezoelectricStatic || this is ModelThermoelasticStatic || this is ModelStructurePhaseFieldStatic)
                    {
                        switch (StructureDimension)
                        {
                            case Dimension.Beam:
                                //rLocal[4]
                                //0  1     2  3     4  
                                //0  1  2  3  4  5  6  7
                                //u  v  T  u  v  T
                                n = tArrayGlobal[j + (int)(j)];
                                break;
                            case Dimension.Plane:
                                //rLocal[8]
                                //0  1     2  3     4  5     6  7
                                //0  1  2  3  4  5  6  7  8  9  10  11
                                //u  v  T  u  v  T
                                n = tArrayGlobal[j + (int)(j / 2)];
                                break;
                            case Dimension.Solid:
                                //rLocal[12]
                                //0  1  2     3  4  5     6  7  8       9   10  11
                                //0  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15
                                //u  v  w  T
                                n = tArrayGlobal[j + (int)(j / 3)];
                                break;
                        }
                    }
                    else
                    {
                        n = tArrayGlobal[j];
                    }
                    if (n != -1)
                    {
                        rGlobal[n] += rLocal[j];
                    }
                }
                //////////////////////////////////////////////////////////////////////
                //if (this is ModelThermoelasticStatic && load is HeatFlux)
                //{
                //    switch (StructureDimension)
                //    {
                //        case Dimension.Plane:
                //            for (int j = 0; j < rLocal.Count; j++)
                //            {
                //                int n = tArrayGlobal[3 * j + 2];
                //                if (n != -1)
                //                {
                //                    rGlobal[n] += rLocal[j];
                //                }
                //            }
                //            break;
                //        case Dimension.Solid:
                //            for (int j = 0; j < rLocal.Count; j++)
                //            {
                //                int n = tArrayGlobal[4 * j + 3];
                //                if (n != -1)
                //                {
                //                    rGlobal[n] += rLocal[j];
                //                }
                //            }
                //            break;
                //    }
                //}
                //else if (this is ModelPiezoelectricStatic)
                //{
                //    switch (StructureDimension)
                //    {
                //        case Dimension.Plane:
                //            for (int j = 0; j < rLocal.Count; j++)
                //            {
                //                int n = tArrayGlobal[j + (int)(j / 2)];
                //                if (n != -1)
                //                {
                //                    rGlobal[n] += rLocal[j];
                //                }
                //            }
                //            break;
                //        case Dimension.Solid:
                //            for (int j = 0; j < rLocal.Count; j++)
                //            {
                //                int n = tArrayGlobal[j + (int)(j / 3)];
                //                if (n != -1)
                //                {
                //                    rGlobal[n] += rLocal[j];
                //                }
                //            }
                //            break;
                //    }
                //}
                //else
                //{
                //    for (int j = 0; j < rLocal.Count; j++)
                //    {
                //        int n = tArrayGlobal[j];
                //        if (n != -1)
                //        {
                //            rGlobal[n] += rLocal[j];
                //        }
                //    }
                //}
                ////////////////////////////////////////////////////////////////
            }
        }

        /// <summary>
        /// Set interface of two patch
        /// </summary>
        /// <param name="masterPatch">master patch</param>
        /// <param name="slavePatch">slave patch</param>
        /// <param name="indexMasterObjectInterface">index edge of master patch</param>
        /// <param name="indexSlaveObjectInterface">index edge of slave patch</param>
        public void SetInterfaceBetweenTwoPatches(int indexMasterPatch, int indexSlavePatch, int indexMasterObjectInterface, int indexSlaveObjectInterface, int[] DOFUnCoupling, bool is1_1Constraint, bool isUseC1Constraint, ref List<ControlPoint> listControlPointCoincidented)
        {
            CouplingTwoPatches con = new CouplingTwoPatches(listPatch[indexMasterPatch], listPatch[indexSlavePatch], indexMasterObjectInterface, indexSlaveObjectInterface, DOFUnCoupling, is1_1Constraint, isUseC1Constraint, ref listControlPointCoincidented, IsMergeControlpoints1_1Constraint);
            if (listIndexConnectionPatch == null)
                listIndexConnectionPatch = new List<int[]>();
            if (listConnection == null)
                listConnection = new List<ICoupling>();
            listConnection.Add(con);
            listIndexConnectionPatch.Add(new int[] { indexMasterPatch, indexSlavePatch });
            if (listPatch[indexMasterPatch] is AbstractPatch2D && listPatch[indexSlavePatch] is AbstractPatch2D)
            {
                AbstractPatch2D masterPatch = (AbstractPatch2D)listPatch[indexMasterPatch];
                AbstractPatch2D slavePatch = (AbstractPatch2D)listPatch[indexSlavePatch];
                var cpsMaster = masterPatch.SelectEndPatchControlPoints(indexMasterObjectInterface);
                var cpsSlave = slavePatch.SelectEndPatchControlPoints(indexSlaveObjectInterface);
                //////Coupling meet points
                if (!IsMergeControlpoints1_1Constraint)
                {
                    if (cpsMaster[0].IsCoincident(cpsSlave[0]) && cpsMaster[cpsMaster.Length - 1].IsCoincident(cpsSlave[cpsSlave.Length - 1]))
                    {
                        SetCouplingTwoControlPoint(cpsMaster[0], cpsSlave[0], indexMasterPatch, indexSlavePatch);
                        SetCouplingTwoControlPoint(cpsMaster[cpsMaster.Length - 1], cpsSlave[cpsSlave.Length - 1], indexMasterPatch, indexSlavePatch);
                    }
                    else if (cpsMaster[0].IsCoincident(cpsSlave[cpsSlave.Length - 1]) && cpsMaster[cpsMaster.Length - 1].IsCoincident(cpsSlave[0]))
                    {
                        SetCouplingTwoControlPoint(cpsMaster[0], cpsSlave[cpsSlave.Length - 1], indexMasterPatch, indexSlavePatch);
                        SetCouplingTwoControlPoint(cpsMaster[cpsMaster.Length - 1], cpsSlave[0], indexMasterPatch, indexSlavePatch);
                    }
                }
            }
        }
        private void SetCouplingTwoControlPoint(ControlPoint masterCp, ControlPoint slaveCp, int indexMasterPatch, int indexSlavePatch)
        {
            if (listControlPointCoupling == null)
            {
                listControlPointCoupling = new List<List<ControlPoint>>();
            }
            if (listIndexPatchCoupling == null)
            {
                listIndexPatchCoupling = new List<List<int>>();
            }
            bool isFoundCoincidentControlPoint = false;
            for (int i = 0; i < listControlPointCoupling.Count; i++)
            {
                if (listControlPointCoupling[i][0].IsCoincident(masterCp))
                {
                    bool isFoundMasterCp = false;
                    bool isFoundSlaveCp = false;
                    for (int j = 0; j < listControlPointCoupling[i].Count; j++)
                    {
                        if (ReferenceEquals(listControlPointCoupling[i][j], masterCp))
                            isFoundMasterCp = true;
                        if (ReferenceEquals(listControlPointCoupling[i][j], slaveCp))
                            isFoundSlaveCp = true;
                    }
                    if (!isFoundMasterCp)
                    {
                        listControlPointCoupling[i].Add(masterCp);
                        listIndexPatchCoupling[i].Add(indexMasterPatch);
                    }

                    if (!isFoundSlaveCp)
                    {
                        listControlPointCoupling[i].Add(slaveCp);
                        listIndexPatchCoupling[i].Add(indexSlavePatch);
                    }
                    isFoundCoincidentControlPoint = true;
                    break;
                }
            }
            if (!isFoundCoincidentControlPoint)
            {
                List<ControlPoint> l = new List<ControlPoint>();
                l.Add(masterCp);
                l.Add(slaveCp);
                List<int> ll = new List<int>();
                ll.Add(indexMasterPatch);
                ll.Add(indexSlavePatch);
                listControlPointCoupling.Add(l);
                listIndexPatchCoupling.Add(ll);
            }
        }
        private void DetachedCoincidentUncoupling()
        {
            if (listIndexPatchCoupling != null)
                for (int i = listIndexPatchCoupling.Count - 1; i > 0; i--)
                {
                    if (listIndexPatchCoupling[i].Count >= 4)
                    {
                        List<List<int>> tempIndex = new List<List<int>>();
                        List<int> tempListIndex = listIndexPatchCoupling[i];
                        for (int j = 0; j < tempListIndex.Count; j++)
                        {
                            bool isNew = true;
                            for (int k = 0; k < tempIndex.Count; k++)
                            {
                                bool isBreak = false;
                                for (int kk = 0; kk < tempIndex[k].Count; kk++)
                                {
                                    if (IsHasCouplingInList(tempListIndex[j], tempIndex[k][kk]))
                                    {
                                        tempIndex[k].Add(tempListIndex[j]);
                                        isNew = false;
                                        isBreak = true;
                                        break;
                                    }
                                }
                                if (isBreak)
                                    break;
                            }
                            if (isNew)
                            {
                                List<int> newList = new List<int>();
                                newList.Add(tempListIndex[j]);
                                tempIndex.Add(newList);
                            }
                        }

                        //Find any connect two groups
                        //for (int j = tempIndex.Count - 1; j >= 1; j--)
                        //{
                        //  for (int k = tempIndex[j].Count - 1; k >= 0; k--)
                        //  {

                        //  }
                        //  tempIndex
                        //}

                        if (tempIndex.Count > 1)
                        {
                            //List<List<ControlPoint>> tempCps = new List<List<ControlPoint>>();
                            List<ControlPoint> tempListCps = listControlPointCoupling[i];
                            listControlPointCoupling.RemoveAt(i);
                            listIndexPatchCoupling.RemoveAt(i);
                            //tempListIndex
                            for (int j = 0; j < tempIndex.Count; j++)
                            {
                                List<ControlPoint> lis = new List<ControlPoint>();
                                for (int k = 0; k < tempIndex[j].Count; k++)
                                {
                                    int idx = tempListIndex.IndexOf(tempIndex[j][k]);
                                    lis.Add(tempListCps[idx]);
                                }
                                listControlPointCoupling.Add(lis);
                            }
                        }
                    }
                }
        }
        private bool IsHasCouplingInList(int index1, int index2)
        {
            foreach (int[] a in listIndexConnectionPatch)
            {
                if ((index1 == a[0] && index2 == a[1]) || (index1 == a[1] && index2 == a[0]))
                    return true;
            }
            return false;
        }

        //public void SetAutomaticInterfaceAllPatch(bool Is1_1Constraint, bool IsUseC1Constraint)
        //{
        //  var auto = new GenerateAutomaticConnectionInterfacePatches(this);
        //  auto.GenerateConnectionInterfaceAuto(Is1_1Constraint, IsUseC1Constraint);
        //}
        List<ControlPoint> listControlPointCoincidented;
        public void SetAutomaticInterfaceAllPatch(bool Is1_1Constraint = true, bool IsUseC1Constraint = false, int master = -1, int directionMaster = -1, params int[] unCouplingPatches)
        {
            var auto = new GenerateAutomaticConnectionInterfacePatches(this);
            List<int[]> listConnectionPatch = null;
            string filenameSaveConnectionPatch = PathProject + "restore_connection_patch.msgpack";
            string directory = Path.GetDirectoryName(filenameSaveConnectionPatch);
            if (!Directory.Exists(directory) || !File.Exists(filenameSaveConnectionPatch))
            {
                IntPair[] pairs = null;
                if (unCouplingPatches != null)
                {
                    if (unCouplingPatches.Length % 2 == 0)
                    {
                        pairs = new IntPair[unCouplingPatches.Length / 2];
                        int c = 0;
                        for (int i = 0; i < unCouplingPatches.Length; i = i + 2)
                        {
                            pairs[c++] = new IntPair(unCouplingPatches[i], unCouplingPatches[i + 1]);
                        }
                    }
                    else
                    {
                        throw new ArgumentException("Uncoupling Patches must be in pair");
                    }
                }

                listConnectionPatch = auto.GenerateInformationConnectionInterfaceAuto(Is1_1Constraint, pairs);
                if (listConnectionPatch.Count > 0)
                {
                    if (master != -1 && directionMaster != -1)
                    {
                        listConnectionPatch = IOIGA.OptimizePair(listConnectionPatch, master, directionMaster);
                    }
                    //else
                    //{
                    //  listConnectionPatch = IOIGA.ReArrangePair(listConnectionPatch);
                    //}  
                    int[,] arrayConnectionPatch = new int[listConnectionPatch.Count, listConnectionPatch[0].Length];
                    for (int i = 0; i < listConnectionPatch.Count; i++)
                    {
                        for (int j = 0; j < listConnectionPatch[0].Length; j++)
                        {
                            arrayConnectionPatch[i, j] = listConnectionPatch[i][j];
                        }
                    }
                    //IOIGA.WriteMessagePackArrayData(filenameSaveConnectionPatch, listConnectionPatch);
                    IOIGA.WriteArrayMessagePackFormatter(filenameSaveConnectionPatch, arrayConnectionPatch);
                }
            }
            else
            {
                int[,] temp = IOIGA.ReadArrayIntMessagePackFormatter(filenameSaveConnectionPatch);
                listConnectionPatch = new List<int[]>();
                for (int i = 0; i < temp.GetLength(0); i++)
                {
                    int[] subarray = new int[temp.GetLength(1)];
                    for (int j = 0; j < temp.GetLength(1); j++)
                    {
                        subarray[j] = temp[i, j];
                    }
                    listConnectionPatch.Add(subarray);
                }
            }
            listControlPointCoincidented = new List<ControlPoint>();
            auto.GenerateConnectionInterfaceAuto(listConnectionPatch, Is1_1Constraint, IsUseC1Constraint, ref listControlPointCoincidented);
        }

        //public void MergeMatchedControlPoint()
        //{
        //  string filenameSaveConnectionPatch = PathProject + "restore_connection_patch.msgpack";
        //  string directory = Path.GetDirectoryName(filenameSaveConnectionPatch);
        //  if (Directory.Exists(directory) && File.Exists(filenameSaveConnectionPatch))
        //  {
        //    //listConnection
        //  }
        //}
        public string PathProject
        { get; set; }

        public string NameProject
        { get; set; }

        internal void SetUGlobal(DoubleVector uGlobal)
        {
            for (int i = 0; i < listPatch.Count; i++)
            {
                AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
                ControlPoint[] cps = pa.GetAllControlPoints();
                if (!IsParallelProcesing)
                {
                    for (int j = 0; j < cps.Length; j++)
                    {
                        ControlPoint currentCps = cps[j];
                        for (int k = 0; k < pa.GetCountField(); k++)
                        {
                            int tArray = currentCps.GetTArrayGlobal()[k];
                            if (tArray != -1)
                            {
                                currentCps.SetULocal(k, uGlobal[tArray]);
                                currentCps.AddResult(listComputeResult[k], uGlobal[tArray]);
                            }
                        }
                    }
                }
                else
                {
                    var degreeOfParallelism = NumberOfCPUs;
                    var tasks = new Task[degreeOfParallelism];
                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = cps.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int j = cps.Length * taskNumberCopy / degreeOfParallelism; j < max; j++)
                                 {
                                     ControlPoint currentCps = cps[j];
                                     for (int k = 0; k < pa.GetCountField(); k++)
                                     {
                                         int tArray = currentCps.GetTArrayGlobal()[k];
                                         if (tArray != -1)
                                         {
                                             currentCps.SetULocal(k, uGlobal[tArray]);
                                             currentCps.AddResult(listComputeResult[k], uGlobal[tArray]);
                                         }
                                     }
                                 }
                             });
                    }
                    Task.WaitAll(tasks);

                    //Parallel.For(0, cps.Length, j =>
                    //{
                    //  ControlPoint currentCps = cps[j];
                    //  for (int k = 0; k < pa.GetCountField(); k++)
                    //  {
                    //    int tArray = currentCps.GetTArrayGlobal()[k];
                    //    if (tArray != -1)
                    //    {
                    //      currentCps.SetULocal(k, uGlobal[tArray]);
                    //      currentCps.AddResult(listComputeResult[k], uGlobal[tArray]);
                    //    }
                    //  }
                    //});
                }
            }
        }

        /// <summary>
        /// Draw initial geometry
        /// </summary>
        /// <param name="viewer">viewer form</param>
        /// <param name="resolution">resolution on a element</param>
        public void DrawGeometry(ViewerForm viewer, int resolution, double scale = 1.0)
        {
            switch (StructureDimension)
            {
                case Dimension.Beam:
                    foreach (AbstractPatch1D patch in listPatch)
                    {
                        NURBSCurve curve = null;
                        curve = ((PatchStructureBeam)patch).GetDeformationCurve(scale);
                        //}
                        curve.Draw(viewer);
                    }
                    break;
                case Dimension.Plate:
                case Dimension.Plane:
                    foreach (AbstractPatch2D patch in listPatch)
                    {
                        NURBSSurface surface = null;
                        //if (TypeModel == TypeModelProblem.Structural || TypeModel == TypeModelProblem.PhaseField || TypeModel == TypeModelProblem.StructuralThermal || TypeModel == TypeModelProblem.Piezoelectric)
                        //{
                        surface = ((PatchStructure2D)patch).GetDeformationSurface(scale);
                        //}
                        surface.Draw(viewer);
                    }
                    break;
                case Dimension.Solid:
                    foreach (AbstractPatch3D patch in listPatch)
                    {
                        NURBSVolume vol = null;
                        //if (TypeModel == TypeModelProblem.Structural)
                        //{
                        vol = ((PatchStructure3D)patch).GetDeformationVolume(scale);
                        //}
                        vol.Draw(viewer);
                    }
                    break;
            }
        }

        public void ExtrapolationMaterialProperty(MaterialPropertyName name)
        {
            for (int iji = 0; iji < listPatch.Count; iji++)
            {
                AbstractPatch patch = listPatch[iji];
                int nen = patch.GetCountLocalBasisFunctions(0); //(p + 1) * (q + 1);
                int nnp1 = patch.GetCountGlobalBasisFunctions(0);
                SparseMatrixBuilder<double> AsBuilder = new SparseMatrixBuilder<double>();
                DoubleCsrSparseMatrix AsSparse = null;// ComputeAsCsrSparseMatrix(patch);
                DoubleMatrix As = null;// new DoubleMatrix(nnp1, nnp1);
                if (!IsSparseData)
                    As = new DoubleMatrix(nnp1, nnp1);
                else
                    AsSparse = ComputeAsCsrSparseMatrix(patch);
                DoubleVector Bs = new DoubleVector(nnp1);
                for (int iel = 0; iel < patch.CountElements(); iel++)
                {
                    int[] econn = new int[nen];

                    for (int k = 0; k < nen; k++)
                    {
                        econn[k] = patch.GetIEN(0, iel, k);
                    }

                    DoubleMatrix Ae = new DoubleMatrix(nen, nen);
                    DoubleVector Be = new DoubleVector(nen);
                    if (patch is AbstractPatch2D)
                    {
                        AbstractElement2D elem = (AbstractElement2D)patch.GetElement(iel);
                        elem.ComputeMaterialPropertyValueAtGaussPoint(name);

                        ////change order///////
                        for (int j = 0; j < elem.GetNumberOfGaussPointOnEachDirection(); j++)
                        {
                            for (int i = 0; i < elem.GetNumberOfGaussPointOnEachDirection(); i++)
                            {
                                GaussPoints gpsij = elem.GetGaussPoint(i, j);
                                DoubleVector Nii = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);//new DoubleVector(nen);

                                Ae += MatrixFunctions.OuterProduct(Nii, Nii);
                                Be += (double)gpsij.GetValue(DataInGausspoint.MaterialPropertyValue) * Nii;
                            }
                        }
                    }
                    else if (patch is AbstractPatch3D)
                    {
                        AbstractElement3D elem = (AbstractElement3D)patch.GetElement(iel);
                        elem.ComputeMaterialPropertyValueAtGaussPoint(name);

                        for (int i = 0; i < elem.GetNumberOfGaussPointOnEachDirection(); i++)
                        {
                            for (int j = 0; j < elem.GetNumberOfGaussPointOnEachDirection(); j++)
                            {
                                for (int k = 0; k < elem.GetNumberOfGaussPointOnEachDirection(); k++)
                                {
                                    GaussPoints gpsijk = elem.GetGaussPoint(i, j, k);
                                    DoubleVector Nii = (DoubleVector)gpsijk.GetValue(DataInGausspoint.Ni); //new DoubleVector(nen);
                                    Ae += MatrixFunctions.OuterProduct(Nii, Nii);
                                    Be += (double)gpsijk.GetValue(DataInGausspoint.MaterialPropertyValue) * Nii;
                                }
                            }
                        }
                    }
                    for (int jj = 0; jj < nen; jj++)
                    {
                        for (int ii = 0; ii < nen; ii++)
                        {
                            if (!IsSparseData)
                                As[econn[ii], econn[jj]] += Ae[ii, jj];
                            else
                                AsSparse[econn[ii], econn[jj]] += Ae[ii, jj];
                        }
                        Bs[econn[jj]] += Be[jj];
                    }
                }
                DoubleVector Ps = (!IsSparseData) ? MatrixFunctions.Solve(As, Bs) : MatrixFunctions.Solve(AsSparse, Bs);
                if (patch is AbstractPatch2D)
                {
                    for (int ii = 0; ii < nnp1; ii++)
                    {
                        int ik = patch.GetINC(0, ii, 0);
                        int jk = patch.GetINC(0, ii, 1);
                        ((NURBSSurface)patch.GetGeometry(0)).ControlPoints[ik, jk].MaterialPropertyValue = Ps[ii];
                    }
                }
                else if (patch is AbstractPatch3D)
                {
                    for (int ii = 0; ii < nnp1; ii++)
                    {
                        int ik = patch.GetINC(0, ii, 0);
                        int jk = patch.GetINC(0, ii, 1);
                        int kk = patch.GetINC(0, ii, 2);

                        ((NURBSVolume)patch.GetGeometry(0)).ControlPoints[ik, jk, kk].MaterialPropertyValue = Ps[ii];
                    }
                }
            }
        }

        public void DrawMaterialDistributionAtGaussPoints(ViewerForm viewer, MaterialPropertyName name)
        {
            ExtrapolationMaterialProperty(name);

            List<double> x = new List<double>();
            List<double> y = new List<double>();
            List<double> z = new List<double>();
            List<double> val = new List<double>();
            PointSetColor ps;
            if (StructureDimension == Dimension.Plane)
            {
                foreach (AbstractPatch2D patch in listPatch)
                {
                    var surface = patch.GetSurface();
                    surface.isDrawControlPoint = false;
                    surface.isDrawControlNet = false;
                    surface.isColorfulFace = false;
                    surface.isDrawSurface = false;
                    surface.isDrawCurve = false;
                    surface.isDrawKnot = false;
                    surface.colorCurve = Color.Black;
                    surface.widthCurve = 1;
                    surface.resolution1 = 5;
                    surface.resolution2 = 5;
                    //surface.opacity = 0.2;
                    surface.Draw(viewer);
                    patch.ComputeDataDrawMaterialDistribution(ref x, ref y, ref z, ref val);
                }
            }
            else if (StructureDimension == Dimension.Solid)
            {
                foreach (AbstractPatch3D patch in listPatch)
                {
                    var vol = patch.GetVolume();
                    vol.isDrawGrid = true;
                    vol.isDrawVolume = false;
                    vol.isDrawControlPoint = false;
                    vol.isDrawControlNet = false;
                    vol.isColorfulFace = false;
                    vol.isDrawCurve = false;
                    vol.isDrawKnot = false;
                    vol.colorGrid = Color.Black;
                    vol.widthGrid = 1;
                    vol.resolution1 = 5;
                    vol.resolution2 = 5;
                    vol.resolution3 = 5;
                    //surface.opacity = 0.2;
                    vol.Draw(viewer);
                    patch.ComputeDataDrawMaterialDistribution(ref x, ref y, ref z, ref val);
                }
            }
            ps = new PointSetColor(x.ToArray(), y.ToArray(), z.ToArray(), val.ToArray());
            ps.SetPointSize(10);
            ps.SetOpacity(0.6);
            viewer.AddObject3D(ps);
            viewer.SetColormapBarVisible(ps, name.ToString());
        }
        public void DrawMaterialDistributionOnPatch(ViewerForm viewer, int resolution, MaterialPropertyName name)
        {
            ExtrapolationMaterialProperty(name);

            if (StructureDimension == Dimension.Plane)
            {
                Contours2D[] contours = new Contours2D[listPatch.Count];
                double minVal = 1e100;
                double maxVal = -1e100;

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    double[,][] data = ((AbstractPatch2D)listPatch[lp]).GetSurface().GetDataOnSurface(resolution, resolution);
                    double[,] xData = new double[data.GetLength(0), data.GetLength(1)];
                    double[,] yData = new double[data.GetLength(0), data.GetLength(1)];
                    double[,] zData = new double[data.GetLength(0), data.GetLength(1)];

                    for (int j = 0; j < data.GetLength(1); j++)
                    {
                        for (int i = 0; i < data.GetLength(0); i++)
                        {
                            xData[i, j] = data[i, j][0];
                            yData[i, j] = data[i, j][1];
                            zData[i, j] = data[i, j][2];
                        }
                    }
                    contours[lp] = new Contours2D(xData, yData, zData);
                    contours[lp].SetScalarValue(((AbstractPatch2D)listPatch[lp]).ComputeMaterialPropertyValue(resolution));
                    contours[lp].SetOpacity(10);
                    contours[lp].ColorType = colorType;

                    var val = contours[lp].GetValue();
                    var miVal = val.Min();
                    var maVal = val.Max();
                    if (miVal < minVal)
                        minVal = miVal;
                    if (maVal > maxVal)
                        maxVal = maVal;
                }

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    contours[lp].SetMinMaxValue(minVal, maxVal);
                    contours[lp].Update(minVal, maxVal);
                    viewer.AddObject3D(contours[lp]);
                }
                viewer.SetColormapBarVisible(contours[0], name.ToString());
            }
            else if (StructureDimension == Dimension.Solid)
            {
                Contours3D[] contours = new Contours3D[listPatch.Count];
                double minVal = 1e100;
                double maxVal = -1e100;

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    double[,,][] data = ((AbstractPatch3D)listPatch[lp]).GetVolume().GetDataOnVolume(resolution, resolution, resolution);

                    double[,,] xData = new double[data.GetLength(0), data.GetLength(1), data.GetLength(2)];
                    double[,,] yData = new double[data.GetLength(0), data.GetLength(1), data.GetLength(2)];
                    double[,,] zData = new double[data.GetLength(0), data.GetLength(1), data.GetLength(2)];

                    for (int k = 0; k < data.GetLength(2); k++)
                        for (int j = 0; j < data.GetLength(1); j++)
                            for (int i = 0; i < data.GetLength(0); i++)
                            {
                                xData[i, j, k] = data[i, j, k][0];
                                yData[i, j, k] = data[i, j, k][1];
                                zData[i, j, k] = data[i, j, k][2];
                            }
                    contours[lp] = new Contours3D(xData, yData, zData);
                    contours[lp].SetScalarValue(((AbstractPatch3D)listPatch[lp]).ComputeMaterialPropertyValue(resolution));
                    contours[lp].SetOpacity(10);
                    contours[lp].ColorType = colorType;

                    var val = contours[lp].GetValue();
                    var miVal = val.Min();
                    var maVal = val.Max();
                    if (miVal < minVal)
                        minVal = miVal;
                    if (maVal > maxVal)
                        maxVal = maVal;
                }

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    contours[lp].SetMinMaxValue(minVal, maxVal);
                    contours[lp].Update(minVal, maxVal);
                    viewer.AddObject3D(contours[lp]);
                }
                viewer.SetColormapBarVisible(contours[0], name.ToString());
            }
        }
        private void ExtrapolationValueGausspoint(DataInGausspoint name)
        {
            for (int iji = 0; iji < listPatch.Count; iji++)
            {
                AbstractPatch patch = listPatch[iji];
                int nen = patch.GetCountLocalBasisFunctions(0); //(p + 1) * (q + 1);
                int nnp1 = patch.GetCountGlobalBasisFunctions(0);
                SparseMatrixBuilder<double> AsBuilder = new SparseMatrixBuilder<double>();
                DoubleCsrSparseMatrix AsSparse = null;// ComputeAsCsrSparseMatrix(patch);
                DoubleMatrix As = null;// new DoubleMatrix(nnp1, nnp1);
                if (!IsSparseData)
                    As = new DoubleMatrix(nnp1, nnp1);
                else
                    AsSparse = ComputeAsCsrSparseMatrix(patch);
                DoubleVector Bs = new DoubleVector(nnp1);
                for (int iel = 0; iel < patch.CountElements(); iel++)
                {
                    int[] econn = new int[nen];

                    for (int k = 0; k < nen; k++)
                    {
                        econn[k] = ((AbstractPatch2D)patch).GetIEN(iel, k);
                    }

                    DoubleMatrix Ae = new DoubleMatrix(nen, nen);
                    DoubleVector Be = new DoubleVector(nen);
                    if (listPatch[iji] is AbstractPatch2D)
                    {
                        AbstractElement2D elem = (AbstractElement2D)listPatch[iji].GetElement(iel);
                        elem.ComputeDrawValueAtGaussPoint(name);

                        ////change order///////
                        for (int j = 0; j < elem.GetNumberOfGaussPointOnEachDirection(); j++)
                        {
                            for (int i = 0; i < elem.GetNumberOfGaussPointOnEachDirection(); i++)
                            {
                                GaussPoints gpsij = elem.GetGaussPoint(i, j);
                                DoubleVector Nii = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);//new DoubleVector(nen);
                                Ae += MatrixFunctions.OuterProduct(Nii, Nii);
                                Be += (double)gpsij.GetValue(DataInGausspoint.DrawValue) * Nii;
                            }
                        }
                    }
                    else if (patch is AbstractPatch3D)
                    {
                        AbstractElement3D elem = (AbstractElement3D)patch.GetElement(iel);
                        elem.ComputeDrawValueAtGaussPoint(name);

                        for (int i = 0; i < elem.GetNumberOfGaussPointOnEachDirection(); i++)
                        {
                            for (int j = 0; j < elem.GetNumberOfGaussPointOnEachDirection(); j++)
                            {
                                for (int k = 0; k < elem.GetNumberOfGaussPointOnEachDirection(); k++)
                                {
                                    GaussPoints gpsijk = elem.GetGaussPoint(i, j, k);
                                    DoubleVector Nii = (DoubleVector)gpsijk.GetValue(DataInGausspoint.Ni); //new DoubleVector(nen);
                                    Ae += MatrixFunctions.OuterProduct(Nii, Nii);
                                    Be += (double)gpsijk.GetValue(DataInGausspoint.DrawValue) * Nii;
                                }
                            }
                        }
                    }
                    for (int jj = 0; jj < nen; jj++)
                    {
                        for (int ii = 0; ii < nen; ii++)
                        {
                            if (!IsSparseData)
                                As[econn[ii], econn[jj]] += Ae[ii, jj];
                            else
                                AsSparse[econn[ii], econn[jj]] += Ae[ii, jj];
                        }
                        Bs[econn[jj]] += Be[jj];
                    }
                }
                DoubleVector Ps = (!IsSparseData) ? MatrixFunctions.Solve(As, Bs) : MatrixFunctions.Solve(AsSparse, Bs);
                if (patch is AbstractPatch2D)
                {
                    for (int ii = 0; ii < listPatch[iji].GetCountGlobalBasisFunctions(0); ii++)
                    {
                        int ik = listPatch[iji].GetINC(0, ii, 0);
                        int jk = listPatch[iji].GetINC(0, ii, 1);
                        ((NURBSSurface)listPatch[iji].GetGeometry(0)).ControlPoints[ik, jk].DrawValue = Ps[ii];
                    }
                }
                else if (patch is AbstractPatch3D)
                {
                    for (int ii = 0; ii < listPatch[iji].GetCountGlobalBasisFunctions(0); ii++)
                    {
                        int ik = listPatch[iji].GetINC(0, ii, 0);
                        int jk = listPatch[iji].GetINC(0, ii, 1);
                        int kk = listPatch[iji].GetINC(0, ii, 2);

                        ((NURBSVolume)listPatch[iji].GetGeometry(0)).ControlPoints[ik, jk, kk].DrawValue = Ps[ii];
                    }
                }
            }
        }
        public void DrawValueDistributionAtGaussPoints(ViewerForm viewer, DataInGausspoint name)
        {
            ExtrapolationValueGausspoint(name);

            List<double> x = new List<double>();
            List<double> y = new List<double>();
            List<double> z = new List<double>();
            List<double> val = new List<double>();
            PointSetColor ps;
            if (StructureDimension == Dimension.Plane)
            {
                foreach (AbstractPatch2D patch in listPatch)
                {
                    var surface = patch.GetSurface();
                    surface.isDrawControlPoint = false;
                    surface.isDrawControlNet = false;
                    surface.isColorfulFace = false;
                    surface.isDrawSurface = false;
                    surface.isDrawCurve = false;
                    surface.isDrawKnot = false;
                    surface.colorCurve = Color.Black;
                    surface.widthCurve = 1;
                    surface.resolution1 = 5;
                    surface.resolution2 = 5;
                    //surface.opacity = 0.2;
                    surface.Draw(viewer);
                    patch.ComputeDataDrawValueDistribution(ref x, ref y, ref z, ref val);
                }
            }
            else if (StructureDimension == Dimension.Solid)
            {
                foreach (AbstractPatch3D patch in listPatch)
                {
                    var vol = patch.GetVolume();
                    vol.isDrawGrid = true;
                    vol.isDrawVolume = false;
                    vol.isDrawControlPoint = false;
                    vol.isDrawControlNet = false;
                    vol.isColorfulFace = false;
                    vol.isDrawCurve = false;
                    vol.isDrawKnot = false;
                    vol.colorGrid = Color.Black;
                    vol.widthGrid = 1;
                    vol.resolution1 = 5;
                    vol.resolution2 = 5;
                    vol.resolution3 = 5;
                    //surface.opacity = 0.2;
                    vol.Draw(viewer);
                    patch.ComputeDataDrawValueDistribution(ref x, ref y, ref z, ref val);
                }
            }
            ps = new PointSetColor(x.ToArray(), y.ToArray(), z.ToArray(), val.ToArray());
            ps.SetPointSize(10);
            ps.SetOpacity(0.6);
            viewer.AddObject3D(ps);
            viewer.SetColormapBarVisible(ps, name.ToString());
        }

        public void AddComputeResult(params Result[] results)
        {
            foreach (Result re in results)
            {
                listComputeResult.Add(re);
            }
        }
        public DoubleVector GetUGlobal(bool full = false)
        {
            DoubleVector uGlobal = null;
            if (!full)
            {
                uGlobal = new DoubleVector(countDOF);
                for (int i = 0; i < listPatch.Count; i++)
                {
                    AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
                    ControlPoint[] cps = pa.GetAllControlPoints();
                    for (int j = 0; j < cps.Length; j++)
                    {
                        ControlPoint currentCps = cps[j];
                        for (int k = 0; k < pa.GetCountField(); k++)
                        {
                            int tArray = currentCps.GetTArrayGlobal()[k];
                            if (tArray != -1)
                                uGlobal[tArray] = currentCps.GetULocal(k); ;
                        }
                    }
                }
            }
            else
            {
                int d = listPatch[0].GetCountField(0);

                ///////////////////////////////////////
                ///single patch///////////////////////
                //////////////////////////////////////
                int n = listPatch[0].GetCountGlobalBasisFunctions(0);
                uGlobal = new DoubleVector(d * n);
                //////////////////////////////////////
                for (int i = 0; i < listPatch.Count; i++)
                {
                    AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
                    ControlPoint[] cps = pa.GetAllControlPoints();
                    for (int j = 0; j < cps.Length; j++)
                    {
                        ControlPoint currentCps = cps[j];
                        for (int k = 0; k < d; k++)
                        {
                            uGlobal[k + j * d] = currentCps.GetULocal(k);
                        }
                    }
                }
            }
            return uGlobal;
        }
        public void WriteResultDataOnPatch(int[] resolution, params Result[] listResult)
        {
            List<PatchDataResult> listPatchDataResult = new List<PatchDataResult>();

            double[,] listMaxMinGlobalResult = InitialMaxMin(listResult.Length);
            if (StructureDimension == Dimension.Beam)
            {
                for (int i = 0; i < CountPatch(); i++)
                {
                    PatchDataResult dataPatch = new PatchDataResult(i);
                    NURBSCurve curve = ((AbstractPatch1D)listPatch[i]).GetCurve();
                    double[][] data1D = curve.GetDataOnCurve(resolution[0]);
                    dataPatch.SetResolution(data1D.Length);
                    double[] x = new double[data1D.Length];
                    double[] y = new double[data1D.Length];
                    double[] z = new double[data1D.Length];
                    int c = 0;
                    for (int ii = 0; ii < data1D.Length; ii++)
                    {
                        x[c] = data1D[ii][0];
                        y[c] = data1D[ii][1];
                        z[c] = data1D[ii][2];
                        c++;
                    }

                    dataPatch.InputResultData(Result.X, x);
                    dataPatch.InputResultData(Result.Y, y);
                    dataPatch.InputResultData(Result.Z, z);
                    for (int j = 0; j < listResult.Length; j++)
                    {
                        double[] val = listPatch[i].CalculateResult(listResult[j], resolution[0]);
                        double valMin = val.Min();
                        double valMax = val.Max();
                        listMaxMinGlobalResult[0, j] = Math.Min(listMaxMinGlobalResult[0, j], valMin);
                        listMaxMinGlobalResult[1, j] = Math.Max(listMaxMinGlobalResult[1, j], valMax);
                        dataPatch.InputResultData(listResult[j], val);
                    }
                    listPatchDataResult.Add(dataPatch);
                    string nameFile = PathProject + "data_on_patch_" + i.ToString() + ".msgpack";
                    listPatchDataResult[i].WriteData(nameFile);
                }
            }
            else if (StructureDimension == Dimension.Plane || StructureDimension == Dimension.Plate)
            {
                for (int i = 0; i < CountPatch(); i++)
                {
                    PatchDataResult dataPatch = new PatchDataResult(i);
                    NURBSSurface surface = ((AbstractPatch2D)listPatch[i]).GetSurface();
                    double[,][] data2D = surface.GetDataOnSurface(resolution[0], resolution[1]);
                    dataPatch.SetResolution(data2D.GetLength(0), data2D.GetLength(1));
                    double[] x = new double[data2D.Length];
                    double[] y = new double[data2D.Length];
                    double[] z = new double[data2D.Length];
                    int c = 0;
                    for (int ii = 0; ii < data2D.GetLength(0); ii++)
                    {
                        for (int jj = 0; jj < data2D.GetLength(1); jj++)
                        {
                            x[c] = data2D[ii, jj][0];
                            y[c] = data2D[ii, jj][1];
                            z[c] = data2D[ii, jj][2];
                            c++;
                        }
                    }

                    dataPatch.InputResultData(Result.X, x);
                    dataPatch.InputResultData(Result.Y, y);
                    dataPatch.InputResultData(Result.Z, z);
                    for (int j = 0; j < listResult.Length; j++)
                    {
                        double[] val = listPatch[i].CalculateResult(listResult[j], resolution[0]);
                        double valMin = val.Min();
                        double valMax = val.Max();
                        listMaxMinGlobalResult[0, j] = Math.Min(listMaxMinGlobalResult[0, j], valMin);
                        listMaxMinGlobalResult[1, j] = Math.Max(listMaxMinGlobalResult[1, j], valMax);
                        dataPatch.InputResultData(listResult[j], val);
                    }
                    listPatchDataResult.Add(dataPatch);
                    string nameFile = PathProject + "data_on_patch_" + i.ToString() + ".msgpack";
                    listPatchDataResult[i].WriteData(nameFile);
                }
            }
            string nameMaxMinFile = PathProject + "max_min.msgpack";
            IOIGA.WriteArrayMessagePackFormatter(nameMaxMinFile, listMaxMinGlobalResult);

            string nameResultOutputFile = PathProject + "list_result_output.msgpack";
            IOIGA.WriteMessagePackData(nameResultOutputFile, listResult);
        }
        private double[,] InitialMaxMin(int countCol)
        {
            double[,] list = new double[2, countCol];
            for (int i = 0; i < countCol; i++)
            {
                list[0, i] = 1e8;
                list[1, i] = -1e8;
            }
            return list;
        }


        /// <summary>
        /// Draw result with specify results
        /// </summary>
        /// <param name="re">kind of result UX, UY, USUM,...</param>
        /// <param name="viewer">viewer form</param>
        /// <param name="resolution">resolution on a element</param>
        /// <param name="isDeformation">draw in deformation or not</param>
        /// <param name="scale">scale of deformation</param>
        public void DrawResult(Result re, ViewerForm viewer, int resolution, bool isDeformation, double scale, Orientation orient = Orientation.Horizontal, ColorType colorType = ColorType.Jet)
        {
            double minVal = 1e100;
            double maxVal = -1e100;
            if (StructureDimension == Dimension.Beam)
            {
                Contours2D[] contours = new Contours2D[listPatch.Count];
                //int[] res = { 60, 20 };
                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    NURBSCurve curve = null;

                    double[][] data = null;

                    if (!isDeformation)
                    {
                        curve = ((AbstractPatch1D)listPatch[lp]).GetCurve();
                    }
                    else
                    {
                        curve = ((AbstractPatch1D)listPatch[lp]).GetDeformationCurve(scale);
                    }

                    //data = surface.GetDataOnSurface(res[lp], res[lp]);
                    data = curve.GetDataOnCurve(resolution);
                    double[,] xData = new double[data.GetLength(0), 1];
                    double[,] yData = new double[data.GetLength(0), 1];
                    double[,] zData = new double[data.GetLength(0), 1];

                    for (int i = 0; i < data.GetLength(0); i++)
                    {
                        xData[i, 0] = data[i][0];
                        yData[i, 0] = data[i][1];
                        zData[i, 0] = data[i][2];
                    }
                    contours[lp] = new Contours2D(xData, yData, zData);
                    //contours[lp].SetScalarValue(listPatch[lp].CalculateResult(re, res[lp]));
                    contours[lp].SetScalarValue(listPatch[lp].CalculateResult(re, resolution));
                    contours[lp].SetOpacity(10);
                    contours[lp].ColorType = colorType;
                    var val = contours[lp].GetValue();
                    var miVal = val.Min();
                    var maVal = val.Max();
                    if (miVal < minVal)
                        minVal = miVal;
                    if (maVal > maxVal)
                        maxVal = maVal;
                }

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    contours[lp].SetMinMaxValue(minVal, maxVal);
                    contours[lp].Update(minVal, maxVal);
                    viewer.AddObject3D(contours[lp]);
                }
                viewer.SetColormapBarVisible(contours[0], re.ToString(), orient);
            }
            else if (StructureDimension == Dimension.Plane || StructureDimension == Dimension.Plate)
            {
                Contours2D[] contours = new Contours2D[listPatch.Count];
                //int[] res = { 60, 20 };
                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    NURBSSurface surface = null;

                    double[,][] data = null;

                    if (!isDeformation)
                    {
                        surface = ((AbstractPatch2D)listPatch[lp]).GetSurface();
                    }
                    else
                    {
                        surface = ((AbstractPatch2D)listPatch[lp]).GetDeformationSurface(scale);
                    }

                    //data = surface.GetDataOnSurface(res[lp], res[lp]);
                    data = surface.GetDataOnSurface(resolution, resolution);
                    double[,] xData = new double[data.GetLength(0), data.GetLength(1)];
                    double[,] yData = new double[data.GetLength(0), data.GetLength(1)];
                    double[,] zData = new double[data.GetLength(0), data.GetLength(1)];

                    for (int i = 0; i < data.GetLength(0); i++)
                        for (int j = 0; j < data.GetLength(1); j++)
                        {
                            xData[i, j] = data[i, j][0];
                            yData[i, j] = data[i, j][1];
                            zData[i, j] = data[i, j][2];
                        }
                    contours[lp] = new Contours2D(xData, yData, zData);
                    //contours[lp].SetScalarValue(listPatch[lp].CalculateResult(re, res[lp]));
                    contours[lp].SetScalarValue(listPatch[lp].CalculateResult(re, resolution));
                    contours[lp].SetOpacity(10);
                    contours[lp].ColorType = colorType;
                    var val = contours[lp].GetValue();
                    var miVal = val.Min();
                    var maVal = val.Max();
                    if (miVal < minVal)
                        minVal = miVal;
                    if (maVal > maxVal)
                        maxVal = maVal;
                }

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    contours[lp].SetMinMaxValue(minVal, maxVal);
                    contours[lp].Update(minVal, maxVal);
                    viewer.AddObject3D(contours[lp]);
                }
                viewer.SetColormapBarVisible(contours[0], re.ToString(), orient);
            }
            else if (StructureDimension == Dimension.Solid)
            {
                Contours3D[] contours = new Contours3D[listPatch.Count];
                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    AbstractPatch3D pa = (AbstractPatch3D)listPatch[lp];
                    double[,,][] data = null;
                    if (!isDeformation)
                        data = pa.GetVolume().GetDataOnVolume(resolution, resolution, resolution);
                    else
                        data = pa.GetDeformationVolume(scale).GetDataOnVolume(resolution, resolution, resolution);
                    double[,,] xData = new double[data.GetLength(0), data.GetLength(1), data.GetLength(2)];
                    double[,,] yData = new double[data.GetLength(0), data.GetLength(1), data.GetLength(2)];
                    double[,,] zData = new double[data.GetLength(0), data.GetLength(1), data.GetLength(2)];

                    for (int i = 0; i < data.GetLength(0); i++)
                        for (int j = 0; j < data.GetLength(1); j++)
                            for (int k = 0; k < data.GetLength(2); k++)
                            {
                                xData[i, j, k] = data[i, j, k][0];
                                yData[i, j, k] = data[i, j, k][1];
                                zData[i, j, k] = data[i, j, k][2];
                            }
                    contours[lp] = new Contours3D(xData, yData, zData);
                    contours[lp].SetScalarValue(pa.CalculateResult(re, resolution));
                    contours[lp].ColorType = colorType;
                    //contours[lp].SetOpacity(0.5);

                    var val = contours[lp].GetValue();
                    var miVal = val.Min();
                    var maVal = val.Max();
                    if (miVal < minVal)
                        minVal = miVal;
                    if (maVal > maxVal)
                        maxVal = maVal;
                }
                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    contours[lp].SetMinMaxValue(minVal, maxVal);
                    contours[lp].Update(minVal, maxVal);
                    viewer.AddObject3D(contours[lp]);
                }
                viewer.SetColormapBarVisible(contours[0], re.ToString(), orient);
            }
        }

        public void DrawResult(Result re, ViewerForm viewer, int resolution, bool isDeformation, double scale, double minVal, double maxVal, Orientation orient = Orientation.Horizontal, ColorType colorType = ColorType.Jet)
        {
            if (StructureDimension == Dimension.Plane || StructureDimension == Dimension.Plate)
            {
                Contours2D[] contours = new Contours2D[listPatch.Count];

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    NURBSSurface surface = null;

                    double[,][] data = null;

                    if (!isDeformation)
                    {
                        surface = ((AbstractPatch2D)listPatch[lp]).GetSurface();
                    }
                    else
                    {
                        surface = ((AbstractPatch2D)listPatch[lp]).GetDeformationSurface(scale);
                    }

                    data = surface.GetDataOnSurface(resolution, resolution);
                    double[,] xData = new double[data.GetLength(0), data.GetLength(1)];
                    double[,] yData = new double[data.GetLength(0), data.GetLength(1)];
                    double[,] zData = new double[data.GetLength(0), data.GetLength(1)];

                    for (int i = 0; i < data.GetLength(0); i++)
                        for (int j = 0; j < data.GetLength(1); j++)
                        {
                            xData[i, j] = data[i, j][0];
                            yData[i, j] = data[i, j][1];
                            zData[i, j] = data[i, j][2];
                        }
                    contours[lp] = new Contours2D(xData, yData, zData);
                    contours[lp].SetScalarValue(listPatch[lp].CalculateResult(re, resolution));
                    contours[lp].SetOpacity(10);
                    contours[lp].ColorType = colorType;
                }

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    contours[lp].SetMinMaxValue(minVal, maxVal);
                    contours[lp].Update(minVal, maxVal);
                    viewer.AddObject3D(contours[lp]);
                }
                viewer.SetColormapBarVisible(contours[0], re.ToString(), orient);
            }
            else if (StructureDimension == Dimension.Solid)
            {
                Contours3D[] contours = new Contours3D[listPatch.Count];

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    AbstractPatch3D pa = (AbstractPatch3D)listPatch[lp];
                    double[,,][] data = null;
                    if (!isDeformation)
                        data = pa.GetVolume().GetDataOnVolume(resolution, resolution, resolution);
                    else
                        data = pa.GetDeformationVolume(scale).GetDataOnVolume(resolution, resolution, resolution);
                    double[,,] xData = new double[data.GetLength(0), data.GetLength(1), data.GetLength(2)];
                    double[,,] yData = new double[data.GetLength(0), data.GetLength(1), data.GetLength(2)];
                    double[,,] zData = new double[data.GetLength(0), data.GetLength(1), data.GetLength(2)];

                    for (int i = 0; i < data.GetLength(0); i++)
                        for (int j = 0; j < data.GetLength(1); j++)
                            for (int k = 0; k < data.GetLength(2); k++)
                            {
                                xData[i, j, k] = data[i, j, k][0];
                                yData[i, j, k] = data[i, j, k][1];
                                zData[i, j, k] = data[i, j, k][2];
                            }
                    contours[lp] = new Contours3D(xData, yData, zData);
                    contours[lp].SetScalarValue(pa.CalculateResult(re, resolution));
                    contours[lp].ColorType = colorType;
                    //contours[lp].SetOpacity(0.5);
                }

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    contours[lp].SetMinMaxValue(minVal, maxVal);
                    contours[lp].Update(minVal, maxVal);
                    viewer.AddObject3D(contours[lp]);
                }
                viewer.SetColormapBarVisible(contours[0], re.ToString(), orient);
            }
        }

        public void AddConstraint(params AbstractConstraintValue[] c)
        {
            for (int i = 0; i < c.Length; i++)
                listConstraint.Add(c[i]);
        }

        public void AddInitialConstraint(params AbstractConstraintValue[] c)
        {
            for (int i = 0; i < c.Length; i++)
                listInitialConstraint.Add(c[i]);
        }
        public int CountConstraint()
        { return listConstraint.Count; }


        private void ApplyCouplingInterfaces()
        {
            if (!IsSparseData)
            {
                TDense = new DoubleMatrix(countDOF, countDOF);
                foreach (CouplingTwoPatches c in listConnection)
                {
                    c.MakeCoupling(ref TDense);
                }

                foreach (CouplingTwoPatches c in listConnection)
                {
                    c.MakeCouplingTs1m2(ref TDense);
                }

                ////////////////////////////////////////////////////////////
                //////////////Coupling control point at meeting point ////////
                if (listControlPointCoupling != null && !IsMergeControlpoints1_1Constraint)
                {
                    for (int i = 0; i < listControlPointCoupling.Count; i++)
                    {
                        if (listControlPointCoupling[i].Count > 2)
                        {
                            List<ControlPoint> li = listControlPointCoupling[i];
                            for (int j = 0; j < li.Count - 1; j++)
                            {
                                CouplingTwoControlPoints c = new CouplingTwoControlPoints(li[li.Count - 1], li[j]);//chua xet den bo coupling o dof chi dinh
                                c.MakeCoupling(ref TDense);
                            }
                        }
                    }
                }
                ////


                ////////////////////////////////////////////////////////////
                //#region run_good
                //////////run good//////
                //////Substitute to masterDOF of slaveDOF to have one slaveDOF has many last masterDOFs with multi-level masters
                //// On diagonal must be equal -1
                //for (int i = 0; i < countDOF; i++)////////////////////////////////////
                //{
                //  if (Math.Abs(TDense[i, i]) > 1e-9)
                //  {
                //    double Tii = TDense[i, i];
                //    if (Tii != -1)
                //    {
                //      for (int j = 0; j < countDOF; j++)
                //      {
                //        TDense[i, j] /= -Tii;
                //      }
                //    }
                //    List<int> indexRow = MatrixTool.FindIndexItemNonZeroOnCol(TDense, i);
                //    List<int> indexCol = MatrixTool.FindIndexItemNonZeroOnRow(TDense, i);
                //    for (int j = 0; j < indexRow.Count; j++)
                //    {
                //      if (i != indexRow[j])
                //      {
                //        double f = TDense[indexRow[j], i];

                //        for (int k = 0; k < indexCol.Count; k++)
                //        {
                //          if (i != indexCol[k])
                //          {
                //            TDense[indexRow[j], indexCol[k]] += f * TDense[i, indexCol[k]];
                //          }
                //        }
                //        TDense[indexRow[j], i] = 0;
                //      }
                //    }
                //  }
                //}
                //#endregion
            }
            else
            {
                TSparseBuilder = new SparseMatrixBuilder<double>();
                foreach (CouplingTwoPatches c in listConnection)
                {
                    c.MakeCoupling(ref TSparseBuilder);
                }

                foreach (CouplingTwoPatches c in listConnection)
                {
                    c.MakeCouplingTs1m2(ref TSparseBuilder);
                }

                ////////////////////////////////////////////////////////////
                //////////////Coupling control point at meeting point ////////
                if (listControlPointCoupling != null && !IsMergeControlpoints1_1Constraint)
                {
                    for (int i = 0; i < listControlPointCoupling.Count; i++)
                    {
                        if (listControlPointCoupling[i].Count > 2)
                        {
                            List<ControlPoint> li = listControlPointCoupling[i];
                            for (int j = 0; j < li.Count - 1; j++)
                            {
                                CouplingTwoControlPoints c = new CouplingTwoControlPoints(li[li.Count - 1], li[j]);//chua xet den bo coupling o dof chi dinh
                                c.MakeCoupling(ref TSparseBuilder);
                            }
                        }
                    }
                }
                //  ////////////////////////////////////////////////////////////
                //#region run_good
                //////////run good//////
                //////Substitute to masterDOF of slaveDOF to have one slaveDOF has many last masterDOFs with multi-level masters
                //// On diagonal must be equal -1
                //for (int i = 0; i < countDOF; i++)////////////////////////////////////
                //{
                //  if (TSparseBuilder[i, i] != 0)
                //  {
                //    double Tii = TSparseBuilder[i, i];
                //    if (Tii != -1)
                //    {
                //      for (int j = 0; j < countDOF; j++)
                //      {
                //        TSparseBuilder[i, j] /= -Tii;
                //      }
                //    }
                //    List<int> indexRow = MatrixTool.FindIndexItemNonZeroOnCol(TSparseBuilder, i);
                //    List<int> indexCol = MatrixTool.FindIndexItemNonZeroOnRow(TSparseBuilder, i);
                //    for (int j = 0; j < indexRow.Count; j++)
                //    {
                //      if (i != indexRow[j])
                //      {
                //        double f = TSparseBuilder[indexRow[j], i];

                //        for (int k = 0; k < indexCol.Count; k++)
                //        {
                //          if (i != indexCol[k])
                //          {
                //            TSparseBuilder[indexRow[j], indexCol[k]] += f * TSparseBuilder[i, indexCol[k]];
                //          }
                //        }
                //        TSparseBuilder[indexRow[j], i] = 0;
                //      }
                //    }
                //  }
                //}
                //#endregion
            }
        }

        //public void ApplyBoundaryCondition(ref DoubleMatrix kGlobal, ref DoubleVector rGlobal, double time = 1)
        //{
        //  if (isUseTransformMatrix)
        //  {
        //    kGlobal = MatrixFunctions.Product(T.Transpose(), MatrixFunctions.Product(kGlobal, T));
        //    rGlobal = MatrixFunctions.Product(T.Transpose(), rGlobal);

        //    foreach (ConnectionInterfaceTwoPatches c in listConnection)
        //    {
        //      if (c.GetTSMMatrix() != null)
        //      {
        //        int[] t = c.GetTArrayOnLineMasterSlave(false);
        //        for (int i = 0; i < t.Length; i++)
        //        { kGlobal[t[i], t[i]] = 1; }
        //      }

        //      if (c.GetTS2M1Matrix() != null)
        //      {
        //        int[] t = c.GetTArrayOnLineMasterSlave2(false);
        //        for (int i = 0; i < t.Length; i++)
        //        { kGlobal[t[i], t[i]] = 1; }
        //      }
        //    }
        //  }
        //  foreach (AbstractConstraintValue c in listConstraint)
        //  {
        //    double valueConstraint = 0;
        //    if ((c.GetValueConstraint() is double) || (c.GetValueConstraint() is ConstantFunctionRToR))
        //      valueConstraint = c.GetValueConstraint(time);
        //    ControlPoint[] cpsel = c.GetControlPointsConstrainted();
        //    double[] valueFunctionConstraint = null;
        //    if (c is ConstraintValueEdge2D)
        //    {
        //      valueFunctionConstraint = c.ComputeValueConstraintOnControlPoints(time);
        //    }
        //    else if (c is ConstraintValueFace3D)
        //    {

        //    }

        //    for (int k = 0; k < cpsel.Length; k++)
        //    {
        //      if (!(this is ModelStructurePhaseFieldStatic))
        //      {
        //        if ((c.GetValueConstraint() is FunctionRToR) && (!(c.GetValueConstraint() is ConstantFunctionRToR)))
        //        {
        //          valueConstraint = valueFunctionConstraint[k];
        //        }
        //      }
        //      else
        //      {
        //        valueConstraint = c.GetValueVelocityConstraint(time);
        //        if (!IsApplyBoundaryValue)
        //          valueConstraint = 0;
        //      }
        //      int[] tArrayCp = cpsel[k].GetTArrayGlobal();
        //      int tArray = tArrayCp[c.GetIDOfField()];
        //      if (tArray != -1)
        //      {
        //        DoApplyBoundary(tArray, valueConstraint, ref kGlobal, ref rGlobal);
        //      }
        //    }
        //  }
        //}
        private void ApplyUnchangeBoundaryCondition(double time = 1)
        {
            g = new DoubleVector(countDOF);

            foreach (AbstractConstraintValue c in listConstraint)
            {
                double valueConstraint = 0;
                if ((c.GetValueConstraint() is double) || (c.GetValueConstraint() is ConstantFunctionRToR))
                    valueConstraint = c.GetValueConstraint(time);
                ControlPoint[] cpsel = c.GetControlPointsConstrainted();
                double[] valueFunctionConstraint = null;
                if (c is ConstraintValueEdge2D)
                {
                    valueFunctionConstraint = c.ComputeValueConstraintOnControlPoints(time);
                }
                else if (c is ConstraintValueFace3D)
                {

                }

                Coordination coord = c.GetCoordination();
                //DoubleMatrix mat = MatrixTool.ConvertMathNetMatrixToNMathMatrix(coord.GetTransformationMatrix());
                //DoubleMatrix mat = MatrixTool.ConvertMathNetMatrixToNMathMatrix(coord.GetTransformationMatrix().Inverse());
                DoubleMatrix mat = new DoubleMatrix(coord.GetInverseTransformationMatrixToArray());
                double[] coefConstrain = new double[mat.Cols];//{ mat[c.GetIDOfField(), 0], mat[c.GetIDOfField(), 1], mat[c.GetIDOfField(), 2] };
                for (int i = 0; i < coefConstrain.Length; i++)
                {
                    coefConstrain[i] = mat[c.GetIDOfField(), i];
                }
                for (int k = 0; k < cpsel.Length; k++)
                {
                    if (!(this is ModelStructurePhaseFieldStatic))
                    {
                        if ((c.GetValueConstraint() is FunctionRToR) && (!(c.GetValueConstraint() is ConstantFunctionRToR)))
                        {
                            valueConstraint = valueFunctionConstraint[k];
                        }
                    }
                    else
                    {
                        valueConstraint = c.GetValueVelocityConstraint(time);
                        if (!IsApplyBoundaryValue)
                            valueConstraint = 0;
                    }
                    int[] tArrayCp = cpsel[k].GetTArrayGlobal();
                    //if ((coefConstrain[1] == 0 && coefConstrain[2] == 0) || (coefConstrain[0] == 0 && coefConstrain[1] == 0) || (coefConstrain[0] == 0 && coefConstrain[2] == 0))
                    //{

                    //int tArray = tArrayCp[c.GetIDOfField()];
                    //if (tArray != -1)
                    //{
                    //  List<int> index;
                    //  if (!IsSparseData)
                    //    index = MatrixTool.FindIndexItemNonZeroOnCol(TDense, tArray);
                    //  else
                    //    index = MatrixTool.FindIndexItemNonZeroOnCol(TSparseBuilder, tArray);
                    //  for (int i = 0; i < index.Count; i++)
                    //  {
                    //    if (tArray != index[i])
                    //    {
                    //      if (!IsSparseData)
                    //      {
                    //        g[index[i]] = TDense[index[i], tArray] * valueConstraint;
                    //        TDense[index[i], tArray] = 0;
                    //      }
                    //      else
                    //      {
                    //        g[index[i]] = TSparseBuilder[index[i], tArray] * valueConstraint;
                    //        TSparseBuilder[index[i], tArray] = 0;
                    //      }
                    //    }
                    //  }
                    //  if (!IsSparseData)
                    //  {
                    //    TDense[tArray, tArray] = -1;
                    //  }
                    //  else
                    //  {
                    //    TSparseBuilder[tArray, tArray] = -1;//(dua vao dependent) ////=0;
                    //  }
                    //  g[tArray] = valueConstraint;
                    //}
                    //}
                    //else
                    //{
                    List<int> listDOFConstraint = new List<int>();
                    for (int i = 0; i < 3; i++)
                    {
                        if (Math.Abs(coefConstrain[i]) > 1e-12)
                        {
                            listDOFConstraint.Add(i);
                        }
                    }
                    for (int i = 0; i < listDOFConstraint.Count; i++)
                    {
                        int tArray = tArrayCp[listDOFConstraint[i]];
                        if (tArray != -1)
                        {
                            double coefPivot = coefConstrain[listDOFConstraint[i]];
                            for (int j = 0; j < tArrayCp.Length; j++)
                            {
                                if (tArrayCp[j] != tArray)
                                {
                                    if (!IsSparseData)
                                    {
                                        TDense[tArray, tArrayCp[j]] = Math.Round(TDense[tArray, tArrayCp[j]] + coefConstrain[j] / (-coefPivot), numberOfDecimal);
                                    }
                                    else
                                    {
                                        TSparseBuilder[tArray, tArrayCp[j]] = Math.Round(TSparseBuilder[tArray, tArrayCp[j]] + coefConstrain[j] / (-coefPivot), numberOfDecimal);
                                        //TSparseBuilder[tArray, tArrayCp[j]] += coefConstrain[j] / (-coefPivot);
                                    }
                                }
                            }

                            if (!IsSparseData)
                            {
                                TDense[tArray, tArray] = TDense[tArray, tArray] - 1;
                                //TDense[tArray, tArray] = Math.Round(TDense[tArray, tArray] - 1, 9);
                            }
                            else
                            {
                                TSparseBuilder[tArray, tArray] = TSparseBuilder[tArray, tArray] - 1;
                                //TSparseBuilder[tArray, tArray] = Math.Round(TSparseBuilder[tArray, tArray] - 1, 9);
                            }
                            g[tArray] += -valueConstraint / (-coefPivot);
                        }
                    }
                    //}
                }
            }
        }

        protected void ApplyChangeBoundaryCondition(double time = 1)
        {
            //T = new DoubleMatrix(countDOF, countDOF);
            //g = new DoubleVector(countDOF);

            foreach (AbstractConstraintValue c in listConstraint)
            {
                double valueConstraint = 0;
                if ((c.GetValueConstraint() is double) || (c.GetValueConstraint() is ConstantFunctionRToR))
                    valueConstraint = c.GetValueConstraint(time);
                ControlPoint[] cpsel = c.GetControlPointsConstrainted();
                double[] valueFunctionConstraint = null;
                if (c is ConstraintValueEdge2D)
                {
                    valueFunctionConstraint = c.ComputeValueConstraintOnControlPoints(time);
                }
                else if (c is ConstraintValueFace3D)
                {

                }

                for (int k = 0; k < cpsel.Length; k++)
                {
                    if ((c.GetValueConstraint() is FunctionRToR) && (!(c.GetValueConstraint() is ConstantFunctionRToR)))
                    {
                        valueConstraint = valueFunctionConstraint[k];
                    }
                    int[] tArrayCp = cpsel[k].GetTArrayGlobal();
                    int tArray = tArrayCp[c.GetIDOfField()];
                    if (tArray != -1)
                    {
                        List<int> index;
                        if (!IsSparseData)
                            index = MatrixTool.FindIndexItemNonZeroOnCol(TDense, tArray, AbstractModel.NumberOfCPUs);
                        else
                            index = MatrixTool.FindIndexItemNonZeroOnCol(TSparseBuilder, tArray, AbstractModel.countDOF, AbstractModel.NumberOfCPUs);
                        for (int i = 0; i < index.Count; i++)
                        {
                            if (tArray != index[i])
                            {
                                if (!IsSparseData)
                                {
                                    g[index[i]] = TDense[index[i], tArray] * valueConstraint;
                                    TDense[index[i], tArray] = 0;
                                }
                                else
                                {
                                    g[index[i]] = TSparseBuilder[index[i], tArray] * valueConstraint;
                                    TSparseBuilder[index[i], tArray] = 0;
                                }
                            }
                        }

                        if (!IsSparseData)
                        {
                            TDense[tArray, tArray] = -1;
                        }
                        else
                        {
                            TSparseBuilder[tArray, tArray] = -1;//(dua vao dependent) ////=0;
                        }

                        g[tArray] = valueConstraint;
                        //DoApplyBoundary(tArray, valueConstraint, ref kGlobal, ref rGlobal);
                    }
                }
            }
            //UpdateTMatrixAndTMatrixTranspose();
        }
        //internal void ApplyTAndg(ref SparseMatrixBuilder<double> kGlobal, ref DoubleVector rGlobal, int[] tArrayUnconstraint)
        //{
        //  if (tArrayUnconstraint == null)
        //  {
        //    //TReduce = T;
        //    //gReduce = g;
        //    SparseMatrixBuilder<double> Ttran = MatrixTool.Transpose(T);
        //    //MatrixTool.CorrectSparseMatrixBuilderZeroRows(Ttran, countDOF, countDOF);
        //    //DoubleCsrSparseMatrix TTranspose = new DoubleCsrSparseMatrix(Ttran, countDOF);
        //    rGlobal = MatrixTool.Product(Ttran, countDOF, rGlobal - MatrixTool.Product(kGlobal, countDOF, g));
        //    kGlobal = MatrixTool.Product(Ttran, countDOF, countDOF, MatrixTool.Product(kGlobal, countDOF, countDOF, TMatrix, countDOF, countDOF), countDOF, countDOF);
        //  }
        //  else
        //  {
        //    //var TMatrix = new DoubleCsrSparseMatrix(T, countDOF);
        //    var TReduce = MatrixTool.GetSubMatrix(T, tArrayUnconstraint, tArrayUnconstraint);
        //    var gReduce = MatrixTool.GetSubVector(g, tArrayUnconstraint);
        //    SparseMatrixBuilder<double> TReduceTranspose = MatrixTool.Transpose(TReduce);
        //    kGlobal = MatrixTool.GetSubMatrix(kGlobal, tArrayUnconstraint, tArrayUnconstraint);
        //    rGlobal = MatrixTool.GetSubVector(rGlobal, tArrayUnconstraint);
        //    kGlobal = MatrixTool.Product(TReduceTranspose, countDOF, countDOF, MatrixTool.Product(kGlobal, countDOF, countDOF, TReduce, countDOF, countDOF), countDOF, countDOF);
        //    rGlobal = MatrixTool.Product(TReduceTranspose, countDOF, rGlobal - MatrixTool.Product(kGlobal, countDOF, gReduce));
        //  }
        //  for (int i = 0; i < countDOF; i++)
        //  {
        //    if (kGlobal[i, i] == 0)
        //    {
        //      kGlobal[i, i] = 1;
        //    }
        //  }
        //}

        internal void ApplyTAndg(ref DoubleCsrSparseMatrix kGlobal, ref DoubleVector rGlobal, int[] tArrayUnconstraint)
        {
            if (tArrayUnconstraint == null)
            {
                rGlobal = MatrixFunctions.Product(TSparseTranspose, rGlobal - MatrixFunctions.Product(kGlobal, g));
                kGlobal = MatrixFunctions.Product(TSparseTranspose, MatrixFunctions.Product(kGlobal, TSparse));
            }
            else
            {
                var TReduce = MatrixTool.GetSubMatrix(TSparse, tArrayUnconstraint, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                var gReduce = MatrixTool.GetSubVector(g, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                DoubleCsrSparseMatrix TReduceTranspose = MatrixTool.Transpose(TReduce);
                kGlobal = MatrixTool.GetSubMatrix(kGlobal, tArrayUnconstraint, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                rGlobal = MatrixTool.GetSubVector(rGlobal, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                kGlobal = MatrixFunctions.Product(TReduceTranspose, MatrixFunctions.Product(kGlobal, TReduce));
                rGlobal = MatrixFunctions.Product(TReduceTranspose, rGlobal - MatrixFunctions.Product(kGlobal, gReduce));
            }

            if (!AbstractModel.IsParallelProcesing)
            {
                for (int i = 0; i < kGlobal.Cols; i++)
                {
                    if (kGlobal.Data[i, i] == 0)
                    {
                        kGlobal.Data[i, i] = 1.0;
                    }
                }
            }
            else
            {
                var kTemp = kGlobal;

                var degreeOfParallelism = NumberOfCPUs;
                var tasks = new Task[degreeOfParallelism];
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                         () =>
                         {
                             var max = kTemp.Cols * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = kTemp.Cols * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 if (kTemp.Data[i, i] == 0)
                                 {
                                     kTemp.Data[i, i] = 1.0;
                                 }
                             }
                         });
                }
                Task.WaitAll(tasks);
                kGlobal = kTemp;
            }
        }

        internal void ApplyTAndg(ref DoubleCsrSparseMatrix kGlobal, ref DoubleCsrSparseMatrix mGlobal, int[] tArrayUnconstraint)
        {
            if (tArrayUnconstraint == null)
            {
                mGlobal = MatrixFunctions.Product(TSparseTranspose, MatrixFunctions.Product(mGlobal, TSparse));
                kGlobal = MatrixFunctions.Product(TSparseTranspose, MatrixFunctions.Product(kGlobal, TSparse));
            }
            else
            {
                var TReduce = MatrixTool.GetSubMatrix(TSparse, tArrayUnconstraint, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                var gReduce = MatrixTool.GetSubVector(g, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                DoubleCsrSparseMatrix TReduceTranspose = MatrixTool.Transpose(TReduce);
                kGlobal = MatrixTool.GetSubMatrix(kGlobal, tArrayUnconstraint, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                mGlobal = MatrixTool.GetSubMatrix(mGlobal, tArrayUnconstraint, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                kGlobal = MatrixFunctions.Product(TReduceTranspose, MatrixFunctions.Product(kGlobal, TReduce));
                mGlobal = MatrixFunctions.Product(TReduceTranspose, MatrixFunctions.Product(mGlobal, TReduce));
            }

            if (!AbstractModel.IsParallelProcesing)
            {
                for (int i = 0; i < kGlobal.Cols; i++)
                {
                    if (kGlobal.Data[i, i] == 0)
                    {
                        kGlobal.Data[i, i] = 1.0;
                        mGlobal.Data[i, i] = 1.0;
                    }
                }
            }
            else
            {
                var kTemp = kGlobal;
                var mTemp = mGlobal;

                var degreeOfParallelism = NumberOfCPUs;
                var tasks = new Task[degreeOfParallelism];
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                         () =>
                         {
                             var max = kTemp.Cols * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = kTemp.Cols * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 if (kTemp.Data[i, i] == 0)
                                 {
                                     kTemp.Data[i, i] = 1.0;
                                     mTemp.Data[i, i] = 1.0;
                                 }
                             }
                         });
                }
                Task.WaitAll(tasks);
                kGlobal = kTemp;
                mGlobal = mTemp;
            }
        }
        internal void ApplyTAndg(ref DoubleSymmetricMatrix kGlobal, ref DoubleVector rGlobal, int[] tArrayUnconstraint)
        {
            //if (tArrayUnconstraint == null)
            //{
            rGlobal = MatrixFunctions.Product(TDenseTranspose, rGlobal - MatrixFunctions.Product(kGlobal, g));
            kGlobal = new DoubleSymmetricMatrix(MatrixFunctions.Product(TDenseTranspose, MatrixFunctions.Product(kGlobal.ToGeneralMatrix(), TDense)));
            //}
            //else
            //{
            //  TReduce = MatrixTool.GetSubMatrix(TSparse, tArrayUnconstraint, tArrayUnconstraint);
            //  gReduce = MatrixTool.GetSubVector(g, tArrayUnconstraint);
            //  DoubleMatrix TReduceTranspose = TReduce.Transpose();
            //  kGlobal = MatrixTool.GetSubMatrix(kGlobal, tArrayUnconstraint, tArrayUnconstraint);
            //  rGlobal = MatrixTool.GetSubVector(rGlobal, tArrayUnconstraint);
            //  kGlobal = NMathFunctions.Product(TReduceTranspose, NMathFunctions.Product(kGlobal, TReduce));
            //  rGlobal = NMathFunctions.Product(TReduceTranspose, rGlobal - NMathFunctions.Product(kGlobal, gReduce));
            //}
            if (!AbstractModel.IsParallelProcesing)
            {
                for (int i = 0; i < kGlobal.Cols; i++)
                {
                    if (kGlobal[i, i] == 0)
                    {
                        kGlobal[i, i] = 1.0;
                    }
                }
            }
            else
            {
                var kTemp = kGlobal;

                var degreeOfParallelism = NumberOfCPUs;
                var tasks = new Task[degreeOfParallelism];
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                         () =>
                         {
                             var max = kTemp.Cols * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = kTemp.Cols * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 if (kTemp[i, i] == 0)
                                 {
                                     kTemp[i, i] = 1.0;
                                 }
                             }
                         });
                }
                Task.WaitAll(tasks);
                kGlobal = kTemp;
            }
        }

        internal void ApplyTAndg(ref DoubleMatrix kGlobal, ref DoubleVector rGlobal, int[] tArrayUnconstraint)
        {
            //if (tArrayUnconstraint == null)
            //{
            rGlobal = MatrixFunctions.Product(TDenseTranspose, rGlobal - MatrixFunctions.Product(kGlobal, g));
            kGlobal = MatrixFunctions.Product(TDenseTranspose, MatrixFunctions.Product(kGlobal, TDense));
            //}
            //else
            //{
            //  TReduce = MatrixTool.GetSubMatrix(TSparse, tArrayUnconstraint, tArrayUnconstraint);
            //  gReduce = MatrixTool.GetSubVector(g, tArrayUnconstraint);
            //  DoubleMatrix TReduceTranspose = TReduce.Transpose();
            //  kGlobal = MatrixTool.GetSubMatrix(kGlobal, tArrayUnconstraint, tArrayUnconstraint);
            //  rGlobal = MatrixTool.GetSubVector(rGlobal, tArrayUnconstraint);
            //  kGlobal = NMathFunctions.Product(TReduceTranspose, NMathFunctions.Product(kGlobal, TReduce));
            //  rGlobal = NMathFunctions.Product(TReduceTranspose, rGlobal - NMathFunctions.Product(kGlobal, gReduce));
            //}
            if (!AbstractModel.IsParallelProcesing)
            {
                for (int i = 0; i < kGlobal.Cols; i++)
                {
                    if (kGlobal[i, i] == 0)
                    {
                        kGlobal[i, i] = 1.0;
                    }
                }
            }
            else
            {
                var kTemp = kGlobal;

                var degreeOfParallelism = NumberOfCPUs;
                var tasks = new Task[degreeOfParallelism];
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                         () =>
                         {
                             var max = kTemp.Cols * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = kTemp.Cols * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 if (kTemp[i, i] == 0)
                                 {
                                     kTemp[i, i] = 1.0;
                                 }
                             }
                         });
                }
                Task.WaitAll(tasks);
                kGlobal = kTemp;
            }
        }

        //internal void ApplyTAndg(ref DoubleSymmetricMatrix kGlobal, ref DoubleSymmetricMatrix mGlobal, int[] tArrayUnconstraint)
        //{
        //  if (tArrayUnconstraint == null)
        //  {
        //    mGlobal = new DoubleSymmetricMatrix(MatrixFunctions.Product(TDenseTranspose, MatrixFunctions.Product(mGlobal.ToGeneralMatrix(), TDense)));
        //    kGlobal = new DoubleSymmetricMatrix(MatrixFunctions.Product(TDenseTranspose, MatrixFunctions.Product(kGlobal.ToGeneralMatrix(), TDense)));
        //    if (!AbstractModel.IsParallelProcesing)
        //    {
        //      for (int i = 0; i < kGlobal.Cols; i++)
        //      {
        //        if (kGlobal[i, i] == 0)
        //        {
        //          kGlobal[i, i] = 1.0;
        //          mGlobal[i, i] = 1.0;
        //        }
        //      }
        //    }
        //    else
        //    {
        //      var kTemp = kGlobal;
        //      var mTemp = mGlobal;
        //      var degreeOfParallelism = NumberOfCPUs;
        //      var tasks = new Task[degreeOfParallelism];
        //      for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
        //      {
        //        // capturing taskNumber in lambda wouldn't work correctly
        //        int taskNumberCopy = taskNumber;

        //        tasks[taskNumber] = Task.Factory.StartNew(
        //             () =>
        //             {
        //               var max = kTemp.Cols * (taskNumberCopy + 1) / degreeOfParallelism;
        //               for (int i = kTemp.Cols * taskNumberCopy / degreeOfParallelism; i < max; i++)
        //               {
        //                 if (kTemp[i, i] == 0)
        //                 {
        //                   kTemp[i, i] = 1.0;
        //                   mTemp[i, i] = 1.0;
        //                 }
        //               }
        //             });
        //      }
        //      Task.WaitAll(tasks);
        //      kGlobal = kTemp;
        //      mGlobal = mTemp;
        //    }
        //  }
        //  else
        //  {
        //    DoubleMatrix TReduce = MatrixTool.GetSubMatrix(TDense, tArrayUnconstraint, tArrayUnconstraint);
        //    DoubleMatrix TReduceTranspose = TReduce.Transpose();
        //    kGlobal = MatrixTool.GetSubMatrix(kGlobal, tArrayUnconstraint);
        //    mGlobal = MatrixTool.GetSubMatrix(mGlobal, tArrayUnconstraint);
        //    kGlobal = new DoubleSymmetricMatrix(NMathFunctions.Product(TReduceTranspose, NMathFunctions.Product(kGlobal.ToGeneralMatrix(), TReduce)));
        //    mGlobal = new DoubleSymmetricMatrix(NMathFunctions.Product(TReduceTranspose, NMathFunctions.Product(mGlobal.ToGeneralMatrix(), TReduce)));
        //  }
        //}

        internal void ApplyTAndg(ref DoubleMatrix kGlobal, ref DoubleMatrix mGlobal, int[] tArrayUnconstraint)
        {
            if (tArrayUnconstraint == null)
            {
                mGlobal = MatrixFunctions.Product(TDenseTranspose, MatrixFunctions.Product(mGlobal, TDense));
                kGlobal = MatrixFunctions.Product(TDenseTranspose, MatrixFunctions.Product(kGlobal, TDense));
                if (!AbstractModel.IsParallelProcesing)
                {
                    for (int i = 0; i < kGlobal.Cols; i++)
                    {
                        if (kGlobal[i, i] == 0)
                        {
                            kGlobal[i, i] = 1.0;
                            mGlobal[i, i] = 1.0;
                        }
                    }
                }
                else
                {
                    var kTemp = kGlobal;
                    var mTemp = mGlobal;
                    var degreeOfParallelism = NumberOfCPUs;
                    var tasks = new Task[degreeOfParallelism];
                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = kTemp.Cols * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = kTemp.Cols * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     if (kTemp[i, i] == 0)
                                     {
                                         kTemp[i, i] = 1.0;
                                         mTemp[i, i] = 1.0;
                                     }
                                 }
                             });
                    }
                    Task.WaitAll(tasks);
                    kGlobal = kTemp;
                    mGlobal = mTemp;
                }
            }
            else
            {
                DoubleMatrix TReduce = MatrixTool.GetSubMatrix(TDense, tArrayUnconstraint, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                DoubleMatrix TReduceTranspose = TReduce.Transpose();
                kGlobal = MatrixTool.GetSubMatrix(kGlobal, tArrayUnconstraint, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                mGlobal = MatrixTool.GetSubMatrix(mGlobal, tArrayUnconstraint, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                kGlobal = NMathFunctions.Product(TReduceTranspose, NMathFunctions.Product(kGlobal, TReduce));
                mGlobal = NMathFunctions.Product(TReduceTranspose, NMathFunctions.Product(mGlobal, TReduce));
            }
        }
        internal void ApplyTAndg(ref DoubleSymmetricMatrix kGlobal, int[] tArrayUnconstraint)
        {
            //if (tArrayUnconstraint == null)
            //{
            kGlobal = new DoubleSymmetricMatrix(MatrixFunctions.Product(TDenseTranspose, MatrixFunctions.Product(kGlobal.ToGeneralMatrix(), TDense)));
            //}
            //else
            //{
            //  TReduce = MatrixTool.GetSubMatrix(TSparse, tArrayUnconstraint, tArrayUnconstraint);
            //  gReduce = MatrixTool.GetSubVector(g, tArrayUnconstraint);
            //  DoubleMatrix TReduceTranspose = TReduce.Transpose();
            //  kGlobal = MatrixTool.GetSubMatrix(kGlobal, tArrayUnconstraint, tArrayUnconstraint);
            //  rGlobal = MatrixTool.GetSubVector(rGlobal, tArrayUnconstraint);
            //  kGlobal = NMathFunctions.Product(TReduceTranspose, NMathFunctions.Product(kGlobal, TReduce));
            //  rGlobal = NMathFunctions.Product(TReduceTranspose, rGlobal - NMathFunctions.Product(kGlobal, gReduce));
            //}
            if (!AbstractModel.IsParallelProcesing)
            {
                for (int i = 0; i < kGlobal.Cols; i++)
                {
                    if (kGlobal[i, i] == 0)
                    {
                        kGlobal[i, i] = 1.0;
                    }
                }
            }
            else
            {
                var kTemp = kGlobal;
                var degreeOfParallelism = NumberOfCPUs;
                var tasks = new Task[degreeOfParallelism];
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                         () =>
                         {
                             var max = kTemp.Cols * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = kTemp.Cols * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 if (kTemp[i, i] == 0)
                                 {
                                     kTemp[i, i] = 1.0;
                                 }
                             }
                         });
                }
                Task.WaitAll(tasks);
                kGlobal = kTemp;
            }
        }

        internal void ApplyTAndg(ref DoubleCsrSparseMatrix kGlobal, int[] tArrayUnconstraint)
        {
            if (tArrayUnconstraint == null)
            {
                kGlobal = MatrixFunctions.Product(TSparseTranspose, MatrixFunctions.Product(kGlobal, TSparse));
            }
            else
            {
                var TReduce = MatrixTool.GetSubMatrix(TSparse, tArrayUnconstraint, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                var gReduce = MatrixTool.GetSubVector(g, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                DoubleCsrSparseMatrix TReduceTranspose = MatrixTool.Transpose(TReduce);
                kGlobal = MatrixTool.GetSubMatrix(kGlobal, tArrayUnconstraint, tArrayUnconstraint, AbstractModel.NumberOfCPUs);
                kGlobal = MatrixFunctions.Product(TReduceTranspose, MatrixFunctions.Product(kGlobal, TReduce));
            }

            if (!AbstractModel.IsParallelProcesing)
            {
                for (int i = 0; i < kGlobal.Cols; i++)
                {
                    if (kGlobal.Data[i, i] == 0)
                    {
                        kGlobal.Data[i, i] = 1.0;
                    }
                }
            }
            else
            {
                var kTemp = kGlobal;

                var degreeOfParallelism = NumberOfCPUs;
                var tasks = new Task[degreeOfParallelism];
                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                {
                    // capturing taskNumber in lambda wouldn't work correctly
                    int taskNumberCopy = taskNumber;

                    tasks[taskNumber] = Task.Factory.StartNew(
                         () =>
                         {
                             var max = kTemp.Cols * (taskNumberCopy + 1) / degreeOfParallelism;
                             for (int i = kTemp.Cols * taskNumberCopy / degreeOfParallelism; i < max; i++)
                             {
                                 if (kTemp.Data[i, i] == 0)
                                 {
                                     kTemp.Data[i, i] = 1.0;
                                 }
                             }
                         });
                }
                Task.WaitAll(tasks);
                kGlobal = kTemp;
            }
        }
        protected DoubleMatrix ApplyInitialBoundaryCondition()
        {
            DoubleMatrix uuu = new DoubleMatrix(countDOF, 3);
            if (listInitialConstraint != null)
            {
                if (listInitialConstraint.Count != 0)
                {
                    foreach (AbstractConstraintValue c in listInitialConstraint)
                    {
                        double valueConstraint = c.GetValueConstraint(0);//////////// chua sua dieu kien bien constraint la ham so
                        if (!double.IsNaN(valueConstraint))
                        {
                            ControlPoint[] cpsel = c.GetControlPointsConstrainted();

                            foreach (ControlPoint cp in cpsel)
                            {
                                int[] tArrayCp = cp.GetTArrayGlobal();
                                int tArray = tArrayCp[c.GetIDOfField()];
                                if (tArray != -1)
                                {
                                    uuu[tArray, 0] = valueConstraint;
                                    uuu[tArray, 1] = c.GetValueVelocityConstraint(0);
                                    uuu[tArray, 2] = c.GetValueAccelerationConstraint(0);
                                }
                            }
                        }
                    }
                }
            }
            return uuu;
        }

        private void DoApplyBoundary(int tArray, double valueConstraint, ref DoubleMatrix kGlobal, ref DoubleVector rGlobal)
        {
            int countDOFModel = rGlobal.Length;
            //if (!IsParallelProcesing)
            //{
            for (int iii = 0; iii < countDOFModel; iii++)
            {
                if (iii != tArray)
                {
                    if (kGlobal[iii, tArray] != 0)
                    {
                        rGlobal[iii] -= valueConstraint * kGlobal[iii, tArray];
                        kGlobal[iii, tArray] = 0;
                        kGlobal[tArray, iii] = 0;
                    }
                }
            }
            //}
            //else
            //{
            //    DoubleVector copyVec = rGlobal.DeepenThisCopy();
            //    DoubleMatrix copyMat = kGlobal.DeepenThisCopy();
            //    object monitor = new object();
            //    var degreeOfParallelism = NumberOfCPUs;

            //    var tasks = new Task[degreeOfParallelism];
            //    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
            //    {
            //        // capturing taskNumber in lambda wouldn't work correctly
            //        int taskNumberCopy = taskNumber;

            //        tasks[taskNumber] = Task.Factory.StartNew(
            //            () =>
            //            {
            //                var max = countDOFModel * (taskNumberCopy + 1) / degreeOfParallelism;
            //                for (int i = countDOFModel * taskNumberCopy / degreeOfParallelism; i < max; i++)
            //                {
            //                    if (i != tArray)
            //                    {
            //                        lock (monitor)
            //                        {
            //                            if (copyMat[i, tArray] != 0)
            //                            {
            //                                copyVec[i] -= valueConstraint * copyMat[i, tArray];
            //                                copyMat[i, tArray] = 0;
            //                                copyMat[tArray, i] = 0;
            //                            }
            //                        }
            //                    }
            //                }
            //            });
            //    }
            //    Task.WaitAll(tasks);
            //    rGlobal = copyVec;
            //    kGlobal = copyMat;
            //}
            kGlobal[tArray, tArray] = 1.0;
            rGlobal[tArray] = valueConstraint;
        }

        /// <summary>
        /// For processing to apply boundary condition on K anf F
        /// [dofOnGlobal, value, typeDOFConstrainted]
        /// </summary>
        /// <param name="time"></param>
        /// <returns></returns>
        public double[,] StoreBoundaryCondition(double time = 1)
        {
            List<DoubleVector> listBC = new List<DoubleVector>();
            foreach (AbstractConstraintValue c in listConstraint)
            {
                double valueConstraint = 0;
                if ((c.GetValueConstraint() is double) || (c.GetValueConstraint() is ConstantFunctionRToR))
                    valueConstraint = c.GetValueConstraint(time);
                ControlPoint[] cpsel = c.GetControlPointsConstrainted();
                double[] valueFunctionConstraint = null;
                if (c is ConstraintValueEdge1D)
                {
                    valueFunctionConstraint = c.ComputeValueConstraintOnControlPoints(time);

                }
                else if (c is ConstraintValueEdge2D)
                {
                    valueFunctionConstraint = c.ComputeValueConstraintOnControlPoints(time);

                }
                else if (c is ConstraintValueFace3D)
                {

                }

                for (int k = 0; k < cpsel.Length; k++)
                {
                    if ((c.GetValueConstraint() is FunctionRToR) && (!(c.GetValueConstraint() is ConstantFunctionRToR)))
                    {
                        valueConstraint = valueFunctionConstraint[k];
                    }
                    int[] tArrayCp = cpsel[k].GetTArrayGlobal();
                    int tArray = tArrayCp[c.GetIDOfField()];
                    DoubleVector vBC = new DoubleVector(new double[] { tArray, valueConstraint, c.GetIDOfField() });
                    listBC.Add(vBC);
                }
            }
            if (listBC.Count > 0)
            {
                double[,] matrixBC = new double[listBC.Count, 3];
                if (!AbstractModel.IsParallelProcesing)
                {
                    for (int i = 0; i < listBC.Count; i++)
                    {
                        matrixBC[i, 0] = listBC[i][0];//DOF on global
                        matrixBC[i, 1] = listBC[i][1];//value constraint
                        matrixBC[i, 2] = listBC[i][2];//type of DOF constraint [ux=0, uy=1, uz=2]
                    }
                }
                else
                {
                    var degreeOfParallelism = NumberOfCPUs;
                    var tasks = new Task[degreeOfParallelism];
                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = listBC.Count * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = listBC.Count * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     matrixBC[i, 0] = listBC[i][0];//DOF on global
                                     matrixBC[i, 1] = listBC[i][1];//value constraint
                                     matrixBC[i, 2] = listBC[i][2];//type of DOF constraint [ux=0, uy=1, uz=2]
                                 }
                             });
                    }
                    Task.WaitAll(tasks);
                }
                return matrixBC;
            }
            else { return null; }
        }
        public int[] GetTArrayFromConstraintDof(int IDOfField)
        {
            List<int> tArray = new List<int>();
            for (int i = 0; i < bcMatrix.GetLength(0); i++)
            {
                if (bcMatrix[i, 2] == IDOfField)
                    tArray.Add((int)bcMatrix[i, 0]);
            }
            return tArray.ToArray();
        }
        public void ComputeMaxMinResult(Result re, int resolution, out double minVal, out double maxVal)
        {
            minVal = 1000000000;
            maxVal = -1000000000;
            for (int lp = 0; lp < listPatch.Count; lp++)
            {
                double[] val = listPatch[lp].CalculateResult(re, resolution);
                double min = val.Min();
                double max = val.Max();
                if (minVal > min)
                    minVal = min;
                if (maxVal < max)
                    maxVal = max;
            }
        }

        public void ReadResultByLoadStep(int stepLoad)
        {
            //string[] nameFiles = Directory.GetFiles(PathProject, "data_*.msgpack");
            string fileNameDataTime = PathProject + "data_" + stepLoad + ".msgpack"; //nameFiles[stepLoad];//PathProject + "data_" + stepLoad + ".msgpack";
            DoubleVector u = IOIGA.ReadMessagePackDoubleVector(fileNameDataTime); //DisplacementTime[stepLoad];
            SetUGlobal(u);
        }

        /// <summary>
        /// Get time by step load
        /// </summary>
        /// <param name="step">step load</param>
        /// <returns></returns>
        public double GetTimeByStep(int step)
        {
            //Time = IOIGA.ReadMessagePackData<List<double>>(FileNameTime)[0];
            return Time[step];
        }

        //public int GetNumberOfStep() { return numberOfStepLoad; }
        /// <summary>
        /// Get number of substep in special load step, load step has index 0.
        /// </summary>
        /// <param name="indexStepLoad"></param>
        /// <returns></returns>
        public int GetNumberOfSubstepLoad(int indexStepLoad) { return numberOfSubstepLoad[indexStepLoad]; }
        /// <summary>
        /// Item of array stores number of substep in this load step. Length of the array, which equals 3, means 3 load steps
        /// </summary>
        /// <param name="numberOfSubstepLoad"></param>
        public void SetLoadControlSolver(params int[] numberOfSubstepLoad)
        {
            //this.numberOfStepLoad = numberOfStepLoad;
            this.numberOfSubstepLoad = numberOfSubstepLoad;
        }
        public void SetTransientSolver(double[] deltaT, int[] numberOfSubstepLoad, int[] numberOfStepSave)
        {
            this.deltaT = deltaT;
            this.numberOfSubstepLoad = numberOfSubstepLoad;
            this.numberOfStepSave = numberOfStepSave;
        }
        public void SetIterationSolver(TypeNonlinearSolver type, CriterionConvergence criteria, int maximumOfIteration, double TOL, Norm norm, bool isSkipUnconvergenceStep = false)
        {
            this.type = type;
            convergenceOption = new Convergence(criteria, maximumOfIteration, TOL, norm);
            convergenceOption.IsSkipUnconvergenceStep = isSkipUnconvergenceStep;
        }



        private DoubleCsrSparseMatrix CreateFormOfStiffnessMatrix()
        {
            string filename = PathProject + "sparseMatrixForm.msgpack";
            //Stopwatch sw = new Stopwatch();
            DoubleCsrSparseMatrix ma = null;
            if (File.Exists(filename))
            {
                if (IsCreateSparseMatrixForm)
                {
                    ma = IOIGA.ReadMessagePackDoubleCsrSparseMatrix(filename);
                    if (ma == null)
                        File.Delete(filename);
                }
            }
            if (!File.Exists(filename))
            {
                //sw.Start();

                #region sparsematrixbuilder
                /////////////////////////////////////////////////////////
                SparseMatrixBuilder<double> kGlobalBuilder = new SparseMatrixBuilder<double>();

                if (IsMergeControlpoints1_1Constraint)
                {

                }


                for (int i = 0; i < listPatch.Count; i++)
                {
                    AbstractPatch patch = listPatch[i];
                    int nel = patch.CountElements();
                    if (!IsParallelProcesing)
                    {
                        for (int j = 0; j < nel; j++)
                        {
                            AbstractElement elem = patch.GetElement(j);
                            int[] tArray = elem.GetTArrayGlobal();
                            for (int ii = 0; ii < tArray.Length; ii++)
                            {
                                int enumerateM = tArray[ii];
                                if (enumerateM != -1)
                                    for (int jj = 0; jj < tArray.Length; jj++)
                                    {
                                        int enumerateN = tArray[jj];
                                        if (enumerateN != -1)
                                        {
                                            IntPair key = new IntPair(enumerateM, enumerateN);
                                            if (!kGlobalBuilder.ContainsKey(key))
                                                kGlobalBuilder.Add(key, 0);
                                        }
                                    }
                            }
                        }
                    }
                    else
                    {
                        var degreeOfParallelism = NumberOfCPUs;
                        var tasks = new Task[degreeOfParallelism];
                        object monitor = new object();
                        for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                        {
                            // capturing taskNumber in lambda wouldn't work correctly
                            int taskNumberCopy = taskNumber;

                            tasks[taskNumber] = Task.Factory.StartNew(
                                 () =>
                                 {
                                     var max = nel * (taskNumberCopy + 1) / degreeOfParallelism;
                                     for (int j = nel * taskNumberCopy / degreeOfParallelism; j < max; j++)
                                     {
                                         AbstractElement elem = patch.GetElement(j);
                                         int[] tArray = elem.GetTArrayGlobal();
                                         for (int ii = 0; ii < tArray.Length; ii++)
                                         {
                                             int enumerateM = tArray[ii];
                                             if (enumerateM != -1)
                                                 for (int jj = 0; jj < tArray.Length; jj++)
                                                 {
                                                     int enumerateN = tArray[jj];
                                                     if (enumerateN != -1)
                                                     {
                                                         IntPair key = new IntPair(enumerateM, enumerateN);
                                                         lock (monitor)
                                                         {
                                                             if (!kGlobalBuilder.ContainsKey(key))
                                                                 kGlobalBuilder.Add(key, 0);
                                                         }
                                                     }
                                                 }
                                         }
                                     }
                                 });
                        }
                        Task.WaitAll(tasks);
                    }
                }
                ////////////////////////////////////////
                #endregion
                ma = new DoubleCsrSparseMatrix(kGlobalBuilder, countDOF);
                if (IsCreateSparseMatrixForm)
                {
                    IOIGA.WriteMessagePackDoubleCsrSparseMatrix(filename, ma);
                }
                //sw.Stop();
                //double a = sw.Elapsed.TotalMilliseconds;//3846.2606
            }
            return ma;
        }

        protected DoubleVector SolveStatic(DoubleMatrix kGlobalDense, DoubleCsrSparseMatrix kGlobalSparse, ref DoubleVector rGlobal)
        {
            DoubleVector uGlobal;
            if (IsSparseData)
            {
                ApplyTAndg(ref kGlobalSparse, ref rGlobal, null);

                uGlobal = MatrixFunctions.Solve(kGlobalSparse, rGlobal);
                uGlobal = MatrixFunctions.Product(TSparse, uGlobal) + g;
            }
            else
            {
                //ApplyTAndg(ref kGlobalDense, ref rGlobal, null);
                //uGlobal = MatrixFunctions.Solve(kGlobalDense, rGlobal);
                //uGlobal = MatrixFunctions.Product(TDense, uGlobal) + g;

                int[] tArrayIndependent = this.tArrayIndependent.ToArray();
                int[] tArrayDependent = this.tArrayDependent.ToArray();
                int[] tArrayUnCoupleAndUnconstraint = this.tArrayUncouple.ToArray();
                int[] tArrayIndependentInNew = new int[tArrayIndependent.Length];
                int c = 0;
                for (int i = tArrayUnCoupleAndUnconstraint.Length; i < tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length; i++)
                {
                    tArrayIndependentInNew[c++] = i;
                }

                int[] tArrayUncoupleUnconstraintInNew = new int[tArrayUnCoupleAndUnconstraint.Length];
                for (int i = 0; i < tArrayUnCoupleAndUnconstraint.Length; i++)
                {
                    tArrayUncoupleUnconstraintInNew[i] = i;
                }
                DoubleVector gS = null;
                if (tArrayDependent.Length != 0)
                {
                    gS = MatrixTool.GetSubVector(g, tArrayDependent, AbstractModel.NumberOfCPUs);
                }

                DoubleVector Ru = null;
                DoubleVector Rm = null;
                DoubleVector rhs = new DoubleVector(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);

                DoubleMatrix Kuu = null;
                DoubleMatrix Kus = null;
                DoubleMatrix Kum = null;
                DoubleMatrix Kmm = null;
                if (tArrayUnCoupleAndUnconstraint.Length != 0)
                {
                    Kuu = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnCoupleAndUnconstraint, tArrayUnCoupleAndUnconstraint, AbstractModel.NumberOfCPUs);
                    Ru = MatrixTool.GetSubVector(rGlobal, tArrayUnCoupleAndUnconstraint, AbstractModel.NumberOfCPUs);
                    if (tArrayDependent.Length != 0)
                    {
                        Kus = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnCoupleAndUnconstraint, tArrayDependent, AbstractModel.NumberOfCPUs);
                        Ru += -(NMathFunctions.Product(Kus, gS));
                    }
                }

                DoubleMatrix TsmDense = null;
                if (tArrayDependent.Length != 0 && tArrayIndependent.Length != 0)
                {
                    TsmDense = MatrixTool.GetSubMatrix(TDense, tArrayDependent, tArrayIndependent, AbstractModel.NumberOfCPUs);
                    DoubleMatrix TsmTranspose = TsmDense.Transpose();
                    DoubleMatrix Kms = MatrixTool.GetSubMatrix(kGlobalDense, tArrayIndependent, tArrayDependent, AbstractModel.NumberOfCPUs);
                    DoubleMatrix Kss = MatrixTool.GetSubMatrix(kGlobalDense, tArrayDependent, tArrayDependent, AbstractModel.NumberOfCPUs);
                    DoubleVector Rs = MatrixTool.GetSubVector(rGlobal, tArrayDependent, AbstractModel.NumberOfCPUs);
                    Kmm = MatrixTool.GetSubMatrix(kGlobalDense, tArrayIndependent, tArrayIndependent, AbstractModel.NumberOfCPUs) + NMathFunctions.Product(TsmTranspose, Kms.Transpose() + NMathFunctions.Product(Kss, TsmDense)) + NMathFunctions.Product(Kms, TsmDense);
                    Rm = MatrixTool.GetSubVector(rGlobal, tArrayIndependent, AbstractModel.NumberOfCPUs) + NMathFunctions.Product(TsmTranspose, Rs);
                    if (tArrayUnCoupleAndUnconstraint.Length != 0)
                        Kum = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnCoupleAndUnconstraint, tArrayIndependent, AbstractModel.NumberOfCPUs) + NMathFunctions.Product(Kus, TsmDense);
                    Rm += -(NMathFunctions.Product(Kms, gS) + NMathFunctions.Product(TsmTranspose, NMathFunctions.Product(Kss, gS)));
                }
                DoubleMatrix lhsDense = new DoubleMatrix(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length, tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
                //rhs = new DoubleVector(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
                if (!AbstractModel.IsParallelProcesing)
                {
                    for (int i = 0; i < tArrayUnCoupleAndUnconstraint.Length; i++)
                    {
                        for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
                        {
                            lhsDense[i, j] = Kuu[i, j];
                        }
                        if (Kum != null)
                        {
                            for (int j = 0; j < tArrayIndependent.Length; j++)
                            {
                                lhsDense[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
                            }
                        }
                        rhs[i] = Ru[i];
                    }
                    for (int i = 0; i < tArrayIndependent.Length; i++)
                    {
                        for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
                        {
                            lhsDense[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
                        }

                        for (int j = 0; j < tArrayIndependent.Length; j++)
                        {
                            lhsDense[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
                        }
                        rhs[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
                    }
                }
                else
                {
                    var degreeOfParallelism = NumberOfCPUs;
                    var tasks = new Task[degreeOfParallelism];
                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = tArrayUnCoupleAndUnconstraint.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = tArrayUnCoupleAndUnconstraint.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
                                     {
                                         lhsDense[i, j] = Kuu[i, j];
                                     }
                                     if (Kum != null)
                                     {
                                         for (int j = 0; j < tArrayIndependent.Length; j++)
                                         {
                                             lhsDense[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
                                         }
                                     }
                                     rhs[i] = Ru[i];
                                 }
                             });
                    }
                    Task.WaitAll(tasks);

                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = tArrayIndependent.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = tArrayIndependent.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
                                     {
                                         lhsDense[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
                                     }

                                     for (int j = 0; j < tArrayIndependent.Length; j++)
                                     {
                                         lhsDense[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
                                     }
                                     rhs[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
                                 }
                             });
                    }
                    Task.WaitAll(tasks);
                }

                DoubleVector uLocalS = gS;
                DoubleVector uLocalUM = null;
                DoubleVector uLocalU = null;
                DoubleVector uLocalM = null;
                if (rGlobal.Length != 0)
                {
                    bool a = MatrixTool.IsSymmetricMatrix(lhsDense);
                    uLocalUM = MatrixFunctions.Solve(lhsDense, rhs);
                    uLocalU = MatrixTool.GetSubVector(uLocalUM, tArrayUncoupleUnconstraintInNew, AbstractModel.NumberOfCPUs);
                    uLocalM = MatrixTool.GetSubVector(uLocalUM, tArrayIndependentInNew, AbstractModel.NumberOfCPUs);
                    if (uLocalM != null)
                    {
                        uLocalS += MatrixFunctions.Product(TsmDense, uLocalM);
                    }
                }

                uGlobal = new DoubleVector(countDOF);
                if (!AbstractModel.IsParallelProcesing)
                {
                    for (int i = 0; i < tArrayUncoupleUnconstraintInNew.Length; i++)
                    {
                        uGlobal[tArrayUnCoupleAndUnconstraint[i]] = uLocalU[i];
                    }
                    for (int i = 0; i < tArrayIndependentInNew.Length; i++)
                    {
                        uGlobal[tArrayIndependent[i]] = uLocalM[i];
                    }
                    for (int i = 0; i < tArrayDependent.Length; i++)
                    {
                        uGlobal[tArrayDependent[i]] = uLocalS[i];
                    }
                }
                else
                {
                    var degreeOfParallelism = NumberOfCPUs;
                    var tasks = new Task[degreeOfParallelism];
                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = tArrayUncoupleUnconstraintInNew.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = tArrayUncoupleUnconstraintInNew.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     uGlobal[tArrayUnCoupleAndUnconstraint[i]] = uLocalU[i];
                                 }
                             });
                    }
                    Task.WaitAll(tasks);
                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = tArrayIndependentInNew.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = tArrayIndependentInNew.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     uGlobal[tArrayIndependent[i]] = uLocalM[i];
                                 }
                             });
                    }
                    Task.WaitAll(tasks);

                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = tArrayDependent.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = tArrayDependent.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     uGlobal[tArrayDependent[i]] = uLocalS[i];
                                 }
                             });
                    }
                    Task.WaitAll(tasks);
                }
            }

            return uGlobal;
        }

        protected DoubleVector SolveStatic(DoubleSymmetricMatrix kGlobalDense, DoubleCsrSparseMatrix kGlobalSparse, ref DoubleVector rGlobal)
        {
            DoubleVector uGlobal;
            if (IsSparseData)
            {
                ApplyTAndg(ref kGlobalSparse, ref rGlobal, null);

                uGlobal = MatrixFunctions.Solve(kGlobalSparse, rGlobal);
                uGlobal = MatrixFunctions.Product(TSparse, uGlobal) + g;
            }
            else
            {
                //ApplyTAndg(ref kGlobalDense, ref rGlobal, null);
                //uGlobal = MatrixFunctions.Solve(kGlobalDense, rGlobal);
                //uGlobal = MatrixFunctions.Product(TDense, uGlobal) + g;

                int[] tArrayIndependent = this.tArrayIndependent.ToArray();
                int[] tArrayDependent = this.tArrayDependent.ToArray();
                int[] tArrayUnCoupleAndUnconstraint = this.tArrayUncouple.ToArray();
                int[] tArrayIndependentInNew = new int[tArrayIndependent.Length];
                int c = 0;
                for (int i = tArrayUnCoupleAndUnconstraint.Length; i < tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length; i++)
                {
                    tArrayIndependentInNew[c++] = i;
                }

                int[] tArrayUncoupleUnconstraintInNew = new int[tArrayUnCoupleAndUnconstraint.Length];
                for (int i = 0; i < tArrayUnCoupleAndUnconstraint.Length; i++)
                {
                    tArrayUncoupleUnconstraintInNew[i] = i;
                }
                DoubleVector gS = null;
                if (tArrayDependent.Length != 0)
                {
                    gS = MatrixTool.GetSubVector(g, tArrayDependent, AbstractModel.NumberOfCPUs);
                }

                DoubleVector Ru = null;
                DoubleVector Rm = null;
                DoubleVector rhs = new DoubleVector(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);

                DoubleSymmetricMatrix Kuu = null;
                DoubleMatrix Kus = null;
                DoubleMatrix Kum = null;
                DoubleMatrix Kmm = null;
                if (tArrayUnCoupleAndUnconstraint.Length != 0)
                {
                    Kuu = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnCoupleAndUnconstraint, AbstractModel.NumberOfCPUs);
                    Ru = MatrixTool.GetSubVector(rGlobal, tArrayUnCoupleAndUnconstraint, AbstractModel.NumberOfCPUs);
                    if (tArrayDependent.Length != 0)
                    {
                        Kus = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnCoupleAndUnconstraint, tArrayDependent, AbstractModel.NumberOfCPUs);
                        Ru += -(NMathFunctions.Product(Kus, gS));
                    }
                }

                DoubleMatrix TsmDense = null;
                if (tArrayDependent.Length != 0 && tArrayIndependent.Length != 0)
                {
                    TsmDense = MatrixTool.GetSubMatrix(TDense, tArrayDependent, tArrayIndependent, AbstractModel.NumberOfCPUs);
                    DoubleMatrix TsmTranspose = TsmDense.Transpose();
                    DoubleMatrix Kms = MatrixTool.GetSubMatrix(kGlobalDense, tArrayIndependent, tArrayDependent, AbstractModel.NumberOfCPUs);
                    DoubleMatrix Kss = MatrixTool.GetSubMatrix(kGlobalDense, tArrayDependent, tArrayDependent, AbstractModel.NumberOfCPUs);
                    DoubleVector Rs = MatrixTool.GetSubVector(rGlobal, tArrayDependent, AbstractModel.NumberOfCPUs);
                    Kmm = MatrixTool.GetSubMatrix(kGlobalDense, tArrayIndependent, tArrayIndependent, AbstractModel.NumberOfCPUs) + NMathFunctions.Product(TsmTranspose, Kms.Transpose() + NMathFunctions.Product(Kss, TsmDense)) + NMathFunctions.Product(Kms, TsmDense);
                    Rm = MatrixTool.GetSubVector(rGlobal, tArrayIndependent, AbstractModel.NumberOfCPUs) + NMathFunctions.Product(TsmTranspose, Rs);
                    if (tArrayUnCoupleAndUnconstraint.Length != 0)
                        Kum = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnCoupleAndUnconstraint, tArrayIndependent, AbstractModel.NumberOfCPUs) + NMathFunctions.Product(Kus, TsmDense);
                    Rm += -(NMathFunctions.Product(Kms, gS) + NMathFunctions.Product(TsmTranspose, NMathFunctions.Product(Kss, gS)));
                }
                DoubleMatrix lhsDense = new DoubleMatrix(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length, tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
                //rhs = new DoubleVector(tArrayUnCoupleAndUnconstraint.Length + tArrayIndependent.Length);
                if (!AbstractModel.IsParallelProcesing)
                {
                    for (int i = 0; i < tArrayUnCoupleAndUnconstraint.Length; i++)
                    {
                        for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
                        {
                            lhsDense[i, j] = Kuu[i, j];
                        }
                        if (Kum != null)
                        {
                            for (int j = 0; j < tArrayIndependent.Length; j++)
                            {
                                lhsDense[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
                            }
                        }
                        rhs[i] = Ru[i];
                    }
                    for (int i = 0; i < tArrayIndependent.Length; i++)
                    {
                        for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
                        {
                            lhsDense[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
                        }

                        for (int j = 0; j < tArrayIndependent.Length; j++)
                        {
                            lhsDense[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
                        }
                        rhs[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
                    }
                }
                else
                {
                    var degreeOfParallelism = NumberOfCPUs;
                    var tasks = new Task[degreeOfParallelism];
                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = tArrayUnCoupleAndUnconstraint.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = tArrayUnCoupleAndUnconstraint.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
                                     {
                                         lhsDense[i, j] = Kuu[i, j];
                                     }
                                     if (Kum != null)
                                     {
                                         for (int j = 0; j < tArrayIndependent.Length; j++)
                                         {
                                             lhsDense[i, j + tArrayUnCoupleAndUnconstraint.Length] = Kum[i, j];
                                         }
                                     }
                                     rhs[i] = Ru[i];
                                 }
                             });
                    }
                    Task.WaitAll(tasks);

                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = tArrayIndependent.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = tArrayIndependent.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     for (int j = 0; j < tArrayUnCoupleAndUnconstraint.Length; j++)
                                     {
                                         lhsDense[i + tArrayUnCoupleAndUnconstraint.Length, j] = Kum[j, i];
                                     }

                                     for (int j = 0; j < tArrayIndependent.Length; j++)
                                     {
                                         lhsDense[i + tArrayUnCoupleAndUnconstraint.Length, j + tArrayUnCoupleAndUnconstraint.Length] = Kmm[i, j];
                                     }
                                     rhs[i + tArrayUnCoupleAndUnconstraint.Length] = Rm[i];
                                 }
                             });
                    }
                    Task.WaitAll(tasks);
                }

                DoubleVector uLocalS = gS;
                DoubleVector uLocalUM = null;
                DoubleVector uLocalU = null;
                DoubleVector uLocalM = null;
                if (rGlobal.Length != 0)
                {
                    uLocalUM = MatrixFunctions.Solve(lhsDense, rhs);
                    uLocalU = MatrixTool.GetSubVector(uLocalUM, tArrayUncoupleUnconstraintInNew, AbstractModel.NumberOfCPUs);
                    uLocalM = MatrixTool.GetSubVector(uLocalUM, tArrayIndependentInNew, AbstractModel.NumberOfCPUs);
                    if (uLocalM != null)
                    {
                        uLocalS += MatrixFunctions.Product(TsmDense, uLocalM);
                    }
                }

                uGlobal = new DoubleVector(countDOF);
                if (!AbstractModel.IsParallelProcesing)
                {
                    for (int i = 0; i < tArrayUncoupleUnconstraintInNew.Length; i++)
                    {
                        uGlobal[tArrayUnCoupleAndUnconstraint[i]] = uLocalU[i];
                    }
                    for (int i = 0; i < tArrayIndependentInNew.Length; i++)
                    {
                        uGlobal[tArrayIndependent[i]] = uLocalM[i];
                    }
                    for (int i = 0; i < tArrayDependent.Length; i++)
                    {
                        uGlobal[tArrayDependent[i]] = uLocalS[i];
                    }
                }
                else
                {
                    var degreeOfParallelism = NumberOfCPUs;
                    var tasks = new Task[degreeOfParallelism];
                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = tArrayUncoupleUnconstraintInNew.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = tArrayUncoupleUnconstraintInNew.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     uGlobal[tArrayUnCoupleAndUnconstraint[i]] = uLocalU[i];
                                 }
                             });
                    }
                    Task.WaitAll(tasks);
                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = tArrayIndependentInNew.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = tArrayIndependentInNew.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     uGlobal[tArrayIndependent[i]] = uLocalM[i];
                                 }
                             });
                    }
                    Task.WaitAll(tasks);

                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                             () =>
                             {
                                 var max = tArrayDependent.Length * (taskNumberCopy + 1) / degreeOfParallelism;
                                 for (int i = tArrayDependent.Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                 {
                                     uGlobal[tArrayDependent[i]] = uLocalS[i];
                                 }
                             });
                    }
                    Task.WaitAll(tasks);
                }
            }

            return uGlobal;
        }
        protected void WriteInformationProblem(string nameTypeProblem)
        {
            if (IsWriteLogFile)
            {
                IOIGA.WriteLogFileAndConsole(logFile, true, "-----------------------------------------------------------");
                IOIGA.WriteLogFileAndConsole(logFile, true, nameTypeProblem);
                IOIGA.WriteLogFileAndConsole(logFile, true, "-----------------------------------------------------------");
                IOIGA.WriteLogFileAndConsole(logFile, true, "Information of problem:");
                IOIGA.WriteLogFileAndConsole(logFile, true, "Date: " + DateTime.Now);
                IOIGA.WriteLogFileAndConsole(logFile, true, (IsParallelProcesing) ? "Parallel computing with " + NumberOfCPUs + " CPUs" : "Serial computing");
                IOIGA.WriteLogFileAndConsole(logFile, true, "Size of problem: " + countDOF + " dofs");
                string str = "";
                for (int i = 0; i < typeOfFieldsInMultifield.Count; i++)
                {
                    str += typeOfFieldsInMultifield[i].ToString() + " ";
                }
                IOIGA.WriteLogFileAndConsole(logFile, true, "Type of fields: " + str);
                str = "";
                for (int i = 0; i < listComputeResult.Count; i++)
                {
                    str += listComputeResult[i].ToString() + " ";
                }
                IOIGA.WriteLogFileAndConsole(logFile, true, "Variables: " + str);
            }
        }


    }
}
