using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using DEMSoft.Drawing;
using System.IO;
using System.Linq;
using DEMSoft.Function;
using CenterSpace.NMath.Core;
using MathNet.Numerics.Data.Matlab;
using MathNet.Numerics.LinearAlgebra;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
    /// <summary>
    /// Structure 2D problem (plane stress, plane strain). Default is plane stress
    /// </summary>
    public class ModelStructureModal : AbstractModelStructure
    {
        private double[] eigenValue;
        private DoubleMatrix eigenVector;
        private int numberOfMode;
        /// <summary>
        /// Constructor class
        /// </summary>
        public ModelStructureModal(Dimension structureDimension, string pathProject = "temp", string nameProject = "project")
                     : base(TypeAnalysisModel.Modal, structureDimension, pathProject, nameProject)
        {
            numberOfMode = 6;
        }
        public override void PreProcessing()
        {
            base.PreProcessing();
            tArray = new int[typeOfFieldsInMultifield.Count][];
            tArrayConstraint = new int[typeOfFieldsInMultifield.Count][];
            tArrayUnconstraint = new int[typeOfFieldsInMultifield.Count][];
            if (!IsParallelProcesing)
            {
                tArray[0] = GetTArrayGlobalMechanic();
                if (bcMatrix != null)
                {
                    tArrayUnconstraint[0] = GetTArrayGlobalMechanicUnConstraint();
                    tArrayConstraint[0] = GetTArrayGlobalMechanicConstraint();
                }
                else
                {
                    tArrayUnconstraint[0] = tArray[0];
                }
            }
            else
            {
                if (bcMatrix != null)
                {
                    Task<int[]>[] tasks ={
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalMechanic()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalMechanicUnConstraint()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalMechanicConstraint())};
                    //Block until all tasks complete.
                    Task.WaitAll(tasks);
                    tArray[0] = tasks[0].Result;
                    tArrayUnconstraint[0] = tasks[1].Result;
                    tArrayConstraint[0] = tasks[2].Result;
                }
                else
                {
                    tArrayUnconstraint[0] = GetTArrayGlobalMechanic();
                }
            }
        }
        private int[] GetTArrayGlobalMechanic()
        {
            List<int> listTArrayGlobalMechanic = new List<int>();
            for (int i = 0; i < listPatch.Count; i++)
            {
                AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
                int nen = pa.GetCountGlobalBasisFunctions();
                for (int k = 0; k < nen; k++)
                {
                    for (int j = 0; j < (int)StructureDimension; j++)
                    {
                        int temp = pa.GetIDInGlobal(j, k);
                        if (temp != -1)
                        {
                            if (listTArrayGlobalMechanic.IndexOf(temp) == -1)
                            {
                                listTArrayGlobalMechanic.Add(temp);
                            }
                        }
                    }
                }
            }
            return listTArrayGlobalMechanic.ToArray();
        }
        private int[] GetTArrayGlobalMechanicConstraint()
        {
            List<int> listTArrayGlobalMechanic = new List<int>();
            for (int i = 0; i < listPatch.Count; i++)
            {
                AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
                int nen = pa.GetCountGlobalBasisFunctions();
                for (int k = 0; k < nen; k++)
                {
                    for (int j = 0; j < (int)StructureDimension; j++)
                    {
                        int temp = pa.GetIDInGlobal(j, k);
                        if (temp != -1)
                        {
                            int idxConstraint = MatrixTool.IndexOfCol(bcMatrix, 0, temp, AbstractModel.NumberOfCPUs);
                            if (idxConstraint != -1)
                                if ((listTArrayGlobalMechanic.IndexOf(temp) == -1) && (bcMatrix[idxConstraint, 2] != (int)StructureDimension/* khong phai la bac tu do volt*/))
                                {
                                    listTArrayGlobalMechanic.Add(temp);
                                }
                        }
                    }
                }
            }
            return listTArrayGlobalMechanic.ToArray();
        }

        private int[] GetTArrayGlobalMechanicUnConstraint()
        {
            List<int> listTArrayGlobalMechanic = new List<int>();
            for (int i = 0; i < listPatch.Count; i++)
            {
                AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
                int nen = pa.GetCountGlobalBasisFunctions();
                for (int k = 0; k < nen; k++)
                {
                    for (int j = 0; j < (int)StructureDimension; j++)
                    {
                        int temp = pa.GetIDInGlobal(j, k);
                        if (temp != -1)
                        {
                            int idxConstraint = MatrixTool.IndexOfCol(bcMatrix, 0, temp, AbstractModel.NumberOfCPUs);
                            if ((listTArrayGlobalMechanic.IndexOf(temp) == -1) && (idxConstraint == -1))
                            {
                                listTArrayGlobalMechanic.Add(temp);
                            }
                        }
                    }
                }
            }
            return listTArrayGlobalMechanic.ToArray();
        }

        /// <summary>
        /// Solve problem
        /// </summary>
        public override void Solve()
        {
            WriteInformationProblem("Modal structural module");

            if (listPatch.Count == 0)
                throw new NullReferenceException("Model must be assigned Patch");
            DoubleMatrix kGlobalDense = null; // stiffness
            DoubleMatrix mGlobalDense = null; // "mass"
            int[] tArrayUnCon = null;//tArrayUnconstraint[0];//null;
                                     //tArrayUnconstraint[0] cho ket qua dang rieng tot hon
                                     // Assembly
            if (!IsSparseData)
            {
                AssemblyStiffnessMatrix(out kGlobalDense);
                AssemblyMassMatrix(out mGlobalDense);
                ApplyTAndg(ref kGlobalDense, ref mGlobalDense, tArrayUnCon);
            }
            else
            {
                AssemblyStiffnessMatrix(ref kGlobalSparse);
                AssemblyMassMatrix(ref mGlobalSparse);
                ApplyTAndg(ref kGlobalSparse, ref mGlobalSparse, tArrayUnCon);
            }

            //bool IsMatlabSolver = true;
            if (IsMatlabSolver)
            {
                ////////////////////////////////////////////////////////
                ///////// Matlab ////////////////
                ////////////////////////////////////////////////////////        
                if (matlab == null)
                    ////// Create the MATLAB instance 
                    matlab = new MLApp.MLApp();

                var matrices = new List<MatlabMatrix>();
                matrices.Add(MatlabWriter.Pack(MatrixTool.ConvertNMathMatrixToMathNetMatrix(kGlobalDense), "K"));
                matrices.Add(MatlabWriter.Pack(MatrixTool.ConvertNMathMatrixToMathNetMatrix(mGlobalDense), "M"));

                Matrix<double> bcdof = MatrixTool.ConvertIntArrayToMathNetMatrix(tArrayConstraint[0]);
                for (int i = 0; i < bcdof.RowCount; i++)
                {
                    bcdof[i, 0]++;
                }

                matrices.Add(MatlabWriter.Pack(bcdof, "bcdof"));

                Matrix<double> numMode = Matrix<double>.Build.Dense(1, 1);
                numMode[0, 0] = numberOfMode;
                matrices.Add(MatlabWriter.Pack(numMode, "numberOfMode"));

                if (File.Exists("file.mat"))
                    File.Delete("file.mat");
                MatlabWriter.Store("file.mat", matrices);

                ////////////////////////////////////////////////////////////
                ////////////////// Store without file.mat //////////////////
                ////////////////////////////////////////////////////////////
                matlab.Execute(@"cd '" + Directory.GetCurrentDirectory() + "'");
                matlab.Execute("userpath('" + Directory.GetCurrentDirectory() + "')");
                matlab.Execute("clear all");
                matlab.Execute("load('file.mat')");
                matlab.Execute("global numMode");
                matlab.Execute("numMode=numberOfMode");
                matlab.Execute("[lamda, modeshape]= eigens(K, M, bcdof);");
                matlab.Execute("lamda=sqrt(lamda(1:numberOfMode));");
                matlab.Execute("modeshape=modeshape(:,1:numberOfMode);");
                //matlab.Execute("[lamda, modeshape]= eigens(K, M);");
                eigenValue = new DoubleMatrix(matlab.GetVariable("lamda", "base")).Col(0).ToArray();
                eigenVector = new DoubleMatrix(matlab.GetVariable("modeshape", "base"));
                ///////////////////////////////////////
            }
            else
            {
                eigenValue = new double[numberOfMode];
                eigenVector = new DoubleMatrix(countDOF, numberOfMode);
                bool isSymmetryMatrix = false;
                if (!IsSparseData)
                {
                    isSymmetryMatrix = MatrixTool.IsSymmetricMatrix(kGlobalDense) && MatrixTool.IsSymmetricMatrix(mGlobalDense);
                }
                else
                { isSymmetryMatrix = MatrixTool.IsSymmetricMatrix(kGlobalSparse) && MatrixTool.IsSymmetricMatrix(mGlobalSparse); }

                if (isSymmetryMatrix)
                {
                    ArnoldiEigenvalueOptions opt = new ArnoldiEigenvalueOptions();
                    opt.FindEigenvaluesGreaterThan(0);
                    int numOfConstraint = tArrayDependent.Count;//tArrayDependent.Count;// tArrayConstraint[0].Length;
                    if (tArrayUnCon != null)
                        numOfConstraint = 0;
                    opt.NumEigenvalues = numberOfMode + numOfConstraint;

                    DoubleVector eigenValueTemp = null;
                    DoubleMatrix eigenVectorTemp = null;
                    ArnoldiEigenvalueSolution sol = null;
                    if (!IsSparseData)
                    {
                        sol = ArnoldiEigenvalueSolver.Solve(new DoubleSymmetricMatrix(kGlobalDense), new DoubleSymmetricMatrix(mGlobalDense), opt);
                    }
                    else
                    {
                        DoubleSymCsrSparseMatrix k = MatrixTool.ConvertDoubleCsrSparseMatrixToDoubleSymCsrSparseMatrix(kGlobalSparse);
                        DoubleSymCsrSparseMatrix m = MatrixTool.ConvertDoubleCsrSparseMatrixToDoubleSymCsrSparseMatrix(mGlobalSparse);
                        sol = ArnoldiEigenvalueSolver.Solve(k, m, opt);
                    }
                    eigenValueTemp = sol.Eigenvalues;
                    eigenVectorTemp = sol.Eigenvectors;
                    int i = 0;
                    int c = 0;
                    //if (tArrayUnCon == null)
                    //{
                    while (c < numberOfMode && i < eigenValueTemp.Length)
                    {
                        //Check is constraint dof
                        bool isNonzero = true;
                        if (eigenValueTemp[i] - 1.0 <= 1e-2)//(Math.Abs(eigenVectorTemp[tArrayDependent[0], i]) > 1e-8)
                            isNonzero = false;
                        if (isNonzero)
                        {
                            eigenValue[c] = Math.Sqrt(eigenValueTemp[i]);
                            if (!AbstractModel.IsParallelProcesing)
                            {
                                DoubleVector vectemp = null;
                                if (tArrayUnCon != null)
                                {
                                    DoubleMatrix TReduce = MatrixTool.GetSubMatrix(TDense, tArray[0], tArrayUnconstraint[0], AbstractModel.NumberOfCPUs);
                                    vectemp = MatrixFunctions.Product(TReduce, eigenVectorTemp.Col(i));
                                }
                                for (int j = 0; j < countDOF; j++)
                                {
                                    if (tArrayUnCon == null)
                                    {
                                        eigenVector[j, c] = eigenVectorTemp[j, i];
                                    }
                                    else
                                    {
                                        eigenVector[j, c] = vectemp[j];
                                    }
                                }
                            }
                            else
                            {
                                var degreeOfParallelism = Environment.ProcessorCount;
                                var tasks = new Task[degreeOfParallelism];

                                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                                {
                                    // capturing taskNumber in lambda wouldn't work correctly
                                    int taskNumberCopy = taskNumber;

                                    tasks[taskNumber] = Task.Factory.StartNew(
                                        () =>
                                        {
                                            var max = countDOF * (taskNumberCopy + 1) / degreeOfParallelism;
                                            for (int j = countDOF * taskNumberCopy / degreeOfParallelism; j < max; j++)
                                            {
                                                eigenVector[j, c] = eigenVectorTemp[j, i];
                                            }
                                        });
                                }
                                Task.WaitAll(tasks);
                            }
                            c++;
                        }
                        i++;
                    }
                    ////normalize eigen vector, max value is set equal 1
                    //for (i = 0; i < eigenVector.Cols; i++)
                    //{
                    //  double maxValAbs = Math.Abs(eigenVector[0, i]);
                    //  int indexMaxVal = 0;
                    //  for (int j = 1; j < countDOF; j++)
                    //  {
                    //    if (maxValAbs < Math.Abs(eigenVector[j, i]))
                    //    {
                    //      maxValAbs = Math.Abs(eigenVector[j, i]);
                    //      indexMaxVal = j;
                    //    }
                    //  }
                    //  double maxVal = eigenVector[indexMaxVal, i];
                    //  if (!AbstractModel.IsParallelProcesing)
                    //  {
                    //    for (int j = 0; j < countDOF; j++)
                    //    {
                    //      eigenVector[j, i] = eigenVector[j, i] / maxVal;
                    //    }
                    //  }
                    //  else
                    //  {
                    //    var degreeOfParallelism = Environment.ProcessorCount;
                    //    var tasks = new Task[degreeOfParallelism];

                    //    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    //    {
                    //      // capturing taskNumber in lambda wouldn't work correctly
                    //      int taskNumberCopy = taskNumber;

                    //      tasks[taskNumber] = Task.Factory.StartNew(
                    //          () =>
                    //          {
                    //            var max = countDOF * (taskNumberCopy + 1) / degreeOfParallelism;
                    //            for (int j = countDOF * taskNumberCopy / degreeOfParallelism; j < max; j++)
                    //            {
                    //              eigenVector[j, i] = eigenVector[j, i] / maxVal;
                    //            }
                    //          });
                    //    }
                    //    Task.WaitAll(tasks);
                    //  }
                    //}
                }
                else
                {
                    DoubleMatrix A = MatrixFunctions.Product(MatrixFunctions.Inverse(mGlobalDense), kGlobalDense);
                    var eigDecomp = new DoubleEigDecomp(A);
                    //eigDecomp.FactorNoPreconditioning(A);
                    var eigenValueTemp1 = eigDecomp.EigenValues;
                    var eigenVectorTemp1 = eigDecomp.RightEigenVectors;
                    var eigenValueTempOrder = eigenValueTemp1.OrderBy(item => item.Real);

                    int i = 0;
                    int c = 0;
                    //if (tArrayUnCon == null)
                    //{
                    while (c < numberOfMode && i < eigenValueTempOrder.Count())
                    {
                        //Check is constraint dof
                        bool isNonzero = true;
                        if (eigenValueTempOrder.ElementAt(i).Real - 1.0 <= 1e-2)//(Math.Abs(eigenVectorTemp[tArrayDependent[0], i]) > 1e-8)
                            isNonzero = false;
                        if (isNonzero)
                        {
                            int idx = IndexOfEnumerable(eigenValueTemp1, eigenValueTempOrder.ElementAt(i));
                            eigenValue[c] = Math.Sqrt(eigenValueTemp1.ElementAt(idx).Real);
                            DoubleComplexVector vec = eigenVectorTemp1.Col(idx);
                            if (!AbstractModel.IsParallelProcesing)
                            {
                                //DoubleVector vectemp = null;
                                //if (tArrayUnCon != null)
                                //{
                                //  DoubleMatrix TReduce = MatrixTool.GetSubMatrix(TDense, tArray[0], tArrayUnconstraint[0]);
                                //  vectemp = MatrixFunctions.Product(TReduce, eigenVectorTemp.Col(i));
                                //}
                                for (int j = 0; j < countDOF; j++)
                                {
                                    //if (tArrayUnCon == null)
                                    //{
                                    eigenVector[j, c] = vec[j].Real;
                                    //}
                                    //else
                                    //{
                                    //  eigenVector[j, c] = vectemp[j];
                                    //}
                                }
                            }
                            else
                            {
                                var degreeOfParallelism = Environment.ProcessorCount;
                                var tasks = new Task[degreeOfParallelism];

                                for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                                {
                                    // capturing taskNumber in lambda wouldn't work correctly
                                    int taskNumberCopy = taskNumber;

                                    tasks[taskNumber] = Task.Factory.StartNew(
                                        () =>
                                        {
                                            var max = countDOF * (taskNumberCopy + 1) / degreeOfParallelism;
                                            for (int j = countDOF * taskNumberCopy / degreeOfParallelism; j < max; j++)
                                            {
                                                eigenVector[j, c] = vec[j].Real;
                                            }
                                        });
                                }
                                Task.WaitAll(tasks);
                            }
                            c++;
                        }
                        i++;
                    }
                }
            }

            ////normalize eigen vector, max value is set equal 1
            for (int i = 0; i < eigenVector.Cols; i++)
            {
                double maxValAbs = Math.Abs(eigenVector[0, i]);
                int indexMaxVal = 0;
                for (int j = 1; j < countDOF; j++)
                {
                    if (maxValAbs < Math.Abs(eigenVector[j, i]))
                    {
                        maxValAbs = Math.Abs(eigenVector[j, i]);
                        indexMaxVal = j;
                    }
                }
                double maxVal = eigenVector[indexMaxVal, i];
                if (!AbstractModel.IsParallelProcesing)
                {
                    for (int j = 0; j < countDOF; j++)
                    {
                        eigenVector[j, i] = eigenVector[j, i] / maxVal;
                    }
                }
                else
                {
                    var degreeOfParallelism = Environment.ProcessorCount;
                    var tasks = new Task[degreeOfParallelism];

                    for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
                    {
                        // capturing taskNumber in lambda wouldn't work correctly
                        int taskNumberCopy = taskNumber;

                        tasks[taskNumber] = Task.Factory.StartNew(
                            () =>
                            {
                                var max = countDOF * (taskNumberCopy + 1) / degreeOfParallelism;
                                for (int j = countDOF * taskNumberCopy / degreeOfParallelism; j < max; j++)
                                {
                                    eigenVector[j, i] = eigenVector[j, i] / maxVal;
                                }
                            });
                    }
                    Task.WaitAll(tasks);
                }
            }
            //IOIGA.WriteExcelResult(eigenValue, PathProject + "eig_" + NameProject + ".xls");
            IOIGA.WriteMessagePackDoubleVector(PathProject + "eig_" + NameProject + ".msgpack", new DoubleVector(eigenValue));
            IOIGA.WriteDoubleMatrixMessagePackFormatter(PathProject + "vec_" + NameProject + ".msgpack", eigenVector);
        }
        public static int IndexOfEnumerable(IEnumerable<DoubleComplex> data, DoubleComplex value)
        {
            for (int i = 0; i < data.Count(); i++)
            {
                if (data.ElementAt(i).Equals(value))
                    return i;
            }
            return -1;
        }
        public override void PostProcessing()
        {
            //eigenValue = IOIGA.ReadExcelVectorResult(PathProject + "eig_" + NameProject + ".xls").ToArray();
            eigenValue = IOIGA.ReadMessagePackDoubleVector(PathProject + "eig_" + NameProject + ".msgpack").ToArray();
            eigenVector = IOIGA.ReadDoubleMatrixMessagePackFormatter(PathProject + "vec_" + NameProject + ".msgpack");

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

        public void DrawModeResult(int idMode, ViewerForm viewer, int resolution, double scale)
        {
            if (StructureDimension == Dimension.Plane)
            {
                Contours2D[] contours = new Contours2D[listPatch.Count];
                double minVal = 1000000000;
                double maxVal = -1000000000;

                DoubleVector uGlobal = new DoubleVector(countDOF);
                if (!IsParallelProcesing)
                {
                    for (int i = 0; i < countDOF; i++)
                    {
                        uGlobal[i] = eigenVector[i, idMode];
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
                                var max = countDOF * (taskNumberCopy + 1) / degreeOfParallelism;
                                for (int i = countDOF * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                {
                                    uGlobal[i] = eigenVector[i, idMode];
                                }
                            });
                    }
                    Task.WaitAll(tasks);
                }

                if (IsSparseData)
                {
                    uGlobal = MatrixFunctions.Product(TSparse, uGlobal) + g;
                }
                else
                {
                    uGlobal = MatrixFunctions.Product(TDense, uGlobal) + g;
                }
                //Distribute displacement to control points
                SetUGlobal(uGlobal);

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    AbstractPatch2D pa = (AbstractPatch2D)listPatch[lp];
                    double[,][] data = pa.GetDeformationSurface(scale).GetDataOnSurface(resolution, resolution);
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
                    contours[lp].SetScalarValue(pa.CalculateResult(Result.USUM, resolution));
                    //contours[lp].SetOpacity(10);
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
                viewer.SetColormapBarVisible(contours[0], "Free frequency: " + eigenValue[idMode]);
            }
            else if (StructureDimension == Dimension.Solid)
            {
                Contours3D[] contours = new Contours3D[listPatch.Count];
                double minVal = 1000000000;
                double maxVal = -1000000000;

                DoubleVector uGlobal = new DoubleVector(countDOF);
                if (!IsParallelProcesing)
                {
                    for (int i = 0; i < countDOF; i++)
                    {
                        uGlobal[i] = eigenVector[i, idMode];
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
                                var max = countDOF * (taskNumberCopy + 1) / degreeOfParallelism;
                                for (int i = countDOF * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                {
                                    uGlobal[i] = eigenVector[i, idMode];
                                }
                            });
                    }
                    Task.WaitAll(tasks);
                }
                if (IsSparseData)
                {
                    uGlobal = MatrixFunctions.Product(TSparse, uGlobal) + g;
                }
                else
                {
                    uGlobal = MatrixFunctions.Product(TDense, uGlobal) + g;
                }
                //Distribute displacement to control points
                SetUGlobal(uGlobal);

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    AbstractPatch3D pa = (AbstractPatch3D)listPatch[lp];
                    double[,,][] data = pa.GetDeformationVolume(scale).GetDataOnVolume(resolution, resolution, resolution);
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
                    contours[lp].SetScalarValue(pa.CalculateResult(Result.USUM, resolution));
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
                viewer.SetColormapBarVisible(contours[0], "Free frequency: " + eigenValue[idMode]);
            }
            else if (StructureDimension == Dimension.Plate)
            {
                Contours2D[] contours = new Contours2D[listPatch.Count];
                double minVal = 1000000000;
                double maxVal = -1000000000;

                DoubleVector uGlobal = new DoubleVector(countDOF);
                if (!IsParallelProcesing)
                {
                    for (int i = 0; i < countDOF; i++)
                    {
                        uGlobal[i] = eigenVector[i, idMode];
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
                                var max = countDOF * (taskNumberCopy + 1) / degreeOfParallelism;
                                for (int i = countDOF * taskNumberCopy / degreeOfParallelism; i < max; i++)
                                {
                                    uGlobal[i] = eigenVector[i, idMode];
                                }
                            });
                    }
                    Task.WaitAll(tasks);
                }

                if (IsSparseData)
                {
                    uGlobal = MatrixFunctions.Product(TSparse, uGlobal) + g;
                }
                else
                {
                    uGlobal = MatrixFunctions.Product(TDense, uGlobal) + g;
                }
                //Distribute displacement to control points
                SetUGlobal(uGlobal);

                for (int lp = 0; lp < listPatch.Count; lp++)
                {
                    AbstractPatch2D pa = (AbstractPatch2D)listPatch[lp];
                    FunctionRToR fz = ((PatchStructurePlate)pa).KinematicsFunction;
                    double[,][] data = pa.GetDeformationSurface(scale).GetDataOnSurface(resolution, resolution);
                    double[,] xData = new double[data.GetLength(0), data.GetLength(1)];
                    double[,] yData = new double[data.GetLength(0), data.GetLength(1)];
                    double[,] zData = new double[data.GetLength(0), data.GetLength(1)];

                    for (int i = 0; i < data.GetLength(0); i++)
                        for (int j = 0; j < data.GetLength(1); j++)
                        {
                            xData[i, j] = data[i, j][0];
                            yData[i, j] = data[i, j][1];
                            zData[i, j] = data[i, j][2];
                            //if (TypePlate == TypePlate.FSDT)
                            //{
                            //    xData[i, j] = data[i, j][0];
                            //    yData[i, j] = data[i, j][1];
                            //    zData[i, j] = data[i, j][2];
                            //}
                            //else
                            //{
                            //    xData[i, j] = data[i, j][0] + fz.ValueAt(0) * data[i, j][3];
                            //    yData[i, j] = data[i, j][1] + fz.ValueAt(0) * data[i, j][4];
                            //    zData[i, j] = data[i, j][2];
                            //}
                        }
                    contours[lp] = new Contours2D(xData, yData, zData);
                    contours[lp].SetScalarValue(pa.CalculateResult(Result.UZ, resolution));
                    //contours[lp].SetOpacity(10);

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
                viewer.SetColormapBarVisible(contours[0], "Free frequency: " + eigenValue[idMode]);
            }
        }

        public double[] GetEigenValue()
        { return eigenValue; }

        public double[,] GetEigenVector()
        { return eigenVector.ToArray(); }

        public void SetNumberOfMode(int n)
        { numberOfMode = n; }
    }
}
