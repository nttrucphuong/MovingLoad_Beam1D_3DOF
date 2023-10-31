using System;
using System.Threading.Tasks;
using DEMSoft.Drawing;
using System.IO;
using System.Linq;
using CenterSpace.NMath.Core;
using DEMSoft.Function;
using DEMSoft.Common;
//using MathNet.Numerics.Data.Matlab;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Structure 2D problem (plane stress, plane strain). Default is plane stress
  /// </summary>
  public class ModelPiezoelectricModal : AbstractModelPiezoelectric
  {
    private double[] eigenValue;
    private DoubleMatrix eigenVector;
    private int numberOfMode;
    /// <summary>
    /// Constructor class
    /// </summary>
    public ModelPiezoelectricModal(Dimension structureDimension, string pathProject = "temp", string nameProject = "project")
             : base(TypeAnalysisModel.Modal, structureDimension, pathProject, nameProject)
    {
      numberOfMode = 6;
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
      // Assembly
      if (!IsSparseData)
      {
        AssemblyStiffnessMatrix(out kGlobalDense);
        AssemblyMassMatrix(out mGlobalDense);
        ApplyTAndg(ref kGlobalDense, ref mGlobalDense, null);
      }
      else
      {
        AssemblyStiffnessMatrix(ref kGlobalSparse);
        AssemblyMassMatrix(ref mGlobalSparse);
        ApplyTAndg(ref kGlobalSparse, ref mGlobalSparse, null);
      }
      ///////////////////////////////////////
      ArnoldiEigenvalueOptions opt = new ArnoldiEigenvalueOptions();
      //opt.FindEigenvaluesOfSmallestMagnitude();

      var numOfConstraint = tArrayDependent.Count;// tArrayConstraint[0].Length;
      opt.NumEigenvalues = numberOfMode + numOfConstraint;
      ArnoldiEigenvalueSolution sol = null;
      if (!IsSparseData)
      {
        DoubleMatrix K11Dense = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnconstraint[0], tArrayUnconstraint[0], AbstractModel.NumberOfCPUs);
        DoubleMatrix K33Dense = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnconstraint[1], tArrayUnconstraint[1], AbstractModel.NumberOfCPUs);
        DoubleMatrix K31Dense = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnconstraint[1], tArrayUnconstraint[0], AbstractModel.NumberOfCPUs);
        DoubleMatrix K13Dense = MatrixTool.GetSubMatrix(kGlobalDense, tArrayUnconstraint[0], tArrayUnconstraint[1], AbstractModel.NumberOfCPUs);

        DoubleMatrix knewDense;
        if (K33Dense != null)
        {
          DoubleMatrix K33Inv = MatrixFunctions.PseudoInverse(K33Dense);
          knewDense = K11Dense - MatrixFunctions.Product(K13Dense, MatrixFunctions.Product(K33Inv, K31Dense));
        }
        else
        {
          knewDense = K11Dense;
        }

        DoubleMatrix MuuDense = MatrixTool.GetSubMatrix(mGlobalDense, tArrayUnconstraint[0], tArrayUnconstraint[0], AbstractModel.NumberOfCPUs);
        sol = ArnoldiEigenvalueSolver.Solve(new DoubleSymmetricMatrix(knewDense), new DoubleSymmetricMatrix(MuuDense), opt);
      }
      else
      {
        //DoubleCsrSparseMatrix K11Sparse = MatrixTool.GetSubMatrix(kGlobalSparse, tArrayUnconstraint[0], tArrayUnconstraint[0]);
        //DoubleCsrSparseMatrix K33Sparse = MatrixTool.GetSubMatrix(kGlobalSparse, tArrayUnconstraint[1], tArrayUnconstraint[1]);
        //DoubleCsrSparseMatrix K31Sparse = MatrixTool.GetSubMatrix(kGlobalSparse, tArrayUnconstraint[1], tArrayUnconstraint[0]);
        //DoubleCsrSparseMatrix K13Sparse = MatrixTool.GetSubMatrix(kGlobalSparse, tArrayUnconstraint[0], tArrayUnconstraint[1]);

        //DoubleCsrSparseMatrix knewSparse;
        //if (K33Sparse != null)
        //{
        //  DoubleMatrix K33Inv = MatrixFunctions.PseudoInverse(K33Sparse.ToDenseMatrix());
        //  knewSparse = K11Sparse - MatrixFunctions.Product(K13Sparse.ToDenseMatrix(), MatrixFunctions.Product(K33Inv, K31Sparse.ToDenseMatrix()));
        //}
        //else
        //{
        //  knewSparse = K11Sparse;
        //}

        //DoubleCsrSparseMatrix mNewSparse = MatrixTool.GetSubMatrix(mGlobalSparse, tArrayUnconstraint[0], tArrayUnconstraint[0]);
        //DoubleSymEigDecomp solver = new DoubleSymEigDecomp(NMathFunctions.Product(MatrixFunctions.Inverse(mNewSparse), knewSparse));
        //DoubleSymCsrSparseMatrix k = MatrixTool.ConvertDoubleCsrSparseMatrixToDoubleSymCsrSparseMatrix(knewSparse);
        //DoubleSymCsrSparseMatrix m = MatrixTool.ConvertDoubleCsrSparseMatrixToDoubleSymCsrSparseMatrix(mNewSparse);
        //sol = ArnoldiEigenvalueSolver.Solve(k, m, opt);
      }
      var eigenValueTemp = sol.Eigenvalues;
      var eigenVectorTemp = sol.Eigenvectors;
      eigenValue = new double[numberOfMode];
      eigenVector = new DoubleMatrix(countDOF, numberOfMode);
      int i = 0;
      int c = 0;
      while (c < numberOfMode)
      {
        //Check is constraint dof
        bool isNonzero = true;
        if (Math.Abs(eigenVectorTemp[tArrayDependent[0], i]) > 1e-8)
          isNonzero = false;
        if (isNonzero)
        {
          eigenValue[c] = Math.Sqrt(eigenValueTemp[i]);
          if (!AbstractModel.IsParallelProcesing)
          {
            for (int j = 0; j < countDOF; j++)
            {
              eigenVector[j, c] = eigenVectorTemp[j, i];
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

      IOIGA.WriteExcelResult(eigenValue, PathProject + "eig_" + NameProject + ".xls");
      IOIGA.WriteDoubleMatrixMessagePackFormatter(PathProject + "vec_" + NameProject + ".msgpack", eigenVector);
    }

    public override void PostProcessing()
    {
      eigenValue = IOIGA.ReadExcelVectorResult(PathProject + "eig_" + NameProject + ".xls", AbstractModel.NumberOfCPUs).ToArray();
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
