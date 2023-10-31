using System.Collections.Generic;
using System.Threading.Tasks;
using DEMSoft.NURBS;
using CenterSpace.NMath.Core;
using System;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Abstract of NURBS patch mesh
  /// </summary>
  public abstract class AbstractPatchOneField : AbstractPatch
  {
    public AbstractPatchOneField(int countDimension)
      : base(countDimension)
    {
    }
    /// <summary>
    /// Get surface 
    /// </summary>
    /// <returns></returns>
    public AbstractParametricGeometry GetGeometry()
    {
      return base.GetGeometry(0);
    }

    /// <summary>
    /// Get number of fields
    /// </summary>
    /// <returns></returns>
    public int GetCountField()
    {
      return base.GetCountField(0);
    }

    ///// <summary>
    ///// Apply boundary condition of edge
    ///// </summary>
    ///// <param name="indexEdgeFace"></param>
    ///// <param name="piecewiseLoad">Constraint function coressponse the time</param>
    ///// <param name="constraint">apply values of constraint dof, if no constraint, apply double.NaN</param>
    //public void ApplyBoundaryCondition(int indexEdgeFace, FunctionRToR piecewiseLoad, params double[] constraint)
    //{
    //   int numBound = 0;
    //   if (this is AbstractPatch2D)
    //      numBound = 4;
    //   else if (this is AbstractPatch3D)
    //      numBound = 6;
    //   if (this.constraint == null)
    //      this.constraint = new ConstraintValue[numBound];
    //   this.constraint[indexEdgeFace] = new ConstraintValue(piecewiseLoad, constraint);
    //}

    //public ConstraintValue GetConstraintOnEachEdgeFace(int indexEdgeFace)
    //{
    //   return constraint[indexEdgeFace];
    //}

    //public bool IsConstraintOnEdgeFace()
    //{ return constraint != null; }

    /// <summary>
    /// Create INC array
    /// </summary>
    protected override abstract void CreateINC();

    /// <summary>
    /// Create INC array
    /// </summary>
    protected override abstract void CreateIEN();

    /// <summary>
    /// Create IPN array (index of element)
    /// </summary>
    protected override abstract void CreateIPN();

    /// <summary>
    /// Create ID array
    /// </summary>
    public override abstract int[] EnumerateInPatch();//Enumerate DOF

    /// <summary>
    /// Create ID on Global (Multipatch)
    /// </summary>
    /// <param name="countDof"></param>
    /// <returns></returns>
    public override abstract int EnumerateInGlobalMultiPatch(int countDof);

    /// <summary>
    /// Get corresponding NURBS coordinate with a global basis function number and a parametric direction number
    /// </summary>
    /// <param name="idx">
    /// [0] ux - velocity x, uy - velocity y -----
    /// [1] theta - volume change, p - pressure</param>
    /// <param name="indexGlobalBasisFunction">Index of global basis function in 1D array</param>
    /// <param name="indexCoordinate">a parametric direction number</param>
    /// <returns></returns>
    public int GetINC(int indexGlobalBasisFunction, int indexCoordinate)
    {
      return base.GetINC(0, indexGlobalBasisFunction, indexCoordinate);
    }


    /// <summary>
    /// Find index of global basis function corresponding global coordinate
    /// </summary>
    /// <param name="idx">index of global basis functions in global coordinate</param>
    /// <returns></returns>
    public int FindIndexOfGlobalBasisFunction(params int[] idx)
    {
      return base.FindIndexOfGlobalBasisFunction(0, idx);
    }


    /// <summary>
    /// Get corresponding global basis function number with a local basis function number and an element number
    /// </summary>
    /// <param name="indexElement">index of element</param>
    /// <param name="indexLocalBasisFunction">a local basis function number in element</param>
    /// <returns></returns>
    public int GetIEN(int indexElement, int indexLocalBasisFunction)
    {
      return base.GetIEN(0, indexElement, indexLocalBasisFunction);
    }

    /// <summary>
    /// Get corresponding equation number with a global basis function number and a dof number.
    /// Return -1 then this DOF was constrainted.
    /// </summary>
    /// <param name="indexFields">index of 
    /// s</param>
    /// <param name="indexGlobalBasisFunction">index of a global basis function</param>
    /// <returns></returns>
    public int GetIDInPatch(int indexFields, int indexGlobalBasisFunction)
    {
      return GetEnumerateInPatch(0, indexFields, indexGlobalBasisFunction);
    }

    public int GetIDInGlobal(int indexFields, int indexGlobalBasisFunction)
    {
      return GetEnumerateInGlobal(0, indexFields, indexGlobalBasisFunction);
    }

    public void SetIDGlobal(int idDof, int indexFields, int indexGlobalBasisFunction)
    {
      base.SetEnumerateGlobal(0, idDof, indexFields, indexGlobalBasisFunction);
    }

    /// <summary>
    /// Get corresponding equation number with a dof number, a local basis function number and  an element number
    /// </summary>
    /// <param name="indexLocalBasisFunction">index of local basis function</param>
    /// <param name="indexElement">index of element</param>
    /// <param name="indexDOF">index of DOF</param>
    /// <returns></returns>
    public int GetLM(int indexLocalBasisFunction, int indexElement, int indexDOF)
    {
      return base.GetLM(0, indexLocalBasisFunction, indexElement, indexDOF);
    }

    /// <summary>
    /// Find index of element at coordinate
    /// </summary>
    /// <param name="xi"></param>
    /// <returns></returns>
    public abstract int FindIndexOfElementAt(params double[] xi);

    /// <summary>
    /// Select control points in region
    /// </summary>
    /// <param name="xLower">x-lower</param>
    /// <param name="xUpper">x-upper</param>
    /// <param name="yLower">y-lower</param>
    /// <param name="yUpper">y-upper</param>
    /// <param name="zLower">z-lower</param>
    /// <param name="zUpper">z-upper</param>
    /// <returns></returns>
    public ControlPoint[] SelectControlPointsByRegionBox(double xLower, double xUpper, double yLower, double yUpper, double zLower, double zUpper)
    {
      return base.SelectControlPoints(0, new RegionBox(xLower, xUpper, yLower, yUpper, zLower, zUpper));
    }
    public ControlPoint[] SelectControlPoints(IRegion loc)
    {
      return base.SelectControlPoints(0, loc);
    }
    public abstract ControlPoint[] GetAllControlPoints();

    ///// <summary>
    ///// Get number of patchs
    ///// </summary>
    ///// <returns></returns>
    //public override abstract int GetCountElement();

    /// <summary>
    /// Get number of local basis functions
    /// </summary>
    /// <returns></returns>
    public override abstract int GetCountLocalBasisFunctions(int idx = 0);

    /// <summary>
    /// Get number of global basis functions
    /// </summary>
    /// <returns></returns>
    public override abstract int GetCountGlobalBasisFunctions(int idx = 0);

    //public override abstract void SetUGlobal(DoubleVector uGlobal);

    /// <summary>
    /// Select serial of patchs from mesh by index on each direction
    /// </summary>
    /// <param name="index">[0], [1]: u-direction, [2], [3]: v-direction, [4], [5]: w-direction</param>
    /// <returns></returns>
    public override abstract List<AbstractElement> SelectEndPatchElement(int index);

    /// <summary>
    /// From enumerate number, determine index of Control point on Global and index of field of Control point
    /// </summary>
    /// <param name="enumerate">enumerate number</param>
    /// <param name="indexGlobalCps">index of Global control point on 1 array</param>
    /// <param name="indexField">index of field of control point</param>
    public void FindIndexofIDArray(int enumerate, ref int indexGlobalCps, ref int indexField)
    {
      //if (!AbstractModel.IsParallelProcesing)
      for (int i = 0; i < enumeratePatch[0].GetLength(0); i++)
      {
        for (int j = 0; j < enumeratePatch[0].GetLength(1); j++)
        {
          if (enumeratePatch[0][i, j] == enumerate)
          {
            indexField = i;
            indexGlobalCps = j;
            return;
          }
        }
      }
      //else
      //{
      //    int indexFieldt = -1;
      //    int indexGlobalCpst = -1;
      //    Parallel.For(0, enumeratePatch[0].GetLength(0), (i, state) =>
      //    {
      //        for (int j = 0; j < enumeratePatch[0].GetLength(1); j++)
      //        {
      //            if (enumeratePatch[0][i, j] == enumerate)
      //            {
      //                indexFieldt = i;
      //                indexGlobalCpst = j;
      //                state.Break();
      //            }
      //        }
      //    });
      //    indexField = indexFieldt;
      //    indexGlobalCps = indexGlobalCpst;
      //}
    }

    /// <summary>
    /// Select serial of patchs from mesh by index on each direction
    /// </summary>
    /// <param name="index">index = [0]-[1] : u-direction, [2]-[3] : v-direction</param>
    /// <returns></returns>
    public ControlPoint[] SelectEndPatchControlPoints(int index)
    {
      ControlPoint[] selCps = null;
      if (this is AbstractPatch2D)
      {
        selCps = ((Abstract2DParametricGeometry)GetGeometry()).SelectControlPointInTwoEndLines(index);
      }
      else if (this is AbstractPatch3D)
      {
        ControlPoint[,] temp = ((Abstract3DParametricGeometry)GetGeometry()).SelectControlPointInTwoEndAreas(index);
        selCps = new ControlPoint[temp.GetLength(0) * temp.GetLength(1)];
        int count = 0;
        for (int i = 0; i < temp.GetLength(0); i++)
          for (int j = 0; j < temp.GetLength(1); j++)
            selCps[count++] = temp[i, j];
      }
      return selCps;
    }

    public ControlPoint[] SelectNearEndPatchControlPoints(int index)
    {
      ControlPoint[] selCps = null;
      if (this is AbstractPatch2D)
      {
        selCps = ((Abstract2DParametricGeometry)GetGeometry()).SelectControlPointNearInTwoEndLines(index);
      }
      else if (this is AbstractPatch3D)
      {
        ControlPoint[,] temp = ((Abstract3DParametricGeometry)GetGeometry()).SelectControlPointNearInTwoEndAreas(index);
        selCps = new ControlPoint[temp.GetLength(0) * temp.GetLength(1)];
        int count = 0;
        for (int i = 0; i < temp.GetLength(0); i++)
          for (int j = 0; j < temp.GetLength(1); j++)
            selCps[count++] = temp[i, j];
      }
      return selCps;
    }

    public ControlPoint[] SelectNearEndPatchControlPointsExceptOnBoundary(int index)
    {
      ControlPoint[] selCps = null;
      if (this is AbstractPatch2D)
      {
        ControlPoint[] selCpsTemp = ((Abstract2DParametricGeometry)GetGeometry()).SelectControlPointNearInTwoEndLines(index);
        selCps = new ControlPoint[selCpsTemp.Length - 2];
        int c = 0;
        for (int i = 1; i < selCpsTemp.Length - 1; i++)
        {
          selCps[c++] = selCpsTemp[i];
        }
      }
      else if (this is AbstractPatch3D)
      {
        ControlPoint[,] temp = ((Abstract3DParametricGeometry)GetGeometry()).SelectControlPointNearInTwoEndAreas(index);
        selCps = new ControlPoint[temp.GetLength(0) * temp.GetLength(1)];
        int count = 0;
        for (int i = 0; i < temp.GetLength(0); i++)
          for (int j = 0; j < temp.GetLength(1); j++)
            selCps[count++] = temp[i, j];
      }
      return selCps;
    }

    public override double ComputeSharpCrackSurface()
    {
      double As = 0;
      int nel = CountElements();//(n - p) * (m - q);//number of elements
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int i = 0; i < nel; i++)//Loop through elements
        {
          AbstractElement elem = listElement[i];
          As += elem.ComputeSharpCrackSurface();
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
                  lock (monitor)
                  {
                    As += elem.ComputeSharpCrackSurface();
                  }
                }
              });
        }

        Task.WaitAll(tasks);
      }
      return As;
    }
    public override void ComputeInternalForcePatch(ref DoubleVector residualGlobal)
    {
      int nel = CountElements();//(n - p) * (m - q);//number of elements
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int i = 0; i < nel; i++)//Loop through elements
        {
          AbstractElement elem = listElement[i];
          DoubleVector rel = null;

          elem.ComputeInternalForceElement(out rel);
          int[] tArray = elem.GetTArrayGlobal();
          MatrixTool.AssemblyVector(rel, residualGlobal, tArray);
        }
      }
      else
      {
        DoubleVector rGlobalPatch = residualGlobal;
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
                  DoubleVector rel = null;

                  elem.ComputeInternalForceElement(out rel);
                  int[] tArray = elem.GetTArrayGlobal();
                  lock (monitor)
                  {
                    MatrixTool.AssemblyVector(rel, rGlobalPatch, tArray);
                  }
                }
              });
        }

        Task.WaitAll(tasks);
        residualGlobal = rGlobalPatch;
      }


      //    //DoubleVector Repa = new DoubleVector(countDOF);
      //    /////////////////////////////////////////////////////////////////////////////////
      //    //////// Parallel ///////////////////////////////////////////////////////////////
      //    /////////////////////////////////////////////////////////////////////////////////
      //    //object monitor = new object();
      //    //Parallel.For<DoubleVector>(0, GetCountElement(), new ParallelOptions { MaxDegreeOfParallelism = numberOfCPU },
      //    //    () => new DoubleVector(countDOF),
      //    //    (i, state, subtotal) =>
      //    //    {
      //    //        AbstractStructure2DElement elem = (AbstractStructure2DElement)listElement[i];
      //    //        DoubleVector rel = null;

      //    //        elem.ComputeInternalForceElement(ref rel);
      //    //        int[] tArray = elem.GetTArray();
      //    //        for (int j = 0; j < tArray.Length; j++)
      //    //        {
      //    //            int enumerateM = tArray[j];
      //    //            if (enumerateM != -1)
      //    //                subtotal[enumerateM] += rel[j];
      //    //        }
      //    //        return subtotal;
      //    //    }, (finalResult) => { lock (monitor) Repa += finalResult; }
      //    //);
      //    //Re = Repa;
      //    int countDOFModel = residualGlobal.Length;
      //    object monitor = new object();
      //    DoubleVector rePatch = new DoubleVector(countDOFModel);
      //    Parallel.For<DoubleVector>(0, nel, () => new DoubleVector(countDOFModel),
      //                    (i, state, localResult) =>
      //                    {
      //                        AbstractElementStructure2D elem = (AbstractElementStructure2D)listElement[i];
      //                        DoubleVector rel = null;
      //                        elem.ComputeInternalForceElement(out rel);
      //                        int[] tArray = elem.GetTArrayGlobal();
      //                        MatrixTool.AssemblyVector(rel, localResult, tArray);
      //                        return localResult;
      //                    }, localFinal => { lock (monitor) { rePatch += localFinal; } }
      //        );
      //    residualGlobal = residualGlobal + rePatch;
      //}
    }

    public override void ComputeStiffnessMatrixPatch(ref DoubleMatrix kGlobal)
    {
      int nel = CountElements();//number of elements
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int i = 0; i < nel; i++)//Loop through elements
        {
          ComputeKeAndAssemplyKGlobal(kGlobal, i);
        }
      }
      else
      {
        DoubleMatrix kGlobalPatch = kGlobal;
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
                  DoubleMatrix Ke = ComputeKe(elem);
                  if (!AbstractModel.IsComputeKeParallel)
                  {
                    int[] tArray = elem.GetTArrayGlobal();
                    lock (monitor)
                    {
                      MatrixTool.AssemblyMatrix(Ke, kGlobalPatch, tArray);
                    }
                  }
                }
              });
        }

        Task.WaitAll(tasks);

        if (AbstractModel.IsComputeKeParallel)
        {
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
                    DoubleMatrix Ke = elem.GetStiffnessMatrix();
                    ///////////////////
                    int[] tArray = elem.GetTArrayGlobal();
                    lock (monitor)
                    {
                      MatrixTool.AssemblyMatrix(Ke, kGlobalPatch, tArray);
                    }
                  }
                });
          }
          Task.WaitAll(tasks);
        }
        kGlobal = kGlobalPatch;
      }
    }

    public override void ComputeStiffnessMatrixPatch(ref DoubleCsrSparseMatrix kGlobal)
    {
      int nel = CountElements();//number of elements
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int i = 0; i < nel; i++)//Loop through elements
        {
          ComputeKeAndAssemplyKGlobal(kGlobal, i);
        }
      }
      else
      {
        DoubleCsrSparseMatrix kGlobalPatch = kGlobal;
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
                  DoubleMatrix Ke = ComputeKe(elem);
                  if (!AbstractModel.IsComputeKeParallel)
                  {
                    int[] tArray = elem.GetTArrayGlobal();
                    lock (monitor)
                    {
                      MatrixTool.AssemblyMatrix(Ke, kGlobalPatch, tArray);
                    }
                  }
                }
              });
        }

        Task.WaitAll(tasks);

        if (AbstractModel.IsComputeKeParallel)
        {
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
                    DoubleMatrix Ke = elem.GetStiffnessMatrix();
                    ///////////////////
                    int[] tArray = elem.GetTArrayGlobal();
                    lock (monitor)
                    {
                      MatrixTool.AssemblyMatrix(Ke, kGlobalPatch, tArray);
                    }
                  }
                });
          }
          Task.WaitAll(tasks);
        }
        kGlobal = kGlobalPatch;
      }
    }

    private void ComputeKeAndAssemplyKGlobal(DoubleMatrix kGlobal, int i)
    {
      AbstractElement elem = listElement[i];
      DoubleMatrix Ke = ComputeKe(elem);
      ///////////////////
      int[] tArray = elem.GetTArrayGlobal();
      MatrixTool.AssemblyMatrix(Ke, kGlobal, tArray);
    }

    private void ComputeKGeAndAssemplyKGGlobal(double[,] stressTensor, DoubleMatrix kGGlobal, int i)
    {
      AbstractElement elem = listElement[i];
      DoubleMatrix KGe = ComputeKGe(elem, stressTensor);
      ///////////////////
      int[] tArray = elem.GetTArrayGlobal();
      MatrixTool.AssemblyMatrix(KGe, kGGlobal, tArray);
    }

    private void ComputeKGeAndAssemplyKGGlobal(double[,] stressTensor, DoubleCsrSparseMatrix kGGlobal, int i)
    {
      AbstractElement elem = listElement[i];
      DoubleMatrix KGe = ComputeKGe(elem, stressTensor);
      ///////////////////
      int[] tArray = elem.GetTArrayGlobal();
      MatrixTool.AssemblyMatrix(KGe, kGGlobal, tArray);
    }
    private void ComputeKeAndAssemplyKGlobal(DoubleCsrSparseMatrix kGlobal, int i)
    {
      AbstractElement elem = listElement[i];
      DoubleMatrix Ke = ComputeKe(elem);
      ///////////////////
      int[] tArray = elem.GetTArrayGlobal();
      MatrixTool.AssemblyMatrix(Ke, kGlobal, tArray);
    }

    private DoubleMatrix ComputeKe(AbstractElement elem)
    {
      DoubleMatrix Ke;
      if (AbstractModel.TypeAnalysisModel == TypeAnalysisModel.TopologyOptimization)
      {
        Ke = ((ElementStructureElastic2D)elem).GetStiffnessMatrix();
        double Ee = ((ElementStructureElastic2D)elem).GetCurrentModulus();
        Ke = Ee * Ke;
      }
      else
        elem.ComputeStiffnessMatrixElement(out Ke);
      if (elem is ElementThermal2D)
      {
        if (((ElementThermal2D)elem).GetHeatTransferConvection() != null)
        {
          ((ElementThermal2D)elem).ComputeHeatTransferConvectionStiffnessMatrixElement(ref Ke);
        }
      }
      //////////////// sua lai
      else if (elem is ElementThermal3D)
      {
        if (((ElementThermal3D)elem).GetHeatTransferConvection() != null)
        {
          ((ElementThermal3D)elem).ComputeHeatTransferConvectionStiffnessMatrixElement(ref Ke);
        }
      }

      if (elem is ElementThermoelastic2D)
      {
        if (((ElementThermoelastic2D)elem).GetHeatTransferConvection() != null)
        {
          ((ElementThermoelastic2D)elem).ComputeHeatTransferConvectionStiffnessMatrixElement(ref Ke);
        }
      }
      //////////////// sua lai
      else if (elem is ElementThermoelastic3D)
      {
        if (((ElementThermoelastic3D)elem).GetHeatTransferConvection() != null)
        {
          ((ElementThermoelastic3D)elem).ComputeHeatTransferConvectionStiffnessMatrixElement(ref Ke);
        }
      }
      if (AbstractModel.IsComputeKeParallel)
      {
        elem.SetStiffnessMatrix(Ke);
      }
      return Ke;
    }

    private DoubleMatrix ComputeKGe(AbstractElement elem, double[,] stressTensor)
    {
      elem.ComputeGeometricStiffnessMatrixElement(stressTensor, out DoubleMatrix KGe);
      elem.SetGeometricStiffnessMatrix(KGe);
      return KGe;
    }

    public override void ComputeMassMatrixPatch(ref DoubleCsrSparseMatrix mGlobal)
    {
      int nel = CountElements();//(n - p) * (m - q);//number of elements
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int i = 0; i < nel; i++)//Loop through elements
        {
          AbstractElement elem = listElement[i];
          DoubleMatrix Me = null;

          elem.ComputeMassMatrixElement(out Me);
          int[] tArray = elem.GetTArrayGlobal();
          MatrixTool.AssemblyMatrix(Me, mGlobal, tArray);
        }
      }
      else
      {
        DoubleCsrSparseMatrix mGlobalPatch = mGlobal;
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
                  DoubleMatrix Me = null;

                  elem.ComputeMassMatrixElement(out Me);
                  int[] tArray = elem.GetTArrayGlobal();
                  lock (monitor)
                  {
                    MatrixTool.AssemblyMatrix(Me, mGlobalPatch, tArray);
                  }
                }
              });
        }

        Task.WaitAll(tasks);
      }
    }

    public override void ComputeMassMatrixPatch(ref DoubleMatrix mGlobal)
    {
      int nel = CountElements();//(n - p) * (m - q);//number of elements
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int i = 0; i < nel; i++)//Loop through elements
        {
          AbstractElement elem = listElement[i];
          DoubleMatrix Me = null;

          elem.ComputeMassMatrixElement(out Me);
          int[] tArray = elem.GetTArrayGlobal();
          MatrixTool.AssemblyMatrix(Me, mGlobal, tArray);
        }
      }
      else
      {
        DoubleMatrix mGlobalPatch = mGlobal;
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
                  DoubleMatrix Me = null;

                  elem.ComputeMassMatrixElement(out Me);
                  int[] tArray = elem.GetTArrayGlobal();
                  lock (monitor)
                  {
                    MatrixTool.AssemblyMatrix(Me, mGlobalPatch, tArray);
                  }
                }
              });
        }

        Task.WaitAll(tasks);

        mGlobal = mGlobalPatch;
      }
    }

    public override void ComputeGeometricStiffnessMatrixPatch(double[,] stressTensor, ref DoubleMatrix kGGlobal)
    {
      int nel = CountElements();//number of elements
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int i = 0; i < nel; i++)//Loop through elements
        {
          ComputeKGeAndAssemplyKGGlobal(stressTensor, kGGlobal, i);
        }
      }
      else
      {
        DoubleMatrix kGGlobalPatch = kGGlobal;
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
                  DoubleMatrix KGe = ComputeKGe(elem, stressTensor);
                  if (!AbstractModel.IsComputeKeParallel)
                  {
                    int[] tArray = elem.GetTArrayGlobal();
                    lock (monitor)
                    {
                      MatrixTool.AssemblyMatrix(KGe, kGGlobalPatch, tArray);
                    }
                  }
                }
              });
        }

        Task.WaitAll(tasks);

        if (AbstractModel.IsComputeKeParallel)
        {
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
                    DoubleMatrix Ke = elem.GetGeometricStiffnessMatrix();
                    ///////////////////
                    int[] tArray = elem.GetTArrayGlobal();
                    lock (monitor)
                    {
                      MatrixTool.AssemblyMatrix(Ke, kGGlobalPatch, tArray);
                    }
                  }
                });
          }
          Task.WaitAll(tasks);
        }
        kGGlobal = kGGlobalPatch;
      }
    }

    public override void ComputeGeometricStiffnessMatrixPatch(double[,] stressTensor, ref DoubleCsrSparseMatrix kGGlobal)
    {
      int nel = CountElements();//number of elements
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int i = 0; i < nel; i++)//Loop through elements
        {
          ComputeKGeAndAssemplyKGGlobal(stressTensor, kGGlobal, i);
        }
      }
      else
      {
        DoubleCsrSparseMatrix kGGlobalPatch = kGGlobal;
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
                  DoubleMatrix Ke = ComputeKGe(elem, stressTensor);
                  if (!AbstractModel.IsComputeKeParallel)
                  {
                    int[] tArray = elem.GetTArrayGlobal();
                    lock (monitor)
                    {
                      MatrixTool.AssemblyMatrix(Ke, kGGlobalPatch, tArray);
                    }
                  }
                }
              });
        }

        Task.WaitAll(tasks);

        if (AbstractModel.IsComputeKeParallel)
        {
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
                    DoubleMatrix KGe = elem.GetGeometricStiffnessMatrix();
                    ///////////////////
                    int[] tArray = elem.GetTArrayGlobal();
                    lock (monitor)
                    {
                      MatrixTool.AssemblyMatrix(KGe, kGGlobalPatch, tArray);
                    }
                  }
                });
          }
          Task.WaitAll(tasks);
        }
        kGGlobal = kGGlobalPatch;
      }
    }
  }
}
