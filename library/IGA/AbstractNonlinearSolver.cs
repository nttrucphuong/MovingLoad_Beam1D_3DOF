using System;
using System.Linq;
using CenterSpace.NMath.Core;
using System.Collections.Generic;
using DEMSoft.NURBS;
using DEMSoft.Function;
using System.IO;
using System.Threading.Tasks;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Type of nonlinear solver
  /// </summary>
  /// <see cref="Nam-Ho Kim - Introduction to Nonlinear"/>
  public enum TypeNonlinearSolver
  {
    LoadControlledNewtonRaphson,
    LoadControlledModifiedNewtonRaphson,
    DisplacementControlledNewtonRaphson,
    DisplacementControlledModifiedNewtonRaphson
  }

  public enum TypeSchemeNonlinearSolver
  {
    /// <summary>
    /// Solve all variables in one step in couple-fields
    /// </summary>
    Monolithic,
    /// <summary>
    /// Solve sequence variables in one step in couple-fields.Solve each independence fields in one step for next iterator step.
    /// </summary>
    Staggered,
    /// <summary>
    /// Solve single all sequence variables in one step in couple-fields. Solve multiple independence fields in one step.
    /// </summary>
    SingleStaggered
  }
  /// <summary>
  /// Implement the load controlled Newton-Raphson solver for non-linear problems
  /// </summary>
  internal abstract class AbstractNonlinearSolver
  {
    protected bool bConverged = true;
    protected bool bUpdate = true;
    protected TypeNonlinearSolver type; // 0 for Newton-Raphson, 1 for Modified Newton-Raphson
    protected AbstractModelStructure model;
    protected CheckpointDataProcess checkpointData;
    internal static bool IsSolvePhasefieldStaggered;
    protected bool conditionRestoreCheckpoint;
    internal static double deltaT;
    protected int[][] tArrayConstraintNoneZero;
    internal double prevF, lastPrevF0, lastPrevF1;
    
    /// <summary>
    /// Constructor class
    /// </summary>
    /// <param name="model"></param>
    public AbstractNonlinearSolver(AbstractModelStructure model, TypeNonlinearSolver type)
    {
      this.type = type;
      this.model = model;
      checkpointData = new CheckpointDataProcess(model.PathProject);
      checkpointData.model = model;
      conditionRestoreCheckpoint = !model.IsRunFromInitial && File.Exists(model.PathProject + "restore_parameter.msgpack");
      if (conditionRestoreCheckpoint)
      {
        checkpointData.RestoreLastCheckpointData(model.PathProject);
      }
      switch (type)
      {
        case TypeNonlinearSolver.LoadControlledNewtonRaphson:
        case TypeNonlinearSolver.DisplacementControlledNewtonRaphson:
          // Newton-Raphson: update tangent 						// stiffness at every iteration
          bUpdate = true;
          break;
        case TypeNonlinearSolver.LoadControlledModifiedNewtonRaphson:
        case TypeNonlinearSolver.DisplacementControlledModifiedNewtonRaphson:
          // Modified Newton-Raphson: only 						// update tangent stiffness at the 						// first iteration (iter = 0)
          bUpdate = false;
          break;
      }
      Solve(type);
    }

    /// <summary>
    /// Solve
    /// </summary>
    protected abstract void Solve(TypeNonlinearSolver type);
    protected void ComputeMonitorResults()
    {
      MonitorDataHistorism monitorDataHistorism = model.GetMonitorDataHistorism();
      int countData = monitorDataHistorism.CountData();
      /////////////////////////////////////////////////////////////////////////////
      if (countData > 0)
      {
        List<double[]> data = monitorDataHistorism.MonitorDataStore;
        data.Add(new double[countData]);
        for (int j = 0; j < countData; j++)
        {
          MonitorDataNode nodeData = monitorDataHistorism.GetMonitorDataNode(j);
          double[] xi = nodeData.Xi;
          Result result = nodeData.ResultType;
          AbstractPatch patch = nodeData.Patch;
          switch (result)
          {
            case Result.UX:
            case Result.UY:
            case Result.UZ:
            case Result.PHASEFIELD:
            case Result.VOLT:
            case Result.TEMP:
              if (patch is PatchStructure2D)
              {
                data[data.Count - 1][j] = ((PatchStructure2D)patch).GetApproximateAt(result, xi[0], xi[1]);
              }
              else if (patch is PatchStructure3D)
              {
                data[data.Count - 1][j] = ((PatchStructure3D)patch).GetApproximateAt(result, xi[0], xi[1], xi[2]);
              }
              break;
            case Result.USUM:
              if (patch is PatchStructure2D)
              {
                double[] u = new double[2];
                u[0] = ((PatchStructure2D)patch).GetApproximateAt(Result.UX, xi[0], xi[1]);
                u[1] = ((PatchStructure2D)patch).GetApproximateAt(Result.UY, xi[0], xi[1]);
                data[data.Count - 1][j] = Math.Sqrt(u[0] * u[0] + u[1] + u[1]);
              }
              else if (patch is PatchStructure3D)
              {
                double[] u = new double[3];
                u[0] = ((PatchStructure3D)patch).GetApproximateAt(Result.UX, xi[0], xi[1], xi[2]);
                u[1] = ((PatchStructure3D)patch).GetApproximateAt(Result.UY, xi[0], xi[1], xi[2]);
                u[2] = ((PatchStructure3D)patch).GetApproximateAt(Result.UZ, xi[0], xi[1], xi[2]);
                data[data.Count - 1][j] = Math.Sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
              }
              break;
            case Result.SIGMAXX:
            case Result.SIGMAYY:
            case Result.SIGMAXY:
            case Result.SIGMAZZ:
            case Result.SIGMAYZ:
            case Result.SIGMAXZ:
            case Result.SIGMAEQV:
              if (result != Result.SIGMAEQV)
              {
                model.Extrapolation(result);
              }
              else
              {
                model.Extrapolation(Result.SIGMAXX, Result.SIGMAYY, Result.SIGMAZZ, Result.SIGMAXY, Result.SIGMAYZ, Result.SIGMAXZ);
              }
              if ((patch is PatchStructure3D) || ((patch is PatchStructure2D) && (result != Result.SIGMAYZ && result != Result.SIGMAXZ)))
              {
                data[data.Count - 1][j] = patch.GetApproximateAt(result, xi);
              }
              break;
            case Result.EPSILONXX:
            case Result.EPSILONYY:
            case Result.EPSILONXY:
            case Result.EPSILONZZ:
            case Result.EPSILONYZ:
            case Result.EPSILONXZ:
            case Result.EPSILONEQV:
              if (result != Result.EPSILONEQV)
              {
                model.Extrapolation(result);
              }
              else
              {
                model.Extrapolation(Result.EPSILONXX, Result.EPSILONYY, Result.EPSILONZZ, Result.EPSILONXY, Result.EPSILONYZ, Result.EPSILONXZ);
              }
              if ((patch is PatchStructure3D) || ((patch is PatchStructure2D) && (result != Result.EPSILONYZ && result != Result.EPSILONXZ)))
              {
                data[data.Count - 1][j] = ((PatchStructure2D)patch).GetApproximateAt(result, xi);
              }
              break;
            //case Result.EPSILONPXX:
            //case Result.EPSILONPYY:
            //case Result.EPSILONPXY:
            //case Result.EPSILONPZZ:
            //case Result.EPSILONPYZ:
            //case Result.EPSILONPXZ:
            //	ExtrapolationPlasticStrain();
            //	if (patch is Structure2DPatch)
            //	{
            //		if (result != Result.EPSILONPYZ && result != Result.EPSILONPXZ)
            //			data[totalStep - 1, j] = ((Structure2DPatch)patch).GetStrainAt(xi[0], xi[1], result);
            //	}
            //	else if (patch is Structure3DPatch)
            //	{
            //		data[totalStep - 1, j] = ((Structure3DPatch)patch).GetStrainAt(xi[0], xi[1], xi[2], result);
            //	}
            //	break;
            //case Result.EPSILONPEQV:
            //	ExtrapolationEquivalentPlasticStrain();
            //	if (patch is Structure2DPatch)
            //	{
            //		data[totalStep - 1, j] = ((Structure2DPatch)patch).GetStrainAt(xi[0], xi[1], result);
            //	}
            //	else if (patch is Structure3DPatch)
            //	{
            //		data[totalStep - 1, j] = ((Structure3DPatch)patch).GetStrainAt(xi[0], xi[1], xi[2], result);
            //	}
            //	break;
            //case Result.REACFORCESUM:
            //  //data[data.Count - 1][j] = FR;

            //  AbstractConstraintValue cc = model.listConstraint[nodeData.ConstraintID];
            //  int[] tArrayConstraint = cc.GetTArrayConstrainted();
            //  int[] tArrayConstraintZeroX = model.GetTArrayFromConstraintDof(0).Except(tArrayConstraintNoneZero[0]).ToArray();
            //  int[] tArrayConstraintZeroY = model.GetTArrayFromConstraintDof(1).Except(tArrayConstraintNoneZero[0]).ToArray();
            //  int[] tArrayConstraintZeroZ = model.GetTArrayFromConstraintDof(2).Except(tArrayConstraintNoneZero[0]).ToArray();
            //  //int[] tArrayConstraintZeroX = model.GetTArrayFromConstraintDof(0).Except(tArrayConstraintNoneZero[0]).ToArray();
            //  //int[] tArrayConstraintZeroY = model.GetTArrayFromConstraintDof(1).Except(tArrayConstraintNoneZero[0]).ToArray();
            //  //int[] tArrayConstraintZeroZ = model.GetTArrayFromConstraintDof(2).Except(tArrayConstraintNoneZero[0]).ToArray();
            //  DoubleVector vFRX = MatrixTool.GetSubVector(model.residualForce, tArrayConstraintZeroX);
            //  DoubleVector vFRY = MatrixTool.GetSubVector(model.residualForce, tArrayConstraintZeroY);
            //  DoubleVector vFRZ = MatrixTool.GetSubVector(model.residualForce, tArrayConstraintZeroZ);
            //  double valReactionForceX = vFRX.Sum();
            //  double valReactionForceY = vFRY.Sum();
            //  double valReactionForceZ = vFRZ.Sum();
            //  valReactionForce = Math.Sqrt(valReactionForceX * valReactionForceX + valReactionForceY * valReactionForceY + valReactionForceZ * valReactionForceZ);
            //  data[data.Count - 1][j] = Math.Sqrt(valReactionForceX * valReactionForceX + valReactionForceY * valReactionForceY + valReactionForceZ * valReactionForceZ);
            //  break;
            case Result.REACFORCEX:
            case Result.REACFORCEY:
            case Result.REACFORCEZ:
              AbstractConstraintValue c = model.listConstraint[nodeData.ConstraintID];
              int[] tArrayConstraintZero = c.GetTArrayConstrainted();
              //int indexIDOfField = -1;
              //if (result == Result.REACFORCEX)
              //{
              //  indexIDOfField = model.listComputeResult.IndexOf(Result.UX);
              //}
              //else if (result == Result.REACFORCEY)
              //  indexIDOfField = model.listComputeResult.IndexOf(Result.UY);
              //else
              //  indexIDOfField = model.listComputeResult.IndexOf(Result.UZ);
              //int[] tArrayConstraintZero = model.GetTArrayFromConstraintDof(indexIDOfField).Except(tArrayConstraintNoneZero[0]).ToArray();
              DoubleVector vFR = MatrixTool.GetSubVector(model.residualForce, tArrayConstraintZero, AbstractModel.NumberOfCPUs);
              data[data.Count - 1][j] = vFR.Sum();
              break;
          }
        }
        monitorDataHistorism.UpdateUserDefinedDataMonitorStore();
      }
    }
    protected void UpdateConveredGaussPointsEachField()
    {
      if (AbstractModel.TypeModel == TypeModelProblem.PhaseField)
        for (int i = 0; i < model.CountPatch(); i++)
        {
          AbstractPatch patch = model.GetPatch(i);
          if (patch.TypeOfFieldsInMultifield.IndexOf(TypeFields.PhaseField) != -1)
          {
            int numElement = patch.CountElements();
            if (!AbstractModel.IsParallelProcesing)
            {
              for (int j = 0; j < numElement; j++)
              {
                AbstractElement elem = patch.GetElement(j);
                elem.UpdateGausspointValue();
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
                       var max = numElement * (taskNumberCopy + 1) / degreeOfParallelism;
                       for (int j = numElement * taskNumberCopy / degreeOfParallelism; j < max; j++)
                       {
                         AbstractElement elem = (AbstractElement)patch.GetElement(j);
                         elem.UpdateGausspointValue();
                       }
                     });
              }
              Task.WaitAll(tasks);
            }
          }
        }
    }
    protected void UpdateConvergedGaussPoints()
    {
      //////////////////////////////
      ////// update gausspoint variables ////
      ///////////////////////////////
      if (model.GetPatch(0).Material.TypeMaterialStructure == EngineeringData.TypeMaterialStructure.Plasticity || AbstractModel.TypeModel == TypeModelProblem.PhaseField)
        switch (model.StructureDimension)
        {
          case Dimension.Plane:
            for (int i = 0; i < model.CountPatch(); i++)
            {
              AbstractPatch patch = model.GetPatch(i);
              int numElement = patch.CountElements();
              if (!AbstractModel.IsParallelProcesing)
              {
                for (int j = 0; j < numElement; j++)
                {
                  AbstractElement2D elem = (AbstractElement2D)patch.GetElement(j);
                  int numberGps = elem.GetNumberOfGaussPointOnEachDirection();
                  for (int jj = 0; jj < numberGps; jj++)
                    for (int ii = 0; ii < numberGps; ii++)
                    {
                      elem.GetGaussPoint(ii, jj).UpdateConvergedVariablesLoadStep();
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
                         var max = numElement * (taskNumberCopy + 1) / degreeOfParallelism;
                         for (int j = numElement * taskNumberCopy / degreeOfParallelism; j < max; j++)
                         {
                           AbstractElement2D elem = (AbstractElement2D)patch.GetElement(j);
                           int numberGps = elem.GetNumberOfGaussPointOnEachDirection();
                           for (int jj = 0; jj < numberGps; jj++)
                             for (int ii = 0; ii < numberGps; ii++)
                             {
                               elem.GetGaussPoint(ii, jj).UpdateConvergedVariablesLoadStep();
                             }
                         }
                       });
                }
                Task.WaitAll(tasks);
              }
            }
            break;
          case Dimension.Solid:
            for (int i = 0; i < model.CountPatch(); i++)
            {
              AbstractPatch patch = model.GetPatch(i);
              int numElement = patch.CountElements();
              if (!AbstractModel.IsParallelProcesing)
              {
                for (int j = 0; j < numElement; j++)
                {
                  AbstractElement3D elem = (AbstractElement3D)patch.GetElement(j);
                  int numberGps = elem.GetNumberOfGaussPointOnEachDirection();
                  for (int kk = 0; kk < numberGps; kk++)
                    for (int jj = 0; jj < numberGps; jj++)
                      for (int ii = 0; ii < numberGps; ii++)
                      {
                        elem.GetGaussPoint(ii, jj, kk).UpdateConvergedVariablesLoadStep();
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
                         var max = numElement * (taskNumberCopy + 1) / degreeOfParallelism;
                         for (int j = numElement * taskNumberCopy / degreeOfParallelism; j < max; j++)
                         {
                           AbstractElement3D elem = (AbstractElement3D)patch.GetElement(j);
                           int numberGps = elem.GetNumberOfGaussPointOnEachDirection();
                           for (int kk = 0; kk < numberGps; kk++)
                             for (int jj = 0; jj < numberGps; jj++)
                               for (int ii = 0; ii < numberGps; ii++)
                               {
                                 elem.GetGaussPoint(ii, jj, kk).UpdateConvergedVariablesLoadStep();
                               }
                         }
                       });
                }
                Task.WaitAll(tasks);
              }
            }
            break;
        }
    }

    protected double[] ErrorCalculationIteration(DoubleVector duGlobal, DoubleVector residualForce, DoubleVector internalForceGlobali, DoubleVector internalForceGlobali1, int[] tArray)
    {
      double[] err1 = new double[model.tArray.Length];
      DoubleVector reduceResidualForce = null;
      DoubleVector reduceDu = null;
      int[] noConstraintTArray = null;
      switch (model.convergenceOption.criterionConvergence)
      {

        case CriterionConvergence.DisplacementCriterion:
        case CriterionConvergence.DisplacementCriterion22:
        case CriterionConvergence.DisplacementCriterion3:
        case CriterionConvergence.DisplacementCriterion5:
          reduceDu = MatrixTool.GetSubVector(duGlobal, tArray, AbstractModel.NumberOfCPUs);
          err1[0] = Math.Sqrt(MatrixFunctions.Dot(reduceDu, reduceDu));
          break;
        case CriterionConvergence.DisplacementCriterion4:
          DoubleVector reducedis0 = null;
          DoubleVector reducedis1 = null;
          
          switch (AbstractModel.TypeSchemeNonlinearSolver)
          {
            case TypeSchemeNonlinearSolver.Monolithic:
            case TypeSchemeNonlinearSolver.SingleStaggered:
              int[] noConstraintTArray0 = model.tArray[0].Except(model.tArrayConstraint[0]).ToArray();
              int[] noConstraintTArray1 = model.tArray[1].Except(model.tArrayConstraint[1]).ToArray();
              reducedis0 = MatrixTool.GetSubVector(duGlobal, noConstraintTArray0, AbstractModel.NumberOfCPUs);
              reducedis1 = MatrixTool.GetSubVector(duGlobal, noConstraintTArray1, AbstractModel.NumberOfCPUs);
              err1[0] = Math.Sqrt(MatrixFunctions.Dot(reducedis0, reducedis0));
              err1[1] = Math.Sqrt(MatrixFunctions.Dot(reducedis1, reducedis1));
              break;
            case TypeSchemeNonlinearSolver.Staggered:
              if (!IsSolvePhasefieldStaggered)
              {
                noConstraintTArray = tArray.Except(model.tArrayConstraint[0]).ToArray();
              }
              else
              {
                noConstraintTArray = tArray.Except(model.tArrayConstraint[1]).ToArray();
              }
              reduceDu = MatrixTool.GetSubVector(duGlobal, noConstraintTArray, AbstractModel.NumberOfCPUs);
              err1[0] = Math.Sqrt(MatrixFunctions.Dot(reduceDu, reduceDu));
              break;
          }
          break;
        case CriterionConvergence.ForceCriterion:
        case CriterionConvergence.ResidualCriterion3:
          //DoubleVector reduceResidualForce1 = null;
          //reduceResidualForce1 = MatrixTool.GetSubVector(residualForce, tArray);
          //double aa = reduceResidualForce1.Sum();
          DoubleVector reduceResidualForce0 = null;
          DoubleVector reduceResidualForce1 = null;
          //int[] noConstraintTArray = null;
          DoubleVector dR = internalForceGlobali1 - internalForceGlobali;
          switch (AbstractModel.TypeSchemeNonlinearSolver)
          {
            case TypeSchemeNonlinearSolver.Monolithic:
            case TypeSchemeNonlinearSolver.SingleStaggered:
              //int[] noConstraintTArray0 = model.tArray[0].Except(model.tArrayConstraint[0]).ToArray();
              //int[] noConstraintTArray1 = model.tArray[1].Except(model.tArrayConstraint[1]).ToArray();
              int[] noConstraintTArray0 = model.tArray[0];
              int[] noConstraintTArray1 = model.tArray[1];
              reducedis0 = MatrixTool.GetSubVector(dR, noConstraintTArray0, AbstractModel.NumberOfCPUs);
              reducedis1 = MatrixTool.GetSubVector(dR, noConstraintTArray1, AbstractModel.NumberOfCPUs);
              err1[0] = Math.Sqrt(MatrixFunctions.Dot(reducedis0, reducedis0));
              err1[1] = Math.Sqrt(MatrixFunctions.Dot(reducedis1, reducedis1));
              break;
            case TypeSchemeNonlinearSolver.Staggered:
              noConstraintTArray = tArray;
              reduceResidualForce = MatrixTool.GetSubVector(dR, noConstraintTArray, AbstractModel.NumberOfCPUs);
              //if (!IsSolvePhasefieldStaggered)
              //{
              //}
              //else
              //{
              //  noConstraintTArray = tArray.Except(model.tArrayConstraint[1]).ToArray();
              //  reduceResidualForce = MatrixTool.GetSubVector(residualForce, noConstraintTArray);
              //}
              err1[0] = Math.Sqrt(MatrixFunctions.Dot(reduceResidualForce, reduceResidualForce));
              break;
          }
          break;
        case CriterionConvergence.EnergyCriterion:
          reduceDu = MatrixTool.GetSubVector(duGlobal, tArray, AbstractModel.NumberOfCPUs);
          err1[0] = Math.Abs(MatrixFunctions.Dot(reduceDu, internalForceGlobali1 + internalForceGlobali));
          break;
        case CriterionConvergence.ResidualCriterion:
        case CriterionConvergence.ResidualCriterion2:
          reduceResidualForce = MatrixTool.GetSubVector(residualForce, tArray, AbstractModel.NumberOfCPUs);
          var curF = Math.Sqrt(MatrixFunctions.Dot(reduceResidualForce, reduceResidualForce));
          err1[0] = Math.Abs(curF - prevF);
          prevF = curF;
          break;
        case CriterionConvergence.DisplacementCriterion2:
          reduceDu = MatrixTool.GetSubVector(duGlobal, tArray, AbstractModel.NumberOfCPUs);
          err1[0] = MatrixFunctions.Dot(reduceDu, reduceDu);
          break;
      }
      return err1;
    }
    /// <summary>
    /// Solve for one iteration.
    /// </summary>
    /// <param name="Fstep"></param>
    /// <param name="ui"></param>
    /// <param name="err1"></param>
    /// <param name="bUpdate"></param>
    /// <returns></returns>
    protected abstract DoubleVector SolveOneIterationDense(DoubleMatrix kGlobal, DoubleVector Fstep, DoubleVector ui, ref double[] err1, ref DoubleVector residualForce);

    protected abstract DoubleVector SolveOneIterationSparse(DoubleCsrSparseMatrix kGlobalSparse, DoubleVector Fstep, DoubleVector ui, ref double[] err1, ref DoubleVector residualForce);
    protected void ComputeDisplacementVector(int stepLoad, int[] tArrayContraint, out DoubleVector DUNoneZeroVector, out int[] tArrayContraintNoneZero)
    {
      List<double> DUNone = new List<double>();
      List<int> tArrayNone = new List<int>();
      DUNoneZeroVector = null;
      tArrayContraintNoneZero = null;
      foreach (AbstractConstraintValue c in model.listConstraint)
      {
        double valueConstraint = 0;
        if ((c.GetValueConstraint() is double) || (c.GetValueConstraint() is ConstantFunctionRToR))
          valueConstraint = c.GetValueConstraint(stepLoad + 1) - c.GetValueConstraint(stepLoad);
        ControlPoint[] cpsel = c.GetControlPointsConstrainted();
        double[] valueFunctionConstraint = null;
        if (c is ConstraintValueEdge2D && (c.GetValueConstraint() is FunctionRToR) && (!(c.GetValueConstraint() is ConstantFunctionRToR)))
        {
          valueFunctionConstraint = (new DoubleVector(c.ComputeValueConstraintOnControlPoints(stepLoad + 1)) - new DoubleVector(c.ComputeValueConstraintOnControlPoints(stepLoad))).ToArray();
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
          int indexConstraint = Array.IndexOf(tArrayContraint, tArray);
          if (indexConstraint != -1 && valueConstraint != 0 && tArrayNone.IndexOf(tArray) == -1/*not exist in list*/)
          {
            DUNone.Add(valueConstraint);
            tArrayNone.Add(tArray);
          }
        }
      }
      if (DUNone.Count > 0)
      {
        DUNoneZeroVector = new DoubleVector(DUNone.ToArray());
        tArrayContraintNoneZero = tArrayNone.ToArray();
      }
    }

    
  }

  //public static class CrossThreadExtensions
  //{
  //  public static void PerformSafely(this Form target, System.Action action)
  //  {
  //    if (target.InvokeRequired)
  //    {
  //      target.Invoke(action);
  //    }
  //    else
  //    {
  //      action();
  //    }
  //  }

  //  public static void PerformSafely<T1>(this Control target, Action<T1> action, T1 parameter)
  //  {
  //    if (target.InvokeRequired)
  //    {
  //      target.Invoke(action, parameter);
  //    }
  //    else
  //    {
  //      action(parameter);
  //    }
  //  }

  //  public static void PerformSafely<T1, T2>(this Control target, Action<T1, T2> action, T1 p1, T2 p2)
  //  {
  //    if (target.InvokeRequired)
  //    {
  //      target.Invoke(action, p1, p2);
  //    }
  //    else
  //    {
  //      action(p1, p2);
  //    }
  //  }
  //}
}
