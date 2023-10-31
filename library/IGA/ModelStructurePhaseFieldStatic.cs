using System;
using System.Collections.Generic;
using DEMSoft.NURBS;
using CenterSpace.NMath.Core;
using System.Threading.Tasks;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public enum TypeDegradationEnergy
  {
    /// <summary>
    /// Molnar (2017) 2D and 3D Abaqus implementation of a robust staggered phase-field solution for modeling brittle fracture
    /// </summary>
    Isotropic,
    /// <summary>
    /// Miehe (2010) A phase field model for rate-independent crack propagation: robust algorithmic implementation based on operator splits.
    /// </summary>
    Anisotropic,
    /// <summary>
    /// Ambati(2014) A review on phase-field models of brittle fractur and a new fast hybrid formulation
    /// </summary>
    Hybrid
  }
  public enum TypeOrderPhaseField
  {
    SecondOrder,
    FourthOrder
  }
  public enum TypeModelPhasefield
  {
    AT1,
    /// <summary>
    /// A quadratic degradation function, Miehe model
    /// </summary>
    AT2,
    /// <summary>
    /// A cubic degradation function
    /// </summary>
    Cubic,
    /// <summary>
    /// Linear softening curve Wu(2019)
    /// </summary>
    PFCZM_Linear,
    /// <summary>
    /// Linear softening curve Wu(2019)
    /// </summary>
    PFCZM_BLinear,
    /// <summary>
    /// Exponential softening curve Wu(2019)
    /// </summary>
    PFCZM_Exponential,
    /// <summary>
    /// Hyperbolic softening curve Wu(2019)
    /// </summary>
    PFCZM_Hyperbolic,
    /// <summary>
    /// Cornelissen softening curve Wu(2019)
    /// </summary>
    PFCZM_Cornelissen
  }

  public enum TypeFailureCriterion
  {
    Rankine,
    ModifiedVonMises,
    Wu2017
  }
  /// <summary>
  /// Structure 2D problem (plane stress, plane strain). Default is plane stress
  /// </summary>
  public class ModelStructurePhaseFieldStatic : AbstractModelStructure
  {
    /// <summary>
    /// Length Scale Parameter, Characteristic Crack Width [mm]
    /// </summary>
    public static double l0 = 0.1;
    /// <summary>
    /// Griffith’s or Irwin’s characteristic length lch, Characteristic Crack Width [mm]
    /// </summary>
    //public static double lch = 15;//E0*gC/ft^2  [mm] brittle Wu2018 page 774
    public static TypeModelPhasefield typePhaseFieldModel = TypeModelPhasefield.AT2;
    ///// <summary>
    ///// Strain Energy Release Rate [N/mm]
    ///// </summary>
    //public static double Gc = 1;
    /// <summary>
    /// Positive integer is used in penalty, default is 3
    /// </summary>
    public static int n = 2;
    /// <summary>
    /// Parameter of well-Conditioned System
    /// </summary>
    public static double k = 1e-4;
    /// <summary>
    /// The parameter epsi controls the magnitude of the penalty term,
    /// and should be set to a value that is large enough to sufficiently 
    /// enforce the irreversibility condition, but not too large as to result
    /// in an ill-conditioned system.
    /// </summary>
    public static double epsi = 1e-6;// 100;//2;
    /// <summary>
    /// An isotropic or anisotropic constitutive assumption for the degradation of energy due to fracture Miehe(2010)
    /// </summary>
    public static TypeDegradationEnergy typeDegradationModel = TypeDegradationEnergy.Isotropic;
    public static TypeOrderPhaseField typeOrderPhaseField = TypeOrderPhaseField.SecondOrder;
    public static TypeFailureCriterion typeFailureCriterion = TypeFailureCriterion.Rankine;
    protected int countDOFUnConstraint;

    private bool bConverged = true;
    private List<CrackEdge> cracks;
    private List<Void2D> voids;

    /// <summary>
    /// Constructor class
    /// </summary>
    public ModelStructurePhaseFieldStatic(Dimension structureDimension, string pathProject = "temp", string nameProject = "project")
                 : base(TypeAnalysisModel.Static, structureDimension, pathProject, nameProject)
    {
      if (StructureDimension == Dimension.Plane)
      {
        listComputeResult.Add(Result.PHASEFIELD);
      }
      else if (StructureDimension == Dimension.Solid)
      {
        listComputeResult.Add(Result.PHASEFIELD);
      }
      TypeModel = TypeModelProblem.PhaseField;
      typeOfFieldsInMultifield.Add(TypeFields.PhaseField);
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
        tArray[1] = GetTArrayGlobalPhaseField();
        if (bcMatrix != null)
        {
          tArrayUnconstraint[0] = GetTArrayGlobalMechanicUnConstraint();
          tArrayConstraint[0] = GetTArrayGlobalMechanicConstraint();
          tArrayUnconstraint[1] = GetTArrayGlobalPhaseFieldUnConstraint();
          tArrayConstraint[1] = GetTArrayGlobalPhaseFieldConstraint();
        }
        else
        {
          tArrayUnconstraint[0] = tArray[0];
          tArrayUnconstraint[1] = tArray[1];
        }
      }
      else
      {
        if (bcMatrix != null)
        {
          Task<int[]>[] tasks ={
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalMechanic()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalPhaseField()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalMechanicUnConstraint()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalMechanicConstraint()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalPhaseFieldUnConstraint()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalPhaseFieldConstraint())};
          //Block until all tasks complete.
          Task.WaitAll(tasks);
          tArray[0] = tasks[0].Result;
          tArray[1] = tasks[1].Result;
          tArrayUnconstraint[0] = tasks[2].Result;
          tArrayConstraint[0] = tasks[3].Result;
          tArrayUnconstraint[1] = tasks[4].Result;
          tArrayConstraint[1] = tasks[5].Result;
        }
        else
        {
          Task<int[]>[] tasks ={
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalMechanic()),
            Task<int[]>.Factory.StartNew(() => GetTArrayGlobalPhaseField()) };
          //Block until all tasks complete.
          Task.WaitAll(tasks);
          tArray[0] = tasks[0].Result;
          tArray[1] = tasks[1].Result;
          tArrayUnconstraint[0] = tArray[0];
          tArrayUnconstraint[1] = tArray[1];
        }
      }
      if (cracks != null)
      {
        InitializePreCrack();
      }
      if (voids != null)
      {
        InitializeVoids();
      }
    }
    /// <summary>
    /// Set phase field parameters
    /// </summary>
    /// <param name="l0">Length Scale Parameter, Characteristic Crack Width [mm] (default is 0.1)</param>
    /// <param name="k">Parameter of well-Conditioned System (default is 0.0001)</param>
    /// <param name="eta">A viscosity parameter (default is 1e-6)</param>
    public void SetPhaseFieldParameters(double l0, /*double lch,*/ double k, double eta, TypeModelPhasefield typePhaseFieldModel = TypeModelPhasefield.AT2, TypeFailureCriterion typeFailureCriterion = TypeFailureCriterion.Rankine)
    {
      //ModelStructurePhaseFieldStatic.lch = lch;
      ModelStructurePhaseFieldStatic.typePhaseFieldModel = typePhaseFieldModel;
      ModelStructurePhaseFieldStatic.typeFailureCriterion = typeFailureCriterion;
      ModelStructurePhaseFieldStatic.l0 = l0;
      ModelStructurePhaseFieldStatic.k = k;
      ModelStructurePhaseFieldStatic.epsi = eta;
    }

    /// <summary>
    /// Solve problem
    /// </summary>
    public override void Solve()
    {

      if (listPatch.Count == 0)
        throw new NullReferenceException("Model must be assigned Patch");
      WriteInformationProblem("Static phase-field module");

      switch (type)
      {
        case TypeNonlinearSolver.LoadControlledNewtonRaphson:
        case TypeNonlinearSolver.LoadControlledModifiedNewtonRaphson:
          new SolverLoadControlled(this, type);
          break;
        case TypeNonlinearSolver.DisplacementControlledNewtonRaphson:
        case TypeNonlinearSolver.DisplacementControlledModifiedNewtonRaphson:
          new SolverDisplacementControlled(this, type);
          break;
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

    private int[] GetTArrayGlobalPhaseField()
    {
      List<int> listTArrayGlobalPhaseField = new List<int>();
      for (int i = 0; i < listPatch.Count; i++)
      {
        AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
        if (TypeModel == TypeModelProblem.PhaseField && pa.TypeOfFieldsInMultifield.IndexOf(TypeFields.PhaseField) != -1)
        {
          int nen = pa.GetCountGlobalBasisFunctions();
          for (int k = 0; k < nen; k++)
          {
            int temp = pa.GetIDInGlobal((int)StructureDimension, k);
            if (temp != -1)
            {
              if (listTArrayGlobalPhaseField.IndexOf(temp) == -1)
              {
                listTArrayGlobalPhaseField.Add(temp);
              }
            }
          }
        }
      }
      return listTArrayGlobalPhaseField.ToArray();
    }
    private int[] GetTArrayGlobalPhaseFieldConstraint()
    {
      List<int> listTArrayGlobalPhaseField = new List<int>();
      for (int i = 0; i < listPatch.Count; i++)
      {
        AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
        if (TypeModel == TypeModelProblem.PhaseField && pa.TypeOfFieldsInMultifield.IndexOf(TypeFields.PhaseField) != -1)
        {
          int nen = pa.GetCountGlobalBasisFunctions();
          for (int k = 0; k < nen; k++)
          {
            int temp = pa.GetIDInGlobal((int)StructureDimension, k);
            if (temp != -1)
            {
              int idxConstraint = MatrixTool.IndexOfCol(bcMatrix, 0, temp, AbstractModel.NumberOfCPUs);
              if (idxConstraint != -1)
                if ((listTArrayGlobalPhaseField.IndexOf(temp) == -1) && (bcMatrix[idxConstraint, 2] == (int)StructureDimension/*bac tu do volt*/))
                {
                  listTArrayGlobalPhaseField.Add(temp);
                }
            }
          }
        }
      }
      return listTArrayGlobalPhaseField.ToArray();
    }

    private int[] GetTArrayGlobalPhaseFieldUnConstraint()
    {
      List<int> listTArrayGlobalPhaseField = new List<int>();
      for (int i = 0; i < listPatch.Count; i++)
      {
        AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
        if (TypeModel == TypeModelProblem.PhaseField && pa.TypeOfFieldsInMultifield.IndexOf(TypeFields.PhaseField) != -1)
        {
          int nen = pa.GetCountGlobalBasisFunctions();
          for (int k = 0; k < nen; k++)
          {
            int temp = pa.GetIDInGlobal((int)StructureDimension, k);
            if (temp != -1)
            {
              int idxConstraint = MatrixTool.IndexOfCol(bcMatrix, 0, temp, AbstractModel.NumberOfCPUs);
              if ((listTArrayGlobalPhaseField.IndexOf(temp) == -1) && (idxConstraint == -1))
              {
                listTArrayGlobalPhaseField.Add(temp);
              }
            }
          }
        }
      }
      return listTArrayGlobalPhaseField.ToArray();
    }

    public void AddCrackPath(CrackEdge[] cracks)
    {
      if (this.cracks == null)
      {
        this.cracks = new List<CrackEdge>();
      }
      foreach (var c in cracks)
      {
        this.cracks.Add(c);
      }
    }

    public void AddVoid(Void2D[] voids)
    {
      if (this.voids == null)
      {
        this.voids = new List<Void2D>();
      }
      foreach (var c in voids)
      {
        this.voids.Add(c);
      }
    }

    private void InitializePreCrack()
    {
      //foreach (var c in cracks)
      //{
      //  for (int i = 0; i < c.IndexPatchHasCrossedCrack.Length; i++)
      //  {
      //    AbstractPatch patch = GetPatch(c.IndexPatchHasCrossedCrack[i]);
      //    int numElement = patch.CountElements();
      //    if (!AbstractModel.IsParallelProcesing)
      //    {
      //      for (int j = 0; j < numElement; j++)
      //      {
      //        ElementStructurePhaseField2D elem = (ElementStructurePhaseField2D)patch.GetElement(j);
      //        elem.InitialPreCrackGausspointValue(c);
      //      }
      //    }
      //    else
      //    {
      //      var degreeOfParallelism = AbstractModel.NumberOfCPUs;
      //      var tasks = new Task[degreeOfParallelism];
      //      for (int taskNumber = 0; taskNumber < degreeOfParallelism; taskNumber++)
      //      {
      //        // capturing taskNumber in lambda wouldn't work correctly
      //        int taskNumberCopy = taskNumber;

      //        tasks[taskNumber] = Task.Factory.StartNew(
      //             () =>
      //             {
      //               var max = numElement * (taskNumberCopy + 1) / degreeOfParallelism;
      //               for (int j = numElement * taskNumberCopy / degreeOfParallelism; j < max; j++)
      //               {
      //                 ElementStructurePhaseField2D elem = (ElementStructurePhaseField2D)patch.GetElement(j);
      //                 elem.InitialPreCrackGausspointValue(c);
      //               }
      //             });
      //      }
      //      Task.WaitAll(tasks);
      //    }
      //  }
      //}
      foreach (var c in cracks)
      {
        double[,] boundCrack = c.GetBoundBoxCoverCrackEdge(l0);
        for (int i = 0; i < c.IndexPatchHasCrossedCrack.Length; i++)
        {
          AbstractPatch patch = GetPatch(c.IndexPatchHasCrossedCrack[i]);
          var selElement = patch.SelectElement(boundCrack[0, 0], boundCrack[0, 1], boundCrack[1, 0], boundCrack[1, 1], 0, 0);
          int countElement = selElement.Count;
          if (!AbstractModel.IsParallelProcesing)
          {
            for (int j = 0; j < countElement; j++)
            {
              ElementStructurePhaseField2D elem = (ElementStructurePhaseField2D)selElement[j];
              elem.InitialPreCrackGausspointValue(c);
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
                     var max = countElement * (taskNumberCopy + 1) / degreeOfParallelism;
                     for (int j = countElement * taskNumberCopy / degreeOfParallelism; j < max; j++)
                     {
                       ElementStructurePhaseField2D elem = (ElementStructurePhaseField2D)selElement[j];
                       elem.InitialPreCrackGausspointValue(c);
                     }
                   });
            }
            Task.WaitAll(tasks);
          }
        }
      }
    }

    private void InitializeVoids()
    {
      foreach (var c in voids)
      {
        double[,] boundCrack = c.GetRegionBoxCoverVoid(2 * l0);
        for (int i = 0; i < c.IndexPatchHasCrossedCrack.Length; i++)
        {
          AbstractPatch patch = GetPatch(c.IndexPatchHasCrossedCrack[i]);
          var selElement = patch.SelectElement(boundCrack[0, 0], boundCrack[0, 1], boundCrack[1, 0], boundCrack[1, 1], 0, 0);
          int numElement = selElement.Count;
          if (!AbstractModel.IsParallelProcesing)
          {
            for (int j = 0; j < numElement; j++)
            {
              ElementStructurePhaseField2D elem = (ElementStructurePhaseField2D)selElement[j];
              elem.InitialVoidGausspointValue(c);
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
                       ElementStructurePhaseField2D elem = (ElementStructurePhaseField2D)selElement[j];
                       elem.InitialVoidGausspointValue(c);
                     }
                   });
            }
            Task.WaitAll(tasks);
          }
        }
      }
    }

    public void ExtrapolationInitialPhasefieldGausspoint(ref DoubleVector uGlobal)
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
            elem.ComputeDrawValueAtGaussPoint(DataInGausspoint.currentPhase);

            ////change order///////
            for (int j = 0; j < elem.GetNumberOfGaussPointOnEachDirection(); j++)
            {
              for (int i = 0; i < elem.GetNumberOfGaussPointOnEachDirection(); i++)
              {
                GaussPoints gpsij = elem.GetGaussPoint(i, j);
                DoubleVector Nii = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);//new DoubleVector(nen);

                Ae += MatrixFunctions.OuterProduct(Nii, Nii);
                Be += (double)gpsij.GetValue(DataInGausspoint.currentPhase) * Nii;
              }
            }
          }
          else if (listPatch[iji] is AbstractPatch3D)
          {
            AbstractElement3D elem = (AbstractElement3D)listPatch[iji].GetElement(iel);
            elem.ComputeDrawValueAtGaussPoint(DataInGausspoint.currentPhase);

            for (int i = 0; i < elem.GetNumberOfGaussPointOnEachDirection(); i++)
            {
              for (int j = 0; j < elem.GetNumberOfGaussPointOnEachDirection(); j++)
              {
                for (int k = 0; k < elem.GetNumberOfGaussPointOnEachDirection(); k++)
                {
                  GaussPoints gpsijk = elem.GetGaussPoint(i, j, k);
                  DoubleVector Nii = (DoubleVector)gpsijk.GetValue(DataInGausspoint.Ni);// new DoubleVector(nen);

                  //for (int kk = 0; kk < nen; kk++)
                  //  Nii[kk] = gpsijk.Ni[kk];
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
          for (int ii = 0; ii < nnp1; ii++)
          {
            int ik = patch.GetINC(0, ii, 0);
            int jk = patch.GetINC(0, ii, 1);
            int[] tArrayGlobal = ((NURBSSurface)listPatch[iji].GetGeometry(0)).ControlPoints[ik, jk].GetTArrayGlobal();
            if (tArrayGlobal.Length > 2)
              uGlobal[tArrayGlobal[2]] = Ps[ii];
          }
        }
        else if (patch is AbstractPatch3D)
        {
          for (int ii = 0; ii < nnp1; ii++)
          {
            int ik = listPatch[iji].GetINC(0, ii, 0);
            int jk = listPatch[iji].GetINC(0, ii, 1);
            int kk = listPatch[iji].GetINC(0, ii, 2);
            int[] tArrayGlobal = ((NURBSVolume)listPatch[iji].GetGeometry(0)).ControlPoints[ik, jk, kk].GetTArrayGlobal();
            if (tArrayGlobal.Length > 3)
              uGlobal[tArrayGlobal[3]] = Ps[ii];
          }
        }
      }
    }
    public double ComputeSharpCrackSurface()
    {
      double As = 0;
      for (int i = 0; i < listPatch.Count; i++)
      {
        AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
        As += pa.ComputeSharpCrackSurface();
      }
      return As;
    }
    //protected void SetUGlobalUnconstraint(DoubleVector uGlobal)
    //{
    //    for (int i = 0; i < listPatch.Count; i++)
    //    {
    //        AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
    //        ControlPoint[] cps = pa.GetAllControlPoints();
    //        //if (!IsParallelProcesing)
    //        //{
    //        for (int j = 0; j < cps.Length; j++)
    //        {
    //            ControlPoint currentCps = cps[j];
    //            for (int k = 0; k < (int)StructureDimension; k++)
    //            {
    //                int tArray = currentCps.GetTArrayGlobal()[k];
    //                if (tArray != -1)
    //                {
    //                    int tArrayUnconstraint = Array.IndexOf(tArrayMechanicUnconstraint, tArray);
    //                    if (tArrayUnconstraint != -1)
    //                    {
    //                        currentCps.SetULocal(k, uGlobal[tArrayUnconstraint]);
    //                        currentCps.AddResult(listComputeResult[k], uGlobal[tArrayUnconstraint]);
    //                    }
    //                }
    //            }
    //        }
    //        //}
    //        //else
    //        //{
    //        //    Parallel.For(0, cps.Length - 1, j =>
    //        //    {
    //        //        ControlPoint currentCps = cps[j];
    //        //        for (int k = 0; k < pa.GetCountField(); k++)
    //        //        {
    //        //            int tArray = currentCps.GetTArrayGlobal()[k];
    //        //            if (tArray != -1)
    //        //                currentCps.SetULocal(k, uGlobal[tArray]);
    //        //        }
    //        //    });
    //        //}
    //    }
    //}

    //protected void SetPGlobalUnconstraint(DoubleVector uGlobal)
    //{
    //    for (int i = 0; i < listPatch.Count; i++)
    //    {
    //        AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
    //        ControlPoint[] cps = pa.GetAllControlPoints();
    //        //if (!IsParallelProcesing)
    //        //{
    //        for (int j = 0; j < cps.Length; j++)
    //        {
    //            ControlPoint currentCps = cps[j];
    //            int tArray = currentCps.GetTArrayGlobal()[(int)StructureDimension];
    //            if (tArray != -1)
    //            {
    //                int tArrayUnconstraint = Array.IndexOf(tArrayPhaseFieldUnconstraint, tArray);
    //                if (tArrayUnconstraint != -1)
    //                {
    //                    currentCps.SetULocal((int)StructureDimension, uGlobal[tArrayUnconstraint]);
    //                    currentCps.AddResult(listComputeResult[(int)StructureDimension], uGlobal[tArrayUnconstraint]);
    //                }
    //            }
    //        }
    //        //}
    //        //else
    //        //{
    //        //    Parallel.For(0, cps.Length - 1, j =>
    //        //    {
    //        //        ControlPoint currentCps = cps[j];
    //        //        for (int k = 0; k < pa.GetCountField(); k++)
    //        //        {
    //        //            int tArray = currentCps.GetTArrayGlobal()[k];
    //        //            if (tArray != -1)
    //        //                currentCps.SetULocal(k, uGlobal[tArray]);
    //        //        }
    //        //    });
    //        //}
    //    }
    //}
    //protected DoubleVector GetConstraint(int[] tArrayConstraint)
    //{
    //    if (tArrayConstraint.Length != 0)
    //    {
    //        DoubleVector u0 = new DoubleVector(tArrayConstraint.Length);
    //        for (int i = 0; i < u0.Length; i++)
    //        {
    //            int idx = MatrixTool.IndexOfCol(bcMatrix, 0, tArrayConstraint[i]);
    //            u0[i] = bcMatrix[idx, 1];
    //        }
    //        return u0;
    //    }
    //    else
    //        return null;
    //}
  }
}
