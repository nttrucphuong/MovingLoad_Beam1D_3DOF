using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using DEMSoft.NURBS;
using DEMSoft.Function;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public enum TypeStructure
  { Plane, Plate, Shell }
  /// <summary>
  /// Abstract class to define NURBS surface to NURBS Patch
  /// </summary>
  public abstract class AbstractPatch2D : AbstractPatchOneField
  {
    public double Thickness
    { get; set; }
    public TypeStructure TypeStructure { get; set; }
    public double[,][,] ExtractionOperator
    { get; set; }
    public override void ComputeExtractionOperator()
    {
      ExtractionOperator = GetSurface().ComputeExtractionOperator();
    }

    /// <summary>
    /// Constructor class
    /// </summary>
    /// <param name="surface">Nurbs surface</param>
    /// <param name="numberOfFields">number of fields in model</param>
    public AbstractPatch2D(Abstract2DParametricGeometry surface, TypeStructure typeStructure)
       : base(2)
    {
      this.TypeStructure = typeStructure;
      //connect = new List<ConnectionInterfaceTwoPatches>();
      geometry = new Abstract2DParametricGeometry[] { surface };

      //CreateIPN();
      //CreateINC();
      //CreateIEN();
      //Initialize();
    }

    /// <summary>
    /// Get NURBS surface
    /// </summary>
    /// <returns></returns>
    public NURBSSurface GetSurface()
    {
      return (NURBSSurface)GetGeometry();
    }

    /// <summary>
    /// Create INC array
    /// </summary>
    protected override void CreateINC()
    {
      var basis = GetSurface().Basis;
      int n = basis.GetCountBasisFunction(0);
      int m = basis.GetCountBasisFunction(1);
      int nen = GetCountLocalBasisFunctions();// (p + 1) * (q + 1); // number of local basis functions
      int nnp = GetCountGlobalBasisFunctions();// n * m;//number of global basis functions
      INC = new int[1][,];
      INC[0] = new int[nnp, 2];
      int A = 0;

      for (int j = 0; j < m; j++)
      {
        for (int i = 0; i < n; i++)
        {
          INC[0][A, 0] = i;
          INC[0][A, 1] = j;
          A++;
        }
      }
    }

    /// <summary>
    /// Create INC array
    /// </summary>
    protected override void CreateIEN()
    {
      var basis = GetSurface().Basis;
      int n = basis.GetCountBasisFunction(0);
      int p = basis.GetDegree(0);
      int m = basis.GetCountBasisFunction(1);
      int q = basis.GetDegree(1);
      var knotVectorNoMultiplicity1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var knotVectorNoMultiplicity2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      int numelem1 = knotVectorNoMultiplicity1.Length - 1;
      int numelem2 = knotVectorNoMultiplicity2.Length - 1;
      int nel = numelem1 * numelem2;// GetNumberOfPatchs();//(n - p) * (m - q);//number of elements
      int nen = (p + 1) * (q + 1); //GetNumberOfLocalBasisFunctions(); // number of local basis functions
      IEN = new int[1][,];
      IEN[0] = new int[nen, nel];
      //int e = 0;
      //int A = 0;
      for (int ej = 0; ej < numelem2; ej++)
      {
        for (int ei = 0; ei < numelem1; ei++)
        {
          double mid1 = (knotVectorNoMultiplicity1[ei] + knotVectorNoMultiplicity1[ei + 1]) / 2.0;
          double mid2 = (knotVectorNoMultiplicity2[ej] + knotVectorNoMultiplicity2[ej + 1]) / 2.0;
          int uspan = basis.FindSpan(mid1, 0);
          int vspan = basis.FindSpan(mid2, 1);
          int nume = FindIndexOfElement(ei, ej);
          int b = 0;
          for (int j = 0; j <= q; j++)
          {
            for (int i = 0; i <= p; i++)
            {
              int num = FindIndexOfGlobalBasisFunction(uspan - p + i, vspan - q + j);
              IEN[0][b, nume] = num;
              b++;
            }
          }
        }
      }
    }

    /// <summary>
    /// Create index of patch on direction
    /// </summary>
    protected override void CreateIPN()
    {
      var basis = GetSurface().Basis;
      int n = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity().Length - 1;
      int m = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity().Length - 1;
      int nnp = n * m;
      IPN = new int[nnp, 2];
      int A = 0;
      for (int j = 0; j < m; j++)
      {
        for (int i = 0; i < n; i++)
        {
          IPN[A, 0] = i;
          IPN[A, 1] = j;
          A++;
        }
      }
    }

    public override int FindIndexOfElementAt(params double[] xi)
    {
      int idxi = -1;
      int idxj = -1;
      var kv1 = GetSurface().Basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kv2 = GetSurface().Basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      for (int i = 0; i < kv1.Length; i++)
      {
        if (xi[0] >= kv1[i])
          idxi = i;
        if (xi[0] == kv1[kv1.Length - 1])
          idxi = i - 1;
      }
      for (int i = 0; i < kv2.Length; i++)
      {
        if (xi[1] >= kv2[i])
          idxj = i;
        if (xi[1] == kv2[kv2.Length - 1])
          idxj = i - 1;
      }
      return FindIndexOfElement(idxi, idxj);
    }

    /// <summary>
    /// Get number of patchs
    /// </summary>
    /// <returns></returns>
    public override int CalculateNumberOfElements()
    {
      var basis = GetSurface().Basis;
      int n = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity().Length - 1;
      int m = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity().Length - 1;
      return n * m;
    }

    //public override void SetUGlobal(DoubleVector uGlobal)
    //{
    //    var cps = GetSurface().ControlPoints;
    //    for (int i = 0; i < GetCountGlobalBasisFunctions(); i++)
    //    {
    //        var cp = cps[GetINC(i, 0), GetINC(i, 1)];
    //        cp.SetUGlobal(uGlobal);
    //    }
    //}

    public override List<AbstractElement> SelectEndPatchElement(int index)/////////////////////////////////////////////////////////////////
    {
      int indexCoordinate = index / 2;
      int mod = index % 2;
      NURBSSurface surface = GetSurface();
      var basis = surface.Basis;
      int n = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity().Length - 1;
      int m = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity().Length - 1;
      List<AbstractElement> selPatch = new List<AbstractElement>();
      for (int i = 0; i < n * m; i++)
      {
        if (indexCoordinate == 0)
        {
          if (mod == 0)
          {
            if (GetIPN(i, 1) == 0)
              selPatch.Add(listElement[i]);
          }
          else
          {
            if (GetIPN(i, 1) == m - 1)
              selPatch.Add(listElement[i]);
          }
        }
        else
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
      }
      return selPatch;
    }

    //public ControlPoint[] SelectControlPoints(int index)
    //{
    //   var cps = GetSurface().ControlPoints;
    //   ControlPoint[] selCps = null;
    //   int n = 0;
    //   int indexCoordinate = index / 2;
    //   switch (indexCoordinate)
    //   {
    //      case 0:
    //         n = GetSurface().Basis.GetCountBasisFunction(0);
    //         selCps = new ControlPoint[n];
    //         for (int i = 0; i < n; i++)
    //            selCps[i] = cps[i, index];
    //         break;
    //      case 1:

    //         n = GetSurface().Basis.GetCountBasisFunction(1);
    //         selCps = new ControlPoint[n];
    //         for (int i = 0; i < n; i++)
    //            selCps[i] = cps[index, i];
    //         break;
    //   }
    //   return selCps;
    //}
    public ControlPoint[] SelectControlPoints(int index, int indexCoord)
    {
      var cps = GetSurface().ControlPoints;
      ControlPoint[] selCps = null;
      int n = 0;
      switch (indexCoord)
      {
        case 0:
          n = GetSurface().Basis.GetCountBasisFunction(0);
          selCps = new ControlPoint[n];
          for (int i = 0; i < n; i++)
            selCps[i] = cps[i, index];
          break;
        case 1:
          n = GetSurface().Basis.GetCountBasisFunction(1);
          selCps = new ControlPoint[n];
          for (int i = 0; i < n; i++)
            selCps[i] = cps[index, i];
          break;
      }
      return selCps;
    }
    public override ControlPoint[] GetAllControlPoints()
    {
      var cps = GetSurface().ControlPoints;
      ControlPoint[] allCps = new ControlPoint[cps.GetLength(0) * cps.GetLength(1)];
      int count = 0;
      for (int j = 0; j < cps.GetLength(1); j++)
        for (int i = 0; i < cps.GetLength(0); i++)
        {
          allCps[count++] = cps[i, j];
        }
      return allCps;
    }

    /// <summary>
    /// Create ID array
    /// </summary>
    public override int[] EnumerateInPatch()//Enumerate DOF
    {
      int d = GetCountField();
      var cps = GetSurface().ControlPoints;
      enumeratePatch = new int[1][,];
      enumeratePatch[0] = new int[d, cps.GetLength(0) * cps.GetLength(1)];
      int id = 0;
      for (int i = 0; i < enumeratePatch[0].GetLength(1); i++)
      {
        //int psi = GetINC(i, 0);
        //int eta = GetINC(i, 1);
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
      for (int i = 0; i < cps.GetLength(0); i++)
        for (int j = 0; j < cps.GetLength(1); j++)
        {
          cps[i, j].SetDimension(2);
          cps[i, j].SetNumberOfFields(d);

          int[] tArray = new int[d];
          for (int k = 0; k < d; k++)
            tArray[k] = GetIDInPatch(k, FindIndexOfGlobalBasisFunction(i, j));
          cps[i, j].SetTArray(tArray);

          cps[i, j].Initialize();
        }

      //////////////////////////////////////////////////////////////
      //// Create T-Array on Face///////////////////////////////////
      //////////////////////////////////////////////////////////////
      //foreach (Abstract2DElement elem in listElement)
      //{
      //    Face face = elem.GetFace();
      //    face.EnumerateOnFace();
      //}
      return new int[] { id };
    }

    public override int EnumerateInGlobalMultiPatch(int countDof)
    {
      int d = GetCountField();
      var cps = GetSurface().ControlPoints;
      enumerateGlobal = new int[1][,];
      enumerateGlobal[0] = new int[d, cps.GetLength(0) * cps.GetLength(1)];
      for (int i = 0; i < enumerateGlobal[0].GetLength(1); i++)
      {
        int psi = GetINC(i, 0);
        int eta = GetINC(i, 1);
        var cp = cps[psi, eta];
        ////////////////////////////////////////////////////
        ////// coupling all dof ////////////////
        ////////////////////////////////////////
        //if (cp.GetCoupleControlPoint() == null)
        //{
        //   var constraint = cp.GetConstraints();
        //   for (int j = 0; j < d; j++)
        //   {
        //      if ((constraint == null) || constraint.isFree(j))
        //      {
        //         enumerateGlobal[0][j, i] = countDof;
        //         countDof++;
        //      }
        //      else
        //         enumerateGlobal[0][j, i] = -1;
        //   }
        //}
        //else
        //{
        //   for (int j = 0; j < d; j++)
        //   {
        //      enumerateGlobal[0][j, i] = -2;//is coupled
        //   }
        //}
        ///////////////////////////////////////////
        ///////////////////////////////////////////
        ////// coupling select dof ////////////////
        ///////////////////////////////////////////
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
      for (int j = 0; j < cps.GetLength(1); j++)
        for (int i = 0; i < cps.GetLength(0); i++)
        {
          cps[i, j].SetDimension(2);
          cps[i, j].SetNumberOfFields(d);

          int[] tArrayGlobal = new int[d];
          for (int k = 0; k < d; k++)
            tArrayGlobal[k] = enumerateGlobal[0][k, FindIndexOfGlobalBasisFunction(i, j)];
          cps[i, j].SetTArrayGlobal(tArrayGlobal);
        }

      return countDof;
    }

    public override int GetCountLocalBasisFunctions(int idx = 0)
    {
      NURBSSurface surface = GetSurface();
      int p = surface.Basis.GetDegree(0);
      int q = surface.Basis.GetDegree(1);
      return (p + 1) * (q + 1); // number of local basis functions
    }

    public override int GetCountGlobalBasisFunctions(int idx = 0)
    {
      var basis = GetSurface().Basis;
      int n = basis.GetCountBasisFunction(0);
      int m = basis.GetCountBasisFunction(1);
      return n * m; // number of local basis functions
    }

    //public abstract void ComputeStiffnessMatrixPatch(ref DoubleMatrix kGlobal);
    //public void ComputeDataDrawMaterialDistribution(ref List<double> x, ref List<double> y, ref List<double> z, ref List<double> val)
    //{
    //	for (int i = 0; i < listElement.Count; i++)
    //	{
    //		for (int ii = 0; ii < ((Abstract2DElement)listElement[i]).GetNumberOfGaussPoint(); ii++)
    //			for (int jj = 0; jj < ((Abstract2DElement)listElement[i]).GetNumberOfGaussPoint(); jj++)
    //			{
    //				GaussPoints gps = ((Abstract2DElement)listElement[i]).GetGaussPoint(ii, jj);
    //				double[] pointAt = GetSurface().PointAt(gps.location[0], gps.location[1]);
    //				x.Add(pointAt[0]);
    //				y.Add(pointAt[1]);
    //				z.Add(0);
    //				val.Add(gps.EModulus);
    //			}
    //	}
    //}

    public double ComputeArea()
    {
      return ((Abstract2DParametricGeometry)geometry[0]).ComputeArea();
    }
    public double GetMaterialPropertyValueApproximationAt(double xi, double eta)
    {
      double disp = 0;
      var surface = GetSurface();
      var cps = surface.ControlPoints;
      var basis = (BivariateNURBSBasisFunction)surface.Basis;
      int p = basis.GetDegree(0);
      int uSpan = basis.GetKnotVector(0).FindSpan(xi, p);
      int q = basis.GetDegree(1);
      int vSpan = basis.GetKnotVector(1).FindSpan(eta, q);
      double[,] Nuv = basis.GetValueBivariateBasisFunctions(xi, eta);

      for (int i = 0; i <= p; i++)
        for (int j = 0; j <= q; j++)
        {
          disp += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].MaterialPropertyValue;
        }
      return disp;
    }

    public double[,] ComputeMaterialPropertyValue(int resolution)
    {
      ///////////////////////////////////////////////////////////////////////////////
      ////// Parallel ///////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////
      var surface = GetSurface();
      var data = surface.GetParametricOnSurface(resolution, resolution);
      double[,] val = new double[data[0].Length, data[1].Length];
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int j = 0; j < data[1].Length; j++)
          for (int i = 0; i < data[0].Length; i++)
          {
            val[i, j] = GetMaterialPropertyValueApproximationAt(data[0][i], data[1][j]);
          }
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
                 var max = data[0].Length * (taskNumberCopy + 1) / degreeOfParallelism;
                 for (int i = data[0].Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                 {
                   for (int j = 0; j < data[1].Length; j++)
                   {
                     val[i, j] = GetMaterialPropertyValueApproximationAt(data[0][i], data[1][j]);
                   }
                 }
               });
        }
        Task.WaitAll(tasks);
      }
      return val;
    }

    public override double GetApproximateAt(Result result, params double[] xi)
    {
      if (this is PatchStructurePlate)
      {
        double gzz = ((PatchStructurePlate)this).Getzz;
        FunctionRToR fz = ((PatchStructurePlate)this).KinematicsFunction;
        var typeplate = ((PatchStructurePlate)this).TypePlate;
      }
      double disp = 0;
      var surface = GetSurface();
      var cps = surface.ControlPoints;
      var basis = (BivariateNURBSBasisFunction)surface.Basis;
      int p = basis.GetDegree(0);
      int uSpan = basis.GetKnotVector(0).FindSpan(xi[0], p);
      int q = basis.GetDegree(1);
      int vSpan = basis.GetKnotVector(1).FindSpan(xi[1], q);
      double[,] Nuv = basis.GetValueBivariateBasisFunctions(xi[0], xi[1]);
      double sxx = 0, syy = 0, szz = 0, sxy = 0, sxz = 0, syz = 0;
      for (int j = 0; j <= q; j++)
        for (int i = 0; i <= p; i++)
        {
          switch (result)
          {
            case Result.USUM:
              sxx += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.UX);
              syy += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.UY);
              break;
            case Result.SIGMAEQV:
              sxx += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.SIGMAXX);
              syy += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.SIGMAYY);
              sxy += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.SIGMAXY);
              szz += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.SIGMAZZ);
              break;
            case Result.EPSILONEQV:
              sxx += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.EPSILONXX);
              syy += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.EPSILONYY);
              szz += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.EPSILONZZ);
              sxy += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(Result.EPSILONXY);
              break;
            default:
              disp += Nuv[i, j] * cps[uSpan - p + i, vSpan - q + j].GetResult(result);
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
      double[] param = ((NURBSSurface)geometry[0]).Projection(x[0], x[1], 0);
      return GetApproximateAt(result, param[0], param[1]);
    }

    public override double[] CalculateResult(Result re, int resolution)
    {
      var surface = GetSurface();
      var data = surface.GetParametricOnSurface(resolution, resolution);
      double[] val = new double[data[0].Length * data[1].Length];
      if (!AbstractModel.IsParallelProcesing)
      {
        int count = 0;
        for (int i = 0; i < data[0].Length; i++)
          for (int j = 0; j < data[1].Length; j++)
            val[count++] = GetApproximateAt(re, data[0][i], data[1][j]);
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
                 var max = data[0].Length * (taskNumberCopy + 1) / degreeOfParallelism;
                 for (int i = data[0].Length * taskNumberCopy / degreeOfParallelism; i < max; i++)
                 {
                   for (int j = 0; j < data[1].Length; j++)
                     val[i * data[1].Length + j] = GetApproximateAt(re, data[0][i], data[1][j]);
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
    public NURBSSurface GetDeformationSurface(double scale)
    {
      NURBSSurface surface = GetSurface();
      NURBSSurface surfaceCopy = (NURBSSurface)surface.Clone();
      surfaceCopy.isDrawControlNet = surface.isDrawControlNet;
      surfaceCopy.isDrawControlPoint = surface.isDrawControlPoint;
      surfaceCopy.isDrawCurve = surface.isDrawCurve;
      surfaceCopy.isDrawKnot = surface.isDrawKnot;
      surfaceCopy.isColorfulFace = surface.isColorfulFace;
      surfaceCopy.colorSurface = surface.colorSurface;
      surfaceCopy.colorCurve = surface.colorCurve;
      surfaceCopy.colorControlPoint = surface.colorControlPoint;
      surfaceCopy.colorControlNet = surface.colorControlNet;

      surfaceCopy.resolution1 = surface.resolution1;
      surfaceCopy.resolution2 = surface.resolution2;
      surfaceCopy.opacity = surface.opacity;
      var cps = GetSurface().ControlPoints;
      var cpsCopy = surfaceCopy.ControlPoints;
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int i = 0; i < cpsCopy.GetLength(0); i++)
          for (int j = 0; j < cpsCopy.GetLength(1); j++)
          {
            cpsCopy[i, j].SetCoordinate(0, cpsCopy[i, j].GetCoordinate(0) + scale * cps[i, j].GetResult(Result.UX));
            cpsCopy[i, j].SetCoordinate(1, cpsCopy[i, j].GetCoordinate(1) + scale * cps[i, j].GetResult(Result.UY));
            cpsCopy[i, j].SetCoordinate(2, cpsCopy[i, j].GetCoordinate(2) + scale * cps[i, j].GetResult(Result.UZ));
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
                 var max = cpsCopy.GetLength(0) * (taskNumberCopy + 1) / degreeOfParallelism;
                 for (int i = cpsCopy.GetLength(0) * taskNumberCopy / degreeOfParallelism; i < max; i++)
                 {
                   for (int j = 0; j < cpsCopy.GetLength(1); j++)
                   {
                     cpsCopy[i, j].SetCoordinate(0, cpsCopy[i, j].GetCoordinate(0) + scale * cps[i, j].GetResult(Result.UX));
                     cpsCopy[i, j].SetCoordinate(1, cpsCopy[i, j].GetCoordinate(1) + scale * cps[i, j].GetResult(Result.UY));
                     cpsCopy[i, j].SetCoordinate(2, cpsCopy[i, j].GetCoordinate(2) + scale * cps[i, j].GetResult(Result.UZ));
                   }
                 }
               });
        }
        Task.WaitAll(tasks);
      }
      return surfaceCopy;
    }
  }
}
