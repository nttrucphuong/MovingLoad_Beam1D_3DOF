using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using DEMSoft.Common;
using DEMSoft.NURBS;

namespace DEMSoft.IGA
{
  /// <summary>
  /// Abstract class to define NURBS surface to NURBS Patch
  /// </summary>
  public abstract class AbstractPatch3D : AbstractPatchOneField
  {
    public double[,,][,] ExtractionOperator
    { get; set; }
    public override void ComputeExtractionOperator()
    {
      ExtractionOperator = GetVolume().ComputeExtractionOperator();
    }
    /// <summary>
    /// Constructor class
    /// </summary>
    /// <param name="volume">Nurbs volume</param>
    /// <param name="numberOfFields">number of fields in model</param>
    public AbstractPatch3D(NURBSVolume volume, int numberOfFields)
       : base(3)
    {
      //connect = new List<ConnectionInterfaceTwoPatches>();
      geometry = new NURBSVolume[] { volume };
      //countField = new int[] { numberOfFields };
      //CreateIPN();
      //CreateINC();
      //CreateIEN();
      //Initialize();
    }

    /// <summary>
    /// Get NURBS surface
    /// </summary>
    /// <returns></returns>
    public NURBSVolume GetVolume()
    {
      return (NURBSVolume)GetGeometry();
    }



    /// <summary>
    /// Create INC array
    /// </summary>
    protected override void CreateINC()
    {
      var basis = GetVolume().Basis;
      int n = basis.GetCountBasisFunction(0);
      int m = basis.GetCountBasisFunction(1);
      int l = basis.GetCountBasisFunction(2);
      int nen = GetCountLocalBasisFunctions();// (p + 1) * (q + 1); // number of local basis functions
      int nnp = GetCountGlobalBasisFunctions();// n * m;//number of global basis functions
      INC = new int[1][,];
      INC[0] = new int[nnp, countDimension];
      int A = 0;

      for (int k = 0; k < l; k++)
        for (int j = 0; j < m; j++)
          for (int i = 0; i < n; i++)
          {
            INC[0][A, 0] = i;
            INC[0][A, 1] = j;
            INC[0][A, 2] = k;
            A++;
          }
    }

    /// <summary>
    /// Create INC array
    /// </summary>
    protected override void CreateIEN()
    {
      var basis = GetVolume().Basis;
      int n = basis.GetCountBasisFunction(0);
      int p = basis.GetDegree(0);
      int m = basis.GetCountBasisFunction(1);
      int q = basis.GetDegree(1);
      int l = basis.GetCountBasisFunction(2);
      int r = basis.GetDegree(2);
      var knotVectorNoMultiplicity1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var knotVectorNoMultiplicity2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      var knotVectorNoMultiplicity3 = basis.GetKnotVector(2).GetKnotVectorNoMultiplicity();
      int numelem1 = knotVectorNoMultiplicity1.Length - 1;
      int numelem2 = knotVectorNoMultiplicity2.Length - 1;
      int numelem3 = knotVectorNoMultiplicity3.Length - 1;
      int nel = numelem1 * numelem2 * numelem3;// GetNumberOfPatchs();//(n - p) * (m - q);//number of elements
      int nen = (p + 1) * (q + 1) * (r + 1); //GetNumberOfLocalBasisFunctions(); // number of local basis functions
      IEN = new int[1][,];
      IEN[0] = new int[nen, nel];
      //int e = 0;
      //int A = 0;
      for (int ek = 0; ek < numelem3; ek++)
        for (int ej = 0; ej < numelem2; ej++)
          for (int ei = 0; ei < numelem1; ei++)
          {
            double mid1 = (knotVectorNoMultiplicity1[ei] + knotVectorNoMultiplicity1[ei + 1]) / 2.0;
            double mid2 = (knotVectorNoMultiplicity2[ej] + knotVectorNoMultiplicity2[ej + 1]) / 2.0;
            double mid3 = (knotVectorNoMultiplicity3[ek] + knotVectorNoMultiplicity3[ek + 1]) / 2.0;
            int uspan = basis.FindSpan(mid1, 0);
            int vspan = basis.FindSpan(mid2, 1);
            int wspan = basis.FindSpan(mid3, 2);
            int nume = FindIndexOfElement(ei, ej, ek);
            int b = 0;
            for (int k = 0; k <= r; k++)
              for (int j = 0; j <= q; j++)
                for (int i = 0; i <= p; i++)
                {
                  int num = FindIndexOfGlobalBasisFunction(uspan - p + i, vspan - q + j, wspan - r + k);
                  IEN[0][b, nume] = num;
                  b++;
                }
          }
    }

    /// <summary>
    /// Create ID array
    /// </summary>
    public override int[] EnumerateInPatch()//Enumerate DOF
    {
      int countField = GetCountField();
      var cps = GetVolume().ControlPoints;
      enumeratePatch = new int[1][,];
      enumeratePatch[0] = new int[countField, cps.GetLength(0) * cps.GetLength(1) * cps.GetLength(2)];
      int id = 0;
      for (int i = 0; i < enumeratePatch[0].GetLength(1); i++)
      {
        //int psi = INC[0][i, 0];
        //int eta = INC[0][i, 1];
        //int zeta = INC[0][i, 2];
        for (int j = 0; j < countField; j++)
        {
          enumeratePatch[0][j, i] = id;
          id++;
        }
      }
      countDOF = id;

      for (int i = 0; i < cps.GetLength(0); i++)
        for (int j = 0; j < cps.GetLength(1); j++)
          for (int k = 0; k < cps.GetLength(2); k++)
          {
            cps[i, j, k].SetDimension(3);
            cps[i, j, k].SetNumberOfFields(countField);
            int[] tArray = new int[countField];
            for (int kk = 0; kk < countField; kk++)
              tArray[kk] = GetIDInPatch(kk, FindIndexOfGlobalBasisFunction(i, j, k));
            cps[i, j, k].SetTArray(tArray);
            cps[i, j, k].Initialize();
          }

      //////////////////////////////////////////////////////////////
      //// Create T-Array on Face///////////////////////////////////
      //////////////////////////////////////////////////////////////
      //foreach (Structure3DElement elem in listElement)
      //{
      //    Volume volume = elem.GetVolume();
      //    volume.EnumerateOnVolume();
      //}
      return new int[] { id };
    }

    /// <summary>
    /// Create index of patch on direction
    /// </summary>
    protected override void CreateIPN()
    {
      var basis = GetVolume().Basis;
      int n = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity().Length - 1;
      int m = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity().Length - 1;
      int l = basis.GetKnotVector(2).GetKnotVectorNoMultiplicity().Length - 1;
      int nnp = n * m * l;
      IPN = new int[nnp, countDimension];
      int A = 0;
      for (int k = 0; k < l; k++)
        for (int j = 0; j < m; j++)
          for (int i = 0; i < n; i++)
          {
            IPN[A, 0] = i;
            IPN[A, 1] = j;
            IPN[A, 2] = k;
            A++;
          }
    }

    /// <summary>
    /// Find index of element at coordinate
    /// </summary>
    /// <param name="xi"></param>
    /// <param name="eta"></param>
    /// <param name="zeta"></param>
    /// <returns></returns>
    public override int FindIndexOfElementAt(params double[] xi)
    {
      int idxi = -1;
      int idxj = -1;
      int idxk = -1;
      var kv1 = GetVolume().Basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kv2 = GetVolume().Basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      var kv3 = GetVolume().Basis.GetKnotVector(2).GetKnotVectorNoMultiplicity();
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
      for (int i = 0; i < kv3.Length; i++)
      {
        if (xi[2] >= kv3[i])
          idxk = i;
        if (xi[2] == kv3[kv3.Length - 1])
          idxk = i - 1;
      }
      return FindIndexOfElement(idxi, idxj, idxk);
    }

    /// <summary>
    /// Get number of patchs
    /// </summary>
    /// <returns></returns>
    public override int CalculateNumberOfElements()
    {
      NURBSVolume volume = GetVolume();
      var basis = volume.Basis;
      int n = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity().Length - 1;
      int m = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity().Length - 1;
      int l = basis.GetKnotVector(2).GetKnotVectorNoMultiplicity().Length - 1;
      return n * m * l;
    }

    /// <summary>
    /// Get number of local basis functions
    /// </summary>
    /// <returns></returns>
    public override int GetCountLocalBasisFunctions(int idx = 0)
    {
      NURBSVolume volume = GetVolume();
      int p = volume.Basis.GetDegree(0);
      int q = volume.Basis.GetDegree(1);
      int r = volume.Basis.GetDegree(2);
      return (p + 1) * (q + 1) * (r + 1); // number of local basis functions
    }

    /// <summary>
    /// Get number of global basis functions
    /// </summary>
    /// <returns></returns>
    public override int GetCountGlobalBasisFunctions(int idx = 0)
    {
      var basis = GetVolume().Basis;
      int n = basis.GetCountBasisFunction(0);
      int m = basis.GetCountBasisFunction(1);
      int l = basis.GetCountBasisFunction(2);
      return n * m * l; // number of local basis functions
    }

    //public override void SetUGlobal(DoubleVector uGlobal)
    //{
    //	var cps = GetVolume().ControlPoints;
    //	for (int i = 0; i < GetCountGlobalBasisFunctions(); i++)
    //	{
    //		var cp = cps[GetINC(i, 0), GetINC(i, 1), GetINC(i, 2)];
    //		cp.SetUGlobal(uGlobal);
    //	}
    //}

    public override double GetApproximateAt(Result result, params double[] xi)
    {
      double disp = 0;
      var vol = GetVolume();
      var cps = vol.ControlPoints;
      var basis = (TrivariateNURBSBasisFunction)vol.Basis;
      int p = basis.GetDegree(0);
      int uSpan = basis.GetKnotVector(0).FindSpan(xi[0], p);
      int q = basis.GetDegree(1);
      int vSpan = basis.GetKnotVector(1).FindSpan(xi[1], q);
      int r = basis.GetDegree(2);
      int wSpan = basis.GetKnotVector(2).FindSpan(xi[2], r);
      double[,,] Nuvw = basis.GetValueTrivariateBasisFunctions(xi[0], xi[1], xi[2]);
      double sxx = 0, syy = 0, sxy = 0, szz = 0, syz = 0, sxz = 0;
      for (int i = 0; i <= p; i++)
        for (int j = 0; j <= q; j++)
          for (int k = 0; k <= r; k++)
          {
            switch (result)
            {
              case Result.USUM:
                sxx += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.UX);
                syy += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.UY);
                szz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.UZ);
                break;
              case Result.SIGMAEQV:
                sxx += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.SIGMAXX);
                syy += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.SIGMAYY);
                szz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.SIGMAZZ);
                sxy += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.SIGMAXY);
                syz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.SIGMAYZ);
                sxz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.SIGMAXZ);
                break;
              case Result.EPSILONEQV:
                sxx += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.EPSILONXX);
                syy += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.EPSILONYY);
                szz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.EPSILONZZ);
                sxy += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.EPSILONXY);
                syz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.EPSILONYZ);
                sxz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.EPSILONXZ);
                break;
              default:
                disp += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(result);
                break;
            }
          }
      switch (result)
      {
        case Result.USUM:
          disp = Math.Sqrt(sxx * sxx + syy * syy + szz * szz);
          break;
        case Result.SIGMAEQV:
          disp = Math.Sqrt(1.0 / 2.0 * (Math.Pow(sxx - syy, 2) + Math.Pow(syy - szz, 2) + Math.Pow(sxx - szz, 2) + 6.0 * (Math.Pow(sxy, 2) + Math.Pow(syz, 2) + Math.Pow(sxz, 2))));
          break;
        case Result.EPSILONEQV:
          double xx = 2.0 / 3.0 * sxx - 1.0 / 3.0 * syy - 1.0 / 3.0 * szz;
          double yy = -1.0 / 3.0 * sxx + 2.0 / 3.0 * syy - 1.0 / 3.0 * szz;
          double zz = -1.0 / 3.0 * sxx - 1.0 / 3.0 * syy + 2.0 / 3.0 * szz;
          disp = 2.0 / 3.0 * Math.Sqrt(3.0 / 2.0 * (Math.Pow(xx, 2) + Math.Pow(yy, 2) + Math.Pow(zz, 2)) + 3.0 / 4.0 * (Math.Pow(2 * sxy, 2) + Math.Pow(2 * syz, 2) + Math.Pow(2 * sxz, 2)));
          break;
      }
      return disp;
    }
    public override double GetApproximateOnGlobalAt(Result result, params double[] x)
    {
      double[] param = ((NURBSVolume)geometry[0]).Projection(x[0], x[1], x[2]);
      return GetApproximateAt(result, param[0], param[1], param[2]);
    }
    /// <summary>
    /// Select serial of patchs from mesh by index on each direction
    /// </summary>
    /// <param name="index">[0], [1]: u-direction, [2], [3]: v-direction, [4], [5]: w-direction</param>
    /// <returns></returns>
    public override List<AbstractElement> SelectEndPatchElement(int index)
    {
      //int indexCoordinate = -1;

      //if (index == 0 || index == 1)
      //    indexCoordinate = 0;
      //else if (index == 2 || index == 3)
      //    indexCoordinate = 1;
      //else if (index == 4 || index == 5)
      //    indexCoordinate = 2;

      int indexCoordinate = index / 2;
      int mod = index % 2;

      NURBSVolume volume = GetVolume();
      var basis = volume.Basis;
      int n = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity().Length - 1;
      int m = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity().Length - 1;
      int l = basis.GetKnotVector(2).GetKnotVectorNoMultiplicity().Length - 1;
      List<AbstractElement> selElement = new List<AbstractElement>();
      for (int i = 0; i < n * m * l; i++)
      {
        switch (indexCoordinate)
        {
          case 0:
            if (mod == 0)
            {
              if (GetIPN(i, 2) == 0)
                selElement.Add(listElement[i]);
            }
            else
            {
              if (GetIPN(i, 2) == l - 1)
                selElement.Add(listElement[i]);
            }
            break;
          case 1:
            if (mod == 0)
            {
              if (GetIPN(i, 1) == 0)
                selElement.Add(listElement[i]);
            }
            else
            {
              if (GetIPN(i, 1) == m - 1)
                selElement.Add(listElement[i]);
            }
            break;
          case 2:
            if (mod == 0)
            {
              if (GetIPN(i, 0) == 0)
                selElement.Add(listElement[i]);
            }
            else
            {
              if (GetIPN(i, 0) == n - 1)
                selElement.Add(listElement[i]);
            }
            break;
        }
      }
      return selElement;
    }

    public ControlPoint[,] SelectControlPoints(int index, int indexCoordinate)
    {
      //List<ControlPoint> selCps = new List<ControlPoint>();// ControlPoint[GetVolume().Basis.GetCountBasisFunction(indexCoordinate)];
      //var cps = GetVolume().ControlPoints;
      //for (int i = 0; i < GetVolume().Basis.GetCountBasisFunction(indexCoordinate); i++)
      //{
      //    if (indexCoordinate == 0)
      //        selCps[i] = cps[i, index];
      //    else
      //        selCps[i] = cps[index, i];
      //}
      //return selCps;
      return null;
    }

    /// <summary>
    /// Select serial of patchs from mesh by index on each direction
    /// </summary>
    /// <param name="index">index = 0-first column or row of patch on mesh, 1-last column or row of patch on mesh</param>
    /// <param name="indexCoordinate">index of direction coordinate</param>
    /// <returns></returns>
    public ControlPoint[,] SelectEndPatchControlPoints(int index)
    {
      return GetVolume().SelectControlPointInTwoEndAreas(index);
    }
    public ControlPoint[,] SelectNearEndPatchControlPoints(int index)
    {
      return ((Abstract3DParametricGeometry)GetGeometry()).SelectControlPointNearInTwoEndAreas(index);
    }
    public override ControlPoint[] GetAllControlPoints()
    {
      var cps = GetVolume().ControlPoints;
      ControlPoint[] allCps = new ControlPoint[cps.GetLength(0) * cps.GetLength(1) * cps.GetLength(2)];
      int count = 0;
      for (int k = 0; k < cps.GetLength(2); k++)
        for (int j = 0; j < cps.GetLength(1); j++)
          for (int i = 0; i < cps.GetLength(0); i++)
          {
            allCps[count++] = cps[i, j, k];
          }
      return allCps;
    }

    public override int EnumerateInGlobalMultiPatch(int countDof)
    {
      var cps = GetVolume().ControlPoints;
      int d = GetCountField();
      enumerateGlobal = new int[1][,];
      enumerateGlobal[0] = new int[d, cps.GetLength(0) * cps.GetLength(1) * cps.GetLength(2)];
      for (int i = 0; i < enumeratePatch[0].GetLength(1); i++)
      {
        int psi = INC[0][i, 0];
        int eta = INC[0][i, 1];
        int zeta = INC[0][i, 2];
        var cp = cps[psi, eta, zeta];
        ////////////////////////////////////////////////////////////
        ////////////////// couping all dof //////////////////////////
        //	if (cp.GetCoupleControlPoint() == null)
        //	{
        //		var constraint = cp.GetConstraints();
        //		for (int j = 0; j < countField; j++)
        //		{
        //			if ((constraint == null) || constraint.isFree(j))
        //			{
        //				enumerateGlobal[0][j, i] = countDof;
        //				countDof++;
        //			}
        //			else
        //				enumerateGlobal[0][j, i] = -1;
        //		}
        //	}
        //	else
        //	{
        //		for (int j = 0; j < countField; j++)
        //		{
        //			enumerateGlobal[0][j, i] = -2;//is coupled
        //		}
        //	}
        //////////////////////////////////////////////////////////
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
      }
      for (int k = 0; k < cps.GetLength(2); k++)
        for (int j = 0; j < cps.GetLength(1); j++)
          for (int i = 0; i < cps.GetLength(0); i++)
          {
            cps[i, j, k].SetDimension(3);
            cps[i, j, k].SetNumberOfFields(d);
            int[] tArrayGlobal = new int[d];
            for (int kk = 0; kk < d; kk++)
              tArrayGlobal[kk] = GetIDInGlobal(kk, FindIndexOfGlobalBasisFunction(i, j, k));
            cps[i, j, k].SetTArrayGlobal(tArrayGlobal);
          }

      return countDof;
    }

    public int[] GetCoordinateParameterOnArea(int index)
    {
      int indexCoordinateFace = index / 2;// old version wrong /3
      int[] indexCoordinate = new int[2];
      switch (indexCoordinateFace)
      {
        case 0:
          indexCoordinate[0] = 0;
          indexCoordinate[1] = 1;
          break;
        case 1:
          indexCoordinate[0] = 0;
          indexCoordinate[1] = 2;
          break;
        case 2:
          indexCoordinate[0] = 1;
          indexCoordinate[1] = 2;
          break;
      }
      return indexCoordinate;
    }

    public double ComputeVolume()
    {
      double volume = 0;
      //if (!AbstractModel.IsParallelProcesing)
      for (int i = 0; i < listElement.Count; i++)
      {
        volume += ((AbstractElement3D)listElement[i]).ComupteVolumeOfElement();
      }
      //else
      //{
      //   Parallel.For(0, listElement.Count, i =>
      //   {
      //      volume += ((AbstractElement3D)listElement[i]).ComupteVolumeOfElement();
      //   });
      //}
      return volume;
    }

    public double GetMaterialPropertyValueApproximationAt(double xi, double eta, double zeta)
    {
      double disp = 0;
      var volume = GetVolume();
      var cps = volume.ControlPoints;
      var basis = (TrivariateNURBSBasisFunction)volume.Basis;
      int p = basis.GetDegree(0);
      int uSpan = basis.GetKnotVector(0).FindSpan(xi, p);
      int q = basis.GetDegree(1);
      int vSpan = basis.GetKnotVector(1).FindSpan(eta, q);
      int r = basis.GetDegree(2);
      int wSpan = basis.GetKnotVector(2).FindSpan(zeta, r);
      double[,,] Nuv = basis.GetValueTrivariateBasisFunctions(xi, eta, zeta);

      for (int i = 0; i <= p; i++)
        for (int j = 0; j <= q; j++)
          for (int k = 0; k <= r; k++)
          {
            disp += Nuv[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].MaterialPropertyValue;
          }
      return disp;
    }
    public double[,,] ComputeMaterialPropertyValue(int resolution)
    {
      ///////////////////////////////////////////////////////////////////////////////
      ////// Parallel ///////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////
      var volume = GetVolume();
      var data = volume.GetParametricOnVolume(resolution, resolution, resolution);
      double[,,] val = new double[data[0].Length, data[1].Length, data[2].Length];
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int k = 0; k < data[2].Length; k++)
          for (int j = 0; j < data[1].Length; j++)
            for (int i = 0; i < data[0].Length; i++)
            {
              val[i, j, k] = GetMaterialPropertyValueApproximationAt(data[0][i], data[1][j], data[2][k]);
            }
      }
      else
      {
        //Parallel.For(0, data[2].Length, k =>
        //{
        //  for (int j = 0; j < data[1].Length; j++)
        //    for (int i = 0; i < data[0].Length; i++)
        //    {
        //      val[i, j, k] = GetMaterialPropertyValueApproximationAt(data[0][i], data[1][j], data[2][k]);
        //    }
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
                     for (int k = 0; k < data[2].Length; k++)
                       val[i, j, k] = GetMaterialPropertyValueApproximationAt(data[0][i], data[1][j], data[2][k]);
                 }
               });
        }
        Task.WaitAll(tasks);
      }
      return val;
    }

    public override double[] CalculateResult(Result re, int resolution)
    {
      var volume = GetVolume();
      var data = volume.GetParametricOnVolume(resolution, resolution, resolution);
      double[] val = new double[data[0].Length * data[1].Length * data[2].Length];
      if (!AbstractModel.IsParallelProcesing)
      {
        int count = 0;
        for (int i = 0; i < data[0].Length; i++)
          for (int j = 0; j < data[1].Length; j++)
            for (int k = 0; k < data[2].Length; k++)
              val[count++] = GetApproximateAt(re, data[0][i], data[1][j], data[2][k]);
      }
      else
      {
        //Parallel.For(0, data[0].Length, i =>
        //{
        //  for (int j = 0; j < data[1].Length; j++)
        //    for (int k = 0; k < data[2].Length; k++)
        //      val[i * data[1].Length * data[2].Length + j * data[2].Length + k] = GetApproximateAt(re, data[0][i], data[1][j], data[2][k]);
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
                     for (int k = 0; k < data[2].Length; k++)
                       val[i * data[1].Length * data[2].Length + j * data[2].Length + k] = GetApproximateAt(re, data[0][i], data[1][j], data[2][k]);
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
    public NURBSVolume GetDeformationVolume(double scale)
    {
      NURBSVolume vol = GetVolume();
      NURBSVolume volumeCopy = (NURBSVolume)vol.Clone();
      volumeCopy.isDrawGrid = vol.isDrawGrid;
      volumeCopy.isDrawVolume = vol.isDrawVolume;
      volumeCopy.isDrawControlPoint = vol.isDrawControlPoint;
      volumeCopy.isDrawControlNet = vol.isDrawControlNet;
      volumeCopy.isColorfulFace = vol.isColorfulFace;
      volumeCopy.isDrawCurve = vol.isDrawCurve;
      volumeCopy.isDrawKnot = vol.isDrawKnot;
      volumeCopy.colorVolume = vol.colorVolume;
      volumeCopy.colorGrid = vol.colorGrid;
      volumeCopy.colorCurve = vol.colorCurve;
      volumeCopy.colorControlPoint = vol.colorControlPoint;
      volumeCopy.colorControlNet = vol.colorControlNet;
      volumeCopy.resolution1 = vol.resolution1;
      volumeCopy.resolution2 = vol.resolution2;
      volumeCopy.resolution3 = vol.resolution3;
      volumeCopy.opacity = vol.opacity;

      var cps = GetVolume().ControlPoints;
      var cpsCopy = volumeCopy.ControlPoints;
      if (!AbstractModel.IsParallelProcesing)
      {
        for (int i = 0; i < cpsCopy.GetLength(0); i++)
          for (int j = 0; j < cpsCopy.GetLength(1); j++)
            for (int k = 0; k < cpsCopy.GetLength(2); k++)
            {
              cpsCopy[i, j, k].SetCoordinate(0, cpsCopy[i, j, k].GetCoordinate(0) + scale * cps[i, j, k].GetResult(Result.UX));
              cpsCopy[i, j, k].SetCoordinate(1, cpsCopy[i, j, k].GetCoordinate(1) + scale * cps[i, j, k].GetResult(Result.UY));
              cpsCopy[i, j, k].SetCoordinate(2, cpsCopy[i, j, k].GetCoordinate(2) + scale * cps[i, j, k].GetResult(Result.UZ));
            }
      }
      else
      {
        //Parallel.For(0, cpsCopy.GetLength(0), i =>
        //{
        //  for (int j = 0; j < cpsCopy.GetLength(1); j++)
        //    for (int k = 0; k < cpsCopy.GetLength(2); k++)
        //    {
        //      cpsCopy[i, j, k].SetCoordinate(0, cpsCopy[i, j, k].GetCoordinate(0) + scale * cps[i, j, k].GetResult(Result.UX));
        //      cpsCopy[i, j, k].SetCoordinate(1, cpsCopy[i, j, k].GetCoordinate(1) + scale * cps[i, j, k].GetResult(Result.UY));
        //      cpsCopy[i, j, k].SetCoordinate(2, cpsCopy[i, j, k].GetCoordinate(2) + scale * cps[i, j, k].GetResult(Result.UZ));
        //    }
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
                     for (int k = 0; k < cpsCopy.GetLength(2); k++)
                     {
                       cpsCopy[i, j, k].SetCoordinate(0, cpsCopy[i, j, k].GetCoordinate(0) + scale * cps[i, j, k].GetResult(Result.UX));
                       cpsCopy[i, j, k].SetCoordinate(1, cpsCopy[i, j, k].GetCoordinate(1) + scale * cps[i, j, k].GetResult(Result.UY));
                       cpsCopy[i, j, k].SetCoordinate(2, cpsCopy[i, j, k].GetCoordinate(2) + scale * cps[i, j, k].GetResult(Result.UZ));
                     }
                 }
               });
        }
        Task.WaitAll(tasks);
      }
      return volumeCopy;
    }
  }
}
