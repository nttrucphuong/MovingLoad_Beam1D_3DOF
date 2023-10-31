
using CenterSpace.NMath.Core;
using DEMSoft.NURBS;
using System;

namespace DEMSoft.IGA
{
  public class CouplingTwoControlPoints : ICoupling
  {
    private ControlPoint masterCp;
    private ControlPoint slaveCp;
    private int[] DOFUnCoupling;

    public CouplingTwoControlPoints(ControlPoint masterCp, ControlPoint slaveCp, params int[] DOFUnCoupling)
    {
      this.masterCp = masterCp;
      this.slaveCp = slaveCp;
      this.DOFUnCoupling = DOFUnCoupling;
    }

    public void MakeCoupling(ref DoubleMatrix T)
    {
      int[] tArrayMaster = masterCp.GetTArrayGlobal();
      int[] tArraySlave = slaveCp.GetTArrayGlobal();
      for (int i = 0; i < Math.Min(tArrayMaster.Length, tArraySlave.Length); i++)
      {
        int index = Array.IndexOf(DOFUnCoupling, i);
        if (index == -1)
        {
          if (T[tArraySlave[i], tArraySlave[i]] == 0 || T[tArraySlave[i], tArrayMaster[i]] != 1)
          {
            T[tArraySlave[i], tArraySlave[i]] += -1;
            T[tArraySlave[i], tArrayMaster[i]] += 1;
          }
        }
      }
    }
    public void MakeCoupling(ref SparseMatrixBuilder<double> T)
    {
      int[] tArrayMaster = masterCp.GetTArrayGlobal();
      int[] tArraySlave = slaveCp.GetTArrayGlobal();
      for (int i = 0; i < Math.Min(tArrayMaster.Length, tArraySlave.Length); i++)
      {
        int index = Array.IndexOf(DOFUnCoupling, i);
        if (index == -1)
        {
          if (T[tArraySlave[i], tArraySlave[i]] == 0 || T[tArraySlave[i], tArrayMaster[i]] != 1)
          {
            T[tArraySlave[i], tArraySlave[i]] += -1;
            T[tArraySlave[i], tArrayMaster[i]] += 1;
          }
        }
      }
    }
  }
}
