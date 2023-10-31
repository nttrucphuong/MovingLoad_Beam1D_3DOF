using DEMSoft.Common;
using System.Collections.Generic;
using static DEMSoft.IGA.MonitorDataHistorism;

namespace DEMSoft.IGA
{
  public class MonitorDataNode
  {
    //private double[] xi;
    //private Result storeResult;
    //private AbstractPatch patch;
    public MonitorDataNode(Result re, AbstractPatch patch, double xi, double eta, double zeta)
    {
      Xi = new double[] { xi, eta, zeta };
      ResultType = re;
      Patch = patch;
    }

    public MonitorDataNode(Result re, int constraintID)
    {
      ResultType = re;
      ConstraintID = constraintID;
    }
    //public delegate double DoFormulationOnMonitorData(int[] indexColData);
    private DoFormulationOnMonitorData doFormulation;
    private int[] indexCol;
    public MonitorDataNode(Result re, DoFormulationOnMonitorData doFormulation, params int[] indexCol)
    {
      ResultType = re;
      this.doFormulation = doFormulation;
      this.indexCol = indexCol;
    }
    public double ComputeDelegateFormulation(List<double[]> monitorData)
    {
      List<double> values = new List<double>();
      for (int i = 0; i < indexCol.Length; i++)
      {
        values.Add(monitorData[monitorData.Count - 1][indexCol[i]]);
      }
      return doFormulation(values.ToArray());
    }
    public double[] Xi
    {
      //get { return xi; }
      //set { this.xi = value; }
      get; set;
    }

    public Result ResultType
    {
      //get { return storeResult; }
      //set { storeResult = value; }
      get; set;
    }

    public AbstractPatch Patch
    {
      //get { return patch; }
      //set { patch = value; }
      get; set;
    }
    public int ConstraintID
    { get; set; }
  }
}
