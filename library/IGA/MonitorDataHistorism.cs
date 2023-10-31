using DEMSoft.Common;
using DEMSoft.NURBS;
using DEMSoft.Plot;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZedGraph;

namespace DEMSoft.IGA
{
  public class MonitorDataHistorism
  {
    private List<MonitorDataNode> listMonitorData;
    private GraphRealTime graph;
    public bool IsGraphRealTime
    { get; set; }
    // private string filePath;
    /// <summary>
    /// Store data monitor
    /// </summary>
    public List<double[]> MonitorDataStore
    { get; set; }
    public int[] IndexColOutputMonitor
    { get; set; }
    public MonitorDataHistorism()
    {
      listMonitorData = new List<MonitorDataNode>();
    }
    public int GetLastIndexMonitorDataStore()
    { return MonitorDataStore.Count - 1; }
    public MonitorDataNode GetMonitorDataNode(int index)
    {
      return listMonitorData[index];
    }
    public int[] GetIndexColUserDefined()
    {
      List<int> list = new List<int>();
      for (int i = 0; i < listMonitorData.Count; i++)
      {
        MonitorDataNode node = listMonitorData[i];
        if (node.ResultType == Result.USERDEFINED)
        {
          list.Add(i);
        }
      }

      if (list.Count > 0)
        return list.ToArray();
      else
        return null;
    }
    public void AddMonitorDataNode(Result resultType, AbstractPatch patch, double xi, double eta, double zeta)
    {
      if (resultType != Result.REACFORCEX && resultType != Result.REACFORCEY && resultType != Result.REACFORCEZ)
      {
        listMonitorData.Add(new MonitorDataNode(resultType, patch, xi, eta, zeta));
      }
      else
      {
        throw new ArgumentException("Do not use this method in case of reaction force");
      }
    }

    public void AddMonitorDataNode(Result resultType, int constraintID)
    {
      if (/*resultType == Result.REACFORCESUM || */resultType == Result.REACFORCEX || resultType == Result.REACFORCEY || resultType == Result.REACFORCEZ)
      {
        listMonitorData.Add(new MonitorDataNode(resultType, constraintID));
      }
      else
      {
        throw new ArgumentException("Only use this method in case of reaction force");
      }
    }
    public delegate double DoFormulationOnMonitorData(double[] valueMonitorData);
    public void AddMonitorDataNode(Result resultType, DoFormulationOnMonitorData doFormulation, int[] indexCol)
    {
      if (resultType == Result.USERDEFINED)
      {
        listMonitorData.Add(new MonitorDataNode(resultType, doFormulation, indexCol));
      }
      else
      {
        throw new ArgumentException("Only use this method for user defined");
      }
    }
    //public double ComputeUserDefinedData(int indexUserDefinedData)
    //{
    //  int c = 0;
    //  for (int i = 0; i < listMonitorData.Count; i++)
    //  {
    //    if (listMonitorData[i].ResultType == Result.USERDEFINED && c == indexUserDefinedData)
    //    {
    //      return listMonitorData[i].ComputeDelegateFormulation();
    //    }
    //    else
    //      c++;
    //  }
    //  return 0;
    //}
    public int GetIndexDataByName(Result result)
    {
      for (int i = 0; i < listMonitorData.Count; i++)
      {
        if (listMonitorData[i].ResultType == result)
          return i;
      }
      return -1;
    }
    public int CountData()
    { return listMonitorData.Count; }
    //public GraphRealTime GetGraphRealTimeForm()
    //{ return graph; }
    //internal Thread myThread;
    public void RunGraph(bool conditionRestoreCheckpoint)
    {
      graph = new GraphRealTime();
      //for (int i = 0; i < listOfGraphPaneInfo.Count; i++)
      //{
      //  GraphPaneProperty g = listOfGraphPaneInfo[i];
      //  graph.AddGraphControl(g.nameTab, g.title, g.xLabel, g.yLabel);
      //}
      //for (int i = 0; i < listOfCurveInfo.Count; i++)
      //{
      //  CurveProperty c = listOfCurveInfo[i];
      //  if (c.idZedGraphControl < listOfGraphPaneInfo.Count)
      //    graph.AddCurveList(c.idZedGraphControl, c.label, c.color, c.symbolType);
      //}
      InitializeGraph(conditionRestoreCheckpoint);
    }
    public List<GraphPaneProperty> listOfGraphPaneInfo;
    public List<CurveProperty> listOfCurveInfo;

    public void AddGraphPane(string nameTab, string title, string xLabel, string yLabel)
    {
      if (listOfGraphPaneInfo == null)
        listOfGraphPaneInfo = new List<GraphPaneProperty>();
      GraphPaneProperty g;
      g.nameTab = nameTab;
      g.title = title;
      g.xLabel = xLabel;
      g.yLabel = yLabel;
      listOfGraphPaneInfo.Add(g);
    }
    public void AddCurve(int idZedGraphControl, int indexColDataX, int indexColDataY, string label, Color color)
    {
      if (listOfCurveInfo == null)
        listOfCurveInfo = new List<CurveProperty>();
      CurveProperty c;
      c.idZedGraphControl = idZedGraphControl;
      c.idCurveInControl = -1;
      c.color = color;
      c.symbolType = SymbolType.None;
      c.indexColX = indexColDataX;
      c.indexColY = indexColDataY;
      c.label = label;
      listOfCurveInfo.Add(c);
    }
    //public void InsertPointToCurve(int idZedGraphControl, int idCurve, double x, double y)
    //{
    //  graph.InsertData(idZedGraphControl, idCurve, x, y);
    //}
    //public void InsertPointToCurve(int idZedGraphControl, int idCurve, double[] x, double[] y)
    //{
    //  graph.InsertData(idZedGraphControl, idCurve, x, y);
    //}
    public void UpdateGraph()
    {
      if (listOfCurveInfo != null)
      {
        int countCurve = listOfCurveInfo.Count;
        for (int i = 0; i < countCurve; i++)
        {
          CurveProperty c = listOfCurveInfo[i];
          int count = MonitorDataStore.Count;
          double x = MonitorDataStore[count - 1][c.indexColX];
          double y = MonitorDataStore[count - 1][c.indexColY];
          graph.InsertData(c.idZedGraphControl, c.idCurveInControl, x, y);
        }
      }
    }
    public void InitializeGraph(bool conditionRestoreCheckpoint)
    {
      //int countControl = listOfGraphPaneInfo.Count;
      //for (int i = 0; i < countControl; i++)
      //{
      //  GraphPaneProperty g = listOfGraphPaneInfo[i];
      //  AddGraphPane(g.nameTab, g.title, g.xLabel, g.yLabel);
      //}
      if (listOfGraphPaneInfo != null)
      {
        for (int i = 0; i < listOfGraphPaneInfo.Count; i++)
        {
          GraphPaneProperty g = listOfGraphPaneInfo[i];
          graph.AddGraphControl(g.nameTab, g.title, g.xLabel, g.yLabel);
        }
        int countCurve = listOfCurveInfo.Count;
        for (int i = 0; i < countCurve; i++)
        {
          CurveProperty c = listOfCurveInfo[i];
          int count = MonitorDataStore.Count;
          if (count == 0 || !conditionRestoreCheckpoint)
          {
            graph.AddCurveList(c.idZedGraphControl, c.label, c.color, c.symbolType);
            int numCurve = graph.GetCountCurves(c.idZedGraphControl);
            c.idCurveInControl = numCurve - 1;
            listOfCurveInfo[i] = c;
          }
          else
          {
            if (c.idZedGraphControl < listOfGraphPaneInfo.Count)
            {
              graph.AddCurveList(c.idZedGraphControl, c.label, c.color, c.symbolType);
              int numCurve = graph.GetCountCurves(c.idZedGraphControl);
              c.idCurveInControl = numCurve - 1;
              listOfCurveInfo[i] = c;
            }
            //double[] x = new double[count];
            //double[] y = new double[count];
            //for (int j = 0; j < count; j++)
            //{
            //  x[j] = MonitorDataStore[j][c.indexColX];
            //  y[j] = MonitorDataStore[j][c.indexColY];
            //}
            //graph.AddPoints(c.idZedGraphControl, c.label, x, y, c.color, c.symbolType);
            for (int j = 0; j < count; j++)
            {
              double x = MonitorDataStore[j][c.indexColX];
              double y = MonitorDataStore[j][c.indexColY];
              graph.InsertData(c.idZedGraphControl, c.idCurveInControl, x, y);
            }
          }
          //int numCurve = graph.GetCountCurves(c.idZedGraphControl);
          //c.idCurveInControl = numCurve - 1;
          //listOfCurveInfo[i] = c;
        }
      }
    }
    public void UpdateUserDefinedDataMonitorStore()
    {
      for (int i = 0; i < listMonitorData.Count; i++)
      {
        MonitorDataNode node = listMonitorData[i];
        if (node.ResultType == Result.USERDEFINED)
        {
          MonitorDataStore[MonitorDataStore.Count - 1][i] = listMonitorData[i].ComputeDelegateFormulation(MonitorDataStore);
        }
      }
    }
  }
  public struct GraphPaneProperty
  {
    public string nameTab;
    public string title;
    public string xLabel;
    public string yLabel;
  }

  public struct CurveProperty
  {
    public int idZedGraphControl;
    public int idCurveInControl;
    public int indexColX;
    public int indexColY;
    public string label;
    public Color color;
    public SymbolType symbolType;
  }
}
