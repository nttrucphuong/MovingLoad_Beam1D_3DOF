using System.Collections.Generic;
using System.Linq;
using CenterSpace.NMath.Core;
using DEMSoft.NURBS;

namespace DEMSoft.IGA
{
  class GenerateAutomaticConnectionInterfacePatches
  {
    private AbstractModel model;
    public GenerateAutomaticConnectionInterfacePatches(AbstractModel model)
    {
      this.model = model;
    }

    public void GenerateConnectionInterfaceAuto(bool Is1_1Constraint, bool IsUseC1Constraint, ref List<ControlPoint> listControlPointCoincidented, params IntPair[] unCouplingPatches)
    {
      int numberOfPatch = model.CountPatch();
      for (int i = 0; i < numberOfPatch; i++)
      {
        for (int j = i; j < numberOfPatch; j++)
        {
          if (unCouplingPatches != null && (unCouplingPatches.Contains(new IntPair(i, j)) || unCouplingPatches.Contains(new IntPair(j, i))))
          {
            continue;
          }

          AbstractPatch p1 = model.GetPatch(i);
          AbstractPatch p2 = model.GetPatch(j);
          int[] connection = CheckConnectionInterfaceTwoPatches(p1, p2, Is1_1Constraint);
          if (connection != null)
          {
            model.SetInterfaceBetweenTwoPatches(i, j, connection[0], connection[1], null, Is1_1Constraint, IsUseC1Constraint, ref listControlPointCoincidented);//Chua cai tien bat cap tu dong khong dung 1:1
          }
        }
      }
    }

    public void GenerateConnectionInterfaceAuto(List<int[]> listConnectionPatch, bool Is1_1Constraint, bool IsUseC1Constraint, ref List<ControlPoint> listControlPointCoincidented)
    {
      for (int i = 0; i < listConnectionPatch.Count; i++)
      {
        model.SetInterfaceBetweenTwoPatches(listConnectionPatch[i][0], listConnectionPatch[i][1], listConnectionPatch[i][2], listConnectionPatch[i][3], null, Is1_1Constraint, IsUseC1Constraint, ref listControlPointCoincidented);
      }
    }

    public List<int[]> GenerateInformationConnectionInterfaceAuto(bool Is1_1Constraint, params IntPair[] unCouplingPatches)
    {
      List<int[]> listConnection = new List<int[]>();
      int numberOfPatch = model.CountPatch();
      for (int i = 0; i < numberOfPatch; i++)
      {
        for (int j = i; j < numberOfPatch; j++)
        {
          if (unCouplingPatches != null && (unCouplingPatches.Contains(new IntPair(i, j)) || unCouplingPatches.Contains(new IntPair(j, i))))
          {
            continue;
          }
          AbstractPatch p1 = model.GetPatch(i);
          AbstractPatch p2 = model.GetPatch(j);
          int[] connection = CheckConnectionInterfaceTwoPatches(p1, p2, Is1_1Constraint);
          if (connection != null)
          {
            listConnection.Add(new int[] { i, j, connection[0], connection[1] });
          }
        }
      }
      return listConnection;
    }

    public int[] CheckConnectionInterfaceTwoPatches(AbstractPatch masterPatch, AbstractPatch slavePatch, bool Is1_1Constraint)
    {
      int[] connection = null;

      if (masterPatch is AbstractPatch2D && slavePatch is AbstractPatch2D)
      {
        AbstractPatch2D master2D = (AbstractPatch2D)masterPatch;
        AbstractPatch2D slave2D = (AbstractPatch2D)slavePatch;
        var surfaceMaster = master2D.GetSurface();
        var surfaceSlave = slave2D.GetSurface();
        //for (int indexMaster = 0; indexMaster < 4; indexMaster++)
        //{
        //  for (int indexSlave = 0; indexSlave < 4; indexSlave++)
        //  {
        //    if (masterPatch == slavePatch && indexMaster == indexSlave)
        //      continue;

        //    var cpsMaster = master2D.SelectEndPatchControlPoints(indexMaster);
        //    var cpsSlave = slave2D.SelectEndPatchControlPoints(indexSlave);
        //    int numCpsMaster = cpsMaster.Length;
        //    int numCpsSlave = cpsSlave.Length;

        //    ////////////// coupling 1:1
        //    bool flagContinue = false;
        //    if (Is1_1Constraint)
        //    {
        //      int indexCoordinateMasterLine = indexMaster / 2;
        //      int indexCoordinateSlaveLine = indexSlave / 2;

        //      KnotVector kvMaster = (indexCoordinateMasterLine == 0) ? surfaceMaster.Basis.GetKnotVector(0)
        //          : surfaceMaster.Basis.GetKnotVector(1);

        //      KnotVector kvSlave = (indexCoordinateSlaveLine == 0) ? surfaceSlave.Basis.GetKnotVector(0)
        //          : surfaceSlave.Basis.GetKnotVector(1);

        //      if (!kvMaster.IsMatch(kvSlave))
        //        continue;

        //      if (numCpsMaster != numCpsSlave)
        //        continue;
        //      for (int i = 0; i < numCpsMaster; i++)
        //      {
        //        if ((!cpsMaster[i].IsCoincident(cpsSlave[i])) && (!cpsMaster[i].IsCoincident(cpsSlave[numCpsMaster - 1 - i])))
        //        {
        //          flagContinue = true;
        //          break;
        //        }
        //      }
        //      ////////////////
        //      ////// nonconforming
        //    }
        //    else
        //    {
        //      flagContinue = !((cpsMaster[0].IsCoincident(cpsSlave[0]) && cpsMaster[numCpsMaster - 1].IsCoincident(cpsSlave[numCpsSlave - 1]))
        //      || (cpsMaster[0].IsCoincident(cpsSlave[numCpsSlave - 1]) && cpsMaster[numCpsMaster - 1].IsCoincident(cpsSlave[0])));
        //    }
        //    ///////////////////////
        //    if (flagContinue)
        //      continue;
        //    connection = new int[] { indexMaster, indexSlave };
        //    break;//// maximun 1 connect interface
        //  }
        //}
        if (Is1_1Constraint)
        {
          List<int[]> listTemp = surfaceMaster.GetCompletedMatchInterfaceGeometry(surfaceSlave);
          if (listTemp.Count > 0)
            connection = listTemp[0];
        }
        else
        {
          List<int[]> listTemp = surfaceMaster.GetCompletedCoincidentInterfaceGeometry(surfaceSlave);
          if (listTemp.Count > 0)
            connection = listTemp[0];
        }
      }
      else if (masterPatch is AbstractPatch3D && slavePatch is AbstractPatch3D)
      {
        AbstractPatch3D master3D = (AbstractPatch3D)masterPatch;
        AbstractPatch3D slave3D = (AbstractPatch3D)slavePatch;
        var volumeMaster = master3D.GetVolume();
        var volumeSlave = slave3D.GetVolume();
        //  for (int indexMaster = 0; indexMaster < 6; indexMaster++)
        //  {
        //    for (int indexSlave = 0; indexSlave < 6; indexSlave++)
        //    {
        //      if ((master3D != slave3D) || (indexMaster != indexSlave))
        //      {
        //        //if (Is1_1Constraint)
        //        //{
        //        int[] indexCoordinateMasterArea = master3D.GetCoordinateParameterOnArea(indexMaster);
        //        int[] indexCoordinateSlaveArea = slave3D.GetCoordinateParameterOnArea(indexSlave);

        //        KnotVector[] kvMaster = new KnotVector[2] { volumeMaster.Basis.GetKnotVector(indexCoordinateMasterArea[0])
        //          , volumeMaster.Basis.GetKnotVector(indexCoordinateMasterArea[1])};//0]) };

        //        KnotVector[] kvSlave = new KnotVector[2] { volumeSlave.Basis.GetKnotVector(indexCoordinateSlaveArea[0])
        //          , volumeSlave.Basis.GetKnotVector(indexCoordinateSlaveArea[1]) };

        //        int indexMasterCoincident = -1;
        //        int indexSlaveCoincident = -1;
        //        bool flagContinue = false;
        //        if (Is1_1Constraint)
        //        {
        //          if (!kvMaster[0].IsMatch(kvSlave[0]) && !kvMaster[0].IsMatch(kvSlave[1])
        //&& !kvMaster[1].IsMatch(kvSlave[0]) && !kvMaster[1].IsMatch(kvSlave[1]))
        //            continue;
        //        }
        //        ControlPoint[,] cpsMaster = master3D.SelectEndPatchControlPoints(indexMaster);
        //        ControlPoint[,] cpsSlave = slave3D.SelectEndPatchControlPoints(indexSlave);
        //        int[] numCpsMaster = new int[] { cpsMaster.GetLength(0), cpsMaster.GetLength(1) };
        //        int[] numCpsSlave = new int[] { cpsSlave.GetLength(0), cpsSlave.GetLength(1) };
        //        if (Is1_1Constraint)
        //        {

        //          if (numCpsMaster[0] == numCpsSlave[0])
        //          {
        //            if (numCpsMaster[1] == numCpsSlave[1])
        //            {
        //              indexMasterCoincident = 0;
        //              indexSlaveCoincident = 0;
        //            }
        //            else
        //              continue;
        //          }
        //          else
        //          {
        //            if (numCpsMaster[0] == numCpsSlave[1])
        //            {
        //              if (numCpsMaster[1] == numCpsSlave[0])
        //              {
        //                indexMasterCoincident = 0;
        //                indexSlaveCoincident = 1;
        //              }
        //              else
        //                continue;
        //            }
        //            else
        //              continue;
        //          }

        //          for (int j = 0; j < numCpsMaster[1]; j++)
        //            for (int i = 0; i < numCpsMaster[0]; i++)
        //            {
        //              if (numCpsMaster[0] == numCpsMaster[1])
        //              {
        //                if ((!cpsMaster[i, j].IsCoincident(cpsSlave[i, j]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - i, j])
        //                        && !cpsMaster[i, j].IsCoincident(cpsSlave[i, numCpsSlave[1] - 1 - j]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - i, numCpsSlave[1] - 1 - j]))
        //                && (!cpsMaster[i, j].IsCoincident(cpsSlave[j, i]) && !cpsMaster[i, j].IsCoincident(cpsSlave[j, numCpsSlave[1] - 1 - i])
        //                        && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - j, i]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - j, numCpsSlave[1] - 1 - i])))
        //                {
        //                  flagContinue = true;
        //                  break;
        //                }
        //              }
        //              else
        //              {
        //                if (numCpsMaster[0] == numCpsSlave[0] && numCpsMaster[1] == numCpsSlave[1])
        //                {
        //                  if ((!cpsMaster[i, j].IsCoincident(cpsSlave[i, j]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - i, j])
        //                      && !cpsMaster[i, j].IsCoincident(cpsSlave[i, numCpsSlave[1] - 1 - j]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - i, numCpsSlave[1] - 1 - j])))
        //                  {
        //                    flagContinue = true;
        //                    break;
        //                  }
        //                }
        //                else if (numCpsMaster[0] == numCpsSlave[1] && numCpsMaster[1] == numCpsSlave[0])
        //                {
        //                  if ((!cpsMaster[i, j].IsCoincident(cpsSlave[j, i]) && !cpsMaster[i, j].IsCoincident(cpsSlave[j, numCpsSlave[1] - 1 - i])
        //                          && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - j, i]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - j, numCpsSlave[1] - 1 - i])))
        //                  {
        //                    flagContinue = true;
        //                    break;
        //                  }
        //                }
        //                else
        //                {
        //                  flagContinue = true;
        //                  break;
        //                }
        //              }

        //              //if (numCpsMaster[0] == numCpsSlave[0] && numCpsMaster[1] == numCpsSlave[1])
        //              //{
        //              //    if (!cpsMaster[i, j].IsCoincident(cpsSlave[i, j]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - i, j])
        //              //            && !cpsMaster[i, j].IsCoincident(cpsSlave[i, numCpsSlave[1] - 1 - j]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - i, numCpsSlave[1] - 1 - j]))
        //              //    {
        //              //        flagContinue = true;
        //              //        break;
        //              //    }
        //              //}

        //              //if (numCpsMaster[0] == numCpsSlave[1] && numCpsMaster[1] == numCpsSlave[0] && numCpsMaster[0] != numCpsMaster[1])
        //              //{
        //              //    if (!cpsMaster[i, j].IsCoincident(cpsSlave[j, i]) && !cpsMaster[i, j].IsCoincident(cpsSlave[j, numCpsSlave[1] - 1 - i])
        //              //                && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - j, i]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - j, numCpsSlave[1] - 1 - i]))
        //              //    {
        //              //        flagContinue = true;
        //              //        break;
        //              //    }
        //              //}
        //            }
        //        }
        //        else
        //        {
        //          if (!cpsMaster[0, 0].IsCoincident(cpsSlave[0, 0]) || !cpsMaster[numCpsMaster[0] - 1, 0].IsCoincident(cpsSlave[numCpsSlave[0] - 1, 0])
        //            || !cpsMaster[0, numCpsMaster[1] - 1].IsCoincident(cpsSlave[0, numCpsSlave[1] - 1]))
        //            continue;
        //        }
        //        if (flagContinue)
        //          continue;
        //        connection = new int[] { indexMaster, indexSlave };
        //        break;//// maximun 1 connect interface
        //      }
        //    }
        //  }
        if (Is1_1Constraint)
        {
          List<int[]> listTemp = volumeMaster.GetCompletedMatchInterfaceGeometry(volumeSlave);
          if (listTemp.Count > 0)
            connection = listTemp[0];
        }
        else
        {
          List<int[]> listTemp = volumeMaster.GetCompletedCoincidentInterfaceGeometry(volumeSlave);
          if (listTemp.Count > 0)
            connection = listTemp[0];
        }
      }

      if (connection != null)
        if (masterPatch == slavePatch && connection[0] == connection[1])
          return null;
      return connection;
    }
  }
}
