using System;
using System.Collections.Generic;
using DEMSoft.NURBS;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public class CouplingTwoPatches : ICoupling
  {
    private AbstractPatch masterPatch;
    private AbstractPatch slavePatch;
    private int indexMasterObjectInterface;
    private int indexSlaveObjectInterface;
    //private DoubleMatrix AMaster;
    //private DoubleMatrix ASlave;
    private DoubleMatrix Tsm, Ts2m1, Ts2m2;
    internal DoubleMatrix ASlave1Plus;
    private Dictionary<int, List<int>> dicMaster, dicSlave, dicMaster2, dicSlave2;
    private bool is1_1Constraint;
    private List<int> ListDOFUnCoupling
    { get; set; }

    public int numberOfDecimal = 12;


    public CouplingTwoPatches(AbstractPatch masterPatch, AbstractPatch slavePatch, int indexMaster, int indexSlave, int[] DOFUnCoupling, bool is1_1Constraint, bool isUseC1Constraint, ref List<ControlPoint> listControlPointCoincidented, bool isMergeControlpoints1_1Constraint)
    {
      this.masterPatch = masterPatch;
      this.slavePatch = slavePatch;
      this.indexMasterObjectInterface = indexMaster;
      this.indexSlaveObjectInterface = indexSlave;
      this.is1_1Constraint = is1_1Constraint;
      ListDOFUnCoupling = new List<int>();
      if (DOFUnCoupling != null)
        for (int i = 0; i < DOFUnCoupling.Length; i++)
        {
          ListDOFUnCoupling.Add(DOFUnCoupling[i]);
        }
      if (masterPatch is AbstractPatch2D && slavePatch is AbstractPatch2D)
      {
        AbstractPatch2D master2D = (AbstractPatch2D)masterPatch;
        AbstractPatch2D slave2D = (AbstractPatch2D)slavePatch;
        var surfaceMaster = master2D.GetSurface();
        var surfaceSlave = slave2D.GetSurface();

        KnotVector kvMaster = surfaceMaster.Basis.GetKnotVector(indexMaster / 2);
        KnotVector kvSlave = surfaceSlave.Basis.GetKnotVector(indexSlave / 2);
        var cpsMaster = master2D.SelectEndPatchControlPoints(indexMaster);
        var cpsSlave = slave2D.SelectEndPatchControlPoints(indexSlave);
        //bool isReverse = false;
        if (!(cpsMaster[0].IsCoincident(cpsSlave[0]) && cpsMaster[cpsMaster.Length - 1].IsCoincident(cpsSlave[cpsSlave.Length - 1])))
        {
          throw new ArgumentException("Two patches must same orientation of interface edges, try to use reverse orientation of direction on geometry");
        }

        if (!is1_1Constraint)
        {
          ////// Đổi master thành slave nếu knotvector master nhiều hơn slave
          if (kvMaster.Count > kvSlave.Count)
          {
            this.masterPatch = slavePatch;
            this.slavePatch = masterPatch;
            this.indexMasterObjectInterface = indexSlave;
            this.indexSlaveObjectInterface = indexMaster;
            KnotVector kv = kvMaster;
            kvMaster = kvSlave;
            kvSlave = kv;
          }

          KnotVector kv12 = kvSlave.Substract(kvMaster);
          KnotVector kv21 = kvMaster.Substract(kvSlave);
          DoubleMatrix AMaster2 = null;
          DoubleMatrix ASlave2 = null;
          DoubleMatrix AMaster = ComputeAMatrix(true, isUseC1Constraint, out AMaster2, kv12.GetKnotVector().ToArray());
          DoubleMatrix ASlave = ComputeAMatrix(false, isUseC1Constraint, out ASlave2, kv21.GetKnotVector().ToArray());
          ASlave1Plus = MatrixFunctions.PseudoInverse(ASlave);
          Tsm = MatrixFunctions.Product(ASlave1Plus, AMaster);
          //for (int j = 0; j < Tsm.Cols; j++)
          //{
          //  for (int i = 0; i < Tsm.Rows; i++)
          //  {
          //    Tsm[i, j] = Math.Round(Tsm[i, j], numberOfDecimal);
          //  }
          //}
          DoubleMatrix Csm = null;

          DoubleMatrix ASlave2Plus = null;
          if (isUseC1Constraint)
          {
            ASlave2Plus = MatrixFunctions.PseudoInverse(ASlave2);
            Csm = ComputeCsmMatrix(AMaster, AMaster2, ASlave, ASlave2);
            Ts2m1 = MatrixFunctions.Product(ASlave2Plus, MatrixFunctions.Product(MatrixFunctions.Product(ASlave, ASlave1Plus) + Csm, AMaster));
            Ts2m2 = -MatrixFunctions.Product(ASlave2Plus, MatrixFunctions.Product(Csm, AMaster2));

            //for (int j = 0; j < Ts2m1.Cols; j++)
            //{
            //  for (int i = 0; i < Ts2m1.Rows; i++)
            //  {

            //    Ts2m1[i, j] = Math.Round(Ts2m1[i, j], numberOfDecimal);
            //  }
            //}

            //for (int j = 0; j < Ts2m2.Cols; j++)
            //{
            //  for (int i = 0; i < Ts2m2.Rows; i++)
            //  {

            //    Ts2m2[i, j] = Math.Round(Ts2m2[i, j], numberOfDecimal);
            //  }
            //}
          }

          if (isMergeControlpoints1_1Constraint)
          {
            for (int i = 0; i < cpsMaster.Length; i++)
            {
              if (i == 0 || i == cpsMaster.Length - 1)
              {
                bool isFoundCpInList = false;
                for (int j = 0; j < listControlPointCoincidented.Count; j++)
                {
                  if (cpsMaster[i].IsCoincident(listControlPointCoincidented[j]))
                  {
                    isFoundCpInList = true;
                    cpsMaster[i].SetCoupleControlPoint(listControlPointCoincidented[j]);
                    break;
                  }
                }
                if (!isFoundCpInList)
                {
                  listControlPointCoincidented.Add(cpsMaster[i].Clone());
                  cpsMaster[i].SetCoupleControlPoint(listControlPointCoincidented[listControlPointCoincidented.Count - 1]);
                }
                cpsMaster[i].ListDOFUnCoupling = ListDOFUnCoupling;
              }
            }

            for (int i = 0; i < cpsSlave.Length; i++)
            {
              if (i == 0 || i == cpsSlave.Length - 1)
              {
                for (int j = 0; j < listControlPointCoincidented.Count; j++)
                {
                  if (cpsSlave[i].IsCoincident(listControlPointCoincidented[j]))
                  {
                    cpsSlave[i].SetCoupleControlPoint(listControlPointCoincidented[j]);
                    break;
                  }
                }
                cpsSlave[i].ListDOFUnCoupling = ListDOFUnCoupling;
              }
            }
          }
        }
        else
        {
          if (!kvMaster.IsMatch(kvSlave))
            throw new ArgumentException("Use 1-1 coupling with the same knotvector in two patches coupling");
          //var cpsMaster = master2D.SelectEndPatchControlPoints(indexMaster);
          //var cpsSlave = slave2D.SelectEndPatchControlPoints(indexSlave);
          int numCpsMaster = cpsMaster.Length;
          int numCpsSlave = cpsSlave.Length;
          if (numCpsMaster != numCpsSlave)
            throw new ArgumentException("Number of control point on interface of two edges must equals");
          for (int i = 0; i < numCpsMaster; i++)
          {
            if ((!cpsMaster[i].IsCoincident(cpsSlave[i])) && (!cpsMaster[i].IsCoincident(cpsSlave[numCpsMaster - 1 - i])))
              throw new ArgumentException("Control points on interface of two edges must coincident");
          }

          if (!isMergeControlpoints1_1Constraint)
          {
            for (int i = 0; i < cpsSlave.Length; i++)
            {
              if (cpsMaster[i].IsCoincident(cpsSlave[i]))
                cpsSlave[i].SetCoupleControlPoint(cpsMaster[i]);
              else
                cpsSlave[i].SetCoupleControlPoint(cpsMaster[numCpsMaster - 1 - i]);
              cpsSlave[i].ListDOFUnCoupling = ListDOFUnCoupling;
            }
          }
          else
          {
            for (int i = 0; i < cpsMaster.Length; i++)
            {
              bool isFoundCpInList = false;
              for (int j = 0; j < listControlPointCoincidented.Count; j++)
              {
                if (cpsMaster[i].IsCoincident(listControlPointCoincidented[j]))
                {
                  isFoundCpInList = true;
                  cpsMaster[i].SetCoupleControlPoint(listControlPointCoincidented[j]);
                  break;
                }
              }
              if (!isFoundCpInList)
              {
                listControlPointCoincidented.Add(cpsMaster[i].Clone());
                cpsMaster[i].SetCoupleControlPoint(listControlPointCoincidented[listControlPointCoincidented.Count - 1]);
              }
              cpsMaster[i].ListDOFUnCoupling = ListDOFUnCoupling;
            }

            for (int i = 0; i < cpsSlave.Length; i++)
            {
              for (int j = 0; j < listControlPointCoincidented.Count; j++)
              {
                if (cpsSlave[i].IsCoincident(listControlPointCoincidented[j]))
                {
                  cpsSlave[i].SetCoupleControlPoint(listControlPointCoincidented[j]);
                  break;
                }
              }
              cpsSlave[i].ListDOFUnCoupling = ListDOFUnCoupling;
            }
          }
        }
      }
      else if (masterPatch is AbstractPatch3D && slavePatch is AbstractPatch3D)
      {
        AbstractPatch3D master3D = (AbstractPatch3D)masterPatch;
        AbstractPatch3D slave3D = (AbstractPatch3D)slavePatch;
        var volumeMaster = master3D.GetVolume();
        var volumeSlave = slave3D.GetVolume();
        int[] indexCoordinateMasterArea = master3D.GetCoordinateParameterOnArea(indexMaster);
        int[] indexCoordinateSlaveArea = slave3D.GetCoordinateParameterOnArea(indexSlave);

        KnotVector[] kvMaster = new KnotVector[2] { volumeMaster.Basis.GetKnotVector(indexCoordinateMasterArea[0])
           , volumeMaster.Basis.GetKnotVector(indexCoordinateMasterArea[1])};

        KnotVector[] kvSlave = new KnotVector[2] { volumeSlave.Basis.GetKnotVector(indexCoordinateSlaveArea[0])
           , volumeSlave.Basis.GetKnotVector(indexCoordinateSlaveArea[1]) };
        int indexMasterCoincident = -1;
        int indexSlaveCoincident = -1;

        if (!is1_1Constraint)
        {
          ControlPoint[,] cpsMaster = master3D.SelectEndPatchControlPoints(indexMaster);
          ControlPoint[,] cpsSlave = slave3D.SelectEndPatchControlPoints(indexSlave);
          bool isCheckSameOrientation = false;
          if (cpsMaster[0, 0].IsCoincident(cpsSlave[0, 0]))
          {
            if (cpsMaster[cpsMaster.GetLength(0) - 1, 0].IsCoincident(cpsSlave[cpsSlave.GetLength(0) - 1, 0]))
            {
              if (cpsMaster[0, cpsMaster.GetLength(1) - 1].IsCoincident(cpsSlave[0, cpsSlave.GetLength(1) - 1]))
              {
                isCheckSameOrientation = true;
              }
            }
          }

          if (!(isCheckSameOrientation))
            throw new ArgumentException("Two patches must same orientation of edge on interface faces, try to use reverse orientation of direction on geometry");

          if (kvMaster[0].Count >= kvSlave[0].Count || kvMaster[1].Count >= kvSlave[1].Count)
          {
            this.masterPatch = slavePatch;
            this.slavePatch = masterPatch;
            this.indexMasterObjectInterface = indexSlave;
            this.indexSlaveObjectInterface = indexMaster;
            KnotVector[] kv = kvMaster;
            kvMaster = kvSlave;
            kvSlave = kv;
          }
          KnotVector kv120 = kvSlave[0].Substract(kvMaster[0]);
          KnotVector kv121 = kvSlave[1].Substract(kvMaster[1]);
          KnotVector kv210 = kvMaster[0].Substract(kvSlave[0]);
          KnotVector kv211 = kvMaster[1].Substract(kvSlave[1]);

          DoubleMatrix AMaster2 = null;
          DoubleMatrix ASlave2 = null;
          DoubleMatrix AMaster = ComputeAMatrix(true, isUseC1Constraint, out AMaster2, kv120.GetKnotVector().ToArray(), kv121.GetKnotVector().ToArray());
          DoubleMatrix ASlave = ComputeAMatrix(false, isUseC1Constraint, out ASlave2, kv210.GetKnotVector().ToArray(), kv211.GetKnotVector().ToArray());
          ASlave1Plus = MatrixFunctions.PseudoInverse(ASlave);
          Tsm = MatrixFunctions.Product(ASlave1Plus, AMaster);

          DoubleMatrix Csm = null;

          DoubleMatrix ASlave2Plus = null;
          if (isUseC1Constraint)
          {
            ASlave2Plus = MatrixFunctions.PseudoInverse(ASlave2);
            Csm = ComputeCsmMatrix(AMaster, AMaster2, ASlave, ASlave2);
            Ts2m1 = MatrixFunctions.Product(ASlave2Plus, MatrixFunctions.Product(MatrixFunctions.Product(ASlave, ASlave1Plus) + Csm, AMaster));
            Ts2m2 = -MatrixFunctions.Product(ASlave2Plus, MatrixFunctions.Product(Csm, AMaster2));
          }

          if (isMergeControlpoints1_1Constraint)
          {
            for (int i = 0; i < cpsMaster.GetLength(0); i++)
            {
              for (int k = 0; k < cpsMaster.GetLength(1); k++)
              {
                if ((i == 0 || i == cpsMaster.GetLength(0) - 1) && (k == 0 || k == cpsMaster.GetLength(1) - 1))
                {
                  bool isFoundCpInList = false;
                  for (int j = 0; j < listControlPointCoincidented.Count; j++)
                  {
                    if (cpsMaster[i, k].IsCoincident(listControlPointCoincidented[j]))
                    {
                      isFoundCpInList = true;
                      cpsMaster[i, k].SetCoupleControlPoint(listControlPointCoincidented[j]);
                      break;
                    }
                  }
                  if (!isFoundCpInList)
                  {
                    listControlPointCoincidented.Add(cpsMaster[i, k].Clone());
                    cpsMaster[i, k].SetCoupleControlPoint(listControlPointCoincidented[listControlPointCoincidented.Count - 1]);
                  }
                  cpsMaster[i, k].ListDOFUnCoupling = ListDOFUnCoupling;
                }
              }
            }

            for (int i = 0; i < cpsSlave.GetLength(0); i++)
            {
              for (int k = 0; k < cpsSlave.GetLength(1); k++)
              {
                if ((i == 0 || i == cpsSlave.GetLength(0) - 1) && (k == 0 || k == cpsSlave.GetLength(1) - 1))
                {
                  for (int j = 0; j < listControlPointCoincidented.Count; j++)
                  {
                    if (cpsSlave[i, k].IsCoincident(listControlPointCoincidented[j]))
                    {
                      cpsSlave[i, k].SetCoupleControlPoint(listControlPointCoincidented[j]);
                      break;
                    }
                  }
                  cpsSlave[i, k].ListDOFUnCoupling = ListDOFUnCoupling;
                }
              }
            }
          }

        }
        else
        {
          if (!kvMaster[0].IsMatch(kvSlave[0]) && !kvMaster[0].IsMatch(kvSlave[1])
  && !kvMaster[1].IsMatch(kvSlave[0]) && !kvMaster[1].IsMatch(kvSlave[1]))
            throw new ArgumentException("Use 1-1 coupling with the same knotvector in two patches coupling");
          ControlPoint[,] cpsMaster = master3D.SelectEndPatchControlPoints(indexMaster);
          ControlPoint[,] cpsSlave = slave3D.SelectEndPatchControlPoints(indexSlave);
          int[] numCpsMaster = new int[] { cpsMaster.GetLength(0), cpsMaster.GetLength(1) };
          int[] numCpsSlave = new int[] { cpsSlave.GetLength(0), cpsSlave.GetLength(1) };

          if (numCpsMaster[0] == numCpsSlave[0])
          {
            if (numCpsMaster[1] == numCpsSlave[1])
            {
              indexMasterCoincident = 0;
              indexSlaveCoincident = 0;
            }
            else
              throw new ArgumentException("Number of control point on interface of two areas must equals");
          }
          else
          {
            if (numCpsMaster[0] == numCpsSlave[1])
            {
              if (numCpsMaster[1] == numCpsSlave[0])
              {
                indexMasterCoincident = 0;
                indexSlaveCoincident = 1;
              }
              else
                throw new ArgumentException("Number of control point on interface of two areas must equals");
            }
            else
              throw new ArgumentException("Number of control point on interface of two areas must equals");
          }

          for (int j = 0; j < numCpsMaster[1]; j++)
            for (int i = 0; i < numCpsMaster[0]; i++)
            {
              if (numCpsMaster[0] == numCpsMaster[1])
              {
                if ((!cpsMaster[i, j].IsCoincident(cpsSlave[i, j]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - i, j])
                      && !cpsMaster[i, j].IsCoincident(cpsSlave[i, numCpsSlave[1] - 1 - j]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - i, numCpsSlave[1] - 1 - j]))
                && (!cpsMaster[i, j].IsCoincident(cpsSlave[j, i]) && !cpsMaster[i, j].IsCoincident(cpsSlave[j, numCpsSlave[1] - 1 - i])
                      && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - j, i]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - j, numCpsSlave[1] - 1 - i])))
                  throw new ArgumentException("Control points on interface of two areas must coincident");
              }
              else
              {
                if (numCpsMaster[0] == numCpsSlave[0] && numCpsMaster[1] == numCpsSlave[1])
                {
                  if ((!cpsMaster[i, j].IsCoincident(cpsSlave[i, j]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - i, j])
                             && !cpsMaster[i, j].IsCoincident(cpsSlave[i, numCpsSlave[1] - 1 - j]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - i, numCpsSlave[1] - 1 - j])))
                    throw new ArgumentException("Control points on interface of two areas must coincident");
                }
                else if (numCpsMaster[0] == numCpsSlave[1] && numCpsMaster[1] == numCpsSlave[0])
                {
                  if ((!cpsMaster[i, j].IsCoincident(cpsSlave[j, i]) && !cpsMaster[i, j].IsCoincident(cpsSlave[j, numCpsSlave[1] - 1 - i])
                                && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - j, i]) && !cpsMaster[i, j].IsCoincident(cpsSlave[numCpsSlave[0] - 1 - j, numCpsSlave[1] - 1 - i])))
                    throw new ArgumentException("Control points on interface of two areas must coincident");
                }
                else
                  throw new ArgumentException("Control points on interface of two areas must coincident");
              }
            }
          if (!isMergeControlpoints1_1Constraint)
          {
            for (int j = 0; j < numCpsSlave[1]; j++)
              for (int i = 0; i < numCpsSlave[0]; i++)
              {
                if (numCpsMaster[0] == numCpsMaster[1])
                {
                  if (cpsSlave[i, j].IsCoincident(cpsMaster[i, j]))
                    cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[i, j]);
                  else if (cpsSlave[i, j].IsCoincident(cpsMaster[numCpsMaster[0] - 1 - i, j]))
                    cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[numCpsMaster[0] - 1 - i, j]);
                  else if (cpsSlave[i, j].IsCoincident(cpsMaster[i, numCpsMaster[1] - 1 - j]))
                    cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[i, numCpsMaster[1] - 1 - j]);
                  else if (cpsSlave[i, j].IsCoincident(cpsMaster[numCpsMaster[0] - 1 - i, numCpsMaster[1] - 1 - j]))
                    cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[numCpsMaster[0] - 1 - i, numCpsMaster[1] - 1 - j]);
                  else if (cpsSlave[i, j].IsCoincident(cpsMaster[j, i]))
                    cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[j, i]);
                  else if (cpsSlave[i, j].IsCoincident(cpsMaster[numCpsMaster[0] - 1 - j, i]))
                    cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[numCpsMaster[0] - 1 - j, i]);
                  else if (cpsSlave[i, j].IsCoincident(cpsMaster[j, numCpsMaster[1] - 1 - i]))
                    cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[j, numCpsMaster[1] - 1 - i]);
                  else if (cpsSlave[i, j].IsCoincident(cpsMaster[numCpsMaster[0] - 1 - j, numCpsMaster[1] - 1 - i]))
                    cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[numCpsMaster[0] - 1 - j, numCpsMaster[1] - 1 - i]);
                }
                else
                {
                  if (numCpsMaster[0] == numCpsSlave[0] && numCpsMaster[1] == numCpsSlave[1])
                  {
                    if (cpsSlave[i, j].IsCoincident(cpsMaster[i, j]))
                      cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[i, j]);
                    else if (cpsSlave[i, j].IsCoincident(cpsMaster[numCpsMaster[0] - 1 - i, j]))
                      cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[numCpsMaster[0] - 1 - i, j]);
                    else if (cpsSlave[i, j].IsCoincident(cpsMaster[i, numCpsMaster[1] - 1 - j]))
                      cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[i, numCpsMaster[1] - 1 - j]);
                    else if (cpsSlave[i, j].IsCoincident(cpsMaster[numCpsMaster[0] - 1 - i, numCpsMaster[1] - 1 - j]))
                      cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[numCpsMaster[0] - 1 - i, numCpsMaster[1] - 1 - j]);
                  }
                  else if (numCpsMaster[0] == numCpsSlave[1] && numCpsMaster[1] == numCpsSlave[0])
                  {
                    if (cpsSlave[i, j].IsCoincident(cpsMaster[j, i]))
                      cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[j, i]);
                    else if (cpsSlave[i, j].IsCoincident(cpsMaster[numCpsMaster[0] - 1 - j, i]))
                      cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[numCpsMaster[0] - 1 - j, i]);
                    else if (cpsSlave[i, j].IsCoincident(cpsMaster[j, numCpsMaster[1] - 1 - i]))
                      cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[j, numCpsMaster[1] - 1 - i]);
                    else if (cpsSlave[i, j].IsCoincident(cpsMaster[numCpsMaster[0] - 1 - j, numCpsMaster[1] - 1 - i]))
                      cpsSlave[i, j].SetCoupleControlPoint(cpsMaster[numCpsMaster[0] - 1 - j, numCpsMaster[1] - 1 - i]);
                  }
                }
                cpsSlave[i, j].ListDOFUnCoupling = ListDOFUnCoupling;
              }
          }
          else
          {
            for (int i = 0; i < cpsMaster.GetLength(0); i++)
            {
              for (int k = 0; k < cpsMaster.GetLength(1); k++)
              {
                bool isFoundCpInList = false;
                for (int j = 0; j < listControlPointCoincidented.Count; j++)
                {
                  if (cpsMaster[i, k].IsCoincident(listControlPointCoincidented[j]))
                  {
                    isFoundCpInList = true;
                    cpsMaster[i, k].SetCoupleControlPoint(listControlPointCoincidented[j]);
                    break;
                  }
                }
                if (!isFoundCpInList)
                {
                  listControlPointCoincidented.Add(cpsMaster[i, k].Clone());
                  cpsMaster[i, k].SetCoupleControlPoint(listControlPointCoincidented[listControlPointCoincidented.Count - 1]);
                }
                cpsMaster[i, k].ListDOFUnCoupling = ListDOFUnCoupling;
              }
            }

            for (int i = 0; i < cpsSlave.GetLength(0); i++)
            {
              for (int k = 0; k < cpsSlave.GetLength(1); k++)
              {
                for (int j = 0; j < listControlPointCoincidented.Count; j++)
                {
                  if (cpsSlave[i, k].IsCoincident(listControlPointCoincidented[j]))
                  {
                    cpsSlave[i, k].SetCoupleControlPoint(listControlPointCoincidented[j]);
                    break;
                  }
                }
                cpsSlave[i, k].ListDOFUnCoupling = ListDOFUnCoupling;
              }
            }
          }
        }
      }
      else
        throw new ArgumentException("Master patch and slave patch must same dimension");
    }

    private DoubleMatrix ComputeBMatrix(bool isMaster, double[] xi)
    {
      Abstract2DParametricGeometry surfaceMaster = null;
      int indexCoordinationInterface;
      if (isMaster)
      {
        surfaceMaster = ((AbstractPatch2D)masterPatch).GetSurface();
        indexCoordinationInterface = indexMasterObjectInterface / 2;
      }
      else
      {
        surfaceMaster = ((AbstractPatch2D)slavePatch).GetSurface();
        indexCoordinationInterface = indexSlaveObjectInterface / 2;
      }
      KnotVector kvMaster = surfaceMaster.Basis.GetKnotVector(indexCoordinationInterface).Clone();
      int p = surfaceMaster.Basis.GetDegree(indexCoordinationInterface);
      int numCpsMaster = kvMaster.Count - p - 1;
      int lengthXi = xi.Length;
      if (lengthXi > 0)
      {
        DoubleMatrix[] Tij = new DoubleMatrix[lengthXi];

        for (int jalpha = 0; jalpha < lengthXi; jalpha++)
        {
          Tij[jalpha] = new DoubleMatrix(numCpsMaster + 1 + jalpha, numCpsMaster + jalpha);
          int k = kvMaster.FindK(xi[jalpha]);
          for (int i = 0; i < numCpsMaster + 1 + jalpha; i++)
          {
            double lamda = 0;
            if (0 <= i && i <= k - p)
              lamda = 1.0;
            else if (k - p + 1 <= i && i <= k)
              lamda = (xi[jalpha] - kvMaster[i]) / (kvMaster[i + p] - kvMaster[i]);
            else if (i >= k + 1)
              lamda = 0.0;
            if (i >= 0 && i < numCpsMaster + 1 + jalpha - 1)
              Tij[jalpha][i, i] = lamda;
            if (i > 0 && i < numCpsMaster + 1 + jalpha)
              Tij[jalpha][i, i - 1] = 1.0 - lamda;
          }
          kvMaster.InsertKnot(k + 1, xi[jalpha]);
        }
        DoubleMatrix B = Tij[0];
        for (int i = 1; i < Tij.Length; i++)
        {
          B = MatrixFunctions.Product(Tij[i], B);
        }
        //ControlPoint[] cp = surfaceMaster.SelectControlPointInTwoEndLines(indexSlaveObjectInterface);
        //DoubleVector v1 = new DoubleVector(cp.Length);
        //for (int i = 0; i < cp.Length; i++)
        //{
        //  v1[i] = cp[i][1];
        //}
        //DoubleVector v2 = MatrixFunctions.Product(B, v1);
        return B;
      }
      else
      {
        DoubleMatrix B = new DoubleMatrix(numCpsMaster, numCpsMaster);
        for (int i = 0; i < numCpsMaster; i++)
        {
          B[i, i] = 1.0;
        }
        return B;
      }
    }

    private DoubleMatrix ComputeBMatrix(bool isMaster, double[] xi, double[] eta)
    {
      AbstractPatch3D master3D = null;
      int indexInterface;
      if (isMaster)
      {
        master3D = (AbstractPatch3D)masterPatch;
        indexInterface = indexMasterObjectInterface;
      }
      else
      {
        master3D = (AbstractPatch3D)slavePatch;
        indexInterface = indexSlaveObjectInterface;
      }
      int[] indexCoordinateMasterLine = master3D.GetCoordinateParameterOnArea(indexInterface);
      Abstract3DParametricGeometry volumeMaster = master3D.GetVolume();
      KnotVector[] kvMaster = {
        volumeMaster.Basis.GetKnotVector(indexCoordinateMasterLine[0]).Clone(),
        volumeMaster.Basis.GetKnotVector(indexCoordinateMasterLine[1]).Clone() };
      int[] p = {
        volumeMaster.Basis.GetDegree(indexCoordinateMasterLine[0]),
        volumeMaster.Basis.GetDegree(indexCoordinateMasterLine[1]) };
      int[] numCpsMaster = { kvMaster[0].Count - p[0] - 1, kvMaster[1].Count - p[1] - 1 };
      double[][] xii = { xi, eta };
      int[] lengthXi = { xii[0].Length, xii[1].Length };
      DoubleMatrix[] B = new DoubleMatrix[2];
      if (lengthXi[0] > 0 || lengthXi[1] > 0)
      {
        DoubleMatrix[][] Tij = new DoubleMatrix[2][];
        for (int i = 0; i < 2; i++)
        {
          int lengthXii = lengthXi[i];
          if (lengthXii != 0)
          {
            Tij[i] = lengthXii > 0 ? new DoubleMatrix[lengthXii] : null;
            DoubleMatrix[] Tiji = Tij[i];
            double[] xiii = xii[i];
            KnotVector kvMasteri = kvMaster[i];
            for (int jalpha = 0; jalpha < lengthXii; jalpha++)
            {
              int numCpsMasteri = numCpsMaster[i];
              Tiji[jalpha] = new DoubleMatrix(numCpsMasteri + 1 + jalpha, numCpsMasteri + jalpha);
              DoubleMatrix TijiAlpha = Tiji[jalpha];
              int k = kvMasteri.FindK(xiii[jalpha]);
              for (int j = 0; j < numCpsMasteri + 1 + jalpha; j++)
              {
                double lamda = 0;
                if (0 <= j && j <= k - p[i])
                  lamda = 1.0;
                else if (k - p[i] + 1 <= j && j <= k)
                  lamda = (xiii[jalpha] - kvMasteri[j]) / (kvMasteri[j + p[i]] - kvMasteri[j]);
                else if (j >= k + 1)
                  lamda = 0.0;
                if (j >= 0 && j < numCpsMasteri + 1 + jalpha - 1)
                  TijiAlpha[j, j] = lamda;
                if (j > 0 && j < numCpsMasteri + 1 + jalpha)
                  TijiAlpha[j, j - 1] = 1.0 - lamda;
              }
              kvMaster[i].InsertKnot(k + 1, xiii[jalpha]);
            }
            B[i] = Tiji[0];
            for (int ii = 1; ii < Tiji.Length; ii++)
            {
              B[i] = MatrixFunctions.Product(Tiji[ii], B[i]);
            }
          }
        }
      }
      if (lengthXi[0] <= 0)
      {
        B[0] = new DoubleMatrix(numCpsMaster[0], numCpsMaster[0]);
        for (int i = 0; i < numCpsMaster[0]; i++)
        {
          B[0][i, i] = 1.0;
        }
      }

      if (lengthXi[1] <= 0)
      {
        B[1] = new DoubleMatrix(numCpsMaster[1], numCpsMaster[1]);
        for (int i = 0; i < numCpsMaster[1]; i++)
        {
          B[1][i, i] = 1.0;
        }
      }
      //return B;
      //Stopwatch sw = new Stopwatch();
      //sw.Start();
      DoubleMatrix BB = MatrixTool.OuterProduct(B[0], B[1]);
      //sw.Stop();
      //double time1 = sw.Elapsed.TotalMilliseconds;
      //sw.Reset();
      //sw.Start();
      //DoubleMatrix BB2 = MatrixTool.KroneckerProduct(B[0], B[1]);
      //sw.Stop();
      //double time2 = sw.Elapsed.TotalMilliseconds;
      return BB;
    }
    private DoubleMatrix ComputeCsmMatrix(DoubleMatrix AMaster1, DoubleMatrix AMaster2, DoubleMatrix ASlave1, DoubleMatrix ASlave2)
    {
      //AbstractPatch2D master2D = (AbstractPatch2D)masterPatch;
      //AbstractPatch2D slave2D = (AbstractPatch2D)slavePatch;
      //ControlPoint[] cpsMaster1 = master2D.SelectEndPatchControlPoints(indexMasterObjectInterface);
      //ControlPoint[] cpsMaster2 = master2D.SelectNearEndPatchControlPoints(indexMasterObjectInterface);
      //ControlPoint[] cpsSlave1 = slave2D.SelectEndPatchControlPoints(indexSlaveObjectInterface);
      //ControlPoint[] cpsSlave2 = slave2D.SelectNearEndPatchControlPoints(indexSlaveObjectInterface);
      //double[][] coordMaster1 = new double[3][];
      //double[][] coordMaster2 = new double[3][];
      //double[][] coordSlave1 = new double[3][];
      //double[][] coordSlave2 = new double[3][];
      //DoubleVector[] coordMaster1v = new DoubleVector[3];
      //DoubleVector[] coordMaster2v = new DoubleVector[3];
      //DoubleVector[] coordSlave1v = new DoubleVector[3];
      //DoubleVector[] coordSlave2v = new DoubleVector[3];

      //double[] distanceMaster = new double[AMaster1.Rows];
      //double[] distanceSlave = new double[ASlave1.Rows];
      //for (int i = 0; i < 3; i++)
      //{
      //  coordMaster1[i] = new double[cpsMaster1.Length];
      //  coordMaster2[i] = new double[cpsMaster2.Length];
      //  for (int j = 0; j < cpsMaster1.Length; j++)
      //  {
      //    coordMaster1[i][j] = cpsMaster1[j][i];
      //    coordMaster2[i][j] = cpsMaster2[j][i];
      //  }
      //  coordMaster1v[i] = MatrixFunctions.Product(AMaster1, new DoubleVector(coordMaster1[i]));
      //  coordMaster2v[i] = MatrixFunctions.Product(AMaster2, new DoubleVector(coordMaster2[i]));
      //  for (int j = 0; j < coordMaster1v[i].Length; j++)
      //  {
      //    distanceMaster[j] += Math.Pow(coordMaster1v[i][j] - coordMaster2v[i][j], 2);
      //  }

      //  coordSlave1[i] = new double[cpsSlave1.Length];
      //  coordSlave2[i] = new double[cpsSlave2.Length];
      //  for (int j = 0; j < cpsSlave1.Length; j++)
      //  {
      //    coordSlave1[i][j] = cpsSlave1[j][i];
      //    coordSlave2[i][j] = cpsSlave2[j][i];
      //  }
      //  coordSlave1v[i] = MatrixFunctions.Product(ASlave1, new DoubleVector(coordSlave1[i]));
      //  coordSlave2v[i] = MatrixFunctions.Product(ASlave2, new DoubleVector(coordSlave2[i]));
      //  for (int j = 0; j < coordSlave1v[i].Length; j++)
      //  {
      //    distanceSlave[j] += Math.Pow(coordSlave1v[i][j] - coordSlave2v[i][j], 2);
      //  }
      //}
      //DoubleMatrix Csm = new DoubleMatrix(AMaster1.Rows, AMaster1.Rows);
      //for (int i = 0; i < distanceMaster.Length; i++)
      //{
      //  distanceMaster[i] = Math.Sqrt(distanceMaster[i]);
      //  distanceSlave[i] = Math.Sqrt(distanceSlave[i]);
      //  Csm[i, i] = distanceSlave[i] / distanceMaster[i];
      //}
      DoubleMatrix Csm = null;
      int p = -1;
      KnotVector kvMaster = null;
      KnotVector kvSlave = null;
      int countCpsMaster = 0;
      int countCpsSlave = 0;
      if (masterPatch is AbstractPatch2D)
      {
        AbstractPatch2D master2D = (AbstractPatch2D)masterPatch;
        AbstractPatch2D slave2D = (AbstractPatch2D)slavePatch;
        Abstract2DParametricGeometry masterSurface = master2D.GetSurface();
        Abstract2DParametricGeometry slaveSurface = slave2D.GetSurface();
        int indexMasterCoordinationInterface = indexMasterObjectInterface / 2;
        int indexSlaveCoordinationInterface = indexSlaveObjectInterface / 2;
        //int indexFirstLastMasterInterface = indexMasterObjectInterface % 2;
        //int indexFirstLastSlaveInterface = indexSlaveObjectInterface % 2;
        p = masterSurface.Basis.GetDegree((indexMasterCoordinationInterface == 0) ? 1 : 0);
        kvMaster = masterSurface.Basis.GetKnotVector((indexMasterCoordinationInterface == 0) ? 1 : 0);
        kvSlave = slaveSurface.Basis.GetKnotVector((indexSlaveCoordinationInterface == 0) ? 1 : 0);
        countCpsMaster = masterSurface.Basis.GetCountBasisFunction((indexMasterCoordinationInterface == 0) ? 1 : 0);
        countCpsSlave = slaveSurface.Basis.GetCountBasisFunction((indexSlaveCoordinationInterface == 0) ? 1 : 0);
        //double coefMaster, coefSlave;
        //if (indexFirstLastMasterInterface == 0)
        //{
        //  coefMaster = p / kvMaster[p + 2 - 1];
        //}
        //else
        //{
        //  coefMaster = p / (1.0 - kvMaster[countCpsMaster - 1]);
        //}

        //if (indexFirstLastSlaveInterface == 0)
        //{
        //  coefSlave = p / kvSlave[p + 2 - 1];
        //}
        //else
        //{
        //  coefSlave = p / (1.0 - kvSlave[countCpsSlave - 1]);
        //}
        //Csm = new DoubleMatrix(AMaster1.Rows, AMaster1.Rows);
        //for (int i = 0; i < AMaster1.Rows; i++)
        //{
        //  Csm[i, i] = coefMaster / coefSlave;
        //}
      }
      else if (masterPatch is AbstractPatch3D)
      {
        AbstractPatch3D master3D = (AbstractPatch3D)masterPatch;
        AbstractPatch3D slave3D = (AbstractPatch3D)slavePatch;
        Abstract3DParametricGeometry masterVolume = master3D.GetVolume();
        Abstract3DParametricGeometry slaveVolume = slave3D.GetVolume();
        int[] indexCoordinateMasterSurface = master3D.GetCoordinateParameterOnArea(indexMasterObjectInterface);
        int[] indexCoordinateSlaveSurface = slave3D.GetCoordinateParameterOnArea(indexSlaveObjectInterface);
        int masterNormalDirectionIndex = NormalDirectionIndex(indexCoordinateMasterSurface[0], indexCoordinateMasterSurface[1]);
        int slaveNormalDirectionIndex = NormalDirectionIndex(indexCoordinateSlaveSurface[0], indexCoordinateSlaveSurface[1]);
        p = masterVolume.Basis.GetDegree(masterNormalDirectionIndex);
        kvMaster = masterVolume.Basis.GetKnotVector(masterNormalDirectionIndex);
        kvSlave = slaveVolume.Basis.GetKnotVector(slaveNormalDirectionIndex);
        countCpsMaster = masterVolume.Basis.GetCountBasisFunction(masterNormalDirectionIndex);
        countCpsSlave = slaveVolume.Basis.GetCountBasisFunction(slaveNormalDirectionIndex);
      }
      int indexFirstLastMasterInterface = indexMasterObjectInterface % 2;
      int indexFirstLastSlaveInterface = indexSlaveObjectInterface % 2;
      double coefMaster, coefSlave;
      if (indexFirstLastMasterInterface == 0)
      {
        coefMaster = p / kvMaster[p + 2 - 1];
      }
      else
      {
        coefMaster = p / (1.0 - kvMaster[countCpsMaster - 1]);
      }

      if (indexFirstLastSlaveInterface == 0)
      {
        coefSlave = p / kvSlave[p + 2 - 1];
      }
      else
      {
        coefSlave = p / (1.0 - kvSlave[countCpsSlave - 1]);
      }
      Csm = new DoubleMatrix(AMaster1.Rows, AMaster1.Rows);
      for (int i = 0; i < AMaster1.Rows; i++)
      {
        Csm[i, i] = (coefMaster / coefSlave);
      }
      return Csm;
    }
    private int NormalDirectionIndex(int xi, int eta)
    {
      int n = -1;
      if ((xi == 0 && eta == 1) || (xi == 1 && eta == 0))
        n = 2;
      else if ((xi == 0 && eta == 2) || (xi == 2 && eta == 0))
        n = 1;
      else if ((xi == 1 && eta == 2) || (xi == 2 && eta == 1))
        n = 0;
      return n;
    }
    private DoubleMatrix ComputeAMatrix(bool isMaster, bool isUseC1Constraint, out DoubleMatrix AMatrix2, double[] xi)
    {
      AbstractPatch2D master2D = null;
      int indexInterface;
      if (isMaster)
      {
        master2D = (AbstractPatch2D)masterPatch;
        indexInterface = indexMasterObjectInterface;

      }
      else
      {
        master2D = (AbstractPatch2D)slavePatch;
        indexInterface = indexSlaveObjectInterface;
      }
      int indexCoordinateMasterLine = indexInterface / 2;
      Abstract2DParametricGeometry surfaceMaster = master2D.GetSurface();
      KnotVector kvMaster = surfaceMaster.Basis.GetKnotVector(indexCoordinateMasterLine);
      int p = surfaceMaster.Basis.GetDegree(indexCoordinateMasterLine);
      ControlPoint[] cpsMaster = master2D.SelectEndPatchControlPoints(indexInterface);
      int numCpsMaster = kvMaster.Count - p - 1;
      int numCpsMasterNew = numCpsMaster + xi.Length;
      DoubleMatrix B = ComputeBMatrix(isMaster, xi);
      DoubleVector Wvec = new DoubleVector(numCpsMaster);
      DoubleMatrix W = new DoubleMatrix(numCpsMaster, numCpsMaster);
      ControlPoint[] cpsMaster2 = null;
      DoubleVector Wvec2 = null;
      DoubleMatrix W2 = null;
      if (isUseC1Constraint)
      {
        cpsMaster2 = master2D.SelectNearEndPatchControlPoints(indexInterface);
        Wvec2 = new DoubleVector(numCpsMaster);
        W2 = new DoubleMatrix(numCpsMaster, numCpsMaster);
      }
      for (int i = 0; i < numCpsMaster; i++)
      {
        Wvec[i] = cpsMaster[i][3];
        W[i, i] = cpsMaster[i][3];//weight
        if (isUseC1Constraint)
        {
          Wvec2[i] = cpsMaster2[i][3];
          W2[i, i] = cpsMaster2[i][3];//weight
        }
      }
      DoubleVector WaVec = MatrixFunctions.Product(B, Wvec);
      DoubleMatrix Wa = new DoubleMatrix(numCpsMasterNew, numCpsMasterNew);
      DoubleVector WaVec2 = null;
      DoubleMatrix Wa2 = null;
      if (isUseC1Constraint)
      {
        WaVec2 = MatrixFunctions.Product(B, Wvec2);
        Wa2 = new DoubleMatrix(numCpsMasterNew, numCpsMasterNew);
      }
      for (int i = 0; i < numCpsMasterNew; i++)
      {
        Wa[i, i] = WaVec[i];//weight alpha
        if (isUseC1Constraint)
        {
          Wa2[i, i] = WaVec2[i];//weight alpha
        }
      }
      AMatrix2 = null;
      if (isUseC1Constraint)
      {
        AMatrix2 = MatrixFunctions.Product(MatrixFunctions.Inverse(Wa2), MatrixFunctions.Product(B, W2));
      }
      return MatrixFunctions.Product(MatrixFunctions.Inverse(Wa), MatrixFunctions.Product(B, W));
    }

    private DoubleMatrix ComputeAMatrix(bool isMaster, bool isUseC1Constraint, out DoubleMatrix AMatrix2, double[] xi, double[] eta)
    {
      AbstractPatch3D master3D = null;
      int indexInterface;
      if (isMaster)
      {
        master3D = (AbstractPatch3D)masterPatch;
        indexInterface = indexMasterObjectInterface;
      }
      else
      {
        master3D = (AbstractPatch3D)slavePatch;
        indexInterface = indexSlaveObjectInterface;
      }
      int[] indexCoordinateMasterLine = master3D.GetCoordinateParameterOnArea(indexInterface);
      Abstract3DParametricGeometry volumeMaster = master3D.GetVolume();
      KnotVector[] kvMaster = {volumeMaster.Basis.GetKnotVector(indexCoordinateMasterLine[0]),
        volumeMaster.Basis.GetKnotVector(indexCoordinateMasterLine[1])};

      int[] p = { volumeMaster.Basis.GetDegree(indexCoordinateMasterLine[0]), volumeMaster.Basis.GetDegree(indexCoordinateMasterLine[1]) };
      ControlPoint[,] cpsMaster = master3D.SelectEndPatchControlPoints(indexInterface);
      ControlPoint[,] cpsMaster2 = null;
      int[] numCpsMaster = { kvMaster[0].Count - p[0] - 1, kvMaster[1].Count - p[1] - 1 };
      int[] numCpsMasterNew = { numCpsMaster[0] + xi.Length, numCpsMaster[1] + eta.Length };
      DoubleMatrix B = ComputeBMatrix(isMaster, xi, eta);
      DoubleVector Wvec = new DoubleVector(numCpsMaster[0] * numCpsMaster[1]);
      DoubleMatrix W = new DoubleMatrix(numCpsMaster[0] * numCpsMaster[1], numCpsMaster[0] * numCpsMaster[1]);
      DoubleVector Wvec2 = null;
      DoubleMatrix W2 = null;
      if (isUseC1Constraint)
      {
        cpsMaster2 = master3D.SelectNearEndPatchControlPoints(indexInterface);
        Wvec2 = new DoubleVector(numCpsMaster[0] * numCpsMaster[1]);
        W2 = new DoubleMatrix(numCpsMaster[0] * numCpsMaster[1], numCpsMaster[0] * numCpsMaster[1]);
      }
      int c = 0;
      for (int i = 0; i < numCpsMaster[0]; i++)
      {
        for (int j = 0; j < numCpsMaster[1]; j++)
        {
          Wvec[c] = cpsMaster[i, j][3];
          W[c, c] = cpsMaster[i, j][3];//weight
          if (isUseC1Constraint)
          {
            Wvec2[c] = cpsMaster2[i, j][3];
            W2[c, c] = cpsMaster2[i, j][3];//weight
          }
          c++;
        }
      }
      DoubleVector WaVec = MatrixFunctions.Product(B, Wvec);
      DoubleMatrix Wa = new DoubleMatrix(numCpsMasterNew[0] * numCpsMasterNew[1], numCpsMasterNew[0] * numCpsMasterNew[1]);
      DoubleVector WaVec2 = null;
      DoubleMatrix Wa2 = null;
      if (isUseC1Constraint)
      {
        WaVec2 = MatrixFunctions.Product(B, Wvec2);
        Wa2 = new DoubleMatrix(numCpsMasterNew[0] * numCpsMasterNew[1], numCpsMasterNew[0] * numCpsMasterNew[1]);
      }
      c = 0;
      for (int i = 0; i < numCpsMasterNew[0]; i++)
      {
        for (int j = 0; j < numCpsMasterNew[1]; j++)
        {
          Wa[c, c] = WaVec[c];//weight alpha
          if (isUseC1Constraint)
          {
            Wa2[c, c] = WaVec2[c];//weight alpha
          }
          c++;
        }
      }
      AMatrix2 = null;
      if (isUseC1Constraint)
      {
        AMatrix2 = MatrixFunctions.Product(MatrixFunctions.Inverse(Wa2), MatrixFunctions.Product(B, W2));
      }
      return MatrixFunctions.Product(MatrixFunctions.Inverse(Wa), MatrixFunctions.Product(B, W));
    }

    public int[] GetTArrayOnLineMasterSlave(bool isMaster)
    {
      Dictionary<int, List<int>> dicMaster = new Dictionary<int, List<int>>();

      AbstractPatch2D patchMaster = null;
      int indexInterface;
      if (isMaster)
      {
        patchMaster = ((AbstractPatch2D)masterPatch);
        indexInterface = indexMasterObjectInterface;
      }
      else
      {
        patchMaster = ((AbstractPatch2D)slavePatch);
        indexInterface = indexSlaveObjectInterface;
      }
      Abstract2DParametricGeometry surfaceMaster = patchMaster.GetSurface();
      int countField = patchMaster.GetCountField(0);
      int indexCoordinateMasterLine = indexInterface / 2;
      KnotVector kvMaster = surfaceMaster.Basis.GetKnotVector(indexCoordinateMasterLine).Clone();
      int p = surfaceMaster.Basis.GetDegree(indexCoordinateMasterLine);
      int numCpsMaster = kvMaster.Count - p - 1;

      int[] tArray = new int[(countField - ListDOFUnCoupling.Count) * numCpsMaster];

      var cps = patchMaster.SelectEndPatchControlPoints(indexInterface);
      int c = 0;
      for (int i = 0; i < cps.Length; i++)
      {
        List<int> list = new List<int>();
        for (int d = 0; d < countField; d++)
        {
          if (ListDOFUnCoupling.IndexOf(d) == -1)
          {
            tArray[c++] = cps[i].GetTArrayGlobal()[d];
            list.Add(cps[i].GetTArrayGlobal()[d]);
          }
        }
        dicMaster.Add(i, list);
      }
      if (isMaster)
        this.dicMaster = dicMaster;
      else
        this.dicSlave = dicMaster;
      return tArray;
    }

    public int[] GetTArrayOnLineMasterSlave2(bool isMaster)
    {
      Dictionary<int, List<int>> dicMaster = new Dictionary<int, List<int>>();

      AbstractPatch2D patchMaster = null;
      int indexInterface;
      if (isMaster)
      {
        patchMaster = ((AbstractPatch2D)masterPatch);
        indexInterface = indexMasterObjectInterface;
      }
      else
      {
        patchMaster = ((AbstractPatch2D)slavePatch);
        indexInterface = indexSlaveObjectInterface;
      }
      Abstract2DParametricGeometry surfaceMaster = patchMaster.GetSurface();
      int countDimension = patchMaster.GetCountDimension();
      int countField = patchMaster.GetCountField(0);
      int indexCoordinateMasterLine = indexInterface / 2;
      KnotVector kvMaster = surfaceMaster.Basis.GetKnotVector(indexCoordinateMasterLine).Clone();
      int p = surfaceMaster.Basis.GetDegree(indexCoordinateMasterLine);
      int numCpsMaster = kvMaster.Count - p - 1;
      int countDimensionUnConstraint = 0;
      for (int i = 0; i < countDimension; i++)
      {
        if (ListDOFUnCoupling.IndexOf(i) != -1)
          countDimensionUnConstraint++;
      }
      int[] tArray = new int[(countDimension - countDimensionUnConstraint) * numCpsMaster];

      var cps = patchMaster.SelectNearEndPatchControlPoints(indexInterface);
      int c = 0;
      for (int i = 0; i < cps.Length; i++)
      {
        List<int> list = new List<int>();
        for (int d = 0; d < countDimension; d++)
        {
          if (ListDOFUnCoupling.IndexOf(d) == -1)
          {
            tArray[c++] = cps[i].GetTArrayGlobal()[d];
            list.Add(cps[i].GetTArrayGlobal()[d]);
          }
        }
        dicMaster.Add(i, list);
      }
      if (isMaster)
        this.dicMaster2 = dicMaster;
      else
        this.dicSlave2 = dicMaster;
      return tArray;
    }
    public int[] GetTArrayOnSurfaceMasterSlave(bool isMaster)
    {
      Dictionary<int, List<int>> dicMaster = new Dictionary<int, List<int>>();

      AbstractPatch3D patchMaster = null;
      int indexInterface;
      if (isMaster)
      {
        patchMaster = ((AbstractPatch3D)masterPatch);
        indexInterface = indexMasterObjectInterface;
      }
      else
      {
        patchMaster = ((AbstractPatch3D)slavePatch);
        indexInterface = indexSlaveObjectInterface;
      }
      Abstract3DParametricGeometry volumeMaster = patchMaster.GetVolume();
      int countField = patchMaster.GetCountField(0);
      int[] indexCoordinateMasterLine = patchMaster.GetCoordinateParameterOnArea(indexInterface);
      KnotVector[] kvMaster = {volumeMaster.Basis.GetKnotVector(indexCoordinateMasterLine[0]).Clone(),
        volumeMaster.Basis.GetKnotVector(indexCoordinateMasterLine[1]).Clone()};
      int[] p = { volumeMaster.Basis.GetDegree(indexCoordinateMasterLine[0]), volumeMaster.Basis.GetDegree(indexCoordinateMasterLine[1]) };
      int[] numCpsMaster = { kvMaster[0].Count - p[0] - 1, kvMaster[1].Count - p[1] - 1 };

      int[] tArray = new int[(countField - ListDOFUnCoupling.Count) * numCpsMaster[0] * numCpsMaster[1]];

      var cps = patchMaster.SelectEndPatchControlPoints(indexInterface);
      int c = 0;
      int cc = 0;
      for (int i = 0; i < cps.GetLength(0); i++)
        for (int j = 0; j < cps.GetLength(1); j++)
        {
          List<int> list = new List<int>();
          for (int d = 0; d < countField; d++)
          {
            if (ListDOFUnCoupling.IndexOf(d) == -1)
            {
              tArray[c++] = cps[i, j].GetTArrayGlobal()[d];
              list.Add(cps[i, j].GetTArrayGlobal()[d]);
            }
          }
          dicMaster.Add(cc++, list);
        }
      if (isMaster)
        this.dicMaster = dicMaster;
      else
        this.dicSlave = dicMaster;
      return tArray;
    }
    public int[] GetTArrayOnSurfaceMasterSlave2(bool isMaster)
    {
      Dictionary<int, List<int>> dicMaster = new Dictionary<int, List<int>>();

      AbstractPatch3D patchMaster = null;
      int indexInterface;
      if (isMaster)
      {
        patchMaster = ((AbstractPatch3D)masterPatch);
        indexInterface = indexMasterObjectInterface;
      }
      else
      {
        patchMaster = ((AbstractPatch3D)slavePatch);
        indexInterface = indexSlaveObjectInterface;
      }
      Abstract3DParametricGeometry volumeMaster = patchMaster.GetVolume();
      int countDimension = patchMaster.GetCountDimension();
      int[] indexCoordinateMasterLine = patchMaster.GetCoordinateParameterOnArea(indexInterface);
      KnotVector[] kvMaster = {volumeMaster.Basis.GetKnotVector(indexCoordinateMasterLine[0]).Clone(),
        volumeMaster.Basis.GetKnotVector(indexCoordinateMasterLine[1]).Clone()};
      int[] p = { volumeMaster.Basis.GetDegree(indexCoordinateMasterLine[0]), volumeMaster.Basis.GetDegree(indexCoordinateMasterLine[1]) };
      int[] numCpsMaster = { kvMaster[0].Count - p[0] - 1, kvMaster[1].Count - p[1] - 1 };
      int countDimensionUnConstraint = 0;
      for (int i = 0; i < countDimension; i++)
      {
        if (ListDOFUnCoupling.IndexOf(i) != -1)
          countDimensionUnConstraint++;
      }
      int[] tArray = new int[(countDimension - countDimensionUnConstraint) * numCpsMaster[0] * numCpsMaster[1]];

      var cps = patchMaster.SelectNearEndPatchControlPoints(indexInterface);
      int c = 0;
      int cc = 0;
      for (int i = 0; i < cps.GetLength(0); i++)
        for (int j = 0; j < cps.GetLength(1); j++)
        {
          List<int> list = new List<int>();
          for (int d = 0; d < countDimension; d++)
          {
            if (ListDOFUnCoupling.IndexOf(d) == -1)
            {
              tArray[c++] = cps[i, j].GetTArrayGlobal()[d];
              list.Add(cps[i, j].GetTArrayGlobal()[d]);
            }
          }
          dicMaster.Add(cc++, list);
        }
      if (isMaster)
        this.dicMaster2 = dicMaster;
      else
        this.dicSlave2 = dicMaster;
      return tArray;
    }
    public AbstractPatch GetMasterPatch()
    {
      return masterPatch;
    }

    public AbstractPatch GetSlavePatch()
    {
      return slavePatch;
    }

    public int GetIndexMasterObjectInterface()
    {
      return indexMasterObjectInterface;
    }

    public int GetIndexSlaveObjectInterface()
    {
      return indexSlaveObjectInterface;
    }

    public int GetIndexCoorinateMasterObjectInterface()
    {
      return indexMasterObjectInterface / 2;
    }

    public int GetIndexCoordinateSlaveObjectInterface()
    {
      return indexSlaveObjectInterface / 2;
    }

    public ControlPoint[] GetControlPointMasterObjectInterface()
    {
      return ((AbstractPatchOneField)masterPatch).SelectEndPatchControlPoints(indexMasterObjectInterface);
    }

    public ControlPoint[] GetControlPointSlaveObjectInterface()
    {
      return ((AbstractPatchOneField)slavePatch).SelectEndPatchControlPoints(indexSlaveObjectInterface);
    }
    public DoubleMatrix GetTSMMatrix()
    {
      return Tsm;
    }
    public DoubleMatrix GetTS2M1Matrix()
    {
      return Ts2m1;
    }
    public DoubleMatrix GetTS2M2Matrix()
    {
      return Ts2m2;
    }
    public int FindIndexControlPointFromTArray(bool isMaster, int tArray, out int indexField)
    {
      Dictionary<int, List<int>> dicMaster = null;
      if (isMaster)
        dicMaster = this.dicMaster;
      else
        dicMaster = this.dicSlave;
      for (int i = 0; i < dicMaster.Count; i++)
      {
        List<int> list = dicMaster[i];
        indexField = list.IndexOf(tArray);
        if (list.IndexOf(tArray) != -1)
          return i;
      }
      indexField = -1;
      return -1;
    }

    public int FindIndexControlPointFromTArray2(bool isMaster, int tArray, out int indexField)
    {
      Dictionary<int, List<int>> dicMaster = null;
      if (isMaster)
        dicMaster = this.dicMaster2;
      else
        dicMaster = this.dicSlave2;
      for (int i = 0; i < dicMaster.Count; i++)
      {
        List<int> list = dicMaster[i];
        indexField = list.IndexOf(tArray);
        if (list.IndexOf(tArray) != -1)
          return i;
      }
      indexField = -1;
      return -1;
    }

    //public void MakeCoupling(ref SparseMatrixBuilder<double> T)
    //{
    //  if (Tsm != null)
    //  {
    //    int[] tArrayMaster = null;
    //    int[] tArraySlave = null;
    //    if (masterPatch is AbstractPatch2D)
    //    {
    //      tArrayMaster = GetTArrayOnLineMasterSlave(true);
    //      tArraySlave = GetTArrayOnLineMasterSlave(false);
    //    }
    //    else if (masterPatch is AbstractPatch3D)
    //    {
    //      tArrayMaster = GetTArrayOnSurfaceMasterSlave(true);
    //      tArraySlave = GetTArrayOnSurfaceMasterSlave(false);
    //    }

    //    for (int iSlave = 0; iSlave < tArraySlave.Length; iSlave++)
    //    {
    //      int ti = tArraySlave[iSlave];
    //      T[ti, ti] += -1;
    //      int indexFieldSlave;
    //      int indexCpSlave = FindIndexControlPointFromTArray(false, ti, out indexFieldSlave);
    //      for (int iMaster = 0; iMaster < tArrayMaster.Length; iMaster++)
    //      {
    //        int tj = tArrayMaster[iMaster];
    //        int indexFieldMaster;
    //        int indexCpMaster = FindIndexControlPointFromTArray(true, tj, out indexFieldMaster);
    //        if (indexFieldMaster == indexFieldSlave)
    //        {
    //          T[ti, tj] += Tsm[indexCpSlave, indexCpMaster];
    //        }
    //      }
    //    }

    //    if (Ts2m2 != null)
    //    {
    //      int[] tArrayMaster2 = null;
    //      int[] tArraySlave2 = null;
    //      if (masterPatch is AbstractPatch2D)
    //      {
    //        tArrayMaster2 = GetTArrayOnLineMasterSlave2(true);
    //        tArraySlave2 = GetTArrayOnLineMasterSlave2(false);
    //      }
    //      else if (masterPatch is AbstractPatch3D)
    //      {
    //        tArrayMaster2 = GetTArrayOnSurfaceMasterSlave2(true);
    //        tArraySlave2 = GetTArrayOnSurfaceMasterSlave2(false);
    //      }
    //      for (int iSlave2 = 0; iSlave2 < tArraySlave2.Length; iSlave2++)
    //      {
    //        int ti = tArraySlave2[iSlave2];
    //        T[ti, ti] += -1;
    //        int indexFieldSlave2;
    //        int indexCpSlave2 = FindIndexControlPointFromTArray2(false, ti, out indexFieldSlave2);

    //        for (int iMaster2 = 0; iMaster2 < tArrayMaster2.Length; iMaster2++)
    //        {
    //          int tj = tArrayMaster2[iMaster2];
    //          int indexFieldMaster2;
    //          int indexCpMaster2 = FindIndexControlPointFromTArray2(true, tj, out indexFieldMaster2);
    //          if (indexFieldMaster2 == indexFieldSlave2)
    //          {
    //            T[ti, tj] += Ts2m2[indexCpSlave2, indexCpMaster2];
    //          }
    //        }

    //        for (int iMaster = 0; iMaster < tArrayMaster.Length; iMaster++)
    //        {
    //          int tj = tArrayMaster[iMaster];
    //          int indexFieldMaster1;
    //          int indexCpMaster1 = FindIndexControlPointFromTArray(true, tj, out indexFieldMaster1);
    //          if (indexFieldMaster1 == indexFieldSlave2)
    //          {
    //            T[ti, tj] += Ts2m1[indexCpSlave2, indexCpMaster1];
    //          }
    //        }
    //      }
    //    }
    //  }
    //}

    public void MakeCoupling(ref SparseMatrixBuilder<double> T)
    {
      if (Tsm != null)
      {
        int[] tArrayMaster = null;
        int[] tArraySlave = null;
        if (masterPatch is AbstractPatch2D)
        {
          tArrayMaster = GetTArrayOnLineMasterSlave(true);
          tArraySlave = GetTArrayOnLineMasterSlave(false);
        }
        else if (masterPatch is AbstractPatch3D)
        {
          tArrayMaster = GetTArrayOnSurfaceMasterSlave(true);
          tArraySlave = GetTArrayOnSurfaceMasterSlave(false);
        }
        for (int iSlave = 0; iSlave < tArraySlave.Length; iSlave++)
        {
          int ti = tArraySlave[iSlave];
          int indexFieldSlave;
          int indexCpSlave = FindIndexControlPointFromTArray(false, ti, out indexFieldSlave);
          //if (T[ti, ti] == 0)
          //{
          T[ti, ti] += -1;
          for (int iMaster = 0; iMaster < tArrayMaster.Length; iMaster++)
          {
            int tj = tArrayMaster[iMaster];
            int indexFieldMaster;
            int indexCpMaster = FindIndexControlPointFromTArray(true, tj, out indexFieldMaster);
            if (indexFieldMaster == indexFieldSlave)
            {
              T[ti, tj] += Tsm[indexCpSlave, indexCpMaster];
            }
          }
          //}
          //else if (T[ti, ti] == -1)
          //{
          //  for (int iMaster = 0; iMaster < tArrayMaster.Length; iMaster++)
          //  {
          //    int tj = tArrayMaster[iMaster];
          //    if (T[tj, tj] == 0)
          //    {
          //      int indexFieldMaster;
          //      int indexCpMaster = FindIndexControlPointFromTArray(true, tj, out indexFieldMaster);
          //      T[tj, tj] = -1;
          //      if (indexFieldMaster == indexFieldSlave)
          //      {
          //        T[ti, tj] += Tsm[indexCpSlave, indexCpMaster];
          //      }
          //      break;
          //    }
          //  }
          //}
        }
      }
    }

    public void MakeCouplingTs1m2(ref SparseMatrixBuilder<double> T)
    {
      if (Tsm != null)
      {
        int[] tArrayMaster = null;
        if (masterPatch is AbstractPatch2D)
        {
          tArrayMaster = GetTArrayOnLineMasterSlave(true);
        }
        else if (masterPatch is AbstractPatch3D)
        {
          tArrayMaster = GetTArrayOnSurfaceMasterSlave(true);
        }

        if (Ts2m2 != null)
        {
          int[] tArrayMaster2 = null;
          int[] tArraySlave2 = null;
          if (masterPatch is AbstractPatch2D)
          {
            tArrayMaster2 = GetTArrayOnLineMasterSlave2(true);
            tArraySlave2 = GetTArrayOnLineMasterSlave2(false);
          }
          else if (masterPatch is AbstractPatch3D)
          {
            tArrayMaster2 = GetTArrayOnSurfaceMasterSlave2(true);
            tArraySlave2 = GetTArrayOnSurfaceMasterSlave2(false);
          }
          for (int iSlave2 = 0; iSlave2 < tArraySlave2.Length; iSlave2++)
          {
            int ti = tArraySlave2[iSlave2];
            if (T[ti, ti] == 0)
            {
              T[ti, ti] += -1;
              int indexFieldSlave2;
              int indexCpSlave2 = FindIndexControlPointFromTArray2(false, ti, out indexFieldSlave2);

              for (int iMaster = 0; iMaster < tArrayMaster.Length; iMaster++)
              {
                int tj = tArrayMaster[iMaster];
                int indexFieldMaster1;
                int indexCpMaster1 = FindIndexControlPointFromTArray(true, tj, out indexFieldMaster1);
                if (indexFieldMaster1 == indexFieldSlave2)
                {
                  T[ti, tj] += Ts2m1[indexCpSlave2, indexCpMaster1];
                }
              }

              for (int iMaster2 = 0; iMaster2 < tArrayMaster2.Length; iMaster2++)
              {
                int tj = tArrayMaster2[iMaster2];
                int indexFieldMaster2;
                int indexCpMaster2 = FindIndexControlPointFromTArray2(true, tj, out indexFieldMaster2);
                if (indexFieldMaster2 == indexFieldSlave2)
                {
                  T[ti, tj] += Ts2m2[indexCpSlave2, indexCpMaster2];
                }
              }
            }
            //else if (T[ti, ti] == -1)
            //{
            //  bool a = true;
            //}
          }
        }
      }
    }

    //public void MakeCoupling(ref DoubleMatrix T)
    //{
    //if (Tsm != null)
    //{
    //  int[] tArrayMaster = null;
    //  int[] tArraySlave = null;
    //  if (masterPatch is AbstractPatch2D)
    //  {
    //    tArrayMaster = GetTArrayOnLineMasterSlave(true);
    //    tArraySlave = GetTArrayOnLineMasterSlave(false);
    //  }
    //  else if (masterPatch is AbstractPatch3D)
    //  {
    //    tArrayMaster = GetTArrayOnSurfaceMasterSlave(true);
    //    tArraySlave = GetTArrayOnSurfaceMasterSlave(false);
    //  }
    //  for (int iSlave = 0; iSlave < tArraySlave.Length; iSlave++)
    //  {
    //    int ti = tArraySlave[iSlave];
    //    int indexFieldSlave;
    //    int indexCpSlave = FindIndexControlPointFromTArray(false, ti, out indexFieldSlave);
    //    if (T[ti, ti] == 0)
    //    {
    //      T[ti, ti] += -1;
    //      for (int iMaster = 0; iMaster < tArrayMaster.Length; iMaster++)
    //      {
    //        int tj = tArrayMaster[iMaster];
    //        int indexFieldMaster;
    //        int indexCpMaster = FindIndexControlPointFromTArray(true, tj, out indexFieldMaster);
    //        if (indexFieldMaster == indexFieldSlave)
    //        {
    //          T[ti, tj] += Tsm[indexCpSlave, indexCpMaster];
    //        }
    //      }
    //    }
    //    //else if (T[ti, ti] == -1)
    //    //{
    //    //  for (int iMaster = 0; iMaster < tArrayMaster.Length; iMaster++)
    //    //  {
    //    //    int tj = tArrayMaster[iMaster];
    //    //    if (T[tj, tj] == 0)
    //    //    {
    //    //      int indexFieldMaster;
    //    //      int indexCpMaster = FindIndexControlPointFromTArray(true, tj, out indexFieldMaster);
    //    //      T[tj, tj] = -1;
    //    //      if (indexFieldMaster == indexFieldSlave)
    //    //      {
    //    //        T[ti, tj] += Tsm[indexCpSlave, indexCpMaster];
    //    //      }
    //    //      break;
    //    //    }
    //    //  }
    //    //}
    //  }

    //  if (Ts2m2 != null)
    //  {
    //    int[] tArrayMaster2 = null;
    //    int[] tArraySlave2 = null;
    //    if (masterPatch is AbstractPatch2D)
    //    {
    //      tArrayMaster2 = GetTArrayOnLineMasterSlave2(true);
    //      tArraySlave2 = GetTArrayOnLineMasterSlave2(false);
    //    }
    //    else if (masterPatch is AbstractPatch3D)
    //    {
    //      tArrayMaster2 = GetTArrayOnSurfaceMasterSlave2(true);
    //      tArraySlave2 = GetTArrayOnSurfaceMasterSlave2(false);
    //    }
    //    for (int iSlave2 = 0; iSlave2 < tArraySlave2.Length; iSlave2++)
    //    {
    //      int ti = tArraySlave2[iSlave2];
    //      if (T[ti, ti] == 0)
    //      {
    //        T[ti, ti] += -1;
    //        int indexFieldSlave2;
    //        int indexCpSlave2 = FindIndexControlPointFromTArray2(false, ti, out indexFieldSlave2);

    //        for (int iMaster = 0; iMaster < tArrayMaster.Length; iMaster++)
    //        {
    //          int tj = tArrayMaster[iMaster];
    //          int indexFieldMaster1;
    //          int indexCpMaster1 = FindIndexControlPointFromTArray(true, tj, out indexFieldMaster1);
    //          if (indexFieldMaster1 == indexFieldSlave2)
    //          {
    //            T[ti, tj] += Ts2m1[indexCpSlave2, indexCpMaster1];
    //          }
    //        }

    //        for (int iMaster2 = 0; iMaster2 < tArrayMaster2.Length; iMaster2++)
    //        {
    //          int tj = tArrayMaster2[iMaster2];
    //          int indexFieldMaster2;
    //          int indexCpMaster2 = FindIndexControlPointFromTArray2(true, tj, out indexFieldMaster2);
    //          if (indexFieldMaster2 == indexFieldSlave2)
    //          {
    //            T[ti, tj] += Ts2m2[indexCpSlave2, indexCpMaster2];
    //          }
    //        }
    //      }
    //      else if (T[ti, ti] == -1)
    //      {
    //        bool a = true;
    //      }
    //    }
    //  }
    //}
    //}

    public void MakeCoupling(ref DoubleMatrix T)
    {
      if (Tsm != null)
      {
        int[] tArrayMaster = null;
        int[] tArraySlave = null;
        if (masterPatch is AbstractPatch2D)
        {
          tArrayMaster = GetTArrayOnLineMasterSlave(true);
          tArraySlave = GetTArrayOnLineMasterSlave(false);
        }
        else if (masterPatch is AbstractPatch3D)
        {
          tArrayMaster = GetTArrayOnSurfaceMasterSlave(true);
          tArraySlave = GetTArrayOnSurfaceMasterSlave(false);
        }
        for (int iSlave = 0; iSlave < tArraySlave.Length; iSlave++)
        {
          int ti = tArraySlave[iSlave];
          int indexFieldSlave;
          int indexCpSlave = FindIndexControlPointFromTArray(false, ti, out indexFieldSlave);
          //if (T[ti, ti] == 0)
          //{
          T[ti, ti] += -1;
          for (int iMaster = 0; iMaster < tArrayMaster.Length; iMaster++)
          {
            int tj = tArrayMaster[iMaster];
            int indexFieldMaster;
            int indexCpMaster = FindIndexControlPointFromTArray(true, tj, out indexFieldMaster);
            if (indexFieldMaster == indexFieldSlave)
            {
              T[ti, tj] += Tsm[indexCpSlave, indexCpMaster];
            }
          }
          //}
          //else if (T[ti, ti] == -1)
          //{
          //  for (int iMaster = 0; iMaster < tArrayMaster.Length; iMaster++)
          //  {
          //    int tj = tArrayMaster[iMaster];
          //    if (T[tj, tj] == 0)
          //    {
          //      int indexFieldMaster;
          //      int indexCpMaster = FindIndexControlPointFromTArray(true, tj, out indexFieldMaster);
          //      T[tj, tj] = -1;
          //      if (indexFieldMaster == indexFieldSlave)
          //      {
          //        T[ti, tj] += Tsm[indexCpSlave, indexCpMaster];
          //      }
          //      break;
          //    }
          //  }
          //}
        }
      }
    }

    public void MakeCouplingTs1m2(ref DoubleMatrix T)
    {
      if (Tsm != null)
      {
        int[] tArrayMaster = null;
        if (masterPatch is AbstractPatch2D)
        {
          tArrayMaster = GetTArrayOnLineMasterSlave(true);
        }
        else if (masterPatch is AbstractPatch3D)
        {
          tArrayMaster = GetTArrayOnSurfaceMasterSlave(true);
        }

        if (Ts2m2 != null)
        {
          int[] tArrayMaster2 = null;
          int[] tArraySlave2 = null;
          if (masterPatch is AbstractPatch2D)
          {
            tArrayMaster2 = GetTArrayOnLineMasterSlave2(true);
            tArraySlave2 = GetTArrayOnLineMasterSlave2(false);
          }
          else if (masterPatch is AbstractPatch3D)
          {
            tArrayMaster2 = GetTArrayOnSurfaceMasterSlave2(true);
            tArraySlave2 = GetTArrayOnSurfaceMasterSlave2(false);
          }
          for (int iSlave2 = 0; iSlave2 < tArraySlave2.Length; iSlave2++)
          {
            int ti = tArraySlave2[iSlave2];
            if (T[ti, ti] == 0)
            {
              T[ti, ti] += -1;
              int indexFieldSlave2;
              int indexCpSlave2 = FindIndexControlPointFromTArray2(false, ti, out indexFieldSlave2);

              for (int iMaster = 0; iMaster < tArrayMaster.Length; iMaster++)
              {
                int tj = tArrayMaster[iMaster];
                int indexFieldMaster1;
                int indexCpMaster1 = FindIndexControlPointFromTArray(true, tj, out indexFieldMaster1);
                if (indexFieldMaster1 == indexFieldSlave2)
                {
                  T[ti, tj] += Ts2m1[indexCpSlave2, indexCpMaster1];
                }
              }

              for (int iMaster2 = 0; iMaster2 < tArrayMaster2.Length; iMaster2++)
              {
                int tj = tArrayMaster2[iMaster2];
                int indexFieldMaster2;
                int indexCpMaster2 = FindIndexControlPointFromTArray2(true, tj, out indexFieldMaster2);
                if (indexFieldMaster2 == indexFieldSlave2)
                {
                  T[ti, tj] += Ts2m2[indexCpSlave2, indexCpMaster2];
                }
              }
            }
          }
        }
      }
    }



    public bool Is1_1Constraint
    {
      get { return is1_1Constraint; }
    }
  }
}
