using CenterSpace.NMath.Core;
using DEMSoft.NURBS;

namespace DEMSoft.IGA
{
  public class Face : AbstractMeshPart
  {
    private int index;
    private Edge[] edges;
    private AbstractElement element;
    private Volume volume;
    private int[] p;

    public Face(AbstractElement element)
    {
      this.element = element;
      var patch = element.GetPatch();
      if (patch is AbstractPatchOneField)
      {
        SetNumberOfFields(((AbstractPatchOneField)patch).GetCountField());
        SetDimension(2);
      }
      else
      {
        SetNumberOfFields(patch.GetCountField(0) + patch.GetCountField(1));
        SetDimension(2);
      }
      edges = new Edge[4];
      for (int i = 0; i < edges.Length; i++)
      {
        edges[i] = new Edge(i, this);
      }
    }

    public Face(int index, Volume volume)
    {
      this.index = index;
      this.volume = volume;

      var mesh = volume.GetElement().GetPatch();
      var vo = (NURBSVolume)mesh.GetGeometry(0);
      var basis = vo.Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int p1 = -1;
      int p2 = -1;
      switch (GetIndexCoordinate())
      {
        case 0:
          p1 = p;
          p2 = q;
          break;
        case 1:
          p1 = p;
          p2 = r;
          break;
        case 2:
          p1 = q;
          p2 = r;
          break;
      }
      this.p = new int[] { p1, p2 };
      SetNumberOfFields(volume.GetNumberOfFields());
      SetDimension(3);
    }

    //public int[,] CreateIndexCoordinate()
    //{
    //    int[,] INC = new int[(this.p[0] + 1) * (this.p[1] + 1), 2];
    //    int count = 0;
    //    for (int j = 0; j <= this.p[1]; j++)
    //        for (int i = 0; i <= this.p[0]; i++)
    //        {
    //            INC[count, 0] = i;
    //            INC[count, 1] = j;
    //            count++;
    //        }
    //    return INC;
    //}

    public override int[] GetTArrayGlobal()
    {
      int[] tArray = null;
      var d = GetNumberOfFields();
      if (volume == null)
      {
        if (element.GetPatch() is AbstractPatch2D)
        {
          AbstractPatch2D patch = (AbstractPatch2D)element.GetPatch();
          int nen = patch.GetCountLocalBasisFunctions();
          tArray = new int[d * nen];
          var cps = patch.GetSurface().ControlPoints;
          for (int kk = 0; kk < nen; kk++)
          {
            int ien = patch.GetIEN(0, element.GetID(), kk);
            var tArrayCps = cps[patch.GetINC(0, ien, 0), patch.GetINC(0, ien, 1)].GetTArrayGlobal();
            for (int j = 0; j < d; j++)
            {
              tArray[kk * d + j] = tArrayCps[j];
            }
          }
        }
        else
        {
          //int nen0 = mesh.GetCountLocalBasisFunctions(0);
          //int nen1 = mesh.GetCountLocalBasisFunctions(1);
          //int d0 = mesh.GetCountField(0);
          //int d1 = mesh.GetCountField(1);
          //tArray = new int[d0 * nen0 + d1 * nen1];
          //for (int k = 0; k < mesh.GetCountGeometry(); k++)
          //{
          //    int dk = mesh.GetCountField(k);
          //    int nen = mesh.GetCountLocalBasisFunctions(k);
          //    var cps = ((NURBSSurface)(mesh.GetGeometry(k))).ControlPoints;
          //    for (int kk = 0; kk < nen; kk++)
          //    {
          //        int ien = mesh.GetIEN(k, element.GetID(), kk);
          //        var tArrayCps = cps[mesh.GetINC(k, ien, 0), mesh.GetINC(k, ien, 1)].GetTArray();
          //        for (int j = 0; j < dk; j++)
          //        {
          //            tArray[k * (nen0 * d0) + kk * dk + j] = tArrayCps[j];
          //        }
          //    }
          //}
        }
      }
      else
      {
        int numCps = (p[0] + 1) * (p[1] + 1);
        tArray = new int[d * numCps];
        var cpsFace = GetControlPointsOnFace();

        for (int j = 0; j <= p[1]; j++)
          for (int i = 0; i <= p[0]; i++)
            for (int k = 0; k < d; k++)
            {
              tArray[(j * (p[0] + 1) + i) * d + k] = cpsFace[i, j].GetTArrayGlobal()[k];
            }
      }
      return tArray;
    }

    /// ////////////////////////////// Tien
    /// <summary>
    /// </summary>
    /// <returns></returns>
    public override int[] GetTArray()
    {
      int[] tArray = null;
      var d = GetNumberOfFields();

      int numCps = (p[0] + 1) * (p[1] + 1);
      tArray = new int[d * numCps];
      var cpsFace = GetControlPointsOnFace();

      for (int j = 0; j <= p[1]; j++)
        for (int i = 0; i <= p[0]; i++)
          for (int k = 0; k < d; k++)
          {
            tArray[(j * (p[0] + 1) + i) * d + k] = cpsFace[i, j].GetTArray()[k];
          }
      return tArray;
    }
    ///////////

    public AbstractElement GetElement()
    {
      return element;
    }

    public ControlPoint[,] GetControlPointsOnFace(int idx = 0)
    {
      ControlPoint[,] cpsLocal = null;
      if (volume == null)
      {
        var mesh = element.GetPatch();
        var basis = mesh.GetGeometry(0).Basis;
        int p = basis.GetDegree(0);
        int q = basis.GetDegree(1);
        int idPatch = element.GetID();
        int idxi = mesh.GetIPN(idPatch, 0);
        int idxj = mesh.GetIPN(idPatch, 1);
        var kv1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
        var kv2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
        double mid1 = (kv1[idxi] + kv1[idxi + 1]) / 2.0;
        double mid2 = (kv2[idxj] + kv2[idxj + 1]) / 2.0;
        int spanU = basis.FindSpan(mid1, 0);
        int spanV = basis.FindSpan(mid2, 1);
        cpsLocal = new ControlPoint[p + 1, q + 1];
        for (int i = 0; i < mesh.GetCountLocalBasisFunctions(0); i++)
        {
          int ien = mesh.GetIEN(0, idPatch, i);
          int inc1 = mesh.GetINC(0, ien, 0);
          int inc2 = mesh.GetINC(0, ien, 1);
          cpsLocal[inc1 - spanU + p, inc2 - spanV + q] = ((Abstract2DParametricGeometry)(mesh.GetGeometry(idx))).ControlPoints[inc1, inc2];
        }
      }
      else
      {
        var mesh = volume.GetElement().GetPatch();
        var vo = (NURBSVolume)mesh.GetGeometry(0);
        var basis = vo.Basis;
        int p = basis.GetDegree(0);
        int q = basis.GetDegree(1);
        int r = basis.GetDegree(2);

        cpsLocal = new ControlPoint[this.p[0] + 1, this.p[1] + 1];
        var cpsVolume = volume.GetControlPointsOnVolume();

        for (int j = 0; j <= this.p[1]; j++)
          for (int i = 0; i <= this.p[0]; i++)
          {
            switch (GetIndexCoordinate())
            {
              case 0:
                if (GetIndexFrontBack() == 0)
                  cpsLocal[i, j] = cpsVolume[i, j, 0];
                else
                  cpsLocal[i, j] = cpsVolume[i, j, r];
                break;
              case 1:
                if (GetIndexFrontBack() == 0)
                  cpsLocal[i, j] = cpsVolume[i, 0, j];
                else
                  cpsLocal[i, j] = cpsVolume[i, q, j];
                break;
              case 2:
                if (GetIndexFrontBack() == 0)
                  cpsLocal[i, j] = cpsVolume[0, i, j];
                else
                  cpsLocal[i, j] = cpsVolume[p, i, j];
                break;
            }
          }
      }
      return cpsLocal;
    }

    ///////////////////////////////// Tien
    public int CountDofOnFace(int idx = 0)
    {
      var mesh = volume.GetElement().GetPatch();
      return mesh.GetCountField(0) * ((this.p[0] + 1) * (this.p[1] + 1));
    }
    public int CountControlPointOnFace(int idx = 0)
    {

      return (this.p[0] + 1) * (this.p[1] + 1);
    }
    ////////
    public Edge GetEdge(int index)
    {
      return edges[index];
    }

    public Volume GetVolumeBeAttached()
    {
      return volume;
    }

    public int GetDegree(int index)
    {
      return p[index];
    }

    /// <summary>
    /// Get parametric coordinate on end of face on 2 direction
    /// </summary>
    /// <param name="index">0:first direction, 1:second direction</param>
    /// <returns></returns>
    public double[] GetParametricEndFace(int index)
    {
      int indexDirect = -1;
      AbstractParametricBasisFunction basis = null;
      int idxi = -1;
      AbstractPatch patch = null;
      if (volume != null)
      {
        patch = volume.GetElement().GetPatch();
        var vol = (NURBSVolume)patch.GetGeometry(0);
        basis = vol.Basis;
        switch (GetIndexCoordinate())
        {
          case 0:
            if (GetIndexFrontBack() == 0)
              indexDirect = 0;
            else
              indexDirect = 1;
            break;
          case 1:
            if (GetIndexFrontBack() == 0)
              indexDirect = 0;
            else
              indexDirect = 2;
            break;
          case 2:
            if (GetIndexFrontBack() == 0)
              indexDirect = 1;
            else
              indexDirect = 2;
            break;
        }

        idxi = patch.GetIPN(volume.GetElement().GetID(), indexDirect);
      }
      else
      {
      }
      var kv1 = basis.GetKnotVector(indexDirect).GetKnotVectorNoMultiplicity();
      return new double[] { kv1[idxi], kv1[idxi + 1] };
    }

    public int CountDof()
    {
      if (volume == null)
      {
        var patch = element.GetPatch();
        int d = patch.GetCountField(0);
        int nen = patch.GetCountLocalBasisFunctions(0);
        return nen * d;
      }
      else
      {
        var patch = element.GetPatch();
        int d = patch.GetCountField(0);
        return d * (this.p[0] + 1) * (this.p[1] + 1);
      }
    }

    public int GetIndexFrontBack()
    {
      return index % 2;
    }

    public int GetIndexCoordinate()
    {
      //if (index == 0 || index == 1)
      //    return 0;
      //else if (index == 2 || index == 3)
      //    return 1;
      //else if (index == 4 || index == 5)
      //    return 2;
      //else
      //    return -1;
      return index / 2;
    }

    public double[,] GetBivariateBasisFunctionOnFace(double xi, double eta)
    {
      var mesh = GetElement().GetPatch();
      var sur = (NURBSSurface)mesh.GetGeometry(0);
      var basis = (BivariateNURBSBasisFunction)sur.Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int spanU = basis.FindSpan(xi, 0);
      int spanV = basis.FindSpan(eta, 1);
      double[,] Nij = new double[p + 1, q + 1];
      for (int i = 0; i <= p; i++)
        for (int j = 0; j <= q; j++)
        {
          Nij[i, j] = basis.GetValueBivariateBasisFunction(xi, eta, spanU - p + i, spanV - q + j);//Nijk[i, j, 0];                           
        }
      return Nij;
    }


    public double[,] GetTrivariateBasisFunctionOnFace(double xi, double eta, double zeta)
    {
      var mesh = volume.GetElement().GetPatch();
      var vol = (NURBSVolume)mesh.GetGeometry(0);
      var basis = (TrivariateNURBSBasisFunction)vol.Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int spanU = basis.FindSpan(xi, 0);
      int spanV = basis.FindSpan(eta, 1);
      int spanW = basis.FindSpan(zeta, 2);
      double[,] Nij = new double[this.p[0] + 1, this.p[1] + 1];
      for (int i = 0; i <= this.p[0]; i++)
        for (int j = 0; j <= this.p[1]; j++)
        {
          switch (GetIndexCoordinate())
          {
            case 0:
              if (GetIndexFrontBack() == 0)
                Nij[i, j] = basis.GetValueTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + j, spanW - r + 0);//Nijk[0, i, j];
              else
                Nij[i, j] = basis.GetValueTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + j, spanW - r + r);//Nijk[p, i, j];
              break;
            case 1:
              if (GetIndexFrontBack() == 0)
                Nij[i, j] = basis.GetValueTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + 0, spanW - r + j);//Nijk[i, 0, j];
              else
                Nij[i, j] = basis.GetValueTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + q, spanW - r + j);//Nijk[i, q, j];
              break;
            case 2:
              if (GetIndexFrontBack() == 0)
                Nij[i, j] = basis.GetValueTrivariateBasisFunction(xi, eta, zeta, spanU - p + 0, spanV - q + i, spanW - r + j);//Nijk[i, j, 0];
              else
                Nij[i, j] = basis.GetValueTrivariateBasisFunction(xi, eta, zeta, spanU - p + p, spanV - q + i, spanW - r + j);//Nijk[i, j, r];
              break;
          }
        }
      return Nij;
    }

    private DoubleMatrix JacobianAt(DoubleMatrix dNdxi)
    {
      DoubleMatrix dxdxi = new DoubleMatrix(2, 2);//dimension equals 2 on face
      ControlPoint[,] cps2d = GetControlPointsOnFace();
      int nen = (this.p[0] + 1) * (this.p[1] + 1);
      ControlPoint[] cps1d = new ControlPoint[nen];
      int count = 0;
      for (int j = 0; j <= this.p[1]; j++)
        for (int i = 0; i <= this.p[0]; i++)
        {
          cps1d[count++] = cps2d[i, j];
        }

      for (int j = 0; j < 2; j++)
      {
        for (int k = 0; k < nen; k++)
        {
          var controlPointCoord = cps1d[k].GetCoordinate();

          switch (GetIndexCoordinate())
          {
            case 0:
              dxdxi[0, j] += dNdxi[k, j] * controlPointCoord[0];
              dxdxi[1, j] += dNdxi[k, j] * controlPointCoord[1];
              break;
            case 1:
              dxdxi[0, j] += dNdxi[k, j] * controlPointCoord[0];
              dxdxi[1, j] += dNdxi[k, j] * controlPointCoord[2];
              break;
            case 2:
              dxdxi[0, j] += dNdxi[k, j] * controlPointCoord[1];
              dxdxi[1, j] += dNdxi[k, j] * controlPointCoord[2];
              break;
          }
        }
      }
      return dxdxi;
    }
    private DoubleMatrix Jacobian2At(DoubleMatrix ddNdxi)
    {
      int nen = (this.p[0] + 1) * (this.p[1] + 1);
      DoubleMatrix ddxdxi = new DoubleMatrix(3, 2);
      ControlPoint[,] cps2d = GetControlPointsOnFace();
      ControlPoint[] cps1d = new ControlPoint[nen];
      int count = 0;
      for (int j = 0; j <= this.p[1]; j++)
        for (int i = 0; i <= this.p[0]; i++)
        {
          cps1d[count++] = cps2d[i, j];
        }
      for (int j = 0; j < 3; j++)
      {
        for (int k = 0; k < nen; k++)
        {
          var controlPointCoord = cps1d[k].GetCoordinate();

          switch (GetIndexCoordinate())
          {
            case 0:
              ddxdxi[j, 0] += ddNdxi[k, j] * controlPointCoord[0];
              ddxdxi[j, 1] += ddNdxi[k, j] * controlPointCoord[1];
              break;
            case 1:
              ddxdxi[j, 0] += ddNdxi[k, j] * controlPointCoord[0];
              ddxdxi[j, 1] += ddNdxi[k, j] * controlPointCoord[2];
              break;
            case 2:
              ddxdxi[j, 0] += ddNdxi[k, j] * controlPointCoord[1];
              ddxdxi[j, 1] += ddNdxi[k, j] * controlPointCoord[2];
              break;
          }
        }
      }
      return ddxdxi;
    }
    private DoubleMatrix Jacobian13At(DoubleMatrix d3Ndxi)
    {
      int nen = (this.p[0] + 1) * (this.p[1] + 1);
      DoubleMatrix d3xdxi = new DoubleMatrix(4, 2);
      ControlPoint[,] cps2d = GetControlPointsOnFace();
      ControlPoint[] cps1d = new ControlPoint[nen];
      int count = 0;
      for (int j = 0; j <= this.p[1]; j++)
        for (int i = 0; i <= this.p[0]; i++)
        {
          cps1d[count++] = cps2d[i, j];
        }
      for (int j = 0; j < 4; j++)
      {
        for (int k = 0; k < nen; k++)
        {
          var controlPointCoord = cps1d[k].GetCoordinate();

          switch (GetIndexCoordinate())
          {
            case 0:
              d3xdxi[j, 0] += d3Ndxi[k, j] * controlPointCoord[0];
              d3xdxi[j, 1] += d3Ndxi[k, j] * controlPointCoord[1];
              break;
            case 1:
              d3xdxi[j, 0] += d3Ndxi[k, j] * controlPointCoord[0];
              d3xdxi[j, 1] += d3Ndxi[k, j] * controlPointCoord[2];
              break;
            case 2:
              d3xdxi[j, 0] += d3Ndxi[k, j] * controlPointCoord[1];
              d3xdxi[j, 1] += d3Ndxi[k, j] * controlPointCoord[2];
              break;
          }
        }
      }
      return d3xdxi;
    }

    private DoubleMatrix Jacobian14At(DoubleMatrix d4Ndxi)
    {
      int nen = (this.p[0] + 1) * (this.p[1] + 1);
      DoubleMatrix d4xdxi = new DoubleMatrix(5, 2);
      ControlPoint[,] cps2d = GetControlPointsOnFace();
      ControlPoint[] cps1d = new ControlPoint[nen];
      int count = 0;
      for (int j = 0; j <= this.p[1]; j++)
        for (int i = 0; i <= this.p[0]; i++)
        {
          cps1d[count++] = cps2d[i, j];
        }
      for (int j = 0; j < 5; j++)
      {
        for (int k = 0; k < nen; k++)
        {
          var controlPointCoord = cps1d[k].GetCoordinate();

          switch (GetIndexCoordinate())
          {
            case 0:
              d4xdxi[j, 0] += d4Ndxi[k, j] * controlPointCoord[0];
              d4xdxi[j, 1] += d4Ndxi[k, j] * controlPointCoord[1];
              break;
            case 1:
              d4xdxi[j, 0] += d4Ndxi[k, j] * controlPointCoord[0];
              d4xdxi[j, 1] += d4Ndxi[k, j] * controlPointCoord[2];
              break;
            case 2:
              d4xdxi[j, 0] += d4Ndxi[k, j] * controlPointCoord[1];
              d4xdxi[j, 1] += d4Ndxi[k, j] * controlPointCoord[2];
              break;
          }
        }
      }
      return d4xdxi;
    }

    public DoubleMatrix GetFirstDerivativeTrivariateBasisFunctionOnFace(double xi, double eta, double zeta)
    {
      var mesh = volume.GetElement().GetPatch();
      var vol = (NURBSVolume)mesh.GetGeometry(0);
      var basis = (TrivariateNURBSBasisFunction)vol.Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int spanU = basis.FindSpan(xi, 0);
      int spanV = basis.FindSpan(eta, 1);
      int spanW = basis.FindSpan(zeta, 2);
      //var dNijk = basis.GetDerivativeTrivariateBasisFunctions(xi, eta, zeta, 1);
      DoubleMatrix dNij = new DoubleMatrix((this.p[0] + 1) * (this.p[1] + 1), 2);
      //int count = 0;
      for (int j = 0; j <= this.p[1]; j++)
        for (int i = 0; i <= this.p[0]; i++)
        {
          double[,,] dNdxi = null;
          switch (GetIndexCoordinate())
          {
            case 0:
              if (GetIndexFrontBack() == 0)
                dNdxi = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + j, spanW - r + 0, 1);
              else
                dNdxi = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + j, spanW - r + r, 1);
              dNij[j * (this.p[0] + 1) + i, 0] = dNdxi[1, 0, 0];
              dNij[j * (this.p[0] + 1) + i, 1] = dNdxi[0, 1, 0];
              break;
            case 1:
              if (GetIndexFrontBack() == 0)
                dNdxi = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + 0, spanW - r + j, 1);
              else
                dNdxi = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + q, spanW - r + j, 1);
              dNij[j * (this.p[0] + 1) + i, 0] = dNdxi[1, 0, 0];
              dNij[j * (this.p[0] + 1) + i, 1] = dNdxi[0, 0, 1];
              break;
            case 2:
              if (GetIndexFrontBack() == 0)
                dNdxi = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + 0, spanV - q + i, spanW - r + j, 1);
              else
                dNdxi = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + p, spanV - q + i, spanW - r + j, 1);
              dNij[j * (this.p[0] + 1) + i, 0] = dNdxi[0, 1, 0];
              dNij[j * (this.p[0] + 1) + i, 1] = dNdxi[0, 0, 1];
              break;
          }
          //count++;
        }
      return dNij;
    }

    public DoubleMatrix GetFirstDerivativePhysicalCoordinationOnFace(double xi, double eta, double zeta)
    {
      DoubleMatrix dNdxi = GetFirstDerivativeTrivariateBasisFunctionOnFace(xi, eta, zeta);
      DoubleMatrix J = JacobianAt(dNdxi);
      DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
      DoubleMatrix dNdx = MatrixFunctions.Product(dNdxi, invertJ);
      return dNdx;
    }

    public DoubleMatrix GetSecondDerivativeTrivariateBasisFunctionOnFace(double xi, double eta, double zeta)
    {
      var mesh = volume.GetElement().GetPatch();
      var vol = (NURBSVolume)mesh.GetGeometry(0);
      var basis = (TrivariateNURBSBasisFunction)vol.Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int spanU = basis.FindSpan(xi, 0);
      int spanV = basis.FindSpan(eta, 1);
      int spanW = basis.FindSpan(zeta, 2);
      //var dNijk = basis.GetDerivativeTrivariateBasisFunctions(xi, eta, zeta, 1);
      DoubleMatrix ddNij = new DoubleMatrix((this.p[0] + 1) * (this.p[1] + 1), 3);// column 0: ddXi, column 1: ddEta, column 2: dXidEta
      int count = 0;
      for (int j = 0; j <= this.p[1]; j++)
        for (int i = 0; i <= this.p[0]; i++)
        {
          double[,,] ddNijk = null;
          switch (GetIndexCoordinate())
          {
            case 0:
              if (GetIndexFrontBack() == 0)
                ddNijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + j, spanW - r + 0, 2);
              else
                ddNijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + j, spanW - r + r, 2);
              ddNij[j * (this.p[0] + 1) + i, 0] = ddNijk[2, 0, 0];
              ddNij[j * (this.p[0] + 1) + i, 1] = ddNijk[0, 2, 0];
              ddNij[j * (this.p[0] + 1) + i, 2] = ddNijk[1, 1, 0];
              break;
            case 1:
              if (GetIndexFrontBack() == 0)
                ddNijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + 0, spanW - r + j, 2);
              else
                ddNijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + q, spanW - r + j, 2);
              ddNij[j * (this.p[0] + 1) + i, 0] = ddNijk[2, 0, 0];
              ddNij[j * (this.p[0] + 1) + i, 1] = ddNijk[0, 0, 2];
              ddNij[j * (this.p[0] + 1) + i, 2] = ddNijk[1, 0, 1];
              break;
            case 2:
              if (GetIndexFrontBack() == 0)
                ddNijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + 0, spanV - q + i, spanW - r + j, 2);
              else
                ddNijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + p, spanV - q + i, spanW - r + j, 2);
              ddNij[j * (this.p[0] + 1) + i, 0] = ddNijk[0, 2, 0];
              ddNij[j * (this.p[0] + 1) + i, 1] = ddNijk[0, 0, 2];
              ddNij[j * (this.p[0] + 1) + i, 2] = ddNijk[0, 1, 1];
              break;
          }
          count++;
        }
      return ddNij;
    }
    public DoubleMatrix GetThirdDerivativeTrivariateBasisFunctionOnFace(double xi, double eta, double zeta)
    {
      var mesh = volume.GetElement().GetPatch();
      var vol = (NURBSVolume)mesh.GetGeometry(0);
      var basis = (TrivariateNURBSBasisFunction)vol.Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int spanU = basis.FindSpan(xi, 0);
      int spanV = basis.FindSpan(eta, 1);
      int spanW = basis.FindSpan(zeta, 2);
      DoubleMatrix d3Nij = new DoubleMatrix((this.p[0] + 1) * (this.p[1] + 1), 4);//xxx, yyy, xxy, xyy
      int count = 0;
      for (int j = 0; j <= this.p[1]; j++)
        for (int i = 0; i <= this.p[0]; i++)
        {
          double[,,] d4Nijk = null;
          switch (GetIndexCoordinate())
          {
            case 0:
              if (GetIndexFrontBack() == 0)
                d4Nijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + j, spanW - r + 0, 3);
              else
                d4Nijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + j, spanW - r + r, 3);
              d3Nij[j * (this.p[0] + 1) + i, 0] = d4Nijk[3, 0, 0];
              d3Nij[j * (this.p[0] + 1) + i, 1] = d4Nijk[0, 3, 0];
              d3Nij[j * (this.p[0] + 1) + i, 2] = d4Nijk[2, 1, 0];
              d3Nij[j * (this.p[0] + 1) + i, 3] = d4Nijk[1, 2, 0];
              break;
            case 1:
              if (GetIndexFrontBack() == 0)
                d4Nijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + 0, spanW - r + j, 3);
              else
                d4Nijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + q, spanW - r + j, 3);
              d3Nij[j * (this.p[0] + 1) + i, 0] = d4Nijk[3, 0, 0];
              d3Nij[j * (this.p[0] + 1) + i, 1] = d4Nijk[0, 0, 3];
              d3Nij[j * (this.p[0] + 1) + i, 2] = d4Nijk[2, 0, 1];
              d3Nij[j * (this.p[0] + 1) + i, 3] = d4Nijk[1, 0, 2];
              break;
            case 2:
              if (GetIndexFrontBack() == 0)
                d4Nijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + 0, spanV - q + i, spanW - r + j, 3);
              else
                d4Nijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + p, spanV - q + i, spanW - r + j, 3);
              d3Nij[j * (this.p[0] + 1) + i, 0] = d4Nijk[0, 3, 0];
              d3Nij[j * (this.p[0] + 1) + i, 1] = d4Nijk[0, 0, 3];
              d3Nij[j * (this.p[0] + 1) + i, 2] = d4Nijk[0, 2, 1];
              d3Nij[j * (this.p[0] + 1) + i, 3] = d4Nijk[0, 1, 2];
              break;
          }
          count++;
        }
      return d3Nij;
    }
    public DoubleMatrix GetFourthDerivativeTrivariateBasisFunctionOnFace(double xi, double eta, double zeta)
    {
      var mesh = volume.GetElement().GetPatch();
      var vol = (NURBSVolume)mesh.GetGeometry(0);
      var basis = (TrivariateNURBSBasisFunction)vol.Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int spanU = basis.FindSpan(xi, 0);
      int spanV = basis.FindSpan(eta, 1);
      int spanW = basis.FindSpan(zeta, 2);
      //var dNijk = basis.GetDerivativeTrivariateBasisFunctions(xi, eta, zeta, 1);
      DoubleMatrix d4Nij = new DoubleMatrix((this.p[0] + 1) * (this.p[1] + 1), 5);//xxxx, yyyy, xxxy, xxyy, xyyy
      int count = 0;
      for (int j = 0; j <= this.p[1]; j++)
        for (int i = 0; i <= this.p[0]; i++)
        {
          double[,,] d4Nijk = null;
          switch (GetIndexCoordinate())
          {
            case 0:
              if (GetIndexFrontBack() == 0)
                d4Nijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + j, spanW - r + 0, 4);
              else
                d4Nijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + j, spanW - r + r, 4);
              d4Nij[j * (this.p[0] + 1) + i, 0] = d4Nijk[4, 0, 0];
              d4Nij[j * (this.p[0] + 1) + i, 1] = d4Nijk[0, 4, 0];
              d4Nij[j * (this.p[0] + 1) + i, 2] = d4Nijk[3, 1, 0];
              d4Nij[j * (this.p[0] + 1) + i, 3] = d4Nijk[2, 2, 0];
              d4Nij[j * (this.p[0] + 1) + i, 4] = d4Nijk[1, 3, 0];
              break;
            case 1:
              if (GetIndexFrontBack() == 0)
                d4Nijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + 0, spanW - r + j, 4);
              else
                d4Nijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + q, spanW - r + j, 4);
              d4Nij[j * (this.p[0] + 1) + i, 0] = d4Nijk[4, 0, 0];
              d4Nij[j * (this.p[0] + 1) + i, 1] = d4Nijk[0, 0, 4];
              d4Nij[j * (this.p[0] + 1) + i, 2] = d4Nijk[3, 0, 1];
              d4Nij[j * (this.p[0] + 1) + i, 3] = d4Nijk[2, 0, 2];
              d4Nij[j * (this.p[0] + 1) + i, 4] = d4Nijk[1, 0, 3];
              break;
            case 2:
              if (GetIndexFrontBack() == 0)
                d4Nijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + 0, spanV - q + i, spanW - r + j, 4);
              else
                d4Nijk = basis.GetDerivativeTrivariateBasisFunction(xi, eta, zeta, spanU - p + p, spanV - q + i, spanW - r + j, 4);
              d4Nij[j * (this.p[0] + 1) + i, 0] = d4Nijk[0, 4, 0];
              d4Nij[j * (this.p[0] + 1) + i, 1] = d4Nijk[0, 0, 4];
              d4Nij[j * (this.p[0] + 1) + i, 2] = d4Nijk[0, 3, 1];
              d4Nij[j * (this.p[0] + 1) + i, 3] = d4Nijk[0, 2, 2];
              d4Nij[j * (this.p[0] + 1) + i, 4] = d4Nijk[0, 1, 3];
              break;
          }
          count++;
        }
      return d4Nij;
    }
    public DoubleMatrix GetSecondDerivativePhysicalCoordinationOnFace(double xi, double eta, double zeta)
    {
      DoubleMatrix dNdxi = GetFirstDerivativeTrivariateBasisFunctionOnFace(xi, eta, zeta);//[nen,2]
      DoubleMatrix dxdxi = JacobianAt(dNdxi);//[2,2]
      DoubleMatrix dNdx = MatrixFunctions.Product(dNdxi, MatrixFunctions.Inverse(dxdxi));//[nen,2]*[2,2]-->[nen,2]
      DoubleMatrix ddNdxi = GetSecondDerivativeTrivariateBasisFunctionOnFace(xi, eta, zeta);//[nen,3]
      DoubleMatrix jac2 = new DoubleMatrix(3, 3);
      jac2[0, 0] = dxdxi[0, 0] * dxdxi[0, 0];
      jac2[0, 1] = dxdxi[1, 0] * dxdxi[1, 0];
      jac2[0, 2] = 2 * dxdxi[0, 0] * dxdxi[1, 0];
      jac2[1, 0] = dxdxi[0, 1] * dxdxi[0, 1];
      jac2[1, 1] = dxdxi[1, 1] * dxdxi[1, 1];
      jac2[1, 2] = 2 * dxdxi[0, 1] * dxdxi[1, 1];
      jac2[2, 0] = dxdxi[0, 0] * dxdxi[0, 1];
      jac2[2, 1] = dxdxi[1, 0] * dxdxi[1, 1];
      jac2[2, 2] = dxdxi[1, 0] * dxdxi[0, 1] + dxdxi[0, 0] * dxdxi[1, 1];
      DoubleMatrix jac1 = Jacobian2At(ddNdxi/*[nen,3]*/);//[3,2]
      DoubleMatrix term1 = MatrixFunctions.Product(dNdx, jac1.Transpose());//[nen,2]*[2,3]-->[nen,3]
      DoubleMatrix term2 = ddNdxi - term1;//[nen,3]-[nen,3]-->[nen,3]
      DoubleMatrix ddNdX = MatrixFunctions.Product(MatrixFunctions.Inverse(jac2), term2.Transpose());//[3,3]*[3,nen]-->[3,nen]
      return ddNdX.Transpose();//[nen,3]
    }

    public DoubleMatrix GetFourthDerivativePhysicalCoordinationOnFace(double xi, double eta, double zeta)
    {
      DoubleMatrix dNdxi = GetFirstDerivativeTrivariateBasisFunctionOnFace(xi, eta, zeta);//[nen,2]
      DoubleMatrix dxdxi = JacobianAt(dNdxi);//[2,2]
      DoubleMatrix invertJ = MatrixFunctions.Inverse(dxdxi);
      DoubleMatrix dNdx = MatrixFunctions.Product(dNdxi, invertJ);//[nen,2]*[2,2]-->[nen,2]
      DoubleMatrix ddNdxi = GetSecondDerivativeTrivariateBasisFunctionOnFace(xi, eta, zeta);//[nen,3]
      DoubleMatrix jac2 = new DoubleMatrix(3, 3);
      jac2[0, 0] = dxdxi[0, 0] * dxdxi[0, 0];
      jac2[0, 1] = dxdxi[1, 0] * dxdxi[1, 0];
      jac2[0, 2] = 2 * dxdxi[0, 0] * dxdxi[1, 0];
      jac2[1, 0] = dxdxi[0, 1] * dxdxi[0, 1];
      jac2[1, 1] = dxdxi[1, 1] * dxdxi[1, 1];
      jac2[1, 2] = 2 * dxdxi[0, 1] * dxdxi[1, 1];
      jac2[2, 0] = dxdxi[0, 0] * dxdxi[0, 1];
      jac2[2, 1] = dxdxi[1, 0] * dxdxi[1, 1];
      jac2[2, 2] = dxdxi[1, 0] * dxdxi[0, 1] + dxdxi[0, 0] * dxdxi[1, 1];
      DoubleMatrix ddxdxi = Jacobian2At(ddNdxi/*[nen,3]*/).Transpose();//[3,2]
      DoubleMatrix term1 = MatrixFunctions.Product(dNdx, ddxdxi);//[nen,2]*[2,3]-->[nen,3]
      DoubleMatrix term2 = ddNdxi - term1;//[nen,3]-[nen,3]-->[nen,3]
      DoubleMatrix ddNdX = MatrixFunctions.Product(MatrixFunctions.Inverse(jac2), term2.Transpose()).Transpose();//[3,3]*[3,nen]-->[3,nen]

      DoubleMatrix d3Ndxi = GetThirdDerivativeTrivariateBasisFunctionOnFace(xi, eta, zeta);//[nen,4]
      DoubleMatrix d4Ndxi = GetFourthDerivativeTrivariateBasisFunctionOnFace(xi, eta, zeta);//[nen,5]
      DoubleMatrix d3xdxi = Jacobian13At(d3Ndxi).Transpose();
      DoubleMatrix term13 = MatrixFunctions.Product(dNdx, d3xdxi/*.Transpose()*/);
      DoubleMatrix jac23 = new DoubleMatrix(4, 3);

      jac23[0, 0] = 3 * dxdxi[0, 0] * ddxdxi[0, 0];
      jac23[0, 1] = 3 * dxdxi[1, 0] * ddxdxi[1, 0];
      jac23[0, 2] = 3 * (dxdxi[0, 0] * ddxdxi[1, 0] + dxdxi[1, 0] * ddxdxi[0, 0]);

      jac23[1, 0] = 3 * dxdxi[0, 1] * ddxdxi[0, 1];
      jac23[1, 1] = 3 * dxdxi[1, 1] * ddxdxi[1, 1];
      jac23[1, 2] = 3 * (dxdxi[0, 1] * ddxdxi[1, 1] + dxdxi[1, 1] * ddxdxi[0, 1]);

      jac23[2, 0] = dxdxi[0, 1] * ddxdxi[0, 0] + 2 * dxdxi[0, 0] * ddxdxi[0, 2];
      jac23[2, 1] = dxdxi[1, 1] * ddxdxi[1, 0] + 2 * dxdxi[1, 0] * ddxdxi[1, 2];
      jac23[2, 2] = 2 * (dxdxi[0, 0] * ddxdxi[1, 2] + dxdxi[1, 0] * ddxdxi[0, 2]) + dxdxi[0, 1] * ddxdxi[1, 0] + dxdxi[1, 1] * ddxdxi[0, 0];

      jac23[3, 0] = dxdxi[0, 0] * ddxdxi[0, 1] + 2 * dxdxi[0, 1] * ddxdxi[0, 2];
      jac23[3, 1] = dxdxi[1, 0] * ddxdxi[1, 1] + 2 * dxdxi[1, 1] * ddxdxi[1, 2];
      jac23[3, 2] = 2 * (dxdxi[0, 1] * ddxdxi[1, 2] + dxdxi[1, 1] * ddxdxi[0, 2]) + dxdxi[0, 0] * ddxdxi[1, 1] + dxdxi[1, 0] * ddxdxi[0, 1];

      DoubleMatrix term23 = MatrixFunctions.Product(ddNdX, jac23.Transpose());
      DoubleMatrix term3 = d3Ndxi - term13 - term23;

      DoubleMatrix jac13 = new DoubleMatrix(4, 4);
      jac13[0, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0];
      jac13[0, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
      jac13[0, 2] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0];
      jac13[0, 3] = 3 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0];

      jac13[1, 0] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1];
      jac13[1, 1] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
      jac13[1, 2] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1];
      jac13[1, 3] = 3 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1];

      jac13[2, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1];
      jac13[2, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1];
      jac13[2, 2] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0];
      jac13[2, 3] = dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1];

      jac13[3, 0] = dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1];
      jac13[3, 1] = dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1];
      jac13[3, 2] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1];
      jac13[3, 3] = dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1];

      DoubleMatrix d3NdX = MatrixFunctions.Product(MatrixFunctions.Inverse(jac13), term3.Transpose()).Transpose();

      DoubleMatrix jac14 = new DoubleMatrix(5, 5);
      jac14[0, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0];
      jac14[0, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
      jac14[0, 2] = 4 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0];
      jac14[0, 3] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0];
      jac14[0, 4] = 4 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];

      jac14[1, 0] = dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1];
      jac14[1, 1] = dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
      jac14[1, 2] = 4 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1];
      jac14[1, 3] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1];
      jac14[1, 4] = 4 * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];

      jac14[2, 0] = dxdxi[0, 1] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0];
      jac14[2, 1] = dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];
      jac14[2, 2] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] + dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1];
      jac14[2, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] + 3 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1];
      jac14[2, 4] = 3 * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 0] * dxdxi[1, 0] + dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 0];

      jac14[3, 0] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1];
      jac14[3, 1] = dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1];
      jac14[3, 2] = 2 * dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0];
      jac14[3, 3] = dxdxi[0, 0] * dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0];
      jac14[3, 4] = 2 * dxdxi[0, 0] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 0] * dxdxi[1, 1];

      jac14[4, 0] = dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1];
      jac14[4, 1] = dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];
      jac14[4, 2] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 1] + dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0];
      jac14[4, 3] = 3 * dxdxi[0, 0] * dxdxi[0, 1] * dxdxi[1, 1] * dxdxi[1, 1] + 3 * dxdxi[0, 1] * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1];
      jac14[4, 4] = 3 * dxdxi[0, 1] * dxdxi[1, 0] * dxdxi[1, 1] * dxdxi[1, 1] + dxdxi[0, 0] * dxdxi[1, 1] * dxdxi[1, 1] * dxdxi[1, 1];

      DoubleMatrix d4xdxi = Jacobian14At(d4Ndxi).Transpose();
      DoubleMatrix term14 = MatrixFunctions.Product(dNdx, d4xdxi/*.Transpose()*/);

      DoubleMatrix jac24 = new DoubleMatrix(5, 3);

      jac24[0, 0] = 3 * ddxdxi[0, 0] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * d3xdxi[0, 0];
      jac24[0, 1] = 3 * ddxdxi[1, 0] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * d3xdxi[1, 0];
      jac24[0, 2] = 4 * dxdxi[0, 0] * d3xdxi[1, 0] + 4 * dxdxi[1, 0] * d3xdxi[0, 0] + 6 * ddxdxi[0, 0] * ddxdxi[1, 0];

      jac24[1, 0] = 3 * ddxdxi[0, 1] * ddxdxi[0, 1] + 4 * dxdxi[0, 1] * d3xdxi[0, 1];
      jac24[1, 1] = 3 * ddxdxi[1, 1] * ddxdxi[1, 1] + 4 * dxdxi[1, 1] * d3xdxi[1, 1];
      jac24[1, 2] = 4 * dxdxi[0, 1] * d3xdxi[1, 1] + 4 * dxdxi[1, 1] * d3xdxi[0, 1] + 6 * ddxdxi[0, 1] * ddxdxi[1, 1];

      jac24[2, 0] = dxdxi[0, 1] * d3xdxi[0, 0] + 3 * ddxdxi[0, 0] * ddxdxi[0, 2] + 3 * dxdxi[0, 0] * d3xdxi[0, 2];
      jac24[2, 1] = dxdxi[1, 1] * d3xdxi[1, 0] + 3 * ddxdxi[1, 0] * ddxdxi[1, 2] + 3 * dxdxi[1, 0] * d3xdxi[1, 2];
      jac24[2, 2] = 3 * dxdxi[0, 0] * d3xdxi[1, 2] + 3 * dxdxi[1, 0] * d3xdxi[0, 2] + 3 * ddxdxi[0, 2] * ddxdxi[1, 0] + 3 * ddxdxi[0, 0] * ddxdxi[1, 2] + dxdxi[0, 1] * d3xdxi[1, 0] + dxdxi[1, 1] * d3xdxi[0, 0];

      jac24[3, 0] = ddxdxi[0, 0] * ddxdxi[0, 1] + 2 * ddxdxi[0, 2] * ddxdxi[0, 2] + 2 * dxdxi[0, 0] * d3xdxi[0, 3] + 2 * dxdxi[0, 1] * d3xdxi[0, 2];
      jac24[3, 1] = ddxdxi[1, 0] * ddxdxi[1, 1] + 2 * ddxdxi[1, 2] * ddxdxi[1, 2] + 2 * dxdxi[1, 0] * d3xdxi[1, 3] + 2 * dxdxi[1, 1] * d3xdxi[1, 2];
      jac24[3, 2] = 2 * dxdxi[0, 0] * d3xdxi[1, 3] + 2 * dxdxi[0, 1] * d3xdxi[1, 2] + 2 * dxdxi[1, 0] * d3xdxi[0, 3] + 2 * dxdxi[1, 1] * d3xdxi[0, 2] + 4 * ddxdxi[0, 2] * ddxdxi[1, 2] + ddxdxi[0, 0] * ddxdxi[1, 1] + ddxdxi[0, 1] * ddxdxi[1, 0];

      jac24[4, 0] = dxdxi[0, 0] * d3xdxi[0, 1] + 3 * ddxdxi[0, 1] * ddxdxi[0, 2] + 3 * dxdxi[0, 1] * d3xdxi[0, 3];
      jac24[4, 1] = dxdxi[1, 0] * d3xdxi[1, 1] + 3 * ddxdxi[1, 1] * ddxdxi[1, 2] + 3 * dxdxi[1, 1] * d3xdxi[1, 3];
      jac24[4, 2] = 3 * dxdxi[0, 1] * d3xdxi[1, 3] + 3 * dxdxi[1, 1] * d3xdxi[0, 3] + 3 * ddxdxi[0, 1] * ddxdxi[1, 2] + 3 * ddxdxi[0, 2] * ddxdxi[1, 1] + dxdxi[0, 0] * d3xdxi[1, 1] + dxdxi[1, 0] * d3xdxi[0, 1];

      DoubleMatrix term24 = MatrixFunctions.Product(ddNdX, jac24.Transpose());

      DoubleMatrix jac34 = new DoubleMatrix(5, 4);

      jac34[0, 0] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 0];
      jac34[0, 1] = 6 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 0];
      jac34[0, 2] = 6 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 0] + 12 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 0];
      jac34[0, 3] = 6 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 0] + 12 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[1, 0];

      jac34[1, 0] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 1];
      jac34[1, 1] = 6 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 1];
      jac34[1, 2] = 6 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 1] + 12 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 1];
      jac34[1, 3] = 6 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 1] + 12 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[1, 1];

      jac34[2, 0] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 2] + 3 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[0, 0];
      jac34[2, 1] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 2] + 3 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[1, 0];
      jac34[2, 2] = 3 * dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 2] + 3 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[1, 0] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 2] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 0] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 0];
      jac34[2, 3] = 3 * dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 2] + 3 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[0, 0] + 6 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[1, 2] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[1, 0] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[1, 0];

      jac34[3, 0] = dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[0, 2] + dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[0, 1];
      jac34[3, 1] = dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[1, 2] + dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[1, 1];
      jac34[3, 2] = dxdxi[0, 0] * dxdxi[0, 0] * ddxdxi[1, 1] + dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 0] + 2 * dxdxi[0, 0] * dxdxi[1, 0] * ddxdxi[0, 1] + 2 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 0] + 4 * dxdxi[0, 0] * dxdxi[0, 1] * ddxdxi[1, 2] + 4 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 2] + 4 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 2];
      jac34[3, 3] = dxdxi[1, 0] * dxdxi[1, 0] * ddxdxi[0, 1] + dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 0] + 2 * dxdxi[1, 0] * dxdxi[0, 0] * ddxdxi[1, 1] + 2 * dxdxi[1, 1] * dxdxi[0, 1] * ddxdxi[1, 0] + 4 * dxdxi[1, 0] * dxdxi[1, 1] * ddxdxi[0, 2] + 4 * dxdxi[1, 0] * dxdxi[0, 1] * ddxdxi[1, 2] + 4 * dxdxi[1, 1] * dxdxi[0, 0] * ddxdxi[1, 2];

      jac34[4, 0] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[0, 2] + 3 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[0, 1];
      jac34[4, 1] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[1, 2] + 3 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[1, 1];
      jac34[4, 2] = 3 * dxdxi[0, 1] * dxdxi[0, 1] * ddxdxi[1, 2] + 3 * dxdxi[0, 1] * dxdxi[0, 0] * ddxdxi[1, 1] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[0, 2] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[0, 1] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[0, 1];
      jac34[4, 3] = 3 * dxdxi[1, 1] * dxdxi[1, 1] * ddxdxi[0, 2] + 3 * dxdxi[1, 1] * dxdxi[1, 0] * ddxdxi[0, 1] + 6 * dxdxi[0, 1] * dxdxi[1, 1] * ddxdxi[1, 2] + 3 * dxdxi[0, 1] * dxdxi[1, 0] * ddxdxi[1, 1] + 3 * dxdxi[0, 0] * dxdxi[1, 1] * ddxdxi[1, 1];

      DoubleMatrix term34 = MatrixFunctions.Product(d3NdX, jac34.Transpose());

      DoubleMatrix term4 = d4Ndxi - term14 - term24 - term34;
      DoubleMatrix d4NdX = MatrixFunctions.Product(MatrixFunctions.Inverse(jac14), term4.Transpose()).Transpose();

      return d4NdX;//[nen,3]
    }

    public double[] GetCentriodCoordinate()
    {
      double[] centriod = new double[3];
      var d = GetNumberOfFields();
      if (volume == null)
      {
        if (element.GetPatch() is AbstractPatch2D)
        {
          int id = element.GetID();
          AbstractPatch2D patch = (AbstractPatch2D)element.GetPatch();
          NURBSSurface surface = patch.GetSurface();
          BivariateNURBSBasisFunction basis = (BivariateNURBSBasisFunction)((NURBSSurface)(patch.GetGeometry())).Basis;
          int p = basis.GetDegree(0);
          int q = basis.GetDegree(1);
          int nen = (p + 1) * (q + 1); // number of local basis functions
          var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
          var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
          int idx1 = patch.GetIPN(id, 0);
          int idx2 = patch.GetIPN(id, 1);

          double[] p1 = surface.PointAt(kvNoMulticiply1[idx1], kvNoMulticiply2[idx2]);
          double[] p2 = surface.PointAt(kvNoMulticiply1[idx1 + 1], kvNoMulticiply2[idx2]);
          double[] p3 = surface.PointAt(kvNoMulticiply1[idx1 + 1], kvNoMulticiply2[idx2 + 1]);
          double[] p4 = surface.PointAt(kvNoMulticiply1[idx1], kvNoMulticiply2[idx2 + 1]);

          centriod[0] = (p1[0] + p2[0] + p3[0] + p4[0]) / 4.0;
          centriod[1] = (p1[1] + p2[1] + p3[1] + p4[1]) / 4.0;
          centriod[2] = (p1[2] + p2[2] + p3[2] + p4[2]) / 4.0;
        }
      }
      return centriod;
    }
  }
}
