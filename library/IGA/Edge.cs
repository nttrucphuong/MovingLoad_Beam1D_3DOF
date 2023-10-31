using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms.VisualStyles;
using DEMSoft.NURBS;

namespace DEMSoft.IGA
{
    public class Edge : AbstractMeshPart
    {
        private int index;
        private Face face;
        private AbstractElement element;
        public Edge(AbstractElement element)
        {
            this.element = element;
            var patch = element.GetPatch();
            if (patch is AbstractPatchOneField)
            {
                SetNumberOfFields(((AbstractPatchOneField)patch).GetCountField());
                SetDimension(1);
            }
            else
            {
                SetNumberOfFields(patch.GetCountField(0) + patch.GetCountField(1));
                SetDimension(1);
            }
        }

        /// <summary>
        /// Constructor class
        /// </summary>
        /// <param name="index">0 - first edge on u direction, 1 - end edge on u direction, 2 - first edge on v direction, 3 - end edge on v direction </param>
        /// <param name="face">Face to be attach edge</param>
        public Edge(int index, Face face)
        {
            this.index = index;
            this.face = face;
            SetNumberOfFields(face.GetNumberOfFields());
            SetDimension(face.CountDimension());
        }

        public int GetDegree(int idx)
        {
            var mesh = face.GetElement().GetPatch();
            var surface = (NURBSSurface)mesh.GetGeometry(idx);
            var basis = surface.Basis;
            return basis.GetDegree(index / 2);
        }

        public int GetIndexFrontBack()
        { return index % 2; }

        public Face GetFace()
        { return face; }

        public ControlPoint[] GetControlPointsOnEdge(int idx = 0)
        {
            var mesh = face.GetElement().GetPatch();
            var surface = (NURBSSurface)mesh.GetGeometry(idx);
            var basis = surface.Basis;
            int p = basis.GetDegree(0);
            int q = basis.GetDegree(1);

            int numCps = basis.GetDegree(GetIndexCoordinate()) + 1;
            ControlPoint[] cps = new ControlPoint[numCps];
            var cpsPatch = face.GetControlPointsOnFace();
            int indexFrontBack = index % 2;
            for (int i = 0; i < numCps; i++)
            {
                switch (GetIndexCoordinate())
                {
                    case 0:
                        if (indexFrontBack == 0)
                            cps[i] = cpsPatch[i, 0];
                        else
                            cps[i] = cpsPatch[i, q];
                        break;
                    case 1:
                        if (indexFrontBack == 0)
                            cps[i] = cpsPatch[0, i];
                        else
                            cps[i] = cpsPatch[p, i];
                        break;

                }
            }
            return cps;
        }

        public int CountDofOnEdge(int idx = 0)
        {
            var mesh = face.GetElement().GetPatch();
            var surface = (NURBSSurface)mesh.GetGeometry(idx);
            var basis = surface.Basis;

            return mesh.GetCountField(0) * (basis.GetDegree(GetIndexCoordinate()) + 1);
        }

        //////////////// Tien
        public int CountControlPointOnEdge(int idx = 0)
        {
            var mesh = face.GetElement().GetPatch();
            var surface = (NURBSSurface)mesh.GetGeometry(idx);
            var basis = surface.Basis;

            return basis.GetDegree(GetIndexCoordinate()) + 1;
        }
        //////////////
        //private int[] CreateTArrayOnEdge(int idx = 0)
        //{
        //    var d = GetNumberOfFields();
        //    var mesh = face.GetElement().GetPatch();


        //    var surface = (NURBSSurface)mesh.GetGeometry(idx);
        //    var basis = surface.Basis;
        //    int p = basis.GetDegree(0);
        //    int q = basis.GetDegree(1);
        //    int indexCoordinate = index / GetNumberOfFields();

        //    int numCps = basis.GetDegree(indexCoordinate) + 1;
        //    int[] tArray = new int[d * numCps];
        //    var cpsEdge = GetControlPointsOnEdge(idx);

        //    for (int i = 0; i < numCps; i++)
        //    {
        //        for (int j = 0; j < d; j++)
        //        {
        //            tArray[i * d + j] = cpsEdge[i].GetTArray()[j];
        //        }
        //    }
        //    return tArray;
        //}

        //private int[] CreateTArrayGlobalOnEdge()
        //{
        //    var mesh = face.GetElement().GetPatch();
        //    int d = face.GetElement().GetPatch().GetCountField(0);
        //    var surface = (NURBSSurface)mesh.GetGeometry(0);
        //    var basis = surface.Basis;
        //    int p = basis.GetDegree(0);
        //    int q = basis.GetDegree(1);
        //    int numCps = basis.GetDegree(GetIndexCoordinate()) + 1;
        //    int[] tArray = new int[d * numCps];
        //    var cpsEdge = GetControlPointsOnEdge(0);

        //    for (int i = 0; i < numCps; i++)
        //    {
        //        for (int j = 0; j < d; j++)
        //        {
        //            tArray[i * d + j] = cpsEdge[i].GetTArrayGlobal()[j];
        //        }
        //    }
        //    return tArray;
        //}

        public override int[] GetTArrayGlobal()
        {
            int[] tArray = null;
            if (face != null)
            {
                var mesh = face.GetElement().GetPatch();
                int d = face.GetElement().GetPatch().GetCountField(0);
                var surface = (NURBSSurface)mesh.GetGeometry(0);
                var basis = surface.Basis;
                int p = basis.GetDegree(0);
                int q = basis.GetDegree(1);
                int numCps = basis.GetDegree(GetIndexCoordinate()) + 1;
                tArray = new int[d * numCps];
                var cpsEdge = GetControlPointsOnEdge(0);

                for (int i = 0; i < numCps; i++)
                {
                    for (int j = 0; j < d; j++)
                    {
                        tArray[i * d + j] = cpsEdge[i].GetTArrayGlobal()[j];
                    }
                }
            }
            else
            {
                if (element.GetPatch() is AbstractPatch1D)
                {
                    AbstractPatch1D patch = (AbstractPatch1D)element.GetPatch();
                    int nen = patch.GetCountLocalBasisFunctions();
                    int d = GetNumberOfFields();
                    tArray = new int[d * nen];
                    var cps = patch.GetCurve().ControlPoints;
                    for (int kk = 0; kk < nen; kk++)
                    {
                        int ien = patch.GetIEN(0, element.GetID(), kk);
                        var tArrayCps = cps[patch.GetINC(0, ien, 0)].GetTArrayGlobal();
                        for (int j = 0; j < d; j++)
                        {
                            tArray[kk * d + j] = tArrayCps[j];
                        }
                    }
                }
            }
            return tArray;
        }

        public override int[] GetTArray()
        {
            var mesh = face.GetElement().GetPatch();
            int d = face.GetElement().GetPatch().GetCountField(0);
            var surface = (NURBSSurface)mesh.GetGeometry(0);
            var basis = surface.Basis;
            int p = basis.GetDegree(0);
            int numCps = basis.GetDegree(GetIndexCoordinate()) + 1;
            int[] tArray = new int[d * numCps];
            var cpsEdge = GetControlPointsOnEdge(0);

            for (int i = 0; i < numCps; i++)
            {
                for (int j = 0; j < d; j++)
                {
                    tArray[i * d + j] = cpsEdge[i].GetTArray()[j];
                }
            }
            return tArray;
        }

        public double[] GetParametricEndEdge(int idx = 0)
        {
            var mesh = face.GetElement().GetPatch();
            var surface = (NURBSSurface)mesh.GetGeometry(idx);
            var basis = surface.Basis;


            int idxi = mesh.GetIPN(face.GetElement().GetID(), GetIndexCoordinate());
            var kv1 = basis.GetKnotVector(GetIndexCoordinate()).GetKnotVectorNoMultiplicity();
            return new double[] { kv1[idxi], kv1[idxi + 1] };
        }

        public double[] GetBivariateBasisFunctionOnEdge(int idx, double xi, double eta)
        {
            var patch = face.GetElement().GetPatch();
            var surface = (NURBSSurface)patch.GetGeometry(idx);
            var basis = (BivariateNURBSBasisFunction)surface.Basis;
            int p = basis.GetDegree(0);
            int q = basis.GetDegree(1);

            int numCps = basis.GetDegree(GetIndexCoordinate()) + 1;

            var Nij = basis.GetValueBivariateBasisFunctions(xi, eta);
            double[] Ni = new double[numCps];
            int indexFrontBack = index % 2;
            for (int i = 0; i < numCps; i++)
            {
                switch (GetIndexCoordinate())
                {
                    case 0:
                        if (indexFrontBack == 0)
                            Ni[i] = Nij[i, 0];
                        else
                            Ni[i] = Nij[i, q];
                        break;
                    case 1:
                        if (indexFrontBack == 0)
                            Ni[i] = Nij[0, i];
                        else
                            Ni[i] = Nij[p, i];
                        break;
                }
            }
            return Ni;
        }

        public double[] GetDerivativeBivariateBasisFunctionOnEdge(int idx, double xi, double eta)
        {
            var mesh = face.GetElement().GetPatch();
            var surface = (NURBSSurface)mesh.GetGeometry(idx);
            var basis = (BivariateNURBSBasisFunction)surface.Basis;
            int p = basis.GetDegree(0);
            int q = basis.GetDegree(1);
            int numCps = basis.GetDegree(GetIndexCoordinate()) + 1;
            var dNij = basis.GetDerivativeBivariateBasisFunctions(xi, eta, 1);
            double[] dNi = new double[numCps];
            int indexFrontBack = index % 2;
            for (int i = 0; i < numCps; i++)
            {
                switch (GetIndexCoordinate())
                {
                    case 0:
                        if (indexFrontBack == 0)
                            dNi[i] = dNij[i, 0][1, 0];
                        else
                            dNi[i] = dNij[i, q][1, 0];
                        break;
                    case 1:
                        if (indexFrontBack == 0)
                            dNi[i] = dNij[0, i][0, 1];
                        else
                            dNi[i] = dNij[p, i][0, 1];
                        break;
                }
            }
            return dNi;
        }

        public int GetIndexCoordinate()
        {
            return index / 2;
        }
    }
}
