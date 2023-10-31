using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using DEMSoft.NURBS;

namespace DEMSoft.IGA
{
   /// <summary>
   /// Geometry of one element. Volume is geometry of 3D element
   /// </summary>
   public class Volume : AbstractMeshPart
   {
      /// <summary>
      /// On one coordinate of volume, we have two coordinate [xi(i), xi(i+1)], 0 - xi(i); 1 - xi(i+1)
      /// </summary>
      private int index;
      /// <summary>
      /// Six face of volume
      /// </summary>
      private Face[] face;
      /// <summary>
      /// Element which volume correspond
      /// </summary>
      private AbstractElement3D element;

      /// <summary>
      /// Constructor class
      /// </summary>
      /// <param name="element">3D element</param>
      public Volume(AbstractElement3D element)
      {
         this.element = element;
         var mesh = (AbstractPatchOneField)element.GetPatch();
         SetNumberOfFields(mesh.GetCountField());

         face = new Face[6];
         for (int i = 0; i < face.Length; i++)
         {
            face[i] = new Face(i, this);
         }
      }

      public override int[] GetTArrayGlobal()
      {
         var d = GetNumberOfFields();
         var patch = (AbstractPatchOneField)element.GetPatch();
         int nen = patch.GetCountLocalBasisFunctions();
         int[] tArray = new int[d * nen];
         var cps = ((NURBSVolume)(patch.GetGeometry())).ControlPoints;

         for (int k = 0; k < nen; k++)
         {
            int ien = patch.GetIEN(element.GetID(), k);
            var tArrayCps = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(ien, 2)].GetTArrayGlobal();
            for (int j = 0; j < d; j++)
            {
               tArray[k * d + j] = tArrayCps[j];
            }
         }
         return tArray;
      }


      //private int[] CreateTArrayOnVolume()
      //{
      //    var d = GetNumberOfFields();
      //    var patch = (AbstractPatchOne)element.GetPatch();
      //    int nen = patch.GetCountLocalBasisFunctions();
      //    int[] tArray = new int[d * nen];
      //    var cps = ((NURBSVolume)(patch.GetGeometry())).ControlPoints;

      //    for (int k = 0; k < nen; k++)
      //    {
      //        int ien = patch.GetIEN(element.GetID(), k);
      //        var tArrayCps = cps[patch.GetINC(ien, 0), patch.GetINC(ien, 1), patch.GetINC(ien, 2)].GetTArray();
      //        for (int j = 0; j < d; j++)
      //        {
      //            tArray[k * d + j] = tArrayCps[j];
      //        }
      //    }
      //    return tArray;
      //}

      /// <summary>
      /// Get element which volume coresspond
      /// </summary>
      /// <returns>3D element</returns>
      public AbstractElement3D GetElement()
      {
         return (AbstractElement3D)element;
      }

      /// <summary>
      /// Get control points on volume
      /// </summary>
      /// <returns>control points</returns>
      public ControlPoint[,,] GetControlPointsOnVolume()
      {
         var mesh = (AbstractPatchOneField)element.GetPatch();
         SetNumberOfFields(mesh.GetCountField());

         var basis = mesh.GetGeometry().Basis;
         int p = basis.GetDegree(0);
         int q = basis.GetDegree(1);
         int r = basis.GetDegree(2);
         int idPatch = element.GetID();
         int idxi = mesh.GetIPN(idPatch, 0);
         int idxj = mesh.GetIPN(idPatch, 1);
         int idxk = mesh.GetIPN(idPatch, 2);
         var kv1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
         var kv2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
         var kv3 = basis.GetKnotVector(2).GetKnotVectorNoMultiplicity();
         double mid1 = (kv1[idxi] + kv1[idxi + 1]) / 2.0;
         double mid2 = (kv2[idxj] + kv2[idxj + 1]) / 2.0;
         double mid3 = (kv3[idxk] + kv3[idxk + 1]) / 2.0;
         int spanU = basis.FindSpan(mid1, 0);
         int spanV = basis.FindSpan(mid2, 1);
         int spanW = basis.FindSpan(mid3, 2);

         var cpsLocal = new ControlPoint[p + 1, q + 1, r + 1];
         for (int i = 0; i < mesh.GetCountLocalBasisFunctions(); i++)
         {
            int ien = mesh.GetIEN(idPatch, i);
            int inc1 = mesh.GetINC(ien, 0);
            int inc2 = mesh.GetINC(ien, 1);
            int inc3 = mesh.GetINC(ien, 2);
            cpsLocal[inc1 - spanU + p, inc2 - spanV + q, inc3 - spanW + r] = ((Abstract3DParametricGeometry)(mesh.GetGeometry())).ControlPoints[inc1, inc2, inc3];
         }
         return cpsLocal;
      }

      /// <summary>
      /// Get face of volume
      /// </summary>
      /// <returns></returns>
      public Face GetFace(int index)
      {
         return face[index];
      }

      /// <summary>
      /// Get degree of volume on each direction
      /// </summary>
      /// <param name="indexDirection">Index of direction on three direction, 0 - xi, 1 - eta, 2 - zeta</param>
      /// <returns></returns>
      public int GetDegree(int indexDirection)
      {
         var mesh = (AbstractPatchOneField)element.GetPatch();
         var surface = (NURBSVolume)mesh.GetGeometry();
         var basis = surface.Basis;
         return basis.GetDegree(indexDirection);
      }

      public int[] GetCoordinateParameterOnFace(int index)
      {
         int indexCoordinateFace = face[index].GetIndexCoordinate();
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

      public int CountDof()
      {
         var patch = element.GetPatch();
         int d = patch.GetCountField(0);
         int nen = patch.GetCountLocalBasisFunctions(0);
         return nen * d;
      }
      public int GetIndexFrontBack()
      {
         return index % 2;
      }

      public int GetIndexCoordinate()
      {
         if (index == 0 || index == 1)
            return 2;
         else if (index == 2 || index == 3)
            return 1;
         else if (index == 4 || index == 5)
            return 0;
         else
            return -1;
      }
      public double[,,] GetTrivariateBasisFunctionOnVolume(double xi, double eta, double zeta)
      {
         var patch = GetElement().GetPatch();
         var vol = (NURBSVolume)patch.GetGeometry(0);
         var basis = (TrivariateNURBSBasisFunction)vol.Basis;
         int p = basis.GetDegree(0);
         int q = basis.GetDegree(1);
         int r = basis.GetDegree(2);
         int spanU = basis.FindSpan(xi, 0);
         int spanV = basis.FindSpan(eta, 1);
         int spanW = basis.FindSpan(zeta, 2);
         double[,,] Nijk = new double[p + 1, q + 1, r + 1];
         for (int i = 0; i <= p; i++)
            for (int j = 0; j <= q; j++)
               for (int k = 0; k <= r; k++)
               {
                  Nijk[i, j, k] = basis.GetValueTrivariateBasisFunction(xi, eta, zeta, spanU - p + i, spanV - q + j, spanW - r + k);//Nijk[i, j, 0];                           
               }
         return Nijk;
      }

      public double[] GetCentriodCoordinate()
      {
         double[] centriod = new double[3];
         var d = GetNumberOfFields();
         int id = element.GetID();
         AbstractPatch3D patch = (AbstractPatch3D)element.GetPatch();
         NURBSVolume vol = patch.GetVolume();
         TrivariateNURBSBasisFunction basis = (TrivariateNURBSBasisFunction)vol.Basis;
         int p = basis.GetDegree(0);
         int q = basis.GetDegree(1);
         int r = basis.GetDegree(2);
         int nen = (p + 1) * (q + 1) * (r + 1); // number of local basis functions
         var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
         var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
         var kvNoMulticiply3 = basis.GetKnotVector(2).GetKnotVectorNoMultiplicity();
         int idx1 = patch.GetIPN(id, 0);
         int idx2 = patch.GetIPN(id, 1);
         int idx3 = patch.GetIPN(id, 2);

         double[] p1 = vol.PointAt(kvNoMulticiply1[idx1], kvNoMulticiply2[idx2], kvNoMulticiply3[idx3]);
         double[] p2 = vol.PointAt(kvNoMulticiply1[idx1 + 1], kvNoMulticiply2[idx2], kvNoMulticiply3[idx3]);
         double[] p3 = vol.PointAt(kvNoMulticiply1[idx1 + 1], kvNoMulticiply2[idx2 + 1], kvNoMulticiply3[idx3]);
         double[] p4 = vol.PointAt(kvNoMulticiply1[idx1], kvNoMulticiply2[idx2 + 1], kvNoMulticiply3[idx3]);
         double[] p5 = vol.PointAt(kvNoMulticiply1[idx1], kvNoMulticiply2[idx2], kvNoMulticiply3[idx3 + 1]);
         double[] p6 = vol.PointAt(kvNoMulticiply1[idx1 + 1], kvNoMulticiply2[idx2], kvNoMulticiply3[idx3 + 1]);
         double[] p7 = vol.PointAt(kvNoMulticiply1[idx1 + 1], kvNoMulticiply2[idx2 + 1], kvNoMulticiply3[idx3 + 1]);
         double[] p8 = vol.PointAt(kvNoMulticiply1[idx1], kvNoMulticiply2[idx2 + 1], kvNoMulticiply3[idx3 + 1]);


         centriod[0] = (p1[0] + p2[0] + p3[0] + p4[0] + p5[0] + p6[0] + p7[0] + p8[0]) / 8.0;
         centriod[1] = (p1[1] + p2[1] + p3[1] + p4[1] + p5[1] + p6[1] + p7[1] + p8[1]) / 8.0;
         centriod[2] = (p1[1] + p2[1] + p3[1] + p4[1] + p5[2] + p6[2] + p7[2] + p8[2]) / 8.0;
         return centriod;
      }
   }
}
