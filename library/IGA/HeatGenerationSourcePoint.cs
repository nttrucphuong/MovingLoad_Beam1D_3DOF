using System.Linq;
using DEMSoft.NURBS;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
    public class HeatGenerationSourcePoint : HeatGenerationSourceBody
    {
        private double[] xi;
        public HeatGenerationSourcePoint(AbstractPatch patch, double G, bool isOnPhysicsCoordination, params double[] x)
            : base(patch, G)
        {
            if (isOnPhysicsCoordination)
            {
                if (patch is AbstractPatch2D)
                {
                    NURBSSurface surface = ((AbstractPatch2D)patch).GetSurface();
                    xi = surface.Projection(x[0], x[1], 0);
                }
                else if (patch is AbstractPatch3D)
                {
                    NURBSVolume volume = ((AbstractPatch3D)patch).GetVolume();
                    xi = volume.Projection(x[0], x[1], x[2]);
                }
            }
            else
            {
                xi = x;
            }
        }

        public override void ComputeLocalLoadVector(ref DoubleVector rGlobal)
        {
            AbstractPatch patch = GetPatch();
            double G = GetValueHeatGenerationSource();
            int numElement = ((AbstractPatchOneField)patch).FindIndexOfElementAt(xi);
            if (patch is AbstractPatch2D)
            {
                ElementThermal2D element = (ElementThermal2D)patch.GetElement(numElement);
                BivariateNURBSBasisFunction basic = (BivariateNURBSBasisFunction)(((AbstractPatch2D)patch).GetSurface().Basis);
                Face face = element.GetFace();
                int countDOF = face.CountDof();
                DoubleVector re = new DoubleVector(countDOF);

                double[,] Nij = face.GetBivariateBasisFunctionOnFace(xi[0], xi[1]);
                int count = 0;
                for (int j = 0; j < Nij.GetLength(1); j++)
                    for (int i = 0; i < Nij.GetLength(0); i++)
                    {
                        re[count++] = G * Nij[i, j];
                    }
                int[] tArrayGlobal = face.GetTArrayGlobal();
                for (int i = 0; i < re.Length; i++)
                {
                    rGlobal[tArrayGlobal[i]] += re[i];
                }
            }
            else if (patch is AbstractPatch3D)
            {
                ElementThermal3D element = (ElementThermal3D)patch.GetElement(numElement);
                TrivariateNURBSBasisFunction basic = (TrivariateNURBSBasisFunction)(((AbstractPatch3D)patch).GetVolume().Basis);
                Volume volume = element.GetVolume();
                int countDOF = volume.CountDof();
                DoubleVector re = new DoubleVector(countDOF);

                double[,,] Nijk = volume.GetTrivariateBasisFunctionOnVolume(xi[0], xi[1], xi[2]);
                int count = 0;
                for (int k = 0; k < Nijk.GetLength(2); k++)
                    for (int j = 0; j < Nijk.GetLength(1); j++)
                        for (int i = 0; i < Nijk.GetLength(0); i++)
                        {
                            re[count++] = G * Nijk[i, j, k];
                        }
                int[] tArrayGlobal = volume.GetTArrayGlobal();
                for (int i = 0; i < re.Length; i++)
                {
                    rGlobal[tArrayGlobal[i]] += re[i];
                }
            }
        }
    }
}
