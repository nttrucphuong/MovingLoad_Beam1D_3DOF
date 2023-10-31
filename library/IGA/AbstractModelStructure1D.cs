using CenterSpace.NMath.Core;
using DEMSoft.Common;
using DEMSoft.Function;

namespace DEMSoft.IGA
{
    public enum Structure1DState
    {
        Truss, Beam
    };
    /// <summary>
    /// Structure 1D problem (plane stress, plane strain). Default is plane stress
    /// </summary>
    public abstract class AbstractModelStructure1D : AbstractModel
    {
        /// <summary>
        /// Constructor class
        /// </summary>
        /// <param name="problem">Define state of stress, plane stress or plane strain</param>
        public AbstractModelStructure1D(TypeAnalysisModel type, Dimension structureDimension, string pathProject, string nameProject)
                                  : base(TypeModelProblem.Structural, type, structureDimension, pathProject, nameProject)
        {
            listComputeResult.Add(Result.UX);
            listComputeResult.Add(Result.UY);
            listComputeResult.Add(Result.THETAZ);
            typeOfFieldsInMultifield.Add(TypeFields.Structural);
        }
        public TypeBeam TypeBeam
        { get; set; }
        public void SetKinematicsFunction(FunctionRToR fz, int indexPatch)
        {
            ((PatchStructureBeam)listPatch[indexPatch]).KinematicsFunction = fz;
        }
        public void Givezz(double zz, int indexPatch)
        {
            ((PatchStructureBeam)listPatch[indexPatch]).Getzz = zz;
        }
        public void SetTypeVF(TypeVFunction typeVF, int indexPatch)
        {
            ((PatchStructureBeam)listPatch[indexPatch]).GetTypeVF = typeVF;
        }

        public void AssemblyInternalForce(out DoubleVector ResidualGlobal)
        {
            ResidualGlobal = new DoubleVector(countDOF);
            for (int i = 0; i < listPatch.Count; i++)
            {
                ((PatchStructure1D)listPatch[i]).ComputeInternalForcePatch(ref ResidualGlobal);
            }
        }
    }
}
