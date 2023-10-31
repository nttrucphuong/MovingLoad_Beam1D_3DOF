using CenterSpace.NMath.Core;
using DEMSoft.Common;
using DEMSoft.Function;

namespace DEMSoft.IGA
{
    public enum Structure2DState
    {
        PlaneStress, PlaneStrain, Axisymetric
    };
    /// <summary>
    /// Structure 2D problem (plane stress, plane strain). Default is plane stress
    /// </summary>
    public abstract class AbstractModelStructure : AbstractModel, IModelStructure
    {
        public Structure2DState StressState
        { get; set; }
        /// <summary>
        /// Constructor class
        /// </summary>
        /// <param name="problem">Define state of stress, plane stress or plane strain</param>
        public AbstractModelStructure(TypeAnalysisModel type, Dimension structureDimension, string pathProject, string nameProject)
                                  : base(TypeModelProblem.Structural, type, structureDimension, pathProject, nameProject)
        {
            if (StructureDimension == Dimension.Plane)
            {
                listComputeResult.Add(Result.UX);
                listComputeResult.Add(Result.UY);
            }
            else if (StructureDimension == Dimension.Solid)
            {
                listComputeResult.Add(Result.UX);
                listComputeResult.Add(Result.UY);
                listComputeResult.Add(Result.UZ);
            }
            else if (StructureDimension == Dimension.Plate)
            {
                listComputeResult.Add(Result.UX);
                listComputeResult.Add(Result.UY);
                listComputeResult.Add(Result.UZ);
                listComputeResult.Add(Result.THETAX);
                listComputeResult.Add(Result.THETAY);
            }
            else if (StructureDimension == Dimension.Beam)
            {
                listComputeResult.Add(Result.UX);
                listComputeResult.Add(Result.UY);
                listComputeResult.Add(Result.THETAZ);
            }

            typeOfFieldsInMultifield.Add(TypeFields.Structural);
        }
        public TypePlate TypePlate
        { get; set; }
        public TypeBeam TypeBeam
        { get; set; }
        public void SetThicknessPlate(double thickness, int indexPatch)
        {
            ((AbstractPatch2D)listPatch[indexPatch]).Thickness = thickness;
        }
        public void SetThicknessBeam(double thickness, int indexPatch)
        {
            ((AbstractPatch1D)listPatch[indexPatch]).Thickness = thickness;
        }
        public void SetKinematicsFunction(FunctionRToR fz, int indexPatch)
        {
            if (listPatch[indexPatch] is PatchStructurePlate)
                ((PatchStructurePlate)listPatch[indexPatch]).KinematicsFunction = fz;
            else if (listPatch[indexPatch] is PatchStructureBeam)
                ((PatchStructureBeam)listPatch[indexPatch]).KinematicsFunction = fz;
        }
        public void Givezz(double zz, int indexPatch)
        {
            if (listPatch[indexPatch] is PatchStructurePlate)
                ((PatchStructurePlate)listPatch[indexPatch]).Getzz = zz;
            else if (listPatch[indexPatch] is PatchStructureBeam)
                ((PatchStructureBeam)listPatch[indexPatch]).Getzz = zz;
        }
        public void SetTypeVF(TypeVFunction typeVF, int indexPatch)
        {
            if (listPatch[indexPatch] is PatchStructurePlate)
                ((PatchStructurePlate)listPatch[indexPatch]).GetTypeVF = typeVF;
            else if (listPatch[indexPatch] is PatchStructureBeam)
                ((PatchStructureBeam)listPatch[indexPatch]).GetTypeVF = typeVF;
        }

        public void AssemblyInternalForce(out DoubleVector ResidualGlobal)
        {
            ResidualGlobal = new DoubleVector(countDOF);
            for (int i = 0; i < listPatch.Count; i++)
            {
                AbstractPatchOneField pa = (AbstractPatchOneField)listPatch[i];
                switch (StructureDimension)
                {
                    case Dimension.Beam:
                        ((PatchStructure1D)listPatch[i]).ComputeInternalForcePatch(ref ResidualGlobal);
                        break;
                    case Dimension.Plane:
                        ((PatchStructure2D)listPatch[i]).ComputeInternalForcePatch(ref ResidualGlobal);
                        break;
                    case Dimension.Solid:
                        ((PatchStructure3D)listPatch[i]).ComputeInternalForcePatch(ref ResidualGlobal);
                        break;
                }
            }
        }
    }
}
