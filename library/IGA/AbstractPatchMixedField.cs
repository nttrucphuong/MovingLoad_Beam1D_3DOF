//namespace DEMSoft.IGA
//{
//    /// <summary>
//    /// Abstract of NURBS patch mesh
//    /// </summary>
//    public abstract class AbstractPatchMixedField : AbstractPatch
//    {
//        //protected override abstract void CreateINC();

//        //protected override abstract void CreateIEN();

//        //protected override abstract void CreateIPN();

//        //public override abstract int[] EnumerateInPatch();

//        //public override abstract int GetCountElement();

//        //public override abstract int GetCountLocalBasisFunctions(int idx);

//        //public override abstract int GetCountGlobalBasisFunctions(int idx);

//        //public override abstract void SetUGlobal(DoubleVector uGlobal);
        
//        //public override abstract int EnumerateInGlobalMultiPatch(int countDof);

//        public void FindIndexofIDArray(int enumerate, ref int indexGlobalCps, ref int indexField, ref int indexSurface)
//        {
//            for (int k = 0; k < 2; k++)
//                for (int i = 0; i < enumeratePatch[k].GetLength(0); i++)
//                {
//                    for (int j = 0; j < enumeratePatch[k].GetLength(1); j++)
//                    {
//                        if (enumeratePatch[k][i, j] == enumerate)
//                        {
//                            indexField = i;
//                            indexGlobalCps = j;
//                            indexSurface = k;
//                            return;
//                        }
//                    }
//                }
//        }
//    }
//}
