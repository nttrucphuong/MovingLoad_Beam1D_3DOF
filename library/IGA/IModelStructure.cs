namespace DEMSoft.IGA
{
  interface IModelStructure
    {
        //void AssemblyMassMatrix(out SparseMatrixBuilder<double> mGlobal);
        Structure2DState StressState
        { get; set; }
    }
}
