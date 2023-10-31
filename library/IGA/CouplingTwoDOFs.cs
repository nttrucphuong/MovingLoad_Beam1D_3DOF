using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public class CouplingTwoDOFs : ICoupling
  {
    private int masterDOF;
    private int slaveDOF;

    public CouplingTwoDOFs(int masterDOF, int slaveDOF)
    {
      if (masterDOF == slaveDOF)
        throw new InvalidArgumentException("Two degree of freedoms must be difference");
      this.masterDOF = masterDOF;
      this.slaveDOF = slaveDOF;
    }

    public void MakeCoupling(ref SparseMatrixBuilder<double> T)
    {
      T[slaveDOF, slaveDOF] += -1;
      T[slaveDOF, masterDOF] += 1;
    }

    public void MakeCoupling(ref DoubleMatrix T)
    {
      T[slaveDOF, slaveDOF] += -1;
      T[slaveDOF, masterDOF] += 1;
    }
  }
}
