using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public interface ICoupling
  {
    void MakeCoupling(ref SparseMatrixBuilder<double> T);
    void MakeCoupling(ref DoubleMatrix T);
  }
}
