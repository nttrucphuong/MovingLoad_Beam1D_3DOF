using CenterSpace.NMath.Core;
using System.Collections.Generic;

namespace DEMSoft.IGA
{
  public enum DataInGausspoint
  {
    Ni,
    dNdxi,
    ddNdxi,
    dNdX,
    ddNdX,
    d3NdX,
    d4NdX,
    detJ,
    thickness,
    NBernsteini,
    dNBernsteindxi,
    currentStress,
    lastStress,
    currentBackStress,
    lastBackStress,
    currentPlasticStrain,
    lastPlasticStrain,
    currentAlpha,
    lastAlpha,
    currentStrain,
    lastStrain,
    lastDep,
    currentDep,
    EModulus,
    nu,
    Density,
    NonLocalParameter,
    GeneralNonLocalParameter1,
    GeneralNonLocalParameter2,
    ThermalConductivity,
    CoefficientsThermalExpansion,
    MaterialPropertyValue,
    previousPhase,
    currentPhase,
    lastPhase,
    xiEpsilon,
    DrawValue
  }
  public class DataGausspoint
  {
    public DataInGausspoint typeData { get; set; }
    public object value { get; set; }
  }
  //public class DoubleDataGausspoint : DataGausspoint
  //{
  //  public double value { get; set; }
  //}
  //public class DoubleMatrixDataGausspoint : DataGausspoint
  //{
  //  public DoubleMatrix value { get; set; }
  //}
  //public class DoubleVectorDataGausspoint : DataGausspoint
  //{
  //  public DoubleVector value { get; set; }
  //}
  public class GaussPoints : DEMSoft.Function.GaussPoints
  {
    public List<DataGausspoint> listDataGausspoint;

    //public DoubleMatrix dNdxi;
    //public DoubleMatrix dNdX;
    //public DoubleMatrix ddNdxi;
    //public DoubleMatrix ddNdX;
    //public DoubleVector Ni;
    //public double detJ;

    //public DoubleVector NBernsteini;
    //public DoubleMatrix dNBernsteindxi;

    //public DoubleVector currentStress;
    //public DoubleVector lastStress;
    //public DoubleVector currentBackStress;
    //public DoubleVector lastBackStress;
    //public DoubleVector currentPlasticStrain;
    //public DoubleVector lastPlasticStrain;
    //public double currentAlpha;
    //public double lastAlpha;
    //public DoubleVector currentStrain;
    //public DoubleVector lastStrain;

    //public DoubleMatrix lastDep;
    //public DoubleMatrix currentDep;

    //public double EModulus;
    //public double nu;
    //public double ThermalConductivity;
    //public double coefficientsthermalexpansion;
    //public double materialPropertyValue;
    //public double previousPhase;
    //public double currentPhase;
    //public double lastPhase;
    ///// <summary>
    ///// Maximum positive reference energy
    ///// </summary>
    //public double xiEpsilon;
    //public double DrawValue;
    //public enum Variable
    //{
    //  None,
    //  xiEpsilon,
    //  currentPhase,
    //  currentStress,
    //  currentStrain
    //}

    public GaussPoints(double[] location, double w) : base(location, w)
    {
    }

    public void UpdateConvergedVariablesLoadStep()
    {
      if (AbstractModel.TypeModel == TypeModelProblem.PhaseField)
        UpdateConvergedPhaseField();
      else
        UpdateConvergedPlasticity();

    }

    private void UpdateConvergedPlasticity()
    {
      listDataGausspoint.Find(x => x.typeData == DataInGausspoint.lastStress).value = 1.0 * (DoubleVector)listDataGausspoint.Find(x => x.typeData == DataInGausspoint.currentStress).value;
      listDataGausspoint.Find(x => x.typeData == DataInGausspoint.lastPlasticStrain).value = 1.0 * (DoubleVector)listDataGausspoint.Find(x => x.typeData == DataInGausspoint.currentPlasticStrain).value;
      listDataGausspoint.Find(x => x.typeData == DataInGausspoint.lastBackStress).value = 1.0 * (DoubleVector)listDataGausspoint.Find(x => x.typeData == DataInGausspoint.currentBackStress).value;
      listDataGausspoint.Find(x => x.typeData == DataInGausspoint.lastAlpha).value = 1.0 * (double)listDataGausspoint.Find(x => x.typeData == DataInGausspoint.currentAlpha).value;
      listDataGausspoint.Find(x => x.typeData == DataInGausspoint.lastDep).value = 1.0 * (DoubleVector)listDataGausspoint.Find(x => x.typeData == DataInGausspoint.currentDep).value;
      listDataGausspoint.Find(x => x.typeData == DataInGausspoint.lastStrain).value = 1.0 * (DoubleVector)listDataGausspoint.Find(x => x.typeData == DataInGausspoint.currentStrain).value;

      //lastStress = 1.0 * currentStress;
      //lastPlasticStrain = 1.0 * currentPlasticStrain;
      //lastBackStress = 1.0 * currentBackStress;
      //lastAlpha = currentAlpha;
      //lastDep = 1.0 * currentDep;
      //lastStrain = 1.0 * currentStrain;
    }

    private void UpdateConvergedPhaseField()
    {
      //listDataGausspoint.Find(x => x.typeData == DataInGausspoint.lastPhase).value = 1.0 * (double)listDataGausspoint.Find(x => x.typeData == DataInGausspoint.currentPhase).value;
      SetValue(DataInGausspoint.lastPhase, GetValue(DataInGausspoint.currentPhase));
      SetValue(DataInGausspoint.lastStrain, GetValue(DataInGausspoint.currentStrain));
      SetValue(DataInGausspoint.lastStress, GetValue(DataInGausspoint.currentStress));
      ////lastStress = 1.0 * currentStress;
      ////lastStrain = 1.0 * currentStrain;
      ////previousLastPhase = lastPhase;
      //lastPhase = currentPhase;
    }
    /// <summary>
    /// Find item in list of data of gausspoint
    /// </summary>
    /// <param name="type">Type of variable data in gausspoint</param>
    /// <returns></returns>
    public DataGausspoint Find(DataInGausspoint type)
    {
      return listDataGausspoint.Find(x => x.typeData == type);
    }
    public object GetValue(DataInGausspoint type)
    {
      DataGausspoint data = listDataGausspoint.Find(x => x.typeData == type);
      if (data != null)
        return data.value;
      else
        return null;
    }
    public void SetValue(DataInGausspoint type, object value)
    {
      DataGausspoint data = listDataGausspoint.Find(x => x.typeData == type);
      if (data == null)
      {
        data = new DataGausspoint();
        data.typeData = type;
        data.value = value;
        listDataGausspoint.Add(data);
      }
      else
      {
        data.value = value;
      }
    }
  }
}