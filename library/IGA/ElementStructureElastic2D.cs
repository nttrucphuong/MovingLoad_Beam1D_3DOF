using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public class ElementStructureElastic2D : AbstractElementStructure2D
  {
    public ElementStructureElastic2D(AbstractPatch2D patch, int id)
              : base(patch, id)
    {
    }

    public override void ComputeStiffnessMatrixElement(out DoubleMatrix Ke)
    {
      PatchStructure2D patch = (PatchStructure2D)this.patch;
      int d = patch.GetCountField();
      BivariateNURBSBasisFunction basis = (BivariateNURBSBasisFunction)((NURBSSurface)(patch.GetGeometry())).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int nen = (p + 1) * (q + 1); // number of local basis functions
      var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      int idx1 = patch.GetIPN(id, 0);
      int idx2 = patch.GetIPN(id, 1);
      Ke = new DoubleMatrix(d * nen, d * nen);
      DoubleMatrix B = new DoubleMatrix(4, d * nen);
      DoubleMatrix D = null;
      if (!(Material is FGMUserDefinedGradedMaterial))
      {
        D = CreateMaterialMatrix();
      }
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          int c = 0;
          double x1 = 1;
          double a = 1;
          GaussPoints gpsij = gps[i, j];
          if (patch.StateStress == Structure2DState.Axisymetric)
          {
            c = 1;
            x1 = PointAt(gpsij.location[0], gpsij.location[1])[0];
            a = 2 * Math.PI * x1;
          }
          else if (patch.StateStress == Structure2DState.PlaneStress)
          {
            a = patch.Thickness;
            if (a == 0)
              throw new InvalidArgumentException("Plane stress condition must had thickness");
          }
          double detJbar = 1.0 / 4.0
                  * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                  * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          DoubleMatrix dNdx = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.dNdX);
          DoubleVector N = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);
          for (int k = 0; k < nen; k++)
          {
            B[0, 2 * k] = dNdx[k, 0];
            B[2, 2 * k] = c * N[k] / x1;
            B[3, 2 * k] = dNdx[k, 1];
            B[1, 2 * k + 1] = dNdx[k, 1];
            B[3, 2 * k + 1] = dNdx[k, 0];
          }
          if (Material is FGMUserDefinedGradedMaterial)
          {
            D = CreateMaterialMatrix(i, j);
          }
          DoubleMatrix BTDB = MatrixFunctions.TransposeProduct(B, MatrixFunctions.Product(D, B));
          Ke += a * gpsij.weight * detJbar * Math.Abs((double)gpsij.GetValue(DataInGausspoint.detJ)) * BTDB;
        }
      }
    }

    public override DoubleVector StrainAt(params double[] xi)
    {
      PatchStructure2D patch = (PatchStructure2D)this.patch;
      int d = patch.GetCountDimension();
      int nen = patch.GetCountLocalBasisFunctions();
      var cps = ((NURBSSurface)(patch.GetGeometry())).ControlPoints;
      DoubleMatrix B = new DoubleMatrix(4, d * nen);
      DoubleMatrix gradNE = GradBasisFunction(xi[0], xi[1]);
      DoubleMatrix J = JacobianAt(gradNE);
      DoubleMatrix invertJ = MatrixFunctions.Inverse(J);
      DoubleMatrix gradBasis = MatrixFunctions.Product(gradNE, invertJ);
      DoubleVector N = ValueBasisFunction(xi[0], xi[1]);
      int c = 0;
      double x1 = 1;
      if (patch.StateStress == Structure2DState.Axisymetric)
      {
        c = 1;
        x1 = PointAt(xi[0], xi[1])[0];
      }
      for (int k = 0; k < nen; k++)
      {
        B[0, 2 * k] = gradBasis[k, 0];
        B[1, 2 * k + 1] = gradBasis[k, 1];
        B[2, 2 * k] = c * N[k] / x1;
        B[3, 2 * k] = gradBasis[k, 1];
        B[3, 2 * k + 1] = gradBasis[k, 0];
      }
      DoubleVector U = GetDisplacementLocal();
      DoubleVector strainTemp = MatrixFunctions.Product(B, U);
      //strainTemp[3] /= 2.0;
      if (patch.StateStress == Structure2DState.PlaneStress)
      {
        double nu = ComputeParameterProperty(MaterialPropertyName.PoissonRatio, xi[0], xi[1]);
        strainTemp[2] = -nu / (1.0 - nu) * (strainTemp[0] + strainTemp[1]);
      }
      return strainTemp;
    }

    public override DoubleVector StressAt(params double[] xi)
    {
      PatchStructure2D patch = (PatchStructure2D)this.patch;
      DoubleMatrix D = null;
      if (!(Material is FGMUserDefinedGradedMaterial))
        D = CreateMaterialMatrix();
      else
        D = CreateMaterialMatrix(xi[0], xi[1]);
      DoubleVector strain = StrainAt(xi);
      //strain[3] *= 2.0;
      DoubleVector stress = MatrixFunctions.Product(D, strain);
      if (patch.StateStress == Structure2DState.PlaneStrain)
      {
        double nu = ComputeParameterProperty(MaterialPropertyName.PoissonRatio, xi[0], xi[1]);
        stress[2] = nu * (stress[0] + stress[1]);
      }
      return stress;
    }


    //public override void UpdateStressGaussPoint()
    //{
    //	for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
    //		for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
    //		{
    //			GaussPoints gauss = gps[i, j];
    //			double[] locGauss = gauss.location;
    //			gauss.lastStress = StressAt(locGauss[0], locGauss[1]);
    //		}
    //}

    public override void ComputeInternalForceElement(out DoubleVector fi)
    {
      PatchStructure2D patch = (PatchStructure2D)this.patch;
      int d = patch.GetCountField();
      BivariateNURBSBasisFunction basis = (BivariateNURBSBasisFunction)((NURBSSurface)(patch.GetGeometry(0))).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int nen = (p + 1) * (q + 1); // number of local basis functions
      var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      int idx1 = patch.GetIPN(id, 0);
      int idx2 = patch.GetIPN(id, 1);
      fi = new DoubleVector(d * nen);
      //DoubleVector fiU = fiU = new DoubleVector(2 * nen);
      DoubleMatrix Bu = new DoubleMatrix(4, 2 * nen);

      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          GaussPoints gpsij = gps[i, j];
          int c = 0;
          double x1 = 1;
          double a = 1;
          if (patch.StateStress == Structure2DState.Axisymetric)
          {
            c = 1;
            x1 = PointAt(gpsij.location[0], gpsij.location[1])[0];
            a = 2 * Math.PI * x1;
          }
          else if (patch.StateStress == Structure2DState.PlaneStress)
          {
            a = patch.Thickness;
            if (a == 0)
              throw new InvalidArgumentException("Plane stress condition must had thickness");
          }
          double detJbar = 1.0 / 4.0 * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1]) * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          DoubleMatrix dNdx = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.dNdX);
          DoubleVector N = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);
          double detJ = Math.Abs((double)gpsij.GetValue(DataInGausspoint.detJ));
          for (int k = 0; k < nen; k++)
          {
            //xx yy zz xy
            Bu[0, 2 * k] = dNdx[k, 0];
            Bu[2, 2 * k] = c * N[k] / x1;
            Bu[3, 2 * k] = dNdx[k, 1];
            Bu[1, 2 * k + 1] = dNdx[k, 1];
            Bu[3, 2 * k + 1] = dNdx[k, 0];
          }

          DoubleVector stress = StressAt(gpsij.location);
          fi += a * gpsij.weight * detJbar * detJ * NMathFunctions.Product(Bu.Transpose(), stress);
        }
      }

      //for (int i = 0; i < nen; i++)
      //{
      //  for (int ii = 0; ii < 2; ii++)//x y
      //  {
      //    fi[i * 3 + ii] = fiU[i * 2 + ii];
      //  }
      //}
    }

    private double Ee;
    private double xPhys;

    public void SetCurrentModulus(double Ee)
    {
      this.Ee = Ee;
    }

    public double GetCurrentModulus()
    {
      return Ee;
    }

    public void SetDensityFilter(double xPhys)
    {
      this.xPhys = xPhys;
    }

    public double GetDensityFilter()
    {
      return xPhys;
    }

    public override DoubleVector StrainElasticAt(params double[] xi)
    {
      throw new NotImplementedException();
    }

    public override DoubleVector StrainThermoAt(params double[] xi)
    {
      throw new NotImplementedException();
    }


  }
}
