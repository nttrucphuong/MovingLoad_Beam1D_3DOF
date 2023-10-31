using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;
using DEMSoft.Function;
using System.Linq;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
  public class ElementStructurePhaseField2D : AbstractElementStructure2D
  {
    public ElementStructurePhaseField2D(AbstractPatch2D patch, int id)
        : base(patch, id)
    {
    }

    private DoubleMatrix JTensor()
    {
      DoubleMatrix J = new DoubleMatrix(6, 6);
      J[0, 0] = ItemJTensor(1, 1, 1, 1);
      J[0, 1] = J[1, 0] = ItemJTensor(1, 1, 2, 2);
      J[0, 2] = J[2, 0] = ItemJTensor(1, 1, 3, 3);
      J[0, 3] = J[3, 0] = ItemJTensor(1, 1, 2, 3);
      J[0, 4] = J[4, 0] = ItemJTensor(1, 1, 1, 3);
      J[0, 5] = J[5, 0] = ItemJTensor(1, 1, 1, 2);
      J[1, 1] = ItemJTensor(2, 2, 2, 2);
      J[1, 2] = J[2, 1] = ItemJTensor(2, 2, 3, 3);
      J[1, 3] = J[3, 1] = ItemJTensor(2, 2, 2, 3);
      J[1, 4] = J[4, 1] = ItemJTensor(2, 2, 1, 3);
      J[1, 5] = J[5, 1] = ItemJTensor(2, 2, 1, 2);
      J[2, 2] = ItemJTensor(3, 3, 3, 3);
      J[2, 3] = J[3, 2] = ItemJTensor(3, 3, 2, 3);
      J[2, 4] = J[4, 2] = ItemJTensor(3, 3, 1, 3);
      J[2, 5] = J[5, 2] = ItemJTensor(3, 3, 1, 2);
      J[3, 3] = ItemJTensor(2, 3, 3, 2);
      J[3, 4] = J[4, 3] = ItemJTensor(2, 3, 1, 3);
      J[3, 5] = J[5, 3] = ItemJTensor(2, 3, 1, 2);
      J[4, 4] = ItemJTensor(1, 3, 3, 1);
      J[4, 5] = J[5, 4] = ItemJTensor(1, 3, 1, 2);
      J[5, 5] = ItemJTensor(1, 2, 2, 1);
      return J;
    }
    private double ItemJTensor(int i, int j, int k, int l)
    {
      //return 1.0 / 2.0 * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
      //return KroneckerDelta(i, k) * KroneckerDelta(j, l);//Liu 2016
      return KroneckerDelta(i, j) * KroneckerDelta(k, l);//Formular
    }
    protected double KroneckerDelta(int i, int j)
    {
      return (i == j) ? 1.0 : 0.0;
    }

    //protected DoubleMatrix CreateMaterialMatrix1()
    //{
    //  AnisotropicElasticity prob = (AnisotropicElasticity)Material.GetProperty(MaterialPropertyName.AnisotropicElasticity);
    //  if (prob != null)
    //  {
    //    double[,] ce = prob.GetAnisotropicElasticityMatrix();
    //    DoubleMatrix C = new DoubleMatrix(6, 6);
    //    bool isInverse = prob.GetIsInverseMatrix();
    //    for (int i = 0; i < 6; i++)
    //      for (int j = 0; j < 6; j++)
    //      {
    //        C[i, j] = ce[i, j];
    //      }
    //    if (!isInverse)
    //      return C;
    //    else
    //      return NMathFunctions.Inverse(C);
    //  }
    //  else
    //  {
    //    DoubleMatrix C = null;
    //    OrthotropicElasticity prob1 = (OrthotropicElasticity)Material.GetProperty(MaterialPropertyName.OrthotropicElasticity);
    //    if (prob1 == null)
    //    {
    //      C = new DoubleMatrix(6, 6);
    //      //Xem lai thu tu xx yy zz yz xz xy
    //      C[0, 0] = coefficient(1, 1, 1, 1);
    //      C[0, 1] = C[1, 0] = coefficient(1, 1, 2, 2);
    //      C[0, 2] = C[2, 0] = coefficient(1, 1, 3, 3);
    //      C[0, 3] = C[3, 0] = coefficient(1, 1, 2, 3);
    //      C[0, 4] = C[4, 0] = coefficient(1, 1, 1, 3);
    //      C[0, 5] = C[5, 0] = coefficient(1, 1, 1, 2);
    //      C[1, 1] = coefficient(2, 2, 2, 2);
    //      C[1, 2] = C[2, 1] = coefficient(2, 2, 3, 3);
    //      C[1, 3] = C[3, 1] = coefficient(2, 2, 2, 3);
    //      C[1, 4] = C[4, 1] = coefficient(2, 2, 1, 3);
    //      C[1, 5] = C[5, 1] = coefficient(2, 2, 1, 2);
    //      C[2, 2] = coefficient(3, 3, 3, 3);
    //      C[2, 3] = C[3, 2] = coefficient(3, 3, 2, 3);
    //      C[2, 4] = C[4, 2] = coefficient(3, 3, 1, 3);
    //      C[2, 5] = C[5, 2] = coefficient(3, 3, 1, 2);
    //      C[3, 3] = coefficient(2, 3, 3, 2);
    //      C[3, 4] = C[4, 3] = coefficient(2, 3, 1, 3);
    //      C[3, 5] = C[5, 3] = coefficient(2, 3, 1, 2);
    //      C[4, 4] = coefficient(1, 3, 3, 1);
    //      C[4, 5] = C[5, 4] = coefficient(1, 3, 1, 2);
    //      C[5, 5] = coefficient(1, 2, 2, 1);
    //      C = MatrixTool.GetSubMatrix(C, new int[] { 0, 1, 2, 5 }, new int[] { 0, 1, 2, 5 });
    //    }
    //    else
    //    {
    //      C = new DoubleMatrix(prob1.GetOrthotropicElasticityMatrix());
    //    }
    //    return C;
    //  }
    //}

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

      DoubleMatrix Keuu = null;
      DoubleMatrix Kepp = null;
      DoubleMatrix Keup = null;
      DoubleMatrix Kepu = null;
      DoubleMatrix Bu = new DoubleMatrix(4, 2 * nen);
      DoubleMatrix Bp = new DoubleMatrix(3, nen);
      switch (AbstractModel.TypeSchemeNonlinearSolver)
      {
        case TypeSchemeNonlinearSolver.Monolithic:
          Keuu = new DoubleMatrix(2 * nen, 2 * nen);
          Kepp = new DoubleMatrix(nen, nen);
          Keup = new DoubleMatrix(2 * nen, nen);
          Kepu = new DoubleMatrix(nen, 2 * nen);
          break;
        case TypeSchemeNonlinearSolver.SingleStaggered:
          Keuu = new DoubleMatrix(2 * nen, 2 * nen);
          Kepp = new DoubleMatrix(nen, nen);
          break;
        case TypeSchemeNonlinearSolver.Staggered:
          if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
          {
            Keuu = new DoubleMatrix(2 * nen, 2 * nen);
          }
          else
          {
            Kepp = new DoubleMatrix(nen, nen);
          }
          break;
      }
      double l0 = ModelStructurePhaseFieldStatic.l0;
      //double lch = ModelStructurePhaseFieldStatic.lch;
      double k1 = ModelStructurePhaseFieldStatic.k;
      DoubleMatrix BuT = null;
      DoubleMatrix BpT = null;
      double previousPhase = 0;
      double currentPhase = 0;
      DoubleVector ue = null;////////////
      DoubleVector strain = null;
      DoubleMatrix DsigmaDEpsilon = null;

      DoubleVector principleStrain;
      DoubleVector[] eigenVectors;
      DoubleMatrix BuTDBu = null;
      DoubleMatrix kepp = null;
      double xiEpsilon = 0;
      double xiEpsilonPos = 0;
      double xiEpsilonNeg = 0;
      double lamda = 0, mu = 0, gc = 0, ft = 0, lch = 0;
      if (!(Material is FGMUserDefinedGradedMaterial))
      {
        lamda = Material.GetProperty(MaterialPropertyName.LameParameter).GetValueProperty();
        mu = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
        gc = Material.GetProperty(MaterialPropertyName.CriticalEnergyReleaseRate).GetValueProperty();
        ft = Material.GetProperty(MaterialPropertyName.TensileYieldStrength).GetValueProperty();
      }
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
          DoubleMatrix dNdx = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.dNdX);//gps[i, j].dNdX;
          DoubleVector N = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);//gps[i, j].Ni;
          double detJ = Math.Abs((double)gpsij.GetValue(DataInGausspoint.detJ));//Math.Abs(gps[i, j].detJ);
          for (int k = 0; k < nen; k++)
          {
            //xx yy zz xy
            Bu[0, 2 * k] = dNdx[k, 0];
            Bu[2, 2 * k] = c * N[k] / x1;
            Bu[3, 2 * k] = dNdx[k, 1];
            Bu[1, 2 * k + 1] = dNdx[k, 1];
            Bu[3, 2 * k + 1] = dNdx[k, 0];

            //x y z
            Bp[0, k] = dNdx[k, 0];
            Bp[1, k] = dNdx[k, 1];
            Bp[2, k] = c * N[k] / x1;
          }
          if (Material is FGMUserDefinedGradedMaterial)
          {
            lamda = ComputeParameterProperty(MaterialPropertyName.LameParameter, gpsij.location[0], gpsij.location[1]);
            mu = ComputeParameterProperty(MaterialPropertyName.ShearModulus, gpsij.location[0], gpsij.location[1]);
            gc = ComputeParameterProperty(MaterialPropertyName.CriticalEnergyReleaseRate, gpsij.location[0], gpsij.location[1]);
            ft = ComputeParameterProperty(MaterialPropertyName.TensileYieldStrength, gpsij.location[0], gpsij.location[1]);
          }
          double E = mu * (3 * lamda + 2 * mu) / (lamda + mu);
          lch = E * gc / (ft * ft);

          DoubleMatrix NTN = NMathFunctions.OuterProduct(N, N);
          double valDDOmegaD, valCAlpha, valDDAlphaD;
          switch (AbstractModel.TypeSchemeNonlinearSolver)
          {
            case TypeSchemeNonlinearSolver.Monolithic:
              #region monolithicpenalty
              //BuT = Bu.Transpose();
              //BpT = Bp.Transpose();

              //gpsij.SetValue(DataInGausspoint.previousPhase, (double)gpsij.GetValue(DataInGausspoint.currentPhase));
              //gpsij.SetValue(DataInGausspoint.currentPhase, NMathFunctions.Dot(N, base.GetEachVariableLocal(Result.PHASEFIELD)));
              //previousPhase = (double)gpsij.GetValue(DataInGausspoint.previousPhase);
              //currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);

              //if (currentPhase < previousPhase)
              //{
              //  gpsij.SetValue(DataInGausspoint.currentPhase, previousPhase);
              //  currentPhase = previousPhase;
              //}
              //if (currentPhase < 0)
              //  gpsij.SetValue(DataInGausspoint.currentPhase, 0.0);
              //if (currentPhase > 1)
              //  gpsij.SetValue(DataInGausspoint.currentPhase, 1.0);
              //currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);

              //strain = ComputeStrainGauss(patch, Bu, gpsij);

              //gpsij.SetValue(DataInGausspoint.currentStrain, strain);
              //ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
              //DoubleVector stress = ComputeStress(lamda, mu, currentPhase, principleStrain, eigenVectors);
              //gpsij.SetValue(DataInGausspoint.currentStress, stress);
              //DsigmaDEpsilon = CreateMaterialMatrixIsotropicAnisotropic(lamda, mu, k1, l0, lch, currentPhase, principleStrain, eigenVectors);
              //BuTDBu = new DoubleSymmetricMatrix(NMathFunctions.Product(BuT, NMathFunctions.Product(DsigmaDEpsilon, Bu)));
              //Keuu += a * gpsij.weight * detJbar * detJ * BuTDBu;

              //xiEpsilon = (double)gpsij.GetValue(DataInGausspoint.xiEpsilon);
              //double phi0 = (double)gpsij.GetValue(DataInGausspoint.lastPhase);

              //double d_dot = currentPhase - phi0;//phi1 - phi0;
              ////if (d_dot >= 0)
              ////  d_dot = 0;
              ////else
              ////  d_dot = -d_dot;
              //double epsi = ModelStructurePhaseFieldStatic.epsi;
              //double deltaT = AbstractNonlinearSolver.deltaT;
              //int n = ModelStructurePhaseFieldStatic.n;//2
              //if (ModelStructurePhaseFieldStatic.TypeOrder == TypeOrderPhaseField.SecondOrder)
              //{
              //  kepp = gc * l0 * new DoubleSymmetricMatrix(NMathFunctions.Product(BpT, Bp)) + (gc / l0 + 2 * xiEpsilon + epsi * Math.Pow(d_dot, n - 1) / deltaT) * NTN;
              //}
              //else
              //{
              //  DoubleMatrix ddNdX = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
              //  DoubleVector BB = new DoubleVector(nen);
              //  for (int k = 0; k < nen; k++)
              //  {
              //    BB[k] = ddNdX[k, 0] + ddNdX[k, 1];
              //  }

              //  DoubleSymmetricMatrix F4 = gc * Math.Pow(l0, 3) / 16.0 * new DoubleSymmetricMatrix(NMathFunctions.OuterProduct(BB, BB));
              //  kepp = gc * l0 / 2.0 * new DoubleSymmetricMatrix(NMathFunctions.Product(BpT, Bp)) + (gc / l0 + 2 * xiEpsilon + epsi * Math.Pow(d_dot, n - 1) / deltaT) * NTN + F4;
              //}
              //Kepp += a * gpsij.weight * detJbar * detJ * kepp;


              //DoubleVector dSdphi = -2 * (1 - currentPhase) * ComputeDStressDPhi(lamda, mu, k1, currentPhase, principleStrain, eigenVectors);
              //Keup += a * gpsij.weight * detJbar * detJ * NMathFunctions.Product(BuT, NMathFunctions.OuterProduct(dSdphi, N));
              #endregion

              #region monolithic
              BuT = Bu.Transpose();
              BpT = Bp.Transpose();
              //////////////////// da cap nhat sau moi buoc giai
              //gpsij.SetValue(DataInGausspoint.previousPhase, (double)gpsij.GetValue(DataInGausspoint.currentPhase));
              //gpsij.SetValue(DataInGausspoint.currentPhase, NMathFunctions.Dot(N, base.GetEachVariableLocal(Result.PHASEFIELD)));
              //previousPhase = (double)gpsij.GetValue(DataInGausspoint.previousPhase);
              //currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);

              //if (currentPhase < previousPhase)
              //{
              //  gpsij.SetValue(DataInGausspoint.currentPhase, previousPhase);
              //  currentPhase = previousPhase;
              //}
              //if (currentPhase < 0.0)
              //  gpsij.SetValue(DataInGausspoint.currentPhase, 0.0);
              //if (currentPhase > 1.0)
              //  gpsij.SetValue(DataInGausspoint.currentPhase, 1.0);
              //////////////////////////////////////////////////////
              currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);

              //strain = ComputeStrainGauss(patch, Bu, gpsij);

              //gpsij.SetValue(DataInGausspoint.currentStrain, strain);
              strain = (DoubleVector)gpsij.GetValue(DataInGausspoint.currentStrain);
              ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
              DsigmaDEpsilon = CreateMaterialMatrixIsotropicAnisotropic(lamda, mu, k1, l0, lch, currentPhase, principleStrain, eigenVectors);
              BuTDBu = NMathFunctions.Product(BuT, NMathFunctions.Product(DsigmaDEpsilon, Bu));
              Keuu += a * gpsij.weight * detJbar * detJ * BuTDBu;

              xiEpsilon = (double)gpsij.GetValue(DataInGausspoint.xiEpsilon);
              valDDOmegaD = ddOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
              valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
              valDDAlphaD = ddAlphaD(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
              if (ModelStructurePhaseFieldStatic.typeOrderPhaseField == TypeOrderPhaseField.SecondOrder)
              {
                kepp = 2.0 * gc * l0 / valCAlpha * NMathFunctions.Product(BpT, Bp) + (valDDAlphaD * gc / (valCAlpha * l0) + valDDOmegaD * xiEpsilon) * NTN;
              }
              else
              {
                DoubleMatrix ddNdX = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
                DoubleVector BB = new DoubleVector(nen);
                for (int k = 0; k < nen; k++)
                {
                  BB[k] = ddNdX[k, 0] + ddNdX[k, 1];
                }

                //DoubleMatrix F4 = gc * Math.Pow(l0, 3) / 16.0 * NMathFunctions.OuterProduct(BB, BB);
                //kepp = gc * l0 / 2.0 * NMathFunctions.Product(BpT, Bp) + (gc / l0 + 2 * xiEpsilon) * NTN + F4;
                DoubleMatrix F4 = gc * Math.Pow(l0, 3) / (8.0 * valCAlpha) * NMathFunctions.OuterProduct(BB, BB);
                kepp = gc * l0 / valCAlpha * NMathFunctions.Product(BpT, Bp) + (valDDAlphaD * gc / (l0 * valCAlpha) + valDDOmegaD * xiEpsilon) * NTN + F4;
              }
              Kepp += a * gpsij.weight * detJbar * detJ * kepp;

              double valDOmegaD = dOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
              DoubleVector dSdphi = valDOmegaD * ComputeDStressDPhi(lamda, mu, principleStrain, eigenVectors);
              Keup += a * gpsij.weight * detJbar * detJ * NMathFunctions.Product(BuT, NMathFunctions.OuterProduct(dSdphi, N));
              DoubleVector dHDEpsilon = valDOmegaD * ComputeDHDEpsilon(lamda, mu, principleStrain, eigenVectors);
              Kepu += a * gpsij.weight * detJbar * detJ * NMathFunctions.Product(NMathFunctions.OuterProduct(N, dHDEpsilon), Bu);
              #endregion
              break;
            case TypeSchemeNonlinearSolver.SingleStaggered:
              #region singlestaggered
              BuT = Bu.Transpose();
              BpT = Bp.Transpose();
              //gpsij.SetValue(DataInGausspoint.previousPhase, (double)gpsij.GetValue(DataInGausspoint.currentPhase));
              //gpsij.SetValue(DataInGausspoint.currentPhase, NMathFunctions.Dot(N, base.GetEachVariableLocal(Result.PHASEFIELD)));
              //previousPhase = (double)gpsij.GetValue(DataInGausspoint.previousPhase);
              //currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);

              //if (currentPhase < previousPhase)
              //{
              //  gpsij.SetValue(DataInGausspoint.currentPhase, previousPhase);
              //  currentPhase = previousPhase;
              //}
              //if (currentPhase < 0.0)
              //  gpsij.SetValue(DataInGausspoint.currentPhase, 0.0);
              //if (currentPhase > 1.0)
              //  gpsij.SetValue(DataInGausspoint.currentPhase, 1.0);
              currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);

              //strain = ComputeStrainGauss(patch, Bu, gpsij);
              //gpsij.SetValue(DataInGausspoint.currentStrain, strain);
              strain = (DoubleVector)gpsij.GetValue(DataInGausspoint.currentStrain);
              ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
              DsigmaDEpsilon = CreateMaterialMatrixIsotropicAnisotropic(lamda, mu, k1, l0, lch, currentPhase, principleStrain, eigenVectors);
              BuTDBu = NMathFunctions.Product(BuT, NMathFunctions.Product(DsigmaDEpsilon, Bu));
              Keuu += a * gpsij.weight * detJbar * detJ * BuTDBu;

              xiEpsilon = (double)gpsij.GetValue(DataInGausspoint.xiEpsilon);
              valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);//pi
              valDDAlphaD = ddAlphaD(ModelStructurePhaseFieldStatic.typePhaseFieldModel);//-2
              valDDOmegaD = ddOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
              if (ModelStructurePhaseFieldStatic.typeOrderPhaseField == TypeOrderPhaseField.SecondOrder)
              {
                kepp = 2.0 * gc * l0 / valCAlpha * NMathFunctions.Product(BpT, Bp) + (valDDAlphaD * gc / (l0 * valCAlpha) + valDDOmegaD * xiEpsilon) * NTN;
                //kepp = gc * l0 * new DoubleSymmetricMatrix(NMathFunctions.Product(BpT, Bp)) + (gc / l0 + 2 * xiEpsilon) * NTN;
              }
              else
              {
                DoubleMatrix ddNdX = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
                DoubleVector BB = new DoubleVector(nen);
                for (int k = 0; k < nen; k++)
                {
                  BB[k] = ddNdX[k, 0] + ddNdX[k, 1];
                }
                //DoubleMatrix F4 = gc * Math.Pow(l0, 3) / 16.0 * NMathFunctions.OuterProduct(BB, BB);
                //kepp = gc * l0 / 2.0 * NMathFunctions.Product(BpT, Bp) + (gc / l0 + 2 * xiEpsilon) * NTN + F4;
                DoubleMatrix F4 = gc * Math.Pow(l0, 3) / (8.0 * valCAlpha) * NMathFunctions.OuterProduct(BB, BB);////
                kepp = gc * l0 / valCAlpha * NMathFunctions.Product(BpT, Bp) + (valDDAlphaD * gc / (l0 * valCAlpha) + valDDOmegaD * xiEpsilon) * NTN + F4;
              }
              Kepp += a * gpsij.weight * detJbar * detJ * kepp;
              #endregion
              break;
            case TypeSchemeNonlinearSolver.Staggered:
              #region staggered
              if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
              {
                //displacement
                BuT = Bu.Transpose();
                currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);
                //strain = ComputeStrainGauss(patch, Bu, gpsij);
                //gpsij.SetValue(DataInGausspoint.currentStrain, strain);
                strain = (DoubleVector)gpsij.GetValue(DataInGausspoint.currentStrain);
                ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
                DsigmaDEpsilon = CreateMaterialMatrixIsotropicAnisotropic(lamda, mu, k1, l0, lch, currentPhase, principleStrain, eigenVectors);
                BuTDBu = NMathFunctions.Product(BuT, NMathFunctions.Product(DsigmaDEpsilon, Bu));
                Keuu += a * gpsij.weight * detJbar * detJ * BuTDBu;
              }
              else
              {
                BpT = Bp.Transpose();
                xiEpsilon = (double)gpsij.GetValue(DataInGausspoint.xiEpsilon);
                currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);
                valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);//pi
                valDDAlphaD = ddAlphaD(ModelStructurePhaseFieldStatic.typePhaseFieldModel);//-2
                valDDOmegaD = ddOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                if (ModelStructurePhaseFieldStatic.typeOrderPhaseField == TypeOrderPhaseField.SecondOrder)
                {
                  // l0 biểu diễn độ rộng vết nứt
                  //Singh2018
                  //kepp = gc * l0 * new DoubleSymmetricMatrix(NMathFunctions.Product(BpT, Bp)) + (gc / l0 + 2 * xiEpsilon) * NTN;
                  kepp = 2.0 * gc * l0 / valCAlpha * NMathFunctions.Product(BpT, Bp) + (valDDAlphaD * gc / (l0 * valCAlpha) + valDDOmegaD * xiEpsilon) * NTN;
                }
                else
                {
                  DoubleMatrix ddNdX = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
                  DoubleVector BB = new DoubleVector(nen);
                  for (int k = 0; k < nen; k++)
                  {
                    BB[k] = ddNdX[k, 0] + ddNdX[k, 1];
                  }
                  // l0 biểu diễn độ rộng vết nứt
                  //DoubleMatrix F4 = gc * Math.Pow(l0, 3) / 16.0 * NMathFunctions.OuterProduct(BB, BB);
                  //kepp = gc * l0 / 2.0 * NMathFunctions.Product(BpT, Bp) + (gc / l0 + 2 * xiEpsilon) * NTN + F4;
                  DoubleMatrix F4 = gc * Math.Pow(l0, 3) / (8.0 * valCAlpha) * NMathFunctions.OuterProduct(BB, BB);////
                  kepp = gc * l0 / valCAlpha * NMathFunctions.Product(BpT, Bp) + (valDDAlphaD * gc / (l0 * valCAlpha) + valDDOmegaD * xiEpsilon) * NTN + F4;
                }
                Kepp += a * gpsij.weight * detJbar * detJ * kepp;
              }
              #endregion
              break;
          }
        }
      }
      for (int j = 0; j < nen; j++)
        for (int i = 0; i < nen; i++)
        {
          switch (AbstractModel.TypeSchemeNonlinearSolver)
          {
            case TypeSchemeNonlinearSolver.Monolithic:
              for (int jj = 0; jj < 2; jj++)
              {
                for (int ii = 0; ii < 2; ii++)
                {
                  Ke[i * 3 + ii, j * 3 + jj] = Keuu[i * 2 + ii, j * 2 + jj];
                  Ke[i * 3 + ii, j * 3 + 2] = Keup[i * 2 + ii, j];
                }
                Ke[i * 3 + 2, j * 3 + jj] = Kepu[i, j * 2 + jj];
                //Ke[i * 3 + 2, j * 3 + jj] = Ke[j * 3 + jj, i * 3 + 2] = Keup[j * 2 + jj, i];
              }
              Ke[i * 3 + 2, j * 3 + 2] = Kepp[i, j];

              //for (int jj = 0; jj < 2; jj++)
              //{
              //  for (int ii = 0; ii < 2; ii++)
              //  {
              //    Ke[i * 3 + ii, j * 3 + jj] = Keuu[i * 2 + ii, j * 2 + jj];
              //  }
              //  Ke[i * 3 + 2, j * 3 + jj] = Keup[j * 2 + jj, i];
              //}
              //Ke[i * 3 + 2, j * 3 + 2] = Kepp[i, j];
              break;
            case TypeSchemeNonlinearSolver.SingleStaggered:
              for (int jj = 0; jj < 2; jj++)
              {
                for (int ii = 0; ii < 2; ii++)
                {
                  Ke[i * 3 + ii, j * 3 + jj] = Keuu[i * 2 + ii, j * 2 + jj];
                }
              }
              Ke[i * 3 + 2, j * 3 + 2] = Kepp[i, j];
              break;
            case TypeSchemeNonlinearSolver.Staggered:
              if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
              {
                for (int jj = 0; jj < 2; jj++)
                {
                  for (int ii = 0; ii < 2; ii++)
                  {
                    Ke[i * 3 + ii, j * 3 + jj] = Keuu[i * 2 + ii, j * 2 + jj];
                  }
                }
              }
              else
              {
                Ke[i * 3 + 2, j * 3 + 2] = Kepp[i, j];
              }
              break;
          }
        }
    }

    private DoubleVector ComputeStrainGauss(PatchStructure2D patch, DoubleMatrix Bu, GaussPoints gpsij)
    {
      DoubleVector ue = GetDisplacementLocal();////////////
      DoubleVector strain = NMathFunctions.Product(Bu, ue);
      if (patch.StateStress == Structure2DState.PlaneStress)
      {
        double nu = ComputeParameterProperty(MaterialPropertyName.PoissonRatio, gpsij.location[0], gpsij.location[1]);
        strain[2] = -nu / (1.0 - nu) * (strain[0] + strain[1]);
      }
      return strain;
    }

    private DoubleVector ComputeStress(double lamda, double mu, double lch, double phiCurrent, DoubleVector principleStrain, DoubleVector[] eigenVectors)
    {
      double k1 = ModelStructurePhaseFieldStatic.k;
      double l0 = ModelStructurePhaseFieldStatic.l0;
      //double lch = ModelStructurePhaseFieldStatic.lch;//mm

      //double E = mu * (3 * lamda + 2 * mu) / (lamda + mu);
      //double lch = E * gc / (ft * ft);

      DoubleVector stress = null;
      double valOmega = OmegaD(phiCurrent, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel) + k1;
      switch (ModelStructurePhaseFieldStatic.typeDegradationModel)
      {
        case TypeDegradationEnergy.Isotropic:
        case TypeDegradationEnergy.Hybrid:
          DoubleVector stress0 = null;
          ComputeStressIsotropic(principleStrain, eigenVectors, lamda, mu, out stress0);
          stress = valOmega * stress0;
          break;
        case TypeDegradationEnergy.Anisotropic:
          DoubleVector stress0Pos = null;
          DoubleVector stress0Neg = null;
          ComputeStressAnisotropic(principleStrain, eigenVectors, lamda, mu, out stress0Pos, out stress0Neg);
          stress = valOmega * stress0Pos + stress0Neg;
          break;
      }

      return stress;
    }
    private DoubleVector ComputeEfficientStress(double lamda, double mu, DoubleVector principleStrain, DoubleVector[] eigenVectors)
    {
      DoubleVector stress = null;
      switch (ModelStructurePhaseFieldStatic.typeDegradationModel)
      {
        case TypeDegradationEnergy.Isotropic:
        case TypeDegradationEnergy.Hybrid:
          DoubleVector stress0 = null;
          ComputeStressIsotropic(principleStrain, eigenVectors, lamda, mu, out stress0);
          stress = stress0;
          break;
        case TypeDegradationEnergy.Anisotropic:
          DoubleVector stress0Pos = null;
          DoubleVector stress0Neg = null;
          ComputeStressAnisotropic(principleStrain, eigenVectors, lamda, mu, out stress0Pos, out stress0Neg);
          stress = stress0Pos;
          break;
      }

      return stress;
    }
    private DoubleVector ComputeDStressDPhi(double lamda, double mu, DoubleVector principleStrain, DoubleVector[] eigenVectors)
    {
      DoubleVector stress = null;
      switch (ModelStructurePhaseFieldStatic.typeDegradationModel)
      {
        case TypeDegradationEnergy.Isotropic:
        case TypeDegradationEnergy.Hybrid:
          DoubleVector stress0 = null;
          ComputeStressIsotropic(principleStrain, eigenVectors, lamda, mu, out stress0);
          stress = stress0;
          break;
        case TypeDegradationEnergy.Anisotropic:
          DoubleVector stress0Pos = null;
          DoubleVector stress0Neg = null;
          ComputeStressAnisotropic(principleStrain, eigenVectors, lamda, mu, out stress0Pos, out stress0Neg);
          stress = stress0Pos;
          break;
      }

      return stress;
    }

    private DoubleVector ComputeDHDEpsilon(double lamda, double mu, DoubleVector principleStrain, DoubleVector[] eigenVectors)
    {
      DoubleVector stress = null;
      switch (ModelStructurePhaseFieldStatic.typeDegradationModel)
      {
        case TypeDegradationEnergy.Isotropic:
          DoubleVector stress0 = null;
          ComputeStressIsotropic(principleStrain, eigenVectors, lamda, mu, out stress0);
          stress = stress0;
          break;
        case TypeDegradationEnergy.Hybrid:
        case TypeDegradationEnergy.Anisotropic:
          DoubleVector stress0Pos = null;
          DoubleVector stress0Neg = null;
          ComputeStressAnisotropic(principleStrain, eigenVectors, lamda, mu, out stress0Pos, out stress0Neg);
          stress = stress0Pos;
          break;
      }

      return stress;
    }

    private DoubleMatrix CreateMaterialMatrixIsotropicAnisotropic(double lamda, double mu, double k1, double l0, double lch, double phiCurrent, DoubleVector principleStrain, DoubleVector[] eigenVectors)
    {
      DoubleMatrix DsigmaDEpsilon = null;
      switch (ModelStructurePhaseFieldStatic.typeDegradationModel)
      {
        case TypeDegradationEnergy.Isotropic:
        case TypeDegradationEnergy.Hybrid:
          DsigmaDEpsilon = CreateMaterialMatrixIsotropic(phiCurrent, k1, l0, lch, lamda, mu);
          break;
        case TypeDegradationEnergy.Anisotropic:
          DsigmaDEpsilon = CreateMaterialMatrixAnisotropic(principleStrain, eigenVectors, phiCurrent, k1, l0, lch, lamda, mu);
          break;
      }
      return DsigmaDEpsilon;
    }

    public override double ComputeSharpCrackSurface()
    {
      double Ad = 0;
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
      DoubleMatrix Bp = new DoubleMatrix(3, nen);
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
          DoubleMatrix dNdx = (DoubleMatrix)gps[i, j].GetValue(DataInGausspoint.dNdX);//gps[i, j].dNdX;
          DoubleVector N = (DoubleVector)gps[i, j].GetValue(DataInGausspoint.Ni);//gps[i, j].Ni;
          double detJ = Math.Abs((double)gps[i, j].GetValue(DataInGausspoint.detJ));//Math.Abs(gps[i, j].detJ);
          for (int k = 0; k < nen; k++)
          {
            //x y z
            Bp[0, k] = dNdx[k, 0];
            Bp[1, k] = dNdx[k, 1];
            Bp[2, k] = c * N[k] / x1;
          }
          double l0 = ModelStructurePhaseFieldStatic.l0;
          double currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);
          DoubleVector phie = GetEachVariableLocal(Result.PHASEFIELD);
          double valA = 0;
          if (ModelStructurePhaseFieldStatic.typeOrderPhaseField == TypeOrderPhaseField.SecondOrder)
          {
            //// l0 là kích thước của vết nứt
            double valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
            double valAlphaD = alphaD(currentPhase, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
            DoubleVector Bphi = NMathFunctions.Product(Bp, phie);
            valA = 1.0 / valCAlpha * (1.0 / l0 * valAlphaD + l0 * Math.Pow(Bphi.TwoNormSquared(), 2));
          }
          else
          {

          }
          Ad += a * gps[i, j].weight * detJbar * detJ * valA;
        }
      }

      return Ad;
    }
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
      DoubleVector fiU = null;
      DoubleVector fiP = null;
      DoubleMatrix Bu = new DoubleMatrix(4, 2 * nen);
      DoubleMatrix Bp = new DoubleMatrix(3, nen);

      switch (AbstractModel.TypeSchemeNonlinearSolver)
      {
        case TypeSchemeNonlinearSolver.Monolithic:
        case TypeSchemeNonlinearSolver.SingleStaggered:
          fiU = new DoubleVector(2 * nen);
          fiP = new DoubleVector(nen);
          break;
        case TypeSchemeNonlinearSolver.Staggered:
          if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
          {
            fiU = new DoubleVector(2 * nen);
          }
          else
          {
            fiP = new DoubleVector(nen);
          }
          break;
      }

      double lamda = 0, mu = 0, gc = 0, ft = 0, lch = 0;
      if (!(Material is FGMUserDefinedGradedMaterial))
      {
        lamda = Material.GetProperty(MaterialPropertyName.LameParameter).GetValueProperty();
        mu = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
        gc = Material.GetProperty(MaterialPropertyName.CriticalEnergyReleaseRate).GetValueProperty();
        ft = Material.GetProperty(MaterialPropertyName.TensileYieldStrength).GetValueProperty();
      }
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
          DoubleMatrix dNdx = (DoubleMatrix)gps[i, j].GetValue(DataInGausspoint.dNdX);//gps[i, j].dNdX;
          DoubleVector N = (DoubleVector)gps[i, j].GetValue(DataInGausspoint.Ni);//gps[i, j].Ni;
          double detJ = Math.Abs((double)gps[i, j].GetValue(DataInGausspoint.detJ));//Math.Abs(gps[i, j].detJ);
          for (int k = 0; k < nen; k++)
          {
            //xx yy zz xy
            Bu[0, 2 * k] = dNdx[k, 0];
            Bu[2, 2 * k] = c * N[k] / x1;
            Bu[3, 2 * k] = dNdx[k, 1];
            Bu[1, 2 * k + 1] = dNdx[k, 1];
            Bu[3, 2 * k + 1] = dNdx[k, 0];

            //x y z
            Bp[0, k] = dNdx[k, 0];
            Bp[1, k] = dNdx[k, 1];
            Bp[2, k] = c * N[k] / x1;
          }
          if (Material is FGMUserDefinedGradedMaterial)
          {
            lamda = ComputeParameterProperty(MaterialPropertyName.LameParameter, gps[i, j].location[0], gps[i, j].location[1]);
            mu = ComputeParameterProperty(MaterialPropertyName.ShearModulus, gps[i, j].location[0], gps[i, j].location[1]);
            gc = ComputeParameterProperty(MaterialPropertyName.CriticalEnergyReleaseRate, gps[i, j].location[0], gps[i, j].location[1]);
            ft = ComputeParameterProperty(MaterialPropertyName.TensileYieldStrength, gpsij.location[0], gpsij.location[1]);
          }
          double E = mu * (3 * lamda + 2 * mu) / (lamda + mu);
          lch = E * gc / (ft * ft);
          //double lch = ModelStructurePhaseFieldStatic.lch;
          double l0 = ModelStructurePhaseFieldStatic.l0;
          double k1 = ModelStructurePhaseFieldStatic.k;
          double previousPhase;
          double currentPhase;
          DoubleVector stress = null;
          DoubleVector phie = null;
          double xiEpsilon = 0;
          DoubleMatrix F1 = null;
          DoubleVector F3 = null;
          DoubleVector F2 = null;
          DoubleVector ue = null;
          DoubleVector strain = null;
          DoubleVector principleStrain;
          DoubleVector[] eigenVectors;
          double valDOmegaD;
          double valCAlpha;
          double valDAlphaD;

          switch (AbstractModel.TypeSchemeNonlinearSolver)
          {
            case TypeSchemeNonlinearSolver.Monolithic:
            //#region monolithicPanelty
            //gpsij.SetValue(DataInGausspoint.previousPhase, (double)gpsij.GetValue(DataInGausspoint.currentPhase));
            //phie = GetEachVariableLocal(Result.PHASEFIELD);
            //gpsij.SetValue(DataInGausspoint.currentPhase, NMathFunctions.Dot(N, phie));
            //previousPhase = (double)gpsij.GetValue(DataInGausspoint.previousPhase);
            //currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);

            //if (currentPhase < previousPhase)
            //{
            //  gpsij.SetValue(DataInGausspoint.currentPhase, previousPhase);
            //  currentPhase = previousPhase;
            //}
            //if (currentPhase < 0)
            //  gpsij.SetValue(DataInGausspoint.currentPhase, 0.0);
            //if (currentPhase > 1)
            //  gpsij.SetValue(DataInGausspoint.currentPhase, 1.0);
            //currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);
            ////gps[i, j].previousPhase = gps[i, j].currentPhase;

            ////gps[i, j].currentPhase = NMathFunctions.Dot(N, phie);

            ////if (gps[i, j].currentPhase < gps[i, j].previousPhase)
            ////  gps[i, j].currentPhase = gps[i, j].previousPhase;
            ////if (gps[i, j].currentPhase < 0)
            ////  gps[i, j].currentPhase = 0;
            ////if (gps[i, j].currentPhase > 1)
            ////  gps[i, j].currentPhase = 1;
            ////phiCurrent = gps[i, j].currentPhase;
            //strain = ComputeStrainGauss(patch, Bu, gpsij);
            //gpsij.SetValue(DataInGausspoint.currentStrain, strain);
            ////gps[i, j].currentStrain = NMathFunctions.Product(Bu, ue);
            ////gps[i, j].currentStrain[3] /= 2.0;
            ////strain = gps[i, j].currentStrain;
            //ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
            //DoubleVector stress = ComputeStress(lamda, mu, currentPhase, principleStrain, eigenVectors);
            //gpsij.SetValue(DataInGausspoint.currentStress, stress);

            ////stress = gps[i, j].currentStress = ComputeStress(lamda, mu, k1, currentPhase, principleStrain, eigenVectors);
            //fiU += a * gps[i, j].weight * detJbar * detJ * NMathFunctions.Product(Bu.Transpose(), stress);

            //xiEpsilon = (double)gpsij.GetValue(DataInGausspoint.xiEpsilon);
            //double phi0 = (double)gpsij.GetValue(DataInGausspoint.lastPhase);
            ////xiEpsilon = gps[i, j].xiEpsilon;
            ////double phi0 = gps[i, j].lastPhase;
            //double d_dot = currentPhase - phi0;//phi1 - phi0;
            ////if (d_dot >= 0)
            ////  d_dot = 0;
            ////else
            ////  d_dot = -d_dot;
            //double epsi = ModelStructurePhaseFieldStatic.epsi;
            //double deltaT = AbstractNonlinearSolver.deltaT;
            //int n = ModelStructurePhaseFieldStatic.n;//2
            //if (ModelStructurePhaseFieldStatic.TypeOrder == TypeOrderPhaseField.SecondOrder)
            //{
            //  F3 = (gc / l0 * currentPhase - 2 * (1 - currentPhase) * xiEpsilon + epsi / (n * deltaT) * Math.Pow(d_dot, n)) * N;
            //  F1 = gc * l0 * NMathFunctions.Product(Bp.Transpose(), Bp);
            //}
            //else
            //{
            //  F3 = (gc / l0 * currentPhase - 2 * (1 - currentPhase) * xiEpsilon + epsi / (n * deltaT) * Math.Pow(d_dot, n)) * N;
            //  DoubleMatrix ddNdX = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
            //  DoubleVector BB = new DoubleVector(nen);
            //  for (int k = 0; k < nen; k++)
            //  {
            //    BB[k] = ddNdX[k, 0] + ddNdX[k, 1];
            //  }
            //  F1 = gc * l0 / 2.0 * NMathFunctions.Product(Bp.Transpose(), Bp) + gc * Math.Pow(l0, 3) / 16.0 * NMathFunctions.OuterProduct(BB, BB);
            //}
            ////phie = GetEachVariableLocal(Result.PHASEFIELD);
            //F2 = NMathFunctions.Product(F1, phie);
            //fiP += a * gps[i, j].weight * detJbar * detJ * (F2 + F3);
            //#endregion
            //#region monolithic
            //gpsij.SetValue(DataInGausspoint.previousPhase, (double)gpsij.GetValue(DataInGausspoint.currentPhase));
            //phie = GetEachVariableLocal(Result.PHASEFIELD);
            //gpsij.SetValue(DataInGausspoint.currentPhase, NMathFunctions.Dot(N, phie));
            //previousPhase = (double)gpsij.GetValue(DataInGausspoint.previousPhase);
            //currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);

            //if (currentPhase < previousPhase)
            //{
            //  gpsij.SetValue(DataInGausspoint.currentPhase, previousPhase);
            //  currentPhase = previousPhase;
            //}
            //if (currentPhase < 0)
            //  gpsij.SetValue(DataInGausspoint.currentPhase, 0.0);
            //if (currentPhase > 1)
            //  gpsij.SetValue(DataInGausspoint.currentPhase, 1.0);
            //currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);
            //strain = ComputeStrainGauss(patch, Bu, gpsij);
            //gpsij.SetValue(DataInGausspoint.currentStrain, strain);
            //ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
            //DoubleVector stress = ComputeStress(lamda, mu, currentPhase, principleStrain, eigenVectors);
            //gpsij.SetValue(DataInGausspoint.currentStress, stress);

            //fiU += a * gps[i, j].weight * detJbar * detJ * NMathFunctions.Product(Bu.Transpose(), stress);

            //xiEpsilon = (double)gpsij.GetValue(DataInGausspoint.xiEpsilon);
            //if (ModelStructurePhaseFieldStatic.TypeOrder == TypeOrderPhaseField.SecondOrder)
            //{
            //  double valDOmegaD = dOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
            //  double valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
            //  double valDAlphaD = dAlphaD(currentPhase, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
            //  F3 = (valDAlphaD * gc / (l0 * valCAlpha) + valDOmegaD * xiEpsilon) * N;
            //  F1 = 2.0 * gc * l0 / valCAlpha * NMathFunctions.Product(Bp.Transpose(), Bp);

            //  //F3 = (gc / l0 * currentPhase - 2 * (1 - currentPhase) * xiEpsilon) * N;
            //  //F1 = gc * l0 * NMathFunctions.Product(Bp.Transpose(), Bp);
            //}
            //else
            //{
            //  F3 = (gc / l0 * currentPhase - 2 * (1 - currentPhase) * xiEpsilon) * N;
            //  DoubleMatrix ddNdX = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
            //  DoubleVector BB = new DoubleVector(nen);
            //  for (int k = 0; k < nen; k++)
            //  {
            //    BB[k] = ddNdX[k, 0] + ddNdX[k, 1];
            //  }
            //  F1 = gc * l0 / 2.0 * NMathFunctions.Product(Bp.Transpose(), Bp) + gc * Math.Pow(l0, 3) / 16.0 * NMathFunctions.OuterProduct(BB, BB);
            //}
            //F2 = NMathFunctions.Product(F1, phie);
            //fiP += a * gps[i, j].weight * detJbar * detJ * (F2 + F3);
            //#endregion
            //break;
            case TypeSchemeNonlinearSolver.SingleStaggered:////////////thay doi Wu degradation
              #region singlestaggered
              //gpsij.SetValue(DataInGausspoint.previousPhase, (double)gpsij.GetValue(DataInGausspoint.currentPhase));
              phie = GetEachVariableLocal(Result.PHASEFIELD);
              //gpsij.SetValue(DataInGausspoint.currentPhase, NMathFunctions.Dot(N, phie));
              //previousPhase = (double)gpsij.GetValue(DataInGausspoint.previousPhase);
              //currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);

              //if (currentPhase < previousPhase)
              //{
              //  gpsij.SetValue(DataInGausspoint.currentPhase, previousPhase);
              //  currentPhase = previousPhase;
              //}
              //if (currentPhase < 0)
              //  gpsij.SetValue(DataInGausspoint.currentPhase, 0.0);
              //if (currentPhase > 1)
              //  gpsij.SetValue(DataInGausspoint.currentPhase, 1.0);
              currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);

              //strain = ComputeStrainGauss(patch, Bu, gpsij);
              //gpsij.SetValue(DataInGausspoint.currentStrain, strain);
              //strain = (DoubleVector)gpsij.GetValue(DataInGausspoint.currentStrain);
              //ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
              //stress = ComputeStress(lamda, mu, currentPhase, principleStrain, eigenVectors);
              //gpsij.SetValue(DataInGausspoint.currentStress, stress);
              stress = (DoubleVector)gpsij.GetValue(DataInGausspoint.currentStress);
              fiU += a * gps[i, j].weight * detJbar * detJ * NMathFunctions.Product(Bu.Transpose(), stress);
              //double xiEpsilonPos, xiEpsilonNeg;
              xiEpsilon = (double)gpsij.GetValue(DataInGausspoint.xiEpsilon);
              //xiEpsilon = gps[i, j].xiEpsilon;
              valDOmegaD = dOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
              valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
              valDAlphaD = dAlphaD(currentPhase, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
              if (ModelStructurePhaseFieldStatic.typeOrderPhaseField == TypeOrderPhaseField.SecondOrder)
              {
                // l0 là kích thước của vết nứt

                F3 = (valDAlphaD * gc / (l0 * valCAlpha) + valDOmegaD * xiEpsilon) * N;
                F1 = 2.0 * gc * l0 / valCAlpha * NMathFunctions.Product(Bp.Transpose(), Bp);

                //F3 = (gc / l0 * currentPhase - 2 * (1 - currentPhase) * xiEpsilon) * N;
                //F1 = gc * l0 * NMathFunctions.Product(Bp.Transpose(), Bp);
              }
              else
              {
                //Borden(2014)
                //F3 = (gc / l0 * currentPhase - 2 * (1 - currentPhase) * xiEpsilon) * N;
                F3 = (valDAlphaD * gc / (l0 * valCAlpha) + valDOmegaD * xiEpsilon) * N;
                DoubleMatrix ddNdX = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
                DoubleVector BB = new DoubleVector(nen);
                for (int k = 0; k < nen; k++)
                {
                  BB[k] = ddNdX[k, 0] + ddNdX[k, 1];
                }
                //Weinberg (2015) == Borden (2014) divide by 2*l0/gc
                //F1 = gc * l0 / 2.0 * NMathFunctions.Product(Bp.Transpose(), Bp) + gc * Math.Pow(l0, 3) / 16.0 * NMathFunctions.OuterProduct(BB, BB);
                F1 = gc * l0 / valCAlpha * NMathFunctions.Product(Bp.Transpose(), Bp) + gc * Math.Pow(l0, 3) / (8.0 * valCAlpha) * NMathFunctions.OuterProduct(BB, BB);
              }
              F2 = NMathFunctions.Product(F1, phie);
              fiP += a * gps[i, j].weight * detJbar * detJ * (F2 + F3);
              #endregion
              break;
            case TypeSchemeNonlinearSolver.Staggered:
              #region staggered
              if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
              {
                //displacement
                //currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);
                //currentPhase = gps[i, j].currentPhase;
                //strain = ComputeStrainGauss(patch, Bu, gpsij);
                //strain = (DoubleVector)gpsij.GetValue(DataInGausspoint.currentStrain);
                ////gpsij.SetValue(DataInGausspoint.currentStrain, strain);
                //ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
                //stress = ComputeStress(lamda, mu, currentPhase, principleStrain, eigenVectors);
                //gpsij.SetValue(DataInGausspoint.currentStress, stress);
                //stress = gps[i, j].currentStress = ComputeStress(lamda, mu, k1, currentPhase, principleStrain, eigenVectors);
                stress = (DoubleVector)gpsij.GetValue(DataInGausspoint.currentStress);
                fiU += a * gps[i, j].weight * detJbar * detJ * NMathFunctions.Product(Bu.Transpose(), stress);
              }
              else
              {
                //phasefield             
                //gpsij.SetValue(DataInGausspoint.previousPhase, (double)gpsij.GetValue(DataInGausspoint.currentPhase));
                phie = GetEachVariableLocal(Result.PHASEFIELD);
                //gpsij.SetValue(DataInGausspoint.currentPhase, NMathFunctions.Dot(N, phie));
                //previousPhase = (double)gpsij.GetValue(DataInGausspoint.previousPhase);
                //currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);

                //if (currentPhase < previousPhase)
                //{
                //  gpsij.SetValue(DataInGausspoint.currentPhase, previousPhase);
                //  currentPhase = previousPhase;
                //}
                //if (currentPhase < 0)
                //  gpsij.SetValue(DataInGausspoint.currentPhase, 0.0);
                //if (currentPhase > 1)
                //  gpsij.SetValue(DataInGausspoint.currentPhase, 1.0);
                currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);
                //gps[i, j].previousPhase = gps[i, j].currentPhase;
                //phie = GetEachVariableLocal(Result.PHASEFIELD);
                //gps[i, j].currentPhase = NMathFunctions.Dot(N, phie);

                //if (gps[i, j].currentPhase < gps[i, j].previousPhase)
                //  gps[i, j].currentPhase = gps[i, j].previousPhase;
                //if (gps[i, j].currentPhase < 0)
                //  gps[i, j].currentPhase = 0;
                //if (gps[i, j].currentPhase > 1)
                //  gps[i, j].currentPhase = 1;
                //currentPhase = gps[i, j].currentPhase;
                xiEpsilon = (double)gpsij.GetValue(DataInGausspoint.xiEpsilon);
                //xiEpsilon = gps[i, j].xiEpsilon;

                valDOmegaD = dOmegaD(currentPhase, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                valDAlphaD = dAlphaD(currentPhase, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                if (ModelStructurePhaseFieldStatic.typeOrderPhaseField == TypeOrderPhaseField.SecondOrder)
                {
                  //// l0 là kích thước của vết nứt
                  //F3 = (gc / l0 * currentPhase - 2 * (1 - currentPhase) * xiEpsilon) * N;
                  //F1 = gc * l0 * NMathFunctions.Product(Bp.Transpose(), Bp);
                  F3 = (valDAlphaD * gc / (l0 * valCAlpha) + valDOmegaD * xiEpsilon) * N;
                  F1 = 2.0 * gc * l0 / valCAlpha * NMathFunctions.Product(Bp.Transpose(), Bp);
                }
                else
                {
                  //Borden(2014)
                  //F3 = (gc / l0 * currentPhase - 2 * (1 - currentPhase) * xiEpsilon) * N;
                  F3 = (valDAlphaD * gc / (l0 * valCAlpha) + valDOmegaD * xiEpsilon) * N;
                  DoubleMatrix ddNdX = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.ddNdX); //gps[i, j].ddNdX;
                  DoubleVector BB = new DoubleVector(nen);
                  for (int k = 0; k < nen; k++)
                  {
                    BB[k] = ddNdX[k, 0] + ddNdX[k, 1];
                  }
                  //Weinberg (2015) == Borden (2014) divide by 2*l0/gc
                  //F1 = gc * l0 / 2.0 * NMathFunctions.Product(Bp.Transpose(), Bp) + gc * Math.Pow(l0, 3) / 16.0 * NMathFunctions.OuterProduct(BB, BB);
                  F1 = gc * l0 / valCAlpha * NMathFunctions.Product(Bp.Transpose(), Bp) + gc * Math.Pow(l0, 3) / (8.0 * valCAlpha) * NMathFunctions.OuterProduct(BB, BB);
                }
                F2 = NMathFunctions.Product(F1, phie);
                fiP += a * gps[i, j].weight * detJbar * detJ * (F2 + F3);
              }
              #endregion
              break;
          }
        }
      }

      for (int i = 0; i < nen; i++)
      {
        switch (AbstractModel.TypeSchemeNonlinearSolver)
        {
          case TypeSchemeNonlinearSolver.Monolithic:
          case TypeSchemeNonlinearSolver.SingleStaggered:
            for (int ii = 0; ii < 2; ii++)//x y
            {
              fi[i * 3 + ii] = fiU[i * 2 + ii];
            }
            fi[i * 3 + 2] = fiP[i];//phi
            break;
          case TypeSchemeNonlinearSolver.Staggered:
            if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
            {
              for (int ii = 0; ii < 2; ii++)//x y
              {
                fi[i * 3 + ii] = fiU[i * 2 + ii];
              }
            }
            else
            {
              fi[i * 3 + 2] = fiP[i];//phi
            }
            break;
        }
      }
    }
    /// <summary>
    /// The energetic degradation function
    /// </summary>
    private double OmegaD(double d, double l0, double lch, TypeModelPhasefield type)
    {
      double val = 0;
      switch (type)
      {
        case TypeModelPhasefield.AT1:
        case TypeModelPhasefield.AT2:
          val = Math.Pow(1.0 - d, 2);
          break;
        case TypeModelPhasefield.Cubic:
          double s = 1e-4;
          val = s * (Math.Pow(1.0 - d, 3) - Math.Pow(1.0 - d, 2)) + 3 * Math.Pow(1.0 - d, 2) - 2 * Math.Pow(1.0 - d, 3);
          break;
        case TypeModelPhasefield.PFCZM_Linear:
        case TypeModelPhasefield.PFCZM_BLinear:
        case TypeModelPhasefield.PFCZM_Exponential:
        case TypeModelPhasefield.PFCZM_Hyperbolic:
        case TypeModelPhasefield.PFCZM_Cornelissen:
          double a1 = 4.0 * lch / (Math.PI * l0);
          DetermineCoefficientCZM(type, out double p, out double a2, out double a3);
          double fac1 = Math.Pow(1 - d, p);
          double fac2 = fac1 + a1 * d + a1 * a2 * d * d + a1 * a2 * a3 * d * d * d;

          val = fac1 / fac2;
          break;
      }
      return val;
    }
    private void DetermineCoefficientCZM(TypeModelPhasefield type, out double p, out double a2, out double a3)
    {
      switch (type)
      {
        case TypeModelPhasefield.PFCZM_Linear:
          p = 2;
          a2 = -0.5;
          a3 = 0;
          break;
        case TypeModelPhasefield.PFCZM_BLinear:
          p = 2;
          a2 = 0.03687;
          a3 = 20.8343;
          break;
        case TypeModelPhasefield.PFCZM_Exponential:
          p = 2.5;
          a2 = Math.Pow(2.0, 5.0 / 3.0) - 3.0;
          a3 = 0;
          break;
        case TypeModelPhasefield.PFCZM_Hyperbolic:
          p = 4.0;
          a2 = Math.Pow(2.0, 7.0 / 3.0) - 4.5;
          a3 = 0;
          break;
        case TypeModelPhasefield.PFCZM_Cornelissen:
          p = 2.0;
          a2 = 1.3868;
          a3 = 0.6567;//0.9106;
          break;
        default:
          p = 0;
          a2 = 0;
          a3 = 0;
          break;
      }
    }
    private double dOmegaD(double d, double l0, double lch, TypeModelPhasefield type)
    {
      double val = 0;
      switch (type)
      {
        case TypeModelPhasefield.AT1:
        case TypeModelPhasefield.AT2:
          val = -2.0 * (1.0 - d);
          break;
        case TypeModelPhasefield.Cubic:
          double s = 1e-4;
          val = s * (-3 * Math.Pow(1.0 - d, 2) + 2 * (1.0 - d)) - 6 * (1.0 - d) + 6 * Math.Pow(1.0 - d, 2);
          break;
        case TypeModelPhasefield.PFCZM_Linear:
        case TypeModelPhasefield.PFCZM_BLinear:
        case TypeModelPhasefield.PFCZM_Exponential:
        case TypeModelPhasefield.PFCZM_Hyperbolic:
        case TypeModelPhasefield.PFCZM_Cornelissen:
          double a1 = 4.0 * lch / (Math.PI * l0);
          DetermineCoefficientCZM(type, out double p, out double a2, out double a3);
          double fac1 = Math.Pow(1 - d, p);
          double dfac1 = -p * Math.Pow(1 - d, p - 1);

          double fac2 = fac1 + a1 * d + a1 * a2 * d * d + a1 * a2 * a3 * d * d * d;
          double dfac2 = dfac1 + a1 + 2.0 * a1 * a2 * d + 3.0 * a1 * a2 * a3 * d * d;

          val = (dfac1 * fac2 - fac1 * dfac2) / (fac2 * fac2);
          break;
      }
      return val;
    }
    private double ddOmegaD(double d, double l0, double lch, TypeModelPhasefield type)
    {
      double val = 0;
      switch (type)
      {
        case TypeModelPhasefield.AT1:
        case TypeModelPhasefield.AT2:
          val = 2.0;
          break;
        case TypeModelPhasefield.Cubic:
          double s = 1e-4;
          val = s * (6 * (1.0 - d) - 2.0) + 6.0 - 12 * (1.0 - d);
          break;
        case TypeModelPhasefield.PFCZM_Linear:
        case TypeModelPhasefield.PFCZM_BLinear:
        case TypeModelPhasefield.PFCZM_Exponential:
        case TypeModelPhasefield.PFCZM_Hyperbolic:
        case TypeModelPhasefield.PFCZM_Cornelissen:
          double a1 = 4.0 * lch / (Math.PI * l0);
          DetermineCoefficientCZM(type, out double p, out double a2, out double a3);
          double fac1 = Math.Pow(1 - d, p);
          double dfac1 = -p * Math.Pow(1 - d, p - 1);
          double ddfac1 = p * (p - 1) * Math.Pow(1 - d, p - 2);

          double fac2 = fac1 + a1 * d + a1 * a2 * d * d + a1 * a2 * a3 * d * d * d;
          double dfac2 = dfac1 + a1 + 2.0 * a1 * a2 * d + 3.0 * a1 * a2 * a3 * d * d;
          double ddfac2 = ddfac1 + 2.0 * a1 * a2 + 6.0 * a1 * a2 * a3 * d;

          val = ((ddfac1 * fac2 - fac1 * ddfac2) * fac2 - 2 * (dfac1 * fac2 - fac1 * dfac2) * dfac2) / (fac2 * fac2 * fac2);
          break;
      }
      return val;
    }

    private double alphaD(double d, TypeModelPhasefield type)
    {
      double xi = TypeOfAlphaFunction(type);
      return xi * d + (1.0 - xi) * d * d;
    }
    private double dAlphaD(double d, TypeModelPhasefield type)
    {
      double xi = TypeOfAlphaFunction(type);
      return xi + 2.0 * (1.0 - xi) * d;
    }
    private double ddAlphaD(TypeModelPhasefield type)
    {
      double xi = TypeOfAlphaFunction(type);
      return 2.0 * (1.0 - xi);
    }
    private double cAlpha(TypeModelPhasefield type)
    {
      double xi = TypeOfAlphaFunction(type);
      FunctionRToR alpha = new PolynomialFunctionRToR(0, xi, 1 - xi);

      FunctionRToR SqrtAlpha = new PowerFunctionRToR(4.0, 0.5, alpha);
      return SqrtAlpha.Integration(0, 1, 20);
    }
    private static double TypeOfAlphaFunction(TypeModelPhasefield type)
    {
      double xi = 0;
      switch (type)
      {
        case TypeModelPhasefield.AT1:
          xi = 1.0;
          break;
        case TypeModelPhasefield.AT2:
        case TypeModelPhasefield.Cubic:
          xi = 0;
          break;
        case TypeModelPhasefield.PFCZM_Linear:
        case TypeModelPhasefield.PFCZM_BLinear:
        case TypeModelPhasefield.PFCZM_Exponential:
        case TypeModelPhasefield.PFCZM_Hyperbolic:
        case TypeModelPhasefield.PFCZM_Cornelissen:
          xi = 2.0;
          break;
      }
      return xi;
    }

    internal override void UpdateGausspointValue()
    {
      PatchStructure2D patch = (PatchStructure2D)this.patch;
      int d = patch.GetCountField();
      BivariateNURBSBasisFunction basis = (BivariateNURBSBasisFunction)((NURBSSurface)(patch.GetGeometry())).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int nen = (p + 1) * (q + 1); // number of local basis functions
      DoubleMatrix Bu = new DoubleMatrix(4, 2 * nen);
      double lamda = 0, mu = 0, gc = 0, ft = 0;
      if (!(Material is FGMUserDefinedGradedMaterial))
      {
        lamda = Material.GetProperty(MaterialPropertyName.LameParameter).GetValueProperty();
        mu = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
        gc = Material.GetProperty(MaterialPropertyName.CriticalEnergyReleaseRate).GetValueProperty();
        ft = Material.GetProperty(MaterialPropertyName.TensileYieldStrength).GetValueProperty();
      }
      DoubleMatrix dNdx = null;
      DoubleVector N = null;
      DoubleVector ue = null;////////////
      DoubleVector strain = null;
      double previousPhase;
      double currentPhase;
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          int c = 0;
          double x1 = 0;
          GaussPoints gpsij = gps[i, j];
          if (Material is FGMUserDefinedGradedMaterial)
          {
            lamda = ComputeParameterProperty(MaterialPropertyName.LameParameter, gpsij.location[0], gpsij.location[1]);
            mu = ComputeParameterProperty(MaterialPropertyName.ShearModulus, gpsij.location[0], gpsij.location[1]);
            gc = ComputeParameterProperty(MaterialPropertyName.CriticalEnergyReleaseRate, gpsij.location[0], gpsij.location[1]);
            ft = ComputeParameterProperty(MaterialPropertyName.TensileYieldStrength, gpsij.location[0], gpsij.location[1]);
          }
          double E = mu * (3 * lamda + 2 * mu) / (lamda + mu);
          double lch = E * gc / (ft * ft);
          double xiEpsilon = 0;
          double xiEpsilonPos = 0;
          double xiEpsilonNeg = 0;
          DoubleVector principleStrain = null;
          DoubleVector[] eigenVectors = null;
          DoubleVector phie = null;
          double phiCurrent;
          switch (AbstractModel.TypeSchemeNonlinearSolver)
          {
            case TypeSchemeNonlinearSolver.Monolithic:
            case TypeSchemeNonlinearSolver.SingleStaggered:
            case TypeSchemeNonlinearSolver.Staggered:
              ////update after displacement step
              c = 0;
              x1 = 1;
              if (patch.StateStress == Structure2DState.Axisymetric)
              {
                c = 1;
                x1 = PointAt(gpsij.location[0], gpsij.location[1])[0];
              }
              dNdx = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.dNdX);//gps[i, j].dNdX;
              N = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);//gps[i, j].Ni;
              for (int k = 0; k < nen; k++)
              {
                //xx yy zz xy
                Bu[0, 2 * k] = dNdx[k, 0];
                Bu[1, 2 * k + 1] = dNdx[k, 1];
                Bu[2, 2 * k] = c * N[k] / x1;
                Bu[3, 2 * k] = dNdx[k, 1];
                Bu[3, 2 * k + 1] = dNdx[k, 0];
              }

              strain = ComputeStrainGauss(patch, Bu, gpsij);
              gpsij.SetValue(DataInGausspoint.currentStrain, strain);
              ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);

              switch (ModelStructurePhaseFieldStatic.typeDegradationModel)
              {
                case TypeDegradationEnergy.Isotropic:
                  xiEpsilon = ComputeXiEpsilonIsotropic(principleStrain, eigenVectors, lamda, mu);
                  break;
                case TypeDegradationEnergy.Hybrid:
                case TypeDegradationEnergy.Anisotropic:
                  ComputeXiEpsilonAnisotropic(principleStrain, eigenVectors, lamda, mu, out xiEpsilonPos, out xiEpsilonNeg, ft);
                  xiEpsilon = xiEpsilonPos;
                  break;
              }
              if (xiEpsilon > (double)gpsij.GetValue(DataInGausspoint.xiEpsilon))
              {
                gpsij.SetValue(DataInGausspoint.xiEpsilon, xiEpsilon);
              }
              phie = GetEachVariableLocal(Result.PHASEFIELD);
              phiCurrent = NMathFunctions.Dot(N, phie);
              if (phiCurrent < 0.0)
                phiCurrent = 0.0;
              if (phiCurrent > 1.0)
                phiCurrent = 1.0;

              //if (phiCurrent > currentPhase)
              //{
              //  gpsij.SetValue(DataInGausspoint.currentPhase, phiCurrent);
              //}
              gpsij.SetValue(DataInGausspoint.currentPhase, phiCurrent);
              currentPhase = phiCurrent;//(double)gpsij.GetValue(DataInGausspoint.currentPhase);
              DoubleVector stress = ComputeStress(lamda, mu, lch, currentPhase, principleStrain, eigenVectors);
              gpsij.SetValue(DataInGausspoint.currentStress, stress);
              break;
              //case TypeSchemeNonlinearSolver.Staggered:
              //  if (!AbstractNonlinearSolver.IsSolvePhasefieldStaggered)
              //  {
              //    ////update after displacement step
              //    c = 0;
              //    x1 = 1;
              //    if (patch.StateStress == Structure2DState.Axisymetric)
              //    {
              //      c = 1;
              //      x1 = PointAt(gpsij.location[0], gpsij.location[1])[0];
              //    }
              //    dNdx = (DoubleMatrix)gpsij.GetValue(DataInGausspoint.dNdX);//gps[i, j].dNdX;
              //    N = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);//gps[i, j].Ni;
              //    for (int k = 0; k < nen; k++)
              //    {
              //      //xx yy zz xy
              //      Bu[0, 2 * k] = dNdx[k, 0];
              //      Bu[1, 2 * k + 1] = dNdx[k, 1];
              //      Bu[2, 2 * k] = c * N[k] / x1;
              //      Bu[3, 2 * k] = dNdx[k, 1];
              //      Bu[3, 2 * k + 1] = dNdx[k, 0];
              //    }

              //    strain = ComputeStrainGauss(patch, Bu, gpsij);
              //    gpsij.SetValue(DataInGausspoint.currentStrain, strain);
              //    ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
              //    switch (ModelStructurePhaseFieldStatic.TypeDegradationModel)
              //    {
              //      case TypeDegradationEnergy.Isotropic:
              //        xiEpsilon = ComputeXiEpsilonIsotropic(principleStrain, eigenVectors, lamda, mu);
              //        break;
              //      case TypeDegradationEnergy.Hybrid:
              //      case TypeDegradationEnergy.Anisotropic:
              //        ComputeXiEpsilonAnisotropic(principleStrain, eigenVectors, lamda, mu, out xiEpsilonPos, out xiEpsilonNeg);
              //        xiEpsilon = xiEpsilonPos;
              //        break;
              //    }
              //    if (xiEpsilon > (double)gpsij.GetValue(DataInGausspoint.xiEpsilon))
              //    {
              //      gpsij.SetValue(DataInGausspoint.xiEpsilon, xiEpsilon);
              //    }
              //  }
              //  else
              //  {
              //    //update after phasefield step
              //    N = (DoubleVector)gpsij.GetValue(DataInGausspoint.Ni);//gps[i, j].Ni;

              //    gpsij.SetValue(DataInGausspoint.previousPhase, gpsij.GetValue(DataInGausspoint.currentPhase));
              //    phie = GetEachVariableLocal(Result.PHASEFIELD);
              //    phiCurrent = NMathFunctions.Dot(N, phie);
              //    gpsij.SetValue(DataInGausspoint.currentPhase, phiCurrent);
              //    previousPhase = (double)gpsij.GetValue(DataInGausspoint.previousPhase);
              //    currentPhase = (double)gpsij.GetValue(DataInGausspoint.currentPhase);

              //    if (currentPhase < previousPhase)
              //    {
              //      gpsij.SetValue(DataInGausspoint.currentPhase, previousPhase);
              //      currentPhase = previousPhase;
              //    }
              //    if (currentPhase < 0.0)
              //      gpsij.SetValue(DataInGausspoint.currentPhase, 0.0);
              //    if (currentPhase > 1.0)
              //      gpsij.SetValue(DataInGausspoint.currentPhase, 1.0);
              //  }
              //  break;
          }
        }
      }
    }

    internal void InitialPreCrackGausspointValue(CrackEdge crack)
    {
      double lamda = 0, mu = 0, gc = 0, ft = 0;
      if (!(Material is FGMUserDefinedGradedMaterial))
      {
        lamda = Material.GetProperty(MaterialPropertyName.LameParameter).GetValueProperty();
        mu = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
        gc = Material.GetProperty(MaterialPropertyName.CriticalEnergyReleaseRate).GetValueProperty();
        ft = Material.GetProperty(MaterialPropertyName.TensileYieldStrength).GetValueProperty();
      }
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          GaussPoints gpsij = gps[i, j];
          double[] gaussPosition = PointAt(gpsij.location[0], gpsij.location[1]);
          double d = crack.DistanceFromPointToCrack(gaussPosition[0], gaussPosition[1], 0, out double lengthFromXiToCrackTip);
          if (Material is FGMUserDefinedGradedMaterial)
          {
            gc = ComputeParameterProperty(MaterialPropertyName.CriticalEnergyReleaseRate, gpsij.location[0], gpsij.location[1]);
            lamda = ComputeParameterProperty(MaterialPropertyName.LameParameter, gps[i, j].location[0], gps[i, j].location[1]);
            mu = ComputeParameterProperty(MaterialPropertyName.ShearModulus, gps[i, j].location[0], gps[i, j].location[1]);
            ft = ComputeParameterProperty(MaterialPropertyName.TensileYieldStrength, gpsij.location[0], gpsij.location[1]);
          }
          double l0 = ModelStructurePhaseFieldStatic.l0;
          //if (d > l0 / 4.0)
          //{
          //  gpsij.SetValue(DataInGausspoint.currentPhase, Math.Max(0.0, (double)gpsij.GetValue(DataInGausspoint.currentPhase)));
          //}
          //else
          //{
          //  gpsij.SetValue(DataInGausspoint.currentPhase, 1.0);
          //}



          //////////////////// giải p,u extrapolation
          //if (d > 2 * l0)
          //{
          //  gps[i, j].currentPhase = 0;
          //}
          //else
          //{
          //  gps[i, j].currentPhase = Math.Exp(-2 * d / l0);
          //}
          //////////////////// giải p,u không extrapolation hội tụ rất chậm
          ////////////////////// giải u,p không extrapolation, u 29 vòng lặp, p lặp rất nhiều
          //gps[i, j].currentPhase = Math.Exp(-d / l0);
          //////////////////// giải u,p dùng extrapolation, u 29 vòng lặp, p lặp rất nhiều ko hội tụ 1000 lần
          //gps[i, j].SetValue(DataInGausspoint.currentPhase, Math.Exp(-d / l0));///
          //gps[i, j].SetValue(DataInGausspoint.previousPhase, Math.Exp(-d / l0));///
          //gps[i, j].SetValue(DataInGausspoint.lastPhase, Math.Exp(-d / l0));///

          //gps[i, j].currentPhase = Math.Exp(-d / l0);
          //////////////////// giải u,p dùng extrapolation, u 29 vòng lặp, p lặp rất nhiều ko hội tụ 1000 lần
          //gps[i, j].currentPhase = Math.Exp(-2 * d / l0);
          //gps[i, j].SetValue(DataInGausspoint.currentPhase, Math.Exp(-2 * d / l0));///
          //gps[i, j].SetValue(DataInGausspoint.previousPhase, Math.Exp(-2 * d / l0));///
          //gps[i, j].SetValue(DataInGausspoint.lastPhase, Math.Exp(-2 * d / l0));///
          ////////////////////// giải u,p dùng extrapolation, u 29 vòng lặp, p lặp 194 lần, FR=0.02215
          //if (d > l0 / 2.0)
          //  gps[i, j].xiEpsilon = 0;
          //else
          //{
          //  double c = 0.999;
          //  double B = c / (1 - c);
          //  double h = B * gc / (2 * l0) * (1 - 2 * d / l0);
          //  gps[i, j].xiEpsilon = h;
          //}
          //gps[i, j].currentPhase = Math.Exp(-d / l0);

          ////////////////////// giải u,p dùng extrapolation, u 29 vòng lặp, p lặp 234 lần, FR=0.02215
          //if (d > l0 / 2.0)
          //  gps[i, j].xiEpsilon = 0;
          //else
          //{
          //  double c = 0.9999;
          //  double B = c / (1 - c);
          //  double h = B * gc / (2 * l0) * (1 - 2 * d / l0);
          //  gps[i, j].xiEpsilon = h;
          //}
          //gps[i, j].currentPhase = Math.Exp(-d / l0);
          ////////////////////// giải u,p dùng extrapolation, u 20 vòng lặp, p lặp 325 lần, FR=0.02356
          //// những lần lặp sau ra FR đúng với kết quả nhưng lặp p 200 lần
          //if (d > l0 / 2.0)
          //  gps[i, j].xiEpsilon = 0;
          //else
          //{
          //  double c = 0.9999;
          //  double B = c / (1 - c);
          //  double h = B * gc / (2 * l0) * (1 - 2 * d / l0);
          //  gps[i, j].xiEpsilon = h;
          //}
          //gps[i, j].currentPhase = Math.Exp(-2 * d / l0);
          //////////////////// giải u,p không dùng extrapolation, chưa khảo sát
          //if (d > l0 / 2.0)
          //  gps[i, j].xiEpsilon = 0;
          //else
          //{
          //  double c = 0.9999;
          //  double B = c / (1 - c);
          //  double h = B * gc / (2 * l0) * (1 - 2 * d / l0);
          //  gps[i, j].xiEpsilon = h;
          //}
          //////////////////// giải u,p không dùng extrapolation, u 2 vòng lặp, p chưa test
          //if (d > l0 / 2.0)
          //  gps[i, j].xiEpsilon = 0;
          //else
          //{
          //  double c = 0.9999;
          //  double B = c / (1 - c);
          //  double h = B * gc / (2 * l0) * (1 - 2 * d / l0);
          //  double f = lengthFromXiToCrackTipMin;
          //  gps[i, j].xiEpsilon = h * f;
          //}

          ////////////////////// giải u,p không dùng extrapolation, u 2 vòng lặp, p 772 vong lap, FR=0.02428, quá cứng
          //if (d > l0 / 2.0)
          //  gps[i, j].xiEpsilon = 0;
          //else
          //{
          //  double c = 0.999;
          //  double B = c / (1 - c);
          //  double h = B * gc / (2 * l0) * (1 - 2 * d / l0);
          //  double f = lengthFromXiToCrackTipMin;
          //  gps[i, j].xiEpsilon = h * f;
          //}
          ////////////////////// giải p,u không dùng extrapolation, p 972 vòng lặp, u 42 vong lap, FR=0.0141367
          ///////good
          //if (d > l0 / 2.0)
          //  gps[i, j].xiEpsilon = 0;
          //else
          //{
          //  double c = 0.9999;
          //  double B = c / (1 - c);
          //  double h = B * gc / (2 * l0) * (1 - 2 * d / l0);
          //  double f = lengthFromXiToCrackTipMin;
          //  gps[i, j].xiEpsilon = h * f;
          //}
          //////////////////// giải p,u không dùng extrapolation, p 393 vòng lặp, u 30 vong lap, FR=0.01388
          /////////////////////////////////////////////////////////////
          /////////////////////////////////////////////////////////////
          ////////////////////////////////////////////////////////////////
          /////////////////////////////////////////////////////////////
          ///////good////đang chạy cho bài tension
          if (d > l0 / 2.0)
          {
            gpsij.SetValue(DataInGausspoint.xiEpsilon, Math.Max(0.0, (double)gpsij.GetValue(DataInGausspoint.xiEpsilon)));
          }
          else
          {
            double c = 0.9999;
            double B = c / (1 - c);
            double h = B * gc / (2 * l0) * (1 - 2 * d / l0);
            gpsij.SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gpsij.GetValue(DataInGausspoint.xiEpsilon)));
          }

          /////not good cho bê tông////cải tiến cho vật liệu bê tông
          //if (d > l0 / 2.0)
          //{
          //  gpsij.SetValue(DataInGausspoint.xiEpsilon, Math.Max(0.0, (double)gpsij.GetValue(DataInGausspoint.xiEpsilon)));
          //}
          //else
          //{
          //  double E = mu * (3 * lamda + 2 * mu) / (lamda + mu);
          //  double lch = (ft > 0) ? E * gc / (ft * ft) : 0;
          //  double c = 0.9999;
          //  double valDAlphaD = dAlphaD(c, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
          //  double valDOmegaD = dOmegaD(c, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
          //  double valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
          //  double h = -valDAlphaD * gc / (valCAlpha * l0 * valDOmegaD) * (1.0 - 2.0 * d / l0);
          //  //double h = ft * ft / (2 * E) * (1.0 - 2.0 * d / l0);
          //  //h = Math.Max(h, ft * ft / (2 * E));
          //  gpsij.SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gpsij.GetValue(DataInGausspoint.xiEpsilon)));
          //}
          /////////////////////////////////////////////////////////////
          /////////////////////////////////////////////////////////////
          ////////////////////////////////////////////////////////////////
          /////////////////////////////////////////////////////////////
          //if (d > l0 / 2.0)
          //  gps[i, j].xiEpsilon = Math.Max(0, gps[i, j].xiEpsilon);
          //else
          //{
          //  double B = 1000;
          //  double h = B * gc / (2 * l0) * (1 - 2 * d / l0);
          //  gps[i, j].xiEpsilon = Math.Max(h, gps[i, j].xiEpsilon);
          //}
          ///////////////////////////////
          /////////// Hội tụ tạm, kết quả phase sai, FR lớn 0.017, chuẩn 0.014, thứ tự p, u
          //gps[i, j].currentPhase = Math.Exp(-d / l0);
          //gps[i, j].xiEpsilon = gc / (2 * l0) / (1 - gps[i, j].currentPhase);

          //////////////phase 309 vòng lặp, displacement 61 vòng lặp, Fr=0.01763, thứ tự p, u
          //gps[i, j].xiEpsilon = gc / (2 * l0) / (1 - Math.Exp(-d / l0));
          /////////////////////
          ////////////phase 309 vòng lặp, displacement 61 vòng lặp, Fr=0.0189
          //gps[i, j].xiEpsilon = gc / (2 * l0) / (1 - Math.Exp(-2 * d / l0));
          ///////////////////phase 733 vòng lặp, 
          //double c = 0.999;
          //double B = c / (1 - c);
          //double h = B * gc / (2 * l0) / (1 - Math.Exp(-2 * d / l0));
          //double f = lengthFromXiToCrackTipMin;
          //gps[i, j].xiEpsilon = h * f;
          ///////////////////////////////////////////////////////////////////

          //if (d > l0 / 4.0)
          //  gps[i, j].xiEpsilon = 0;
          //else
          //{
          //  double c = 0.999;
          //  double B = c / (1 - c);
          //  //if (ModelStructurePhaseFieldStatic.TypeOrder == TypeOrderPhaseField.SecondOrder)
          //  //{
          //  double h = B * gc / (2 * l0) * (1 - 4.0 * d / l0);
          //  double f = lengthFromXiToCrackTipMin;
          //  gps[i, j].xiEpsilon = h * f;//now
          //  //gps[i, j].xiEpsilon = B * gc / (2 * l0) * Math.Exp(-2 * d / l0);
          //  //double c = Math.Exp(-2 * d / l0);
          //  //double B = c / (1 - c);
          //  ////if (ModelStructurePhaseFieldStatic.TypeOrder == TypeOrderPhaseField.SecondOrder)
          //  ////{
          //  ////gps[i, j].xiEpsilon = B * gc / (2 * l0) * (1 - 2 * d / l0);
          //  //gps[i, j].xiEpsilon = B * gc / (2 * l0);
          //  //}
          //  //else if (ModelStructurePhaseFieldStatic.TypeOrder == TypeOrderPhaseField.FourthOrder)
          //  //{
          //  //  gps[i, j].xiEpsilon = B * gc / l0 / (1 + 2 * d / l0);
          //  //}
          //}
        }
      }
    }

    internal void InitialVoidGausspointValue(Void2D v)
    {
      double lamda = 0, mu = 0, gc = 0, ft = 0;
      if (!(Material is FGMUserDefinedGradedMaterial))
      {
        lamda = Material.GetProperty(MaterialPropertyName.LameParameter).GetValueProperty();
        mu = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
        gc = Material.GetProperty(MaterialPropertyName.CriticalEnergyReleaseRate).GetValueProperty();
        ft = Material.GetProperty(MaterialPropertyName.TensileYieldStrength).GetValueProperty();
      }
      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {

          GaussPoints gpsij = gps[i, j];
          double[] gaussPosition = PointAt(gpsij.location[0], gpsij.location[1]);
          double l0 = ModelStructurePhaseFieldStatic.l0;
          if (Material is FGMUserDefinedGradedMaterial)
          {
            gc = ComputeParameterProperty(MaterialPropertyName.CriticalEnergyReleaseRate, gpsij.location[0], gpsij.location[1]);
            lamda = ComputeParameterProperty(MaterialPropertyName.LameParameter, gps[i, j].location[0], gps[i, j].location[1]);
            mu = ComputeParameterProperty(MaterialPropertyName.ShearModulus, gps[i, j].location[0], gps[i, j].location[1]);
            ft = ComputeParameterProperty(MaterialPropertyName.TensileYieldStrength, gpsij.location[0], gpsij.location[1]);
          }
          double E = mu * (3 * lamda + 2 * mu) / (lamda + mu);
          //double d = 0;
          //if (!v.IsInside(gaussPosition[0], gaussPosition[1], 0))
          //{
          //  d = v.DistanceFromOutsidePointToVoid(gaussPosition[0], gaussPosition[1], 0);

          //  if (Material is FGMUserDefinedGradedMaterial)
          //  {
          //    gc = ComputeParameterProperty(MaterialPropertyName.CriticalEnergyReleaseRate, gps[i, j].location[0], gps[i, j].location[1]);
          //  }
          //  if (d > l0 / 2.0)
          //  {
          //    gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(0.0, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
          //    //gps[i, j].xiEpsilon = Math.Max(0, gps[i, j].xiEpsilon);
          //  }
          //  else
          //  {
          //    double c = 0.9999;
          //    double B = c / (1 - c);
          //    double h = B * gc / (2 * l0) * (1 - 2 * d / l0);
          //    gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
          //    //gps[i, j].xiEpsilon = Math.Max(h, gps[i, j].xiEpsilon);
          //  }
          //  //gps[i, j].SetValue(DataInGausspoint.currentPhase, Math.Exp(-2 * d / l0));
          //}
          //else
          //{
          //  double c = 0.9999;
          //  double B = c / (1 - c);
          //  double h = B * gc / (2 * l0);
          //  gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
          //}

          double d = v.DistanceWithSignFromPointToVoid(gaussPosition[0], gaussPosition[1], 0);
          switch (ModelStructurePhaseFieldStatic.typePhaseFieldModel)
          {
            case TypeModelPhasefield.AT1:
            case TypeModelPhasefield.AT2:

              /////////////////////////////////////////////////////////////////////
              //if (d > l0 / 2.0)
              //{
              //  gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(0.0, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
              //}
              //else if (d <= 0)
              //{
              //  double c = 0.9999;
              //  double B = c / (1 - c);
              //  double h = B * gc / (2 * l0);
              //  gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
              //}
              //else
              //{
              //  double c = 0.9999;
              //  double B = c / (1 - c);
              //  double h = B * gc / (2 * l0) * (1 - 2 * d / l0);
              //  gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
              //}
              /////////////////////////////////////////////////////////////////////
              if (d > 0)
              {
                gpsij.SetValue(DataInGausspoint.xiEpsilon, Math.Max(0.0, (double)gpsij.GetValue(DataInGausspoint.xiEpsilon)));
              }
              else if (d <= -l0 / 2.0)
              {
                double c = 0.9999;
                double B = c / (1 - c);
                double h = B * gc / (2 * l0);
                gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
              }
              else
              {
                double c = 0.9999;
                double B = c / (1 - c);
                double h = B * gc / (2 * l0) * (1 - 2 * d / l0 - 1);
                gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
              }
              /////////////////////////////////////////////////////////////////////
              //if (d > 0)
              //{
              //  gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(0.0, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
              //}
              //else if (d <= -l0)
              //{
              //  double c = 0.9999;
              //  double B = c / (1 - c);
              //  double h = B * gc / (2 * l0);
              //  gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
              //}
              //else
              //{
              //  double c = 0.9999;
              //  double B = c / (1 - c);
              //  double h = B * gc / (2 * l0) * (-2 * d / l0 - 1);
              //  gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
              //}
              break;
            case TypeModelPhasefield.Cubic:
              if (d > 0)
              {
                gpsij.SetValue(DataInGausspoint.xiEpsilon, Math.Max(0.0, (double)gpsij.GetValue(DataInGausspoint.xiEpsilon)));
              }
              else if (d <= -l0 / 2.0)
              {
                double lch = (ft > 0) ? E * gc / (ft * ft) : 0;
                double c = 0.9999;
                double valDAlphaD = dAlphaD(c, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                double valDOmegaD = dOmegaD(c, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                double valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                double h = -valDAlphaD * gc / (valCAlpha * l0 * valDOmegaD);
                gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
              }
              else
              {
                double lch = (ft > 0) ? E * gc / (ft * ft) : 0;
                double c = 0.9999;
                double valDAlphaD = dAlphaD(c, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                double valDOmegaD = dOmegaD(c, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                double valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                double h = -valDAlphaD * gc / (valCAlpha * l0 * valDOmegaD) * (-2 * d / l0);
                gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
              }
              break;
            case TypeModelPhasefield.PFCZM_Linear:
            case TypeModelPhasefield.PFCZM_BLinear:
            case TypeModelPhasefield.PFCZM_Exponential:
            case TypeModelPhasefield.PFCZM_Hyperbolic:
            case TypeModelPhasefield.PFCZM_Cornelissen:
              if (d > 0)
              {
                double h = ft * ft / (2 * E);
                gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
              }
              else if (d <= -l0 / 2.0)
              {
                double lch = (ft > 0) ? E * gc / (ft * ft) : 0;
                double c = 0.9999;
                double valDAlphaD = dAlphaD(c, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                double valDOmegaD = dOmegaD(c, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                double valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                double h = -valDAlphaD * gc / (valCAlpha * l0 * valDOmegaD);
                h = Math.Max(h, ft * ft / (2 * E));
                gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));

              }
              else
              {
                double lch = (ft > 0) ? E * gc / (ft * ft) : 0;
                double c = 0.9999;
                double valDAlphaD = dAlphaD(c, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                double valDOmegaD = dOmegaD(c, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                double valCAlpha = cAlpha(ModelStructurePhaseFieldStatic.typePhaseFieldModel);
                double h = -valDAlphaD * gc / (valCAlpha * l0 * valDOmegaD) * (-2 * d / l0);
                h = Math.Max(h, ft * ft / (2 * E));
                gps[i, j].SetValue(DataInGausspoint.xiEpsilon, Math.Max(h, (double)gps[i, j].GetValue(DataInGausspoint.xiEpsilon)));
              }
              break;
          }
        }
      }
    }
    //private double MinDistance(List<CrackEdge> cracks, double x, double y, out double lengthRatioFromXiToCrackTipMin)
    //{
    //  double d = 1e9;
    //  double lengthFromXiToCrackTipMin = 0;
    //  lengthRatioFromXiToCrackTipMin = 0;
    //  foreach (var item in cracks)
    //  {
    //    double min = item.DistanceFromPointToCrack(x, y, 0, out double lengthFromXiToCrackTip);
    //    if (min < d)
    //    {
    //      d = min;
    //      lengthFromXiToCrackTipMin = lengthFromXiToCrackTip;
    //      lengthRatioFromXiToCrackTipMin = lengthFromXiToCrackTipMin / item.ComputeLength();
    //    }
    //  }
    //  return d;
    //}
    private double ComputeXiEpsilonIsotropic(DoubleVector principleStrain, DoubleVector[] eigenVectors, double lamda, double mu)
    {
      //var e2 = 1.0 / 2.0 * lamda * Math.Pow(principleStrain[0] + principleStrain[1] + principleStrain[2], 2) + mu * (Math.Pow(principleStrain[0], 2) + Math.Pow(principleStrain[1], 2) + Math.Pow(principleStrain[2], 2));
      //var e3 = 1.0 / 2.0 * lamda * Math.Pow(strain[0] + strain[1] + strain[2], 2) + mu * (Math.Pow(strain[0], 2) + Math.Pow(strain[1], 2) + Math.Pow(strain[2], 2) + 2 * Math.Pow(strain[3] / 2.0, 2));
      //DoubleMatrix D = CreateMaterialMatrix();
      //var e1 = 1.0 / 2.0 * NMathFunctions.Dot(strain, NMathFunctions.Product(D, strain));
      //return 1.0 / 2.0 * NMathFunctions.Dot(strain, NMathFunctions.Product(D, strain));
      return 1.0 / 2.0 * lamda * Math.Pow(principleStrain[0] + principleStrain[1] + principleStrain[2], 2) + mu * (Math.Pow(principleStrain[0], 2) + Math.Pow(principleStrain[1], 2) + Math.Pow(principleStrain[2], 2));
      //return 1.0 / 2.0 * lamda * Math.Pow(strain[0] + strain[1] + strain[2], 2) + mu * (Math.Pow(strain[0], 2) + Math.Pow(strain[1], 2) + Math.Pow(strain[2], 2) + 2 * Math.Pow(strain[3], 2) / 2.0);
    }

    private void ComputeXiEpsilonAnisotropic(DoubleVector principleStrain, DoubleVector[] eigenVectors, double lamda, double mu, out double xiEpsilonPositive, out double xiEpsilonNegative, double ft)
    {
      xiEpsilonPositive = xiEpsilonNegative = 0;
      switch (ModelStructurePhaseFieldStatic.typePhaseFieldModel)
      {
        case TypeModelPhasefield.AT1:
        case TypeModelPhasefield.AT2:
        case TypeModelPhasefield.Cubic:
          double trEpsilon = principleStrain[0] + principleStrain[1] + principleStrain[2];

          DoubleMatrix epPos = PositiveRampFunction(principleStrain[0]) * NMathFunctions.OuterProduct(eigenVectors[0], eigenVectors[0])
            + PositiveRampFunction(principleStrain[1]) * NMathFunctions.OuterProduct(eigenVectors[1], eigenVectors[1])
            + PositiveRampFunction(principleStrain[2]) * NMathFunctions.OuterProduct(eigenVectors[2], eigenVectors[2]);
          DoubleMatrix epNeg = NegativeRampFunction(principleStrain[0]) * NMathFunctions.OuterProduct(eigenVectors[0], eigenVectors[0])
            + NegativeRampFunction(principleStrain[1]) * NMathFunctions.OuterProduct(eigenVectors[1], eigenVectors[1])
            + NegativeRampFunction(principleStrain[2]) * NMathFunctions.OuterProduct(eigenVectors[2], eigenVectors[2]);

          DoubleMatrix epPosSquare = NMathFunctions.Product(epPos, epPos);
          DoubleMatrix epNegSquare = NMathFunctions.Product(epNeg, epNeg);

          xiEpsilonPositive = lamda / 2.0 * Math.Pow(PositiveRampFunction(trEpsilon), 2) + mu * (epPosSquare[0, 0] + epPosSquare[1, 1] + epPosSquare[2, 2]);
          xiEpsilonNegative = lamda / 2.0 * Math.Pow(NegativeRampFunction(trEpsilon), 2) + mu * (epNegSquare[0, 0] + epNegSquare[1, 1] + epNegSquare[2, 2]);
          break;
        case TypeModelPhasefield.PFCZM_Linear:
        case TypeModelPhasefield.PFCZM_BLinear:
        case TypeModelPhasefield.PFCZM_Exponential:
        case TypeModelPhasefield.PFCZM_Hyperbolic:
        case TypeModelPhasefield.PFCZM_Cornelissen:
          double E = mu * (3.0 * lamda + 2.0 * mu) / (lamda + mu);
          DoubleVector stress = ComputeEfficientStress(lamda, mu, principleStrain, eigenVectors);
          double smax = FailureCriterion(ModelStructurePhaseFieldStatic.typeFailureCriterion, stress);
          xiEpsilonPositive = Math.Max(0.5 * Math.Pow(ft, 2.0) / E, 0.5 * Math.Pow(smax, 2.0) / E);
          break;
      }
    }

    private double FailureCriterion(TypeFailureCriterion type, DoubleVector stress)
    {
      double sigmaEq = 0;
      switch (type)
      {
        case TypeFailureCriterion.Rankine:
          ComputePrincipleStress(stress, out DoubleVector eigValue, out DoubleVector[] eigVector);
          sigmaEq = PositiveRampFunction(eigValue.Max());
          break;
        case TypeFailureCriterion.ModifiedVonMises:
          double rhoC = 10.0;
          double I1 = stress[0] + stress[1] + stress[2];
          double J2 = 1.0 / 6.0 * (Math.Pow(stress[0] - stress[1], 2) + Math.Pow(stress[1] - stress[2], 2) + Math.Pow(stress[2] - stress[0], 2)) + Math.Pow(stress[3], 2);
          sigmaEq = (rhoC - 1.0) / (2.0 * rhoC) * I1 + 1.0 / (2.0 * rhoC) * Math.Sqrt(Math.Pow(rhoC - 1, 2) * I1 * I1 + 12.0 * rhoC * J2);
          break;
        case TypeFailureCriterion.Wu2017:
          double beta = 10.0 - 1.0;//fc/ft=10.0
          J2 = 1.0 / 6.0 * (Math.Pow(stress[0] - stress[1], 2) + Math.Pow(stress[1] - stress[2], 2) + Math.Pow(stress[2] - stress[0], 2)) + Math.Pow(stress[3], 2);
          ComputePrincipleStress(stress, out DoubleVector eigValue1, out DoubleVector[] eigVector1);
          double sigmaMax = PositiveRampFunction(eigValue1.Max());
          sigmaEq = 1.0 / (1.0 + beta) * (beta * sigmaMax + Math.Sqrt(3.0 * J2));
          break;
      }
      return sigmaEq;
    }

    /// <summary>
    /// Miehe 2010 A phase field model for rate-independent crack propagation: Robust algorithmic implementation based on operator splits
    /// </summary>
    /// <param name="x"></param>
    /// <returns></returns>
    private double PositiveRampFunction(double x)
    {
      //return (Math.Abs(x) + x) / 2.0;
      if (x > 0)
        return x;
      else
        return 0;
    }
    private double NegativeRampFunction(double x)
    {
      //return (Math.Abs(x) - x) / 2.0;
      if (x < 0)
        return x;
      else
        return 0;
    }

    private void ComputeStressIsotropic(DoubleVector principleStrain, DoubleVector[] eigenVectors, double lamda, double mu, out DoubleVector stress0)
    {
      var prinVector1 = eigenVectors[0];
      var prinVector2 = eigenVectors[1];
      var prinVector3 = eigenVectors[2];

      double trEpsilon = principleStrain[0] + principleStrain[1] + principleStrain[2];
      //DoubleMatrix unit = new DoubleMatrix(3, 3);
      //unit[0, 0] = unit[1, 1] = unit[2, 2] = 1.0;
      DoubleMatrix stressTensor = (lamda * trEpsilon + 2 * mu * principleStrain[0]) * NMathFunctions.OuterProduct(prinVector1, prinVector1)
        + (lamda * trEpsilon + 2 * mu * principleStrain[1]) * NMathFunctions.OuterProduct(prinVector2, prinVector2)
        + (lamda * trEpsilon + 2 * mu * principleStrain[2]) * NMathFunctions.OuterProduct(prinVector3, prinVector3);
      stress0 = new DoubleVector(4);
      stress0[0] = stressTensor[0, 0];
      stress0[1] = stressTensor[1, 1];
      stress0[2] = stressTensor[2, 2];
      stress0[3] = stressTensor[0, 1];
      if (((PatchStructure2D)this.patch).StateStress == Structure2DState.PlaneStrain)
      {
        double nu = lamda / (2.0 * (lamda + mu));
        stress0[2] = nu * (stressTensor[0, 0] + stressTensor[1, 1]);
      }
      else if (((PatchStructure2D)this.patch).StateStress == Structure2DState.PlaneStress)
      {
        stress0[2] = 0;
      }
    }
    private DoubleMatrix[] ComputeAnisotropicMa(DoubleVector[] eigenVectors)
    {

      DoubleMatrix[] Ma = new DoubleMatrix[3];
      Ma[0] = MatrixFunctions.OuterProduct(eigenVectors[0], eigenVectors[0]);
      Ma[1] = MatrixFunctions.OuterProduct(eigenVectors[1], eigenVectors[1]);
      Ma[2] = MatrixFunctions.OuterProduct(eigenVectors[2], eigenVectors[2]);
      return Ma;
    }
    private double ComputeAnisotropicQa(DoubleMatrix Ma, int i, int j, int k, int l)
    {
      return Ma[i - 1, j - 1] * Ma[k - 1, l - 1];
    }
    private DoubleMatrix ComputeAnisotropicQaMatrix(DoubleMatrix Ma)
    {
      DoubleMatrix Qa = new DoubleMatrix(6, 6);
      Qa[0, 0] = ComputeAnisotropicQa(Ma, 1, 1, 1, 1);
      Qa[0, 1] = Qa[1, 0] = ComputeAnisotropicQa(Ma, 1, 1, 2, 2);
      Qa[0, 2] = Qa[2, 0] = ComputeAnisotropicQa(Ma, 1, 1, 3, 3);
      Qa[0, 3] = Qa[3, 0] = ComputeAnisotropicQa(Ma, 1, 1, 2, 3);
      Qa[0, 4] = Qa[4, 0] = ComputeAnisotropicQa(Ma, 1, 1, 1, 3);
      Qa[0, 5] = Qa[5, 0] = ComputeAnisotropicQa(Ma, 1, 1, 1, 2);
      Qa[1, 1] = ComputeAnisotropicQa(Ma, 2, 2, 2, 2);
      Qa[1, 2] = Qa[2, 1] = ComputeAnisotropicQa(Ma, 2, 2, 3, 3);
      Qa[1, 3] = Qa[3, 1] = ComputeAnisotropicQa(Ma, 2, 2, 2, 3);
      Qa[1, 4] = Qa[4, 1] = ComputeAnisotropicQa(Ma, 2, 2, 1, 3);
      Qa[1, 5] = Qa[5, 1] = ComputeAnisotropicQa(Ma, 2, 2, 1, 2);
      Qa[2, 2] = ComputeAnisotropicQa(Ma, 3, 3, 3, 3);
      Qa[2, 3] = Qa[3, 2] = ComputeAnisotropicQa(Ma, 3, 3, 2, 3);
      Qa[2, 4] = Qa[4, 2] = ComputeAnisotropicQa(Ma, 3, 3, 1, 3);
      Qa[2, 5] = Qa[5, 2] = ComputeAnisotropicQa(Ma, 3, 3, 1, 2);
      Qa[3, 3] = ComputeAnisotropicQa(Ma, 2, 3, 3, 2);
      Qa[3, 4] = Qa[4, 3] = ComputeAnisotropicQa(Ma, 2, 3, 1, 3);
      Qa[3, 5] = Qa[5, 3] = ComputeAnisotropicQa(Ma, 2, 3, 1, 2);
      Qa[4, 4] = ComputeAnisotropicQa(Ma, 1, 3, 3, 1);
      Qa[4, 5] = Qa[5, 4] = ComputeAnisotropicQa(Ma, 1, 3, 1, 2);
      Qa[5, 5] = ComputeAnisotropicQa(Ma, 1, 2, 2, 1);
      return Qa;
    }
    private double ComputeAnisotropicGab(DoubleMatrix Ma, DoubleMatrix Mb, int i, int j, int k, int l)
    {
      return Ma[i - 1, k - 1] * Mb[j - 1, l - 1] + Ma[i - 1, l - 1] * Mb[j - 1, k - 1];
    }
    private DoubleMatrix ComputeAnisotropicGabMatrix(DoubleMatrix Ma, DoubleMatrix Mb)
    {
      DoubleMatrix Gab = new DoubleMatrix(6, 6);
      Gab[0, 0] = ComputeAnisotropicGab(Ma, Mb, 1, 1, 1, 1);
      Gab[0, 1] = Gab[1, 0] = ComputeAnisotropicGab(Ma, Mb, 1, 1, 2, 2);
      Gab[0, 2] = Gab[2, 0] = ComputeAnisotropicGab(Ma, Mb, 1, 1, 3, 3);
      Gab[0, 3] = Gab[3, 0] = ComputeAnisotropicGab(Ma, Mb, 1, 1, 2, 3);
      Gab[0, 4] = Gab[4, 0] = ComputeAnisotropicGab(Ma, Mb, 1, 1, 1, 3);
      Gab[0, 5] = Gab[5, 0] = ComputeAnisotropicGab(Ma, Mb, 1, 1, 1, 2);
      Gab[1, 1] = ComputeAnisotropicGab(Ma, Mb, 2, 2, 2, 2);
      Gab[1, 2] = Gab[2, 1] = ComputeAnisotropicGab(Ma, Mb, 2, 2, 3, 3);
      Gab[1, 3] = Gab[3, 1] = ComputeAnisotropicGab(Ma, Mb, 2, 2, 2, 3);
      Gab[1, 4] = Gab[4, 1] = ComputeAnisotropicGab(Ma, Mb, 2, 2, 1, 3);
      Gab[1, 5] = Gab[5, 1] = ComputeAnisotropicGab(Ma, Mb, 2, 2, 1, 2);
      Gab[2, 2] = ComputeAnisotropicGab(Ma, Mb, 3, 3, 3, 3);
      Gab[2, 3] = Gab[3, 2] = ComputeAnisotropicGab(Ma, Mb, 3, 3, 2, 3);
      Gab[2, 4] = Gab[4, 2] = ComputeAnisotropicGab(Ma, Mb, 3, 3, 1, 3);
      Gab[2, 5] = Gab[5, 2] = ComputeAnisotropicGab(Ma, Mb, 3, 3, 1, 2);
      Gab[3, 3] = ComputeAnisotropicGab(Ma, Mb, 2, 3, 3, 2);
      Gab[3, 4] = Gab[4, 3] = ComputeAnisotropicGab(Ma, Mb, 2, 3, 1, 3);
      Gab[3, 5] = Gab[5, 3] = ComputeAnisotropicGab(Ma, Mb, 2, 3, 1, 2);
      Gab[4, 4] = ComputeAnisotropicGab(Ma, Mb, 1, 3, 3, 1);
      Gab[4, 5] = Gab[5, 4] = ComputeAnisotropicGab(Ma, Mb, 1, 3, 1, 2);
      Gab[5, 5] = ComputeAnisotropicGab(Ma, Mb, 1, 2, 2, 1);
      return Gab;
    }
    private double HeavisideFunction(double a)
    {
      if (a > 0)
        return 1;
      else if (a == 0)
        return 0.5;
      else
        return 0;
      //if (a > 0)
      //  return 1;
      //else
      //  return 0;
    }
    private void ComputeAnisotropicProjectionTensor(DoubleVector principleStrain, DoubleVector[] eigenVectors, out DoubleMatrix PPositive, out DoubleMatrix PNegative)
    {
      DoubleMatrix[] Ma = ComputeAnisotropicMa(eigenVectors);
      PPositive = new DoubleMatrix(6, 6);
      PNegative = new DoubleMatrix(6, 6);
      //DoubleVector prinTemp = new DoubleVector(3);
      //prinTemp[0] = principleStrain[0];
      //prinTemp[1] = principleStrain[1];
      //prinTemp[2] = principleStrain[2];
      //double delta = 1e-9;
      //if (principleStrain[0] == principleStrain[1])
      //{
      //  prinTemp[0] = prinTemp[0] * (1 + delta);
      //}
      //if (principleStrain[1] == principleStrain[2])
      //{
      //  prinTemp[2] = prinTemp[2] * (1 - delta);
      //}
      //prinTemp[0] = principleStrain[0] * (1 + delta);
      //prinTemp[1] = principleStrain[1] * (1 - delta);
      //prinTemp[2] = principleStrain[2] / ((1 + delta) * (1 - delta));
      for (int a = 0; a < 3; a++)
      {
        DoubleMatrix secondTermPos = new DoubleMatrix(6, 6);
        DoubleMatrix secondTermNeg = new DoubleMatrix(6, 6);
        for (int b = 0; b < 3; b++)
        {
          if (b != a)
          {
            DoubleMatrix temp = 1.0 / 4.0 * (ComputeAnisotropicGabMatrix(Ma[a], Ma[b]) + ComputeAnisotropicGabMatrix(Ma[b], Ma[a]));//current
            double dab = principleStrain[a] - principleStrain[b];
            secondTermPos += HeavisideFunction(dab) * temp;
            secondTermNeg += HeavisideFunction(-dab) * temp;
            /////////////////////////////////////////////////////////////////
            //DoubleMatrix temp = 1.0 / 2.0 * (ComputeAnisotropicGabMatrix(Ma[a], Ma[b]) + ComputeAnisotropicGabMatrix(Ma[b], Ma[a]));//current
            //double dab = principleStrain[a] - principleStrain[b];
            //secondTermPos += HeavisideFunction(dab) * temp;
            //secondTermNeg += HeavisideFunction(-dab) * temp;
            ////////////////////////////////////////////////////////////////////////
            //DoubleMatrix temp = 1.0 / 2.0 * (ComputeAnisotropicGabMatrix(Ma[a], Ma[b]) + ComputeAnisotropicGabMatrix(Ma[b], Ma[a]));//current
            //double dab = principleStrain[a] - principleStrain[b];
            //if (Math.Abs(dab) <= 1e-4)
            //  dab = Math.Abs(dab);
            //secondTermPos += HeavisideFunction(dab) * temp;
            //secondTermNeg += HeavisideFunction(-dab) * temp;
            /////////////////////////////////////////////////////////////////
            //DoubleMatrix temp = 1.0 / 2.0 * (ComputeAnisotropicGabMatrix(Ma[a], Ma[b]) + ComputeAnisotropicGabMatrix(Ma[b], Ma[a]));
            //secondTermPos += (PositiveRampFunction(prinTemp[a]) - PositiveRampFunction(prinTemp[b])) / (prinTemp[a] - prinTemp[b]) * temp;
            //secondTermNeg += (NegativeRampFunction(prinTemp[a]) - NegativeRampFunction(prinTemp[b])) / (prinTemp[a] - prinTemp[b]) * temp;
          }
        }

        DoubleMatrix Qa = ComputeAnisotropicQaMatrix(Ma[a]);
        PPositive += HeavisideFunction(principleStrain[a]) * Qa + secondTermPos;
        PNegative += HeavisideFunction(-principleStrain[a]) * Qa + secondTermNeg;
      }
    }
    private int[] IndexPricipleStrain(DoubleVector principleStrain)
    {
      double maxValue = principleStrain[0];
      int maxIndex = 0;
      for (int i = 1; i < 3; i++)
      {
        if (maxValue < principleStrain[i])
        {
          maxIndex = i;
        }
      }
      int minIndex = -1;
      int midIndex = -1;
      switch (maxIndex)
      {
        case 0:
          if (principleStrain[1] < principleStrain[2])
          {
            midIndex = 2;
            minIndex = 1;
          }
          else
          {
            midIndex = 1;
            minIndex = 2;
          }
          break;
        case 1:
          if (principleStrain[0] < principleStrain[2])
          {
            midIndex = 2;
            minIndex = 0;
          }
          else
          {
            midIndex = 0;
            minIndex = 2;
          }
          break;
        case 2:
          if (principleStrain[0] < principleStrain[1])
          {
            midIndex = 1;
            minIndex = 0;
          }
          else
          {
            midIndex = 0;
            minIndex = 1;
          }
          break;
          //default:
          //  maxIndex = 0;
          //  midIndex = 1;
          //  minIndex = 2;
          //  break;
      }
      return new int[] { maxIndex, midIndex, minIndex };
    }
    private DoubleMatrix CreateMaterialMatrixIsotropic(double phi, double k, double l0, double lch, double lamda, double mu)
    {
      AnisotropicElasticity prob = (AnisotropicElasticity)Material.GetProperty(MaterialPropertyName.AnisotropicElasticity);
      DoubleMatrix C = null;
      if (prob != null)
      {
        double[,] ce = prob.GetAnisotropicElasticityMatrix();
        C = new DoubleMatrix(6, 6);
        bool isInverse = prob.GetIsInverseMatrix();
        for (int i = 0; i < 6; i++)
          for (int j = 0; j < 6; j++)
          {
            C[i, j] = ce[i, j];
          }
        if (isInverse)
          C = NMathFunctions.Inverse(C);
      }
      else
      {
        OrthotropicElasticity prob1 = (OrthotropicElasticity)Material.GetProperty(MaterialPropertyName.OrthotropicElasticity);
        if (prob1 == null)
        {
          C = new DoubleMatrix(6, 6);
          switch (((PatchStructure2D)patch).StateStress)
          {
            case Structure2DState.PlaneStrain:
            case Structure2DState.Axisymetric:
              //Xem lai thu tu xx yy zz yz xz xy
              C[0, 0] = coefficient(1, 1, 1, 1, lamda, mu);
              C[0, 1] = C[1, 0] = coefficient(1, 1, 2, 2, lamda, mu);
              C[0, 2] = C[2, 0] = coefficient(1, 1, 3, 3, lamda, mu);
              C[0, 3] = C[3, 0] = coefficient(1, 1, 2, 3, lamda, mu);
              C[0, 4] = C[4, 0] = coefficient(1, 1, 1, 3, lamda, mu);
              C[0, 5] = C[5, 0] = coefficient(1, 1, 1, 2, lamda, mu);
              C[1, 1] = coefficient(2, 2, 2, 2, lamda, mu);
              C[1, 2] = C[2, 1] = coefficient(2, 2, 3, 3, lamda, mu);
              C[1, 3] = C[3, 1] = coefficient(2, 2, 2, 3, lamda, mu);
              C[1, 4] = C[4, 1] = coefficient(2, 2, 1, 3, lamda, mu);
              C[1, 5] = C[5, 1] = coefficient(2, 2, 1, 2, lamda, mu);
              C[2, 2] = coefficient(3, 3, 3, 3, lamda, mu);
              C[2, 3] = C[3, 2] = coefficient(3, 3, 2, 3, lamda, mu);
              C[2, 4] = C[4, 2] = coefficient(3, 3, 1, 3, lamda, mu);
              C[2, 5] = C[5, 2] = coefficient(3, 3, 1, 2, lamda, mu);
              C[3, 3] = coefficient(2, 3, 3, 2, lamda, mu);
              C[3, 4] = C[4, 3] = coefficient(2, 3, 1, 3, lamda, mu);
              C[3, 5] = C[5, 3] = coefficient(2, 3, 1, 2, lamda, mu);
              C[4, 4] = coefficient(1, 3, 3, 1, lamda, mu);
              C[4, 5] = C[5, 4] = coefficient(1, 3, 1, 2, lamda, mu);
              C[5, 5] = coefficient(1, 2, 2, 1, lamda, mu);
              break;
            case Structure2DState.PlaneStress:
              double eModulus = mu * (3.0 * lamda + 2.0 * mu) / (lamda + mu);
              double nu = lamda / (2 * (lamda + mu));
              double a = eModulus / (1.0 - nu * nu);
              C[0, 0] = a;
              C[0, 1] = a * nu;
              C[1, 0] = a * nu;
              C[1, 1] = a;
              C[5, 5] = a * (1 - nu) / 2;
              break;
          }
          C = MatrixTool.GetSubMatrix(C, new int[] { 0, 1, 2, 5 }, new int[] { 0, 1, 2, 5 }, AbstractModel.NumberOfCPUs);
        }
        else
        {
          C = new DoubleMatrix(prob1.GetOrthotropicElasticityMatrix());
        }
      }
      return (OmegaD(phi, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel) + k) * C;
    }
    private double coefficient(int i, int j, int k, int l, double lamda, double mu)
    {
      //https://en.wikipedia.org/wiki/Linear_elasticity
      //https://en.wikipedia.org/wiki/Hooke%27s_law#Isotropic_materials
      return lamda * KroneckerDelta(i, j) * KroneckerDelta(k, l) + mu * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k));
      //return K * KroneckerDelta(i, j) * KroneckerDelta(k, l) + mu * (KroneckerDelta(i, k) * KroneckerDelta(j, l) + KroneckerDelta(i, l) * KroneckerDelta(j, k) - 2.0 / 3.0 * KroneckerDelta(i, j) * KroneckerDelta(k, l));
    }
    private DoubleMatrix CreateMaterialMatrixAnisotropic(DoubleVector principleStrain, DoubleVector[] eigenVectors, double phi, double k, double l0, double lch, double lamda, double mu)
    {
      //if (principleStrain[0] != 0 || principleStrain[1] != 0 || principleStrain[2] != 0)
      //{
      double trStrain = principleStrain[0] + principleStrain[1] + principleStrain[2];
      DoubleMatrix PPositive = new DoubleMatrix(6, 6);
      DoubleMatrix PNegative = new DoubleMatrix(6, 6);
      ComputeAnisotropicProjectionTensor(principleStrain, eigenVectors, out PPositive, out PNegative);
      DoubleMatrix J = JTensor();
      DoubleMatrix C = (OmegaD(phi, l0, lch, ModelStructurePhaseFieldStatic.typePhaseFieldModel) + k) * (lamda * HeavisideFunction(trStrain) * J + 2 * mu * PPositive)
        + (lamda * HeavisideFunction(-trStrain) * J + 2 * mu * PNegative);
      return MatrixTool.GetSubMatrix(C, new int[] { 0, 1, 2, 5 }, new int[] { 0, 1, 2, 5 }, AbstractModel.NumberOfCPUs);
      //}
      //else
      //{
      //  return (Math.Pow(1 - phi, 2) + k) * CreateMaterialMatrix();
      //}
    }
    private void ComputePrincipleStrains(DoubleVector strain, out DoubleVector eigenValues, out DoubleVector[] eigenVectors)
    {
      DoubleMatrix strainTensor = new DoubleMatrix(3, 3);
      strainTensor[0, 0] = strain[0];
      strainTensor[1, 1] = strain[1];
      strainTensor[2, 2] = strain[2];
      strainTensor[0, 1] = strainTensor[1, 0] = strain[3] / 2.0;
      var eigDecomp = new DoubleSymEigDecomp(strainTensor);
      eigenValues = eigDecomp.EigenValues;
      eigenVectors = new DoubleVector[3];
      //double delta = 1e-9;
      //eigenVectors[0] = eigDecomp.EigenVector(0) * (1 + delta);
      //eigenVectors[1] = eigDecomp.EigenVector(1) * (1 - delta);
      //eigenVectors[2] = eigDecomp.EigenVector(2) / ((1 + delta) * (1 - delta));
      eigenVectors[0] = eigDecomp.EigenVector(0);
      eigenVectors[1] = eigDecomp.EigenVector(1);
      eigenVectors[2] = eigDecomp.EigenVector(2);
    }
    private void ComputePrincipleStress(DoubleVector stress, out DoubleVector eigenValues, out DoubleVector[] eigenVectors)
    {
      DoubleMatrix stressTensor = new DoubleMatrix(3, 3);
      stressTensor[0, 0] = stress[0];
      stressTensor[1, 1] = stress[1];
      stressTensor[2, 2] = stress[2];
      stressTensor[0, 1] = stressTensor[1, 0] = stress[3];
      var eigDecomp = new DoubleSymEigDecomp(stressTensor);
      eigenValues = eigDecomp.EigenValues;
      eigenVectors = new DoubleVector[3];
      //double delta = 1e-9;
      //eigenVectors[0] = eigDecomp.EigenVector(0) * (1 + delta);
      //eigenVectors[1] = eigDecomp.EigenVector(1) * (1 - delta);
      //eigenVectors[2] = eigDecomp.EigenVector(2) / ((1 + delta) * (1 - delta));
      eigenVectors[0] = eigDecomp.EigenVector(0);
      eigenVectors[1] = eigDecomp.EigenVector(1);
      eigenVectors[2] = eigDecomp.EigenVector(2);
    }

    private void ComputeStressAnisotropic(DoubleVector principleStrain, DoubleVector[] eigenVectors, double lamda, double mu, out DoubleVector stress0Pos, out DoubleVector stress0Neg)
    {
      var prinVector1 = eigenVectors[0];
      var prinVector2 = eigenVectors[1];
      var prinVector3 = eigenVectors[2];

      double trEpsilon = principleStrain[0] + principleStrain[1] + principleStrain[2];
      DoubleMatrix unit = new DoubleMatrix(3, 3);
      unit[0, 0] = unit[1, 1] = unit[2, 2] = 1.0;
      //DoubleMatrix stressTensorPos = (lamda * PositiveRampFunction(trEpsilon) + 2 * mu * PositiveRampFunction(principleStrain[0])) * NMathFunctions.OuterProduct(prinVector1, prinVector1)
      //+ (lamda * PositiveRampFunction(trEpsilon) + 2 * mu * PositiveRampFunction(principleStrain[1])) * NMathFunctions.OuterProduct(prinVector2, prinVector2)
      //+ (lamda * PositiveRampFunction(trEpsilon) + 2 * mu * PositiveRampFunction(principleStrain[2])) * NMathFunctions.OuterProduct(prinVector3, prinVector3);
      //DoubleMatrix stressTensorNeg = (lamda * NegativeRampFunction(trEpsilon) + 2 * mu * NegativeRampFunction(principleStrain[0])) * NMathFunctions.OuterProduct(prinVector1, prinVector1)
      //+ (lamda * NegativeRampFunction(trEpsilon) + 2 * mu * NegativeRampFunction(principleStrain[1])) * NMathFunctions.OuterProduct(prinVector2, prinVector2)
      //+ (lamda * NegativeRampFunction(trEpsilon) + 2 * mu * NegativeRampFunction(principleStrain[2])) * NMathFunctions.OuterProduct(prinVector3, prinVector3);

      DoubleMatrix stressTensorPos = lamda * PositiveRampFunction(trEpsilon) * unit +
      2 * mu * (PositiveRampFunction(principleStrain[0]) * NMathFunctions.OuterProduct(prinVector1, prinVector1)
      + PositiveRampFunction(principleStrain[1]) * NMathFunctions.OuterProduct(prinVector2, prinVector2)
      + PositiveRampFunction(principleStrain[2]) * NMathFunctions.OuterProduct(prinVector3, prinVector3));
      DoubleMatrix stressTensorNeg = lamda * NegativeRampFunction(trEpsilon) * unit + 2 * mu * (NegativeRampFunction(principleStrain[0]) * NMathFunctions.OuterProduct(prinVector1, prinVector1)
        + NegativeRampFunction(principleStrain[1]) * NMathFunctions.OuterProduct(prinVector2, prinVector2)
        + NegativeRampFunction(principleStrain[2]) * NMathFunctions.OuterProduct(prinVector3, prinVector3));

      double nu = lamda / (2.0 * (lamda + mu));

      stress0Pos = new DoubleVector(4);
      stress0Neg = new DoubleVector(4);
      stress0Pos[0] = stressTensorPos[0, 0];
      stress0Pos[1] = stressTensorPos[1, 1];
      stress0Pos[2] = stressTensorPos[2, 2];
      stress0Pos[3] = stressTensorPos[0, 1];

      stress0Neg[0] = stressTensorNeg[0, 0];
      stress0Neg[1] = stressTensorNeg[1, 1];
      stress0Neg[2] = stressTensorNeg[2, 2];
      stress0Neg[3] = stressTensorNeg[0, 1];
      if (((PatchStructure2D)this.patch).StateStress == Structure2DState.PlaneStrain)
      {
        stress0Pos[2] = nu * (stressTensorPos[0, 0] + stressTensorPos[1, 1]);
        stress0Neg[2] = nu * (stressTensorNeg[0, 0] + stressTensorNeg[1, 1]);
      }
      else if (((PatchStructure2D)this.patch).StateStress == Structure2DState.PlaneStress)
      {
        stress0Pos[2] = 0;
        stress0Neg[2] = 0;
      }
    }

    public override DoubleVector StressAt(params double[] xi)
    {
      double lamda = ComputeParameterProperty(MaterialPropertyName.LameParameter, xi[0], xi[1]);
      double mu = ComputeParameterProperty(MaterialPropertyName.ShearModulus, xi[0], xi[1]);
      double gc = ComputeParameterProperty(MaterialPropertyName.CriticalEnergyReleaseRate, xi[0], xi[1]);
      double ft = ComputeParameterProperty(MaterialPropertyName.TensileYieldStrength, xi[0], xi[1]);
      double E = mu * (3 * lamda + 2 * mu) / (lamda + mu);
      double lch = E * gc / (ft * ft);
      DoubleVector strain = StrainAt(xi);
      DoubleVector N = ValueBasisFunction(xi[0], xi[1]);
      double phi = NMathFunctions.Dot(N, GetEachVariableLocal(Result.PHASEFIELD));
      DoubleVector principleStrain;
      DoubleVector[] eigenVectors;
      //strain[3] *= 2.0;
      ComputePrincipleStrains(strain, out principleStrain, out eigenVectors);
      return ComputeStress(lamda, mu, lch, phi, principleStrain, eigenVectors);
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
