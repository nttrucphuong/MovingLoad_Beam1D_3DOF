using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
  public class ElementStructurePlastic2D : AbstractElementStructure2D
  {
    public ElementStructurePlastic2D(AbstractPatch2D mesh, int id)
        : base(mesh, id)
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
      DoubleMatrix D = CreateMaterialMatrix();
      DoubleMatrix Idev = new DoubleMatrix(4, 4);
      Idev[0, 0] = 2.0 / 3.0;
      Idev[0, 1] = -1.0 / 3.0;
      Idev[0, 2] = -1.0 / 3.0;
      Idev[1, 0] = -1.0 / 3.0;
      Idev[1, 1] = 2.0 / 3.0;
      Idev[1, 2] = -1.0 / 3.0;
      Idev[2, 0] = -1.0 / 3.0;
      Idev[2, 1] = -1.0 / 3.0;
      Idev[2, 2] = 2.0 / 3.0;
      Idev[3, 3] = 1.0;
      DoubleVector vone = new DoubleVector(4);
      vone[0] = 1.0;
      vone[1] = 1.0;
      vone[2] = 1.0;
      double K = Material.GetProperty(MaterialPropertyName.BulkModulus).GetValueProperty();
      double G = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
      DoubleMatrix Ce = K * MatrixFunctions.OuterProduct(vone, vone) + 2.0 * G * Idev;


      ////doc file ////


      //outputCoffive(table1);

      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {

          double detJbar = 1.0 / 4.0
                  * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                  * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          //DoubleMatrix dNdxi = gps[i, j].dNdxi;
          //DoubleMatrix J = JacobianAt(dNdxi);
          DoubleMatrix dNdx = (DoubleMatrix)gps[i, j].GetValue(DataInGausspoint.dNdX);//dNdxi * invertJ;
          for (int k = 0; k < nen; k++)
          {
            B[0, 2 * k] = dNdx[k, 0];
            B[3, 2 * k] = dNdx[k, 1];
            B[1, 2 * k + 1] = dNdx[k, 1];
            B[3, 2 * k + 1] = dNdx[k, 0];
          }
          // ============================================
          // Update stress using Return-Mapping algorithm
          // ============================================
          // History (converged) values of last load step
          //DoubleVector stress_n = gps[i, j].lastStress;
          DoubleVector eps_pl_n = (DoubleVector)gps[i, j].GetValue(DataInGausspoint.lastPlasticStrain);//lastPlasticStrain;
          DoubleVector backstress_n = (DoubleVector)gps[i, j].GetValue(DataInGausspoint.lastBackStress);//gps[i, j].lastBackStress;
          double alpha_n = (double)gps[i, j].GetValue(DataInGausspoint.lastAlpha);//gps[i, j].lastAlpha;

          DoubleVector ue = GetDisplacementLocal();////////////
          DoubleVector strain = MatrixFunctions.Product(B, ue); // [eps_xx; epx_yy; epx_zz; 2*epx_xy]
          DoubleVector eps_total = strain.DeepenThisCopy();
          eps_total[3] = eps_total[3] / 2.0; // [eps_xx; epx_yy; epx_zz; epx_xy]
                                             // Trial state
          DoubleVector stress_tr = MatrixFunctions.Product(Ce, (eps_total - eps_pl_n));
          // DoubleVector e_n = Idev*eps_total;/////////////////////
          DoubleVector s_tr = MatrixFunctions.Product(Idev, stress_tr);  // deviatoric stress  
                                                                         // DoubleVector s_tr = 2 * G * (e_n - eps_pl_n);/////////////////////
          DoubleVector beta_tr = s_tr - backstress_n; // relative stress
          double norm_beta = Math.Sqrt(beta_tr[0] * beta_tr[0] + beta_tr[1] * beta_tr[1] + beta_tr[2] * beta_tr[2] + 2.0 * beta_tr[3] * beta_tr[3]);
          // Yield condition
          DoubleVector n_vector = beta_tr / norm_beta;
          if (Material.GetProperty(MaterialPropertyName.BilinearIsotropicHardening) != null || Material.GetProperty(MaterialPropertyName.BilinearKinematicHardening) != null || (Material.GetProperty(MaterialPropertyName.BilinearCombinedIsotropicKinematicHardening) != null))
          {
            double sigmaY = Material.GetProperty(MaterialPropertyName.YieldStrength).GetValueProperty(); // initial yield stress
            double H0 = 0;
            double K0 = 0;
            if (Material.GetProperty(MaterialPropertyName.BilinearIsotropicHardening) != null)
            {
              double E = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
              double Et = Material.GetProperty(MaterialPropertyName.TangentModulus).GetValueProperty();
              H0 = E * Et / (E - Et);// for isotropic hardening
            }
            else if (Material.GetProperty(MaterialPropertyName.BilinearKinematicHardening) != null)
            {
              double E = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
              double Et = Material.GetProperty(MaterialPropertyName.TangentModulus).GetValueProperty();
              K0 = E * Et / (E - Et);// for isotropic hardening
            }
            double ftol = sigmaY * 1e-6;//tolerance for yield

            double f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * (sigmaY + H0 * alpha_n);//trial yield function
            if (f_yield_tr < ftol)
            // elastic
            {
              gps[i, j].SetValue(DataInGausspoint.currentStress, stress_tr);
              gps[i, j].SetValue(DataInGausspoint.currentAlpha, alpha_n);
              gps[i, j].SetValue(DataInGausspoint.currentPlasticStrain, eps_pl_n);
              gps[i, j].SetValue(DataInGausspoint.currentBackStress, backstress_n);
              gps[i, j].SetValue(DataInGausspoint.currentDep, D);
              gps[i, j].SetValue(DataInGausspoint.currentStrain, eps_total);

              //gps[i, j].currentStress = stress_tr;
              //gps[i, j].currentAlpha = alpha_n;
              //gps[i, j].currentPlasticStrain = eps_pl_n;
              //gps[i, j].currentBackStress = backstress_n;
              //gps[i, j].currentDep = D;
              //gps[i, j].currentStrain = eps_total;
            }
            else
            {
              double dgamma = f_yield_tr / (2.0 * G + 2.0 * (H0 + K0) / 3.0);

              // Update consistent Elastoplastic material property tensor
              DoubleMatrix Cep_abstr1 = 2.0 * G * (1.0 / (1.0 + (H0 + K0) / (3.0 * G)) - 2.0 * G * dgamma / norm_beta) * MatrixFunctions.OuterProduct(n_vector, n_vector);

              DoubleMatrix Id = Idev.DeepenThisCopy();
              Id[3, 3] = Id[3, 3] / 2.0;
              DoubleMatrix Cep_abstr2 = (4.0 * G * G * dgamma / norm_beta) * Id;
              DoubleMatrix Cep = D - Cep_abstr1 - Cep_abstr2;
              //////////////////////////////////////////////////////////////////////

              // Store
              gps[i, j].SetValue(DataInGausspoint.currentStress, stress_tr - 2.0 * G * dgamma * n_vector);
              gps[i, j].SetValue(DataInGausspoint.currentAlpha, alpha_n + Math.Sqrt(2.0 / 3.0) * dgamma);
              gps[i, j].SetValue(DataInGausspoint.currentPlasticStrain, eps_pl_n + dgamma * n_vector);
              gps[i, j].SetValue(DataInGausspoint.currentBackStress, backstress_n + 2.0 / 3.0 * K0 * dgamma * n_vector);
              gps[i, j].SetValue(DataInGausspoint.currentDep, Cep);
              gps[i, j].SetValue(DataInGausspoint.currentStrain, eps_total);
              //gps[i, j].currentBackStress = backstress_n + 2.0 / 3.0 * K0 * dgamma * n_vector;
              //gps[i, j].currentAlpha = alpha_n + Math.Sqrt(2.0 / 3.0) * dgamma;
              //gps[i, j].currentPlasticStrain = eps_pl_n + dgamma * n_vector;
              //gps[i, j].currentStress = stress_tr - 2.0 * G * dgamma * n_vector;
              //gps[i, j].currentDep = Cep;
              //gps[i, j].currentStrain = eps_total;
            }
            // Assembly
            DoubleMatrix Cij = (DoubleMatrix)gps[i, j].GetValue(DataInGausspoint.currentDep); // material matrix

            DoubleMatrix BTDB = NMathFunctions.TransposeProduct(B, MatrixFunctions.Product(Cij, B));
            double detJ = Math.Abs((double)gps[i, j].GetValue(DataInGausspoint.detJ));
            Ke += gps[i, j].weight * detJbar * detJ * BTDB;
          }



          ///////////////////////////////////////////////////
          //int nm = 3;
          //////////////////////////////////////////////////
          //if (nm == 0) /// multilinear//////
          //{
          //    double f_yield_tr;
          //    if (alpha_n == 0)
          //    {
          //        f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * table1[1, 1];
          //    }
          //    else
          //    {
          //        f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * (LinearInterpolation(alpha_n, Transformtable(table1), 'x'));
          //    }    // vonMises yield function
          //    if (f_yield_tr < 0)
          //    // elastic
          //    {
          //        gps[i, j].currentStress = stress_tr;
          //        gps[i, j].currentBackStress = backstress_n;
          //        gps[i, j].currentAlpha = alpha_n;
          //        gps[i, j].currentPlasticStrain = eps_pl_n;

          //        gps[i, j].Dep = D;
          //    }
          //    else
          //    {
          //        double dgamma = 0;
          //        double gdgama = -Math.Sqrt(2.0 / 3.0) * LinearInterpolation(alpha_n, Transformtable(table1), 'x') + norm_beta;
          //        double dgdgama;
          //        double alphak;
          //        alphak = alpha_n;
          //        double TOL = 10e-5;
          //        int maxIter = 100;
          //        int nIter = 0;

          //        //  double m = 5;
          //        do
          //        {

          //            nIter++;
          //            gdgama = -Math.Sqrt(2.0 / 3.0) * LinearInterpolation(alpha_n, Transformtable(table1), 'x') + norm_beta - 2.0 * G * dgamma;
          //            dgdgama = -2.0 * G * (1.0 + DplasticFun(alphak, Transformtable(table1)) / (3.0 * G));
          //            dgamma -= (gdgama / dgdgama);
          //            alphak = alpha_n + Math.Sqrt(2.0 / 3.0) * dgamma;
          //        }

          //        while ((Math.Abs(gdgama) > TOL) && (nIter < maxIter));


          //        DoubleVector stress = stress_tr - 2.0 * G * dgamma * n_vector;



          //        // Update consistent Elastoplastic material property tensor
          //        DoubleMatrix Cep_abstr1 = 2.0 * G * (1.0 / (1.0 + (DplasticFun(alpha_n, Transformtable(table1))) / (3.0 * G)) - 2.0 * G * dgamma / norm_beta) * ((DoubleMatrix)n_vector.OuterProduct(n_vector));


          //        DoubleMatrix Id = (DoubleMatrix)Idev.Clone();
          //        Id[3, 3] = Id[3, 3] / 2.0;
          //        DoubleMatrix Cep_abstr2 = (4.0 * G * G * dgamma / norm_beta) * Id;
          //        DoubleMatrix Cep = D - Cep_abstr1 - Cep_abstr2;

          //        double alpha = alpha_n + Math.Sqrt(2.0 / 3.0) * dgamma;
          //        DoubleVector eps_pl = eps_pl_n + dgamma * n_vector;

          //        // Store
          //        gps[i, j].currentAlpha = alpha;
          //        gps[i, j].currentPlasticStrain = eps_pl;
          //        gps[i, j].currentStress = stress;
          //        gps[i, j].currentBackStress = backstress_n;
          //        gps[i, j].Dep = Cep;
          //        //gps[i, j].Dep = D;
          //    }

          //    // Assembly
          //    DoubleMatrix Cij = gps[i, j].Dep; // material matrix
          //    //DoubleVector stressS = gps[i, j].lastStress; // stress

          //    DoubleMatrix BTDB = (DoubleMatrix)B.Transpose() * Cij * B;
          //    double detJ = Math.Abs(J.Determinant());
          //    Ke += gps[i, j].weight * detJbar * detJ * BTDB;
          //}

          //if (nm == 1)
          //{
          //    //double f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * (sigmaY + H0 * alpha_n); // vonMises yield function
          //    //H0 = ModunPlastic(table2, 0.001);
          //    double f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * (sigmaY + H0 * alpha_n);
          //    if (f_yield_tr < 0)
          //    // elastic
          //    {
          //        gps[i, j].currentStress = stress_tr;
          //        gps[i, j].currentAlpha = alpha_n;
          //        gps[i, j].currentPlasticStrain = eps_pl_n;
          //        gps[i, j].currentBackStress = backstress_n;
          //        gps[i, j].Dep = D;
          //    }
          //    else
          //    {
          //        double dgamma = f_yield_tr / (2.0 * G + 2.0 * K0 / 3.0 + 2.0 * H0 / 3.0);
          //        DoubleVector stress = stress_tr - 2.0 * G * dgamma * n_vector;



          //        // Update consistent Elastoplastic material property tensor
          //        DoubleMatrix Cep_abstr1 = 2.0 * G * (1.0 / (1.0 + (H0 + K0) / (3.0 * G)) - 2.0 * G * dgamma / norm_beta) * ((DoubleMatrix)n_vector.OuterProduct(n_vector));


          //        DoubleMatrix Id = (DoubleMatrix)Idev.Clone();
          //        Id[3, 3] = Id[3, 3] / 2.0;
          //        DoubleMatrix Cep_abstr2 = (4.0 * G * G * dgamma / norm_beta) * Id;
          //        DoubleMatrix Cep = D - Cep_abstr1 - Cep_abstr2;

          //        double alpha = alpha_n + Math.Sqrt(2.0 / 3.0) * dgamma;
          //        DoubleVector eps_pl = eps_pl_n + dgamma * n_vector;
          //        DoubleVector qCurrent = backstress_n + Math.Sqrt(2.0 / 3.0) * (K0 * dgamma) * n_vector;


          //        // Store
          //        gps[i, j].currentBackStress = qCurrent;
          //        gps[i, j].currentAlpha = alpha;
          //        gps[i, j].currentPlasticStrain = eps_pl;
          //        gps[i, j].currentStress = stress;
          //        gps[i, j].Dep = Cep;
          //    }

          //    // Assembly
          //    DoubleMatrix Cij = gps[i, j].Dep; // material matrix
          //    //DoubleVector stressS = gps[i, j].lastStress; // stress

          //    DoubleMatrix BTDB = (DoubleMatrix)B.Transpose() * Cij * B;
          //    double detJ = Math.Abs(J.Determinant());
          //    Ke += gps[i, j].weight * detJbar * detJ * BTDB;
          //}

          //if (nm == 2)
          //{
          //    // double C_rog = 705.0;
          //    double C_rog = CoefficientRamOsg(table2);//10.7*6.895;//// thông số ramber-osg 
          //    double m = ExponentRamOsg(table2);// 0.2;
          //    sigmaY = sima_offset(table1, 0.001);
          //    double f_yield_tr;
          //    if (alpha_n == 0)
          //    {
          //        //f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * sima_offset(table1,0.002);
          //        f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * (sigmaY + C_rog * Math.Pow(alpha_n, m));
          //    }
          //    else
          //    {
          //        f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * (C_rog * Math.Pow(alpha_n, m)); // vonMises yield function
          //    }
          //    if (f_yield_tr < 0)
          //    // elastic
          //    {
          //        gps[i, j].currentStress = stress_tr;
          //        gps[i, j].currentAlpha = alpha_n;
          //        gps[i, j].currentPlasticStrain = eps_pl_n;
          //        gps[i, j].currentBackStress = backstress_n;
          //        gps[i, j].Dep = D;
          //    }
          //    else
          //    {

          //        if (alpha_n == 0)
          //        {
          //            alpha_n += 0.0005;
          //        }
          //        double dgamma = 0;
          //        double gdgama = -Math.Sqrt(2.0 / 3.0) * (C_rog * Math.Pow(alpha_n, m)) + norm_beta;
          //        double dgdgama;
          //        double alphak;
          //        alphak = alpha_n;
          //        double TOL = 10e-5;
          //        int maxIter = 50;
          //        int nIter = 0;

          //        double abc = 2.0 / 3.0;

          //        do
          //        {

          //            nIter++;
          //            gdgama = -Math.Sqrt(abc) * (C_rog * Math.Pow(alphak, m)) + norm_beta - (2.0 * G * dgamma);

          //            dgdgama = -2.0 * G * (1.0 + (m * C_rog * Math.Pow(alphak, m - 1)) / (3.0 * G));
          //            dgamma -= (gdgama / dgdgama);
          //            alphak = alpha_n + Math.Sqrt(abc) * dgamma;
          //        }

          //        while (Math.Abs(gdgama) > TOL);


          //        DoubleVector stress = stress_tr - 2.0 * G * dgamma * n_vector;

          //        // Update consistent Elastoplastic material property tensor
          //        DoubleMatrix Cep_abstr1 = 2.0 * G * (1.0 / (1.0 + (m * C_rog * Math.Pow(alpha_n, m - 1)) / (3.0 * G)) - 2.0 * G * dgamma / norm_beta) * ((DoubleMatrix)n_vector.OuterProduct(n_vector));

          //        DoubleMatrix Id = (DoubleMatrix)Idev.Clone();
          //        Id[3, 3] = Id[3, 3] / 2.0;
          //        DoubleMatrix Cep_abstr2 = (4.0 * G * G * dgamma / norm_beta) * Id;
          //        DoubleMatrix Cep = D - Cep_abstr1 - Cep_abstr2;

          //        double alpha = alpha_n + Math.Sqrt(2.0 / 3.0) * dgamma;
          //        DoubleVector eps_pl = eps_pl_n + dgamma * n_vector;
          //        DoubleVector qCurrent = backstress_n + Math.Sqrt(2.0 / 3.0) * (K0 * dgamma) * n_vector;
          //        // Store
          //        gps[i, j].currentBackStress = qCurrent;
          //        gps[i, j].currentAlpha = alpha;
          //        gps[i, j].currentPlasticStrain = eps_pl;
          //        gps[i, j].currentStress = stress;
          //        gps[i, j].Dep = Cep;
          //    }

          //    // Assembly
          //    DoubleMatrix Cij = gps[i, j].Dep; // material matrix
          //    //DoubleVector stressS = gps[i, j].lastStress; // stress
          //    DoubleMatrix BTDB = (DoubleMatrix)B.Transpose() * Cij * B;
          //    double detJ = Math.Abs(J.Determinant());
          //    Ke += gps[i, j].weight * detJbar * detJ * BTDB;
          //}
          ////////////////
          //if (nm == 4)
          //{
          //    // double C_rog = 705.0;
          //    double C_rog = CoefficientRamOsg(table2, 0.11);//10.7*6.895;//// thông số ramber-osg 
          //    double m = ExponentRamOsg2(table2);// 0.2;
          //    sigmaY = sima_offset(table1, 0.002);
          //    double f_yield_tr;
          //    //if (alpha_n == 0)
          //    //{
          //    //    f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * (sigmaY);
          //    //}

          //    //else
          //    //  {
          //    f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * (sigmaY + C_rog * Math.Pow(alpha_n, m)); // vonMises yield function
          //    // }

          //    if (f_yield_tr < 0)
          //    // elastic
          //    {
          //        gps[i, j].currentStress = stress_tr;
          //        gps[i, j].currentAlpha = alpha_n;
          //        gps[i, j].currentPlasticStrain = eps_pl_n;
          //        gps[i, j].currentBackStress = backstress_n;
          //        gps[i, j].Dep = D;
          //    }
          //    else
          //    {


          //        double dgamma = 0;
          //        double gdgama = -Math.Sqrt(2.0 / 3.0) * (sigmaY + C_rog * Math.Pow(alpha_n, m)) + norm_beta;
          //        double dgdgama;
          //        double alphak;
          //        alphak = alpha_n;
          //        double TOL = 10e-5;
          //        int maxIter = 50;
          //        int nIter = 0;
          //        if (alpha_n == 0)
          //        {
          //            alpha_n = 0.0005;
          //        }

          //        double abc = 2.0 / 3.0;

          //        do
          //        {

          //            nIter++;
          //            gdgama = -Math.Sqrt(abc) * (sigmaY + C_rog * Math.Pow(alphak, m)) + norm_beta - (2.0 * G * dgamma);

          //            dgdgama = -2.0 * G * (1.0 + (m * C_rog * Math.Pow(alphak, m - 1)) / (3.0 * G));
          //            dgamma -= (gdgama / dgdgama);
          //            alphak = alpha_n + Math.Sqrt(abc) * dgamma;
          //        }

          //        while (Math.Abs(gdgama) > TOL);


          //        DoubleVector stress = stress_tr - 2.0 * G * dgamma * n_vector;



          //        // Update consistent Elastoplastic material property tensor
          //        DoubleMatrix Cep_abstr1 = 2.0 * G * (1.0 / (1.0 + (m * C_rog * Math.Pow(alpha_n, m - 1)) / (3.0 * G)) - 2.0 * G * dgamma / norm_beta) * ((DoubleMatrix)n_vector.OuterProduct(n_vector));

          //        DoubleMatrix Id = (DoubleMatrix)Idev.Clone();
          //        Id[3, 3] = Id[3, 3] / 2.0;
          //        DoubleMatrix Cep_abstr2 = (4.0 * G * G * dgamma / norm_beta) * Id;
          //        DoubleMatrix Cep = D - Cep_abstr1 - Cep_abstr2;

          //        double alpha = alpha_n + Math.Sqrt(2.0 / 3.0) * dgamma;
          //        DoubleVector eps_pl = eps_pl_n + dgamma * n_vector;
          //        DoubleVector qCurrent = backstress_n + Math.Sqrt(2.0 / 3.0) * (K0 * dgamma) * n_vector;
          //        // Store
          //        gps[i, j].currentBackStress = qCurrent;
          //        gps[i, j].currentAlpha = alpha;
          //        gps[i, j].currentPlasticStrain = eps_pl;
          //        gps[i, j].currentStress = stress;
          //        gps[i, j].Dep = Cep;
          //    }

          //    // Assembly
          //    DoubleMatrix Cij = gps[i, j].Dep; // material matrix
          //    //DoubleVector stressS = gps[i, j].lastStress; // stress
          //    DoubleMatrix BTDB = (DoubleMatrix)B.Transpose() * Cij * B;
          //    double detJ = Math.Abs(J.Determinant());
          //    Ke += gps[i, j].weight * detJbar * detJ * BTDB;
          //}
          ///////////////

          //if (nm == 3)
          //{

          //    //// G(alpha)=simaY+(sima_U- simaY)*(1-e^(-deta*alpha))
          //    double sima_U = 991.0;//// thông số EXP
          //    // double sima_U = SigmaU(table2, sima_offset(table1), Emodun(table1));
          //    sigmaY = sima_offset(table1, 0.002);
          //    double deta = ExponBeta(table2, sima_U);
          //    //double deta = 16.93;
          //    double C_ex = sima_U - sigmaY;
          //    double f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * (sigmaY + C_ex * (1.0 - Math.Exp(-(deta * alpha_n)))); // vonMises yield function
          //    if (f_yield_tr < 0)
          //    // elastic
          //    {
          //        gps[i, j].currentStress = stress_tr;
          //        gps[i, j].currentAlpha = alpha_n;
          //        gps[i, j].currentPlasticStrain = eps_pl_n;
          //        gps[i, j].currentBackStress = backstress_n;
          //        gps[i, j].Dep = D;
          //    }
          //    else
          //    {
          //        double dgamma = 0.0;
          //        //double C_ex = sima_U - sigmaY;
          //        double gdgama = -Math.Sqrt(2.0 / 3.0) * (sigmaY + C_ex * (1.0 - Math.Exp(-deta * alpha_n))) + norm_beta;
          //        double dgdgama;
          //        double alphak;
          //        alphak = alpha_n;
          //        double TOL = 10e-5;
          //        int maxIter = 20;
          //        int nIter = 0;
          //        double abc = 2.0 / 3.0;

          //        do
          //        {

          //            nIter++;
          //            gdgama = -Math.Sqrt(2.0 / 3.0) * (sigmaY + C_ex * (1.0 - Math.Exp(-deta * alphak))) + norm_beta - (2.0 * G * dgamma);
          //            dgdgama = -2.0 * G * (1.0 + C_ex * deta * Math.Exp(-(deta * alphak)) / (3.0 * G));
          //            dgamma -= (gdgama / dgdgama);
          //            alphak = alpha_n + Math.Sqrt(abc) * dgamma;
          //        }

          //        while ((Math.Abs(gdgama) > TOL));//&& (nIter < maxIter));


          //        DoubleVector stress = stress_tr - 2.0 * G * dgamma * n_vector;

          //        // Update consistent Elastoplastic material property tensor
          //        //if (alpha_n == 0)
          //        //{
          //        //    alpha_n += 0.00005;
          //        //}

          //        //DoubleMatrix Cep_abstr1 = 2.0 * G * (1.0 / (1.0 + (0 * H0 + m * C_rog * Math.Pow(alpha_n, m - 1)) / (3.0 * G)) - 2.0 * G * dgamma / norm_beta)
          //        DoubleMatrix Cep_abstr1 = 2.0 * G * (1.0 / (1.0 + (C_ex * deta * Math.Exp(-(deta * alpha_n))) / (3.0 * G)) - 2.0 * G * dgamma / norm_beta) * ((DoubleMatrix)n_vector.OuterProduct(n_vector));

          //        DoubleMatrix Id = (DoubleMatrix)Idev.Clone();
          //        Id[3, 3] = Id[3, 3] / 2.0;
          //        DoubleMatrix Cep_abstr2 = (4.0 * G * G * dgamma / norm_beta) * Id;
          //        DoubleMatrix Cep = D - Cep_abstr1 - Cep_abstr2;

          //        double alpha = alpha_n + Math.Sqrt(2.0 / 3.0) * dgamma;
          //        DoubleVector eps_pl = eps_pl_n + dgamma * n_vector;
          //        DoubleVector qCurrent = backstress_n + Math.Sqrt(2.0 / 3.0) * (K0 * dgamma) * n_vector;
          //        // Store
          //        gps[i, j].currentBackStress = qCurrent;
          //        gps[i, j].currentAlpha = alpha;
          //        gps[i, j].currentPlasticStrain = eps_pl;
          //        gps[i, j].currentStress = stress;
          //        gps[i, j].Dep = Cep;
          //    }

          //    // Assembly
          //    DoubleMatrix Cij = gps[i, j].Dep; // material matrix
          //    //DoubleVector stressS = gps[i, j].lastStress; // stress
          //    DoubleMatrix BTDB = (DoubleMatrix)B.Transpose() * Cij * B;
          //    double detJ = Math.Abs(J.Determinant());
          //    Ke += gps[i, j].weight * detJbar * detJ * BTDB;
          //}

        }
      }
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
      DoubleMatrix B = new DoubleMatrix(4, d * nen);

      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          double detJbar = 1.0 / 4.0
                  * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                  * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2]);
          //DoubleMatrix dNdxi = gps[i, j].dNdxi;
          //DoubleMatrix J = JacobianAt(dNdxi);
          //DoubleMatrix invertJ = (DoubleMatrix)J.Inverse();
          DoubleMatrix dNdx = (DoubleMatrix)gps[i, j].GetValue(DataInGausspoint.dNdX);//dNdxi * invertJ;
          for (int k = 0; k < nen; k++)
          {
            B[0, 2 * k] = dNdx[k, 0];
            B[3, 2 * k] = dNdx[k, 1];
            B[1, 2 * k + 1] = dNdx[k, 1];
            B[3, 2 * k + 1] = dNdx[k, 0];
          }

          DoubleVector stressS = (DoubleVector)gps[i, j].GetValue(DataInGausspoint.currentStress); // stress

          double detJ = Math.Abs((double)gps[i, j].GetValue(DataInGausspoint.detJ));
          fi += gps[i, j].weight * detJbar * detJ * MatrixFunctions.Product(B.Transpose(), stressS);
          // Notice:
          // In the iterative solver (i.e SolverLoadControlled or SolverArclengthControlled), when it is converged, the HistoryValuesLast have to be updated by: HistoryValuesLast = HistoryValuesCurrent
        }
      }
    }

    public double[,] Transformtable(double[,] table)
    {
      ///// ham nay dung cho bai toan multiliner////
      int legh = (table.Length / 2 - 1);
      if (legh < 2)
      {
        throw new NotImplementedException("phai nhap nhieu hon hai diem");
      }
      //strain[1] = 0.0513;
      double Emodun = table[1, 1] / table[0, 1];
      double sigma_offset = table[1, 1];
      double[,] str_tressPlas = new double[2, legh];
      //str_tressPlas[0, 0] = 0.0;
      //str_tressPlas[1, 0] = sigma_offset;
      //  str_tressPlas[1, 0] = sigma_offset;

      for (int i = 0; i < legh; i++)
      {
        str_tressPlas[0, i] = table[0, i + 1] - table[1, i + 1] / Emodun;//// chu y ?????
        str_tressPlas[1, i] = table[1, i + 1];
      }

      //FolderBrowserDialog fbd = new FolderBrowserDialog();
      //DialogResult result = fbd.ShowDialog();
      //string path = "D:\\";
      //if (result == DialogResult.OK)
      //{
      //    path = fbd.SelectedPath;
      //}


      //using (System.IO.StreamWriter file = new System.IO.StreamWriter(path + "/stress_strainpl12.txt", true))
      //{
      //    int leng = str_tressPlas.Length / 2;
      //    for (int ii = 0; ii < leng; ii++)
      //    {
      //        double stress_strainpl = str_tressPlas[0, ii];
      //        file.WriteLine(stress_strainpl);
      //    }

      //}

      return str_tressPlas;


    }
    public double[,] Transformtable(double E, double[,] table, double sigma_offset)
    {
      ///// ham nay dung cho bai toan Ramber-Osgood va ham mu////
      int legh = (table.Length / 2) - 1;
      DoubleVector strainPlas = new DoubleVector(legh);
      DoubleVector sigma = new DoubleVector(legh);
      //strain[1] = 0.0513;
      double[,] str_tressPlas = new double[2, legh];
      str_tressPlas[0, 0] = 0.0;
      str_tressPlas[1, 0] = sigma_offset;
      //  str_tressPlas[1, 0] = sigma_offset;

      //  DoubleVector strain_p = new DoubleVector(legh);
      for (int i = 1; i < legh; i++)
      {
        str_tressPlas[0, i] = table[0, i] - sigma_offset / E;//// chu y ?????
        str_tressPlas[1, i] = table[1, i];
      }

      return str_tressPlas;
    }

    public double LinearInterpolation(double value, double[,] table, char opt)
    {
      double result = 1;
      switch (opt)
      {
        case 'x':
          int i = 1;
          while (table[i, 1] < value)
          {
            i++;
          }
          double x1 = table[0, i - 1];
          double y1 = table[1, i - 1];
          double x2 = table[0, i];
          double y2 = table[1, i];
          result = ((y2 - y1) * (value - x1)) / (x2 - x1) + y1;
          break;
        case 'y':
          i = 1;
          while (table[i, 2] < value)
          {
            i++;
          }
          x1 = table[1, i - 1];
          y1 = table[2, i - 1];
          x2 = table[1, i];
          y2 = table[2, i];
          result = ((x2 - x1) * (value - y1)) / (y2 - y1) + x1;
          break;

      }
      return result;
    }

    public double DplasticFun(double value, double[,] table)
    {
      double result = 0;

      int i = 1;
      while (table[i, 1] < value)
      {
        i++;
      }
      double x1 = table[0, i - 1];
      double y1 = table[1, i - 1];
      double x2 = table[0, i];
      double y2 = table[1, i];
      result = (y2 - y1) / (x2 - x1);

      return result;
    }

    //public double sima_offset(double[,] table, double x_offset)
    //{
    //    double ss1 = table[1, 1];
    //    double e1 = table[1, 0];

    //    double x_offstotal = x_offset + e1;
    //    double E = ss1 / e1;
    //    double result = 0;

    //    int i = 1;
    //    while (table[0, i] < x_offstotal)
    //    {
    //        i++;
    //    }
    //    double x1 = table[0, i - 1];
    //    double y1 = table[1, i - 1];
    //    double x2 = table[0, i];
    //    double y2 = table[1, i];

    //    return result = ((y2 - y1) * (x_offstotal - x1)) / (x2 - x1) + y1;
    //}

    public double sima_offset(double[,] table, double x_offset)
    {
      double ss1 = table[1, 1];
      double e1 = table[0, 1];
      double ss2 = table[1, 2];
      double e2 = table[0, 2];
      //double ss1 = table[1, 0];
      //double e1 =  table[0, 0];
      //double ss2 = table[1, 1];
      //double e2 =  table[0, 1];

      double E = ss1 / e1;

      double ep2 = e2 - ss2 / E;
      double epy = x_offset;

      return (epy * ss2 + ss1 * ep2) / (epy + ep2);
    }

    public double ModunPlastic(double[,] table, double x_offset)
    {
      int legh = (table.Length / 2);
      double sumxx = 0.0;
      double sumX = 0.0;
      double sumY = 0.0;
      double sumXY = 0.0;
      double aa;
      //table[0, 0] = 0.001;
      table[0, 0] = x_offset;
      for (int i = 0; i < legh; i++)
      {

        aa = table[0, i];
        //double bb = table[1, i];
        double bb = table[1, i] - 690;
        sumxx += aa * aa;
        sumX += aa;
        sumY += bb;
        sumXY += aa * bb;

      }

      double ModunPlastic = (sumxx - (sumX * sumX) / (legh)) / (sumXY - (sumX * sumY) / (legh));

      return 1.0 / ModunPlastic;
    }

    public double ExponentRamOsg(double[,] table)
    {
      int legh = (table.Length / 2);
      double sumxx = 0.0;
      double sumX = 0.0;
      double sumY = 0.0;
      double sumXY = 0.0;
      double aa;
      table[0, 0] = 0.0001;
      for (int i = 0; i < legh; i++)
      {

        aa = Math.Log(Math.Abs(table[0, i])) / Math.Log(Math.E);
        double bb = Math.Log(Math.Abs(table[1, i])) / Math.Log(Math.E);
        sumxx += aa * aa;
        sumX += aa;
        sumY += bb;
        sumXY += aa * bb;
      }

      double ExponentRamOsg_n = (sumxx - (sumX * sumX) / (legh)) / (sumXY - (sumX * sumY) / (legh));

      return 1.0 / ExponentRamOsg_n;
    }

    public double ExponentRamOsg2(double[,] table)
    {
      double simay = 690.0;
      int legh = (table.Length / 2);
      double sumxx = 0.0;
      double sumX = 0.0;
      double sumY = 0.0;
      double sumXY = 0.0;

      double aa;
      double bb;
      table[0, 0] = 0.00005;
      for (int i = 0; i < legh; i++)
      {
        table[1, 0] = 690.1;
        aa = Math.Log(Math.Abs(table[0, i])) / Math.Log(Math.E);
        //if (table[1,i]==simay)
        //{
        //    bb = Math.Log(Math.Abs(1.0001*table[1, i] - simay)) / Math.Log(Math.E);
        //    i++;
        //}

        bb = Math.Log(Math.Abs(table[1, i] - simay)) / Math.Log(Math.E);
        sumxx += aa * aa;
        sumX += aa;
        sumY += bb;
        sumXY += aa * bb;


      }

      double ExponentRamOsg_n = (sumxx - (sumX * sumX) / (legh)) / (sumXY - (sumX * sumY) / (legh));

      return 1.0 / ExponentRamOsg_n;
    }



    /// <summary>
    /// ////////////////////
    /// </summary>
    /// <param name="table"></param>
    /// <returns></returns>
    /// 



    public double CoefficientRamOsg(double[,] table)
    {
      table[0, 0] = 0.0001;
      int legh = table.Length / 2;
      double aa;
      double sumX = 0.0;
      double sumY = 0.0;

      for (int i = 0; i < legh; i++)
      {


        aa = Math.Log(Math.Abs(table[0, i])) / Math.Log(Math.E);
        double bb = Math.Log(Math.Abs(table[1, i])) / Math.Log(Math.E);
        sumX += aa;
        sumY += bb;

      }

      double lnCoeffRamOsg_n = (sumY - ExponentRamOsg(table) * sumX) / legh;
      double CoefficientRamOsg_n = Math.Exp(lnCoeffRamOsg_n);

      return CoefficientRamOsg_n;
    }


    /// <summary>
    /// ////////////
    /// </summary>
    /// <param name="table"></param>
    /// <returns></returns>
    public double CoefficientRamOsg(double[,] table, double a)
    {
      table[0, 0] = 0.0001;
      double simay = 690.0;
      int legh = table.Length / 2;
      double aa;
      double bb;
      double sumX = 0.0;
      double sumY = 0.0;
      table[1, 0] = 690.1;

      for (int i = 0; i < legh; i++)
      {

        aa = Math.Log(Math.Abs(table[0, i])) / Math.Log(Math.E);
        //if (table[1, i] == simay)
        //{
        //    bb = Math.Log(Math.Abs(1.0001*table[1, i] - simay)) / Math.Log(Math.E);
        //    i++;
        //}

        bb = Math.Log(Math.Abs(table[1, i] - simay)) / Math.Log(Math.E);
        sumX += aa;
        sumY += bb;

      }

      double lnCoeffRamOsg_n = (sumY - ExponentRamOsg2(table) * sumX) / legh;
      double CoefficientRamOsg_n = Math.Exp(lnCoeffRamOsg_n);

      return CoefficientRamOsg_n;
    }

    public double Emodun(double[,] table)
    {
      double Emodun = table[1, 1] / table[0, 1];
      return Emodun;
    }
    public double SigmaU(double[,] table, double sima02, double E)
    {

      double n = 1.0 / ExponentRamOsg(table);
      double e = sima02 / E;
      double aa = (0.2 + 185.0 * e) / (1.0 - 0.0375 * (n - 5));
      double Sigma_u = sima02 / aa;
      return Sigma_u;

    }

    public double ExponBeta(double[,] table, double Sima_U)
    {
      int legh = (table.Length / 2);
      double sumxx = 0.0;
      double sumX = 0.0;
      double sumY = 0.0;
      double sumXY = 0.0;
      double aa;
      double bb;
      table[0, 0] = 0.0001;
      for (int i = 0; i < legh; i++)
      {

        aa = table[0, i];
        ////while (Sima_U > table[1, i])

        bb = Math.Log(Math.Abs(Sima_U - table[1, i])) / Math.Log(Math.E);
        sumxx += aa * aa;
        sumX += aa;
        sumY += bb;
        sumXY += aa * bb;

      }

      double Expon_Beta = (sumxx - (sumX * sumX) / (legh)) / (sumXY - (sumX * sumY) / (legh));

      return -(1.0 / Expon_Beta);
    }


    public void outputCoffive(double[,] table)
    {

      double[,] str_tressPlas = Transformtable(table);
      double CoefficientRamOsg_n = CoefficientRamOsg(str_tressPlas);
      double ExponentRamOsg_n = ExponentRamOsg(str_tressPlas);
      double CoefficientRamOsg_n1 = CoefficientRamOsg(str_tressPlas, 1.0);
      double ExponentRamOsg_n1 = ExponentRamOsg2(str_tressPlas);
      double Hmod = ModunPlastic(str_tressPlas, 0.0001);
      double Emod = Emodun(table);
      double Expo = ExponBeta(str_tressPlas, 991.0);
      double tangetModu = Emod * Hmod / (Emod + Hmod);


      //FolderBrowserDialog fbd = new FolderBrowserDialog();
      //DialogResult result = fbd.ShowDialog();
      //string path = "D:\\";
      //if (result == DialogResult.OK)
      //{
      //    path = fbd.SelectedPath;
      //}


      using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\VS9 X64Bit\Desktop\outputfile\stress_strainplfinal2-6.txt", true))
      {
        file.WriteLine("CoefficientRamOsg_n: " + CoefficientRamOsg_n.ToString());
        file.WriteLine("ExponentRamOsg_n: " + ExponentRamOsg_n.ToString());
        file.WriteLine("CoefficientRamOsg_n1: " + CoefficientRamOsg_n1.ToString());
        file.WriteLine("ExponentRamOsg_n1: " + ExponentRamOsg_n1.ToString());
        file.WriteLine("Hmodun: " + Hmod.ToString());
        file.WriteLine("Emodun: " + Emod.ToString());
        file.WriteLine("Coefficient_Beta: " + Expo.ToString());
        file.WriteLine("TangetModun: " + tangetModu.ToString());
        int leng = str_tressPlas.Length / 2;

        file.Write("strainPlasic: ");
        for (int ii = 0; ii < leng; ii++)
        {
          double stress_strainpl = str_tressPlas[0, ii];
          file.Write(stress_strainpl);
          file.Write(" ");
        }
        file.WriteLine();
        file.Write("stress: ");
        for (int ii = 0; ii < leng; ii++)
        {
          double stress_strainpl = str_tressPlas[1, ii];
          file.Write(stress_strainpl);
          file.Write(" ");
        }

      }

    }
    public double[,] Readfile(string ten)
    {
      double[,] result = null;

      using (System.IO.StreamReader file = new System.IO.StreamReader(ten, true))
      {
        string text = file.ReadLine();
        int d = 0;
        while (text != null)
        {
          string[] temp = text.Split(new Char[] { ' ' });
          if (result == null)
          {
            result = new double[2, temp.Length];
          }
          for (int i = 0; i < temp.Length; i++)
          {
            result[d, i] = double.Parse(temp[i].Trim());

          }
          d++;
          text = file.ReadLine();
        }
        file.Close();
        file.Dispose();

      }
      return result;
    }

    public override DoubleVector StressAt(params double[] xi)
    {
      throw new NotImplementedException();
    }

    public override DoubleVector StrainAt(params double[] xi)
    {
      throw new NotImplementedException();
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
