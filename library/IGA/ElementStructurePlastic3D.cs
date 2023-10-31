using CenterSpace.NMath.Core;
using DEMSoft.EngineeringData;
using DEMSoft.NURBS;
using System;
using System.Collections.Generic;
namespace DEMSoft.IGA
{
  public class ElementStructurePlastic3D : AbstractElementStructure3D
  {
    public ElementStructurePlastic3D(AbstractPatch3D mesh, int id)
        : base(mesh, id)
    {
    }

    private DoubleMatrix CreateMaterialMatrix()
    {
      var eModulus = Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
      var nu = Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
      DoubleMatrix C = new DoubleMatrix(6, 6);
      double a = eModulus / ((1 + nu) * (1 - 2 * nu));
      C[0, 0] = a * (1 - nu);
      C[0, 1] = a * nu;
      C[0, 2] = a * nu;
      C[1, 0] = a * nu;
      C[1, 1] = a * (1 - nu);
      C[1, 2] = a * nu;
      C[2, 0] = a * nu;
      C[2, 1] = a * nu;
      C[2, 2] = a * (1 - nu);
      C[3, 3] = a * (1 - 2 * nu) / 2;
      C[4, 4] = a * (1 - 2 * nu) / 2;
      C[5, 5] = a * (1 - 2 * nu) / 2;
      return C;
    }

    public override void ComputeStiffnessMatrixElement(out DoubleMatrix Ke)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      int d = patch.GetCountField();
      TrivariateNURBSBasisFunction basis = (TrivariateNURBSBasisFunction)((NURBSVolume)(patch.GetGeometry())).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int nen = (p + 1) * (q + 1) * (r + 1); // number of local basis functions
      var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply3 = basis.GetKnotVector(2).GetKnotVectorNoMultiplicity();
      int idx1 = patch.GetIPN(id, 0);
      int idx2 = patch.GetIPN(id, 1);
      int idx3 = patch.GetIPN(id, 2);
      Ke = new DoubleMatrix(d * nen, d * nen);
      DoubleMatrix B = new DoubleMatrix(6, d * nen);
      DoubleMatrix D = CreateMaterialMatrix();
      DoubleMatrix Idev = new DoubleMatrix(6, 6);
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
      Idev[4, 4] = 1.0;
      Idev[5, 5] = 1.0;
      DoubleVector vone = new DoubleVector(6);
      vone[0] = 1.0;
      vone[1] = 1.0;
      vone[2] = 1.0;

      double K = Material.GetProperty(MaterialPropertyName.BulkModulus).GetValueProperty();
      double G = Material.GetProperty(MaterialPropertyName.ShearModulus).GetValueProperty();
      DoubleMatrix Ce = K * MatrixFunctions.OuterProduct(vone, vone) + 2 * G * Idev;

      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {

            double detJbar = 1.0 / 8.0
                    * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                    * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2])
                    * (kvNoMulticiply3[idx3 + 1] - kvNoMulticiply3[idx3]);
            DoubleMatrix dNdxi = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdxi);
            DoubleMatrix J = JacobianAt(dNdxi);
            DoubleMatrix dNdx = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdX);//dNdxi * invertJ;
            for (int kk = 0; kk < nen; kk++)
            {
              B[0, 3 * kk] = dNdx[kk, 0];
              B[1, 3 * kk + 1] = dNdx[kk, 1];
              B[2, 3 * kk + 2] = dNdx[kk, 2];
              B[3, 3 * kk + 1] = dNdx[kk, 2];
              B[3, 3 * kk + 2] = dNdx[kk, 1];
              B[4, 3 * kk] = dNdx[kk, 2];
              B[4, 3 * kk + 2] = dNdx[kk, 0];
              B[5, 3 * kk] = dNdx[kk, 1];
              B[5, 3 * kk + 1] = dNdx[kk, 0];
            }

            //double [,] table3=new double[2,4];
            //table1[0,0]=0.0;
            //table1[0,1]=2.0;
            //table1[0,2]=4.0;
            //table1[0,3]=5.0;
            //table1[1,0]=10.0;
            //table1[1,1]=22.0;
            //table1[1,2]=34.0;
            //table1[1,3]=45.0;
            //string path = "";
            //OpenFileDialog op = new OpenFileDialog();
            //op.ShowDialog();
            //op.Filter = "txt file|*.txt";
            //if (op.ShowDialog() == DialogResult.OK)
            //{
            //    path = op.FileName;
            //}


            // double[,] table1 = Readfile(@"C:\Users\VS9 X64Bit\Desktop\stress_tra\stress_tra.txt");
            // ============================================
            // Update stress using Return-Mapping algorithm
            // ============================================
            // History (converged) values of last load step
            //DoubleVector stress_n = gps[i, j].lastStress;
            DoubleVector eps_pl_n = (DoubleVector)gps[i, j, k].GetValue(DataInGausspoint.lastPlasticStrain);
            DoubleVector backstress_n = (DoubleVector)gps[i, j, k].GetValue(DataInGausspoint.lastBackStress);
            double alpha_n = (double)gps[i, j, k].GetValue(DataInGausspoint.lastAlpha);
            // Strain
            DoubleVector ue = GetDisplacementLocal();////////////
            DoubleVector strain = MatrixFunctions.Product(B, ue); // [eps_xx; epx_yy; epx_zz; 2*eps_xy; 2*epsxz; 2*epsyz]
            DoubleVector eps_total = strain.DeepenThisCopy();
            // [eps_xx; epx_yy;  epx_xy; epx_xz; epx_yz]
            eps_total[3] = eps_total[3] / 2.0;
            eps_total[4] = eps_total[4] / 2.0;
            eps_total[5] = eps_total[5] / 2.0;
            // Trial state
            DoubleVector stress_tr = MatrixFunctions.Product(Ce, (eps_total - eps_pl_n));
            // DoubleVector e_n = Idev*eps_total;/////////////////////
            DoubleVector s_tr = MatrixFunctions.Product(Idev, stress_tr);  // deviatoric stress  
                                                                           /// DoubleVector s_tr = 2 * G * (e_n - eps_pl_n);/////////////////////
            DoubleVector beta_tr = s_tr - backstress_n; // relative stress
            double norm_beta = Math.Sqrt(beta_tr[0] * beta_tr[0] + beta_tr[1] * beta_tr[1] + beta_tr[2] * beta_tr[2] + 2 * beta_tr[3] * beta_tr[3] + 2 * beta_tr[4] * beta_tr[4] + 2 * beta_tr[5] * beta_tr[5]);
            // Yield condition
            DoubleVector n_vector = beta_tr / norm_beta;

            if (Material.GetProperty(MaterialPropertyName.BilinearIsotropicHardening) != null || Material.GetProperty(MaterialPropertyName.BilinearKinematicHardening) != null)
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
              double f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * (sigmaY + H0 * alpha_n); // vonMises yield function
                                                                                              ////double f_yield_tr = norm_beta - (sigmaY + H0 * alpha_n);
              GaussPoints gpsijk = gps[i, j, k];
              if (f_yield_tr < ftol)
              // elastic
              {
                gpsijk.SetValue(DataInGausspoint.currentStress, stress_tr);
                gpsijk.SetValue(DataInGausspoint.currentAlpha, alpha_n);
                gpsijk.SetValue(DataInGausspoint.currentPlasticStrain, eps_pl_n);
                gpsijk.SetValue(DataInGausspoint.currentBackStress, backstress_n);
                gpsijk.SetValue(DataInGausspoint.currentDep, D);
                gpsijk.SetValue(DataInGausspoint.currentStrain, eps_total);

                //gps[i, j, k].currentStress = stress_tr;
                //gps[i, j, k].currentAlpha = alpha_n;
                //gps[i, j, k].currentPlasticStrain = eps_pl_n;
                //gps[i, j, k].currentBackStress = backstress_n;
                //gps[i, j, k].currentDep = D;
                //gps[i, j, k].currentStrain = eps_total;
              }
              else
              {
                // plastic: Update stress by Return-Mapping algorithm
                double dgamma = f_yield_tr / (2.0 * G + 2.0 * (H0 + K0) / 3.0);

                // Update consistent Elastoplastic material property tensor
                DoubleMatrix Cep_abstr1 = 2.0 * G * (1.0 / (1.0 + (H0 + K0) / (3.0 * G)) - 2.0 * G * dgamma / norm_beta) * MatrixFunctions.OuterProduct(n_vector, n_vector);

                DoubleMatrix Id = Idev.DeepenThisCopy();
                Id[3, 3] = Id[3, 3] / 2.0;
                Id[4, 4] = Id[4, 4] / 2.0;
                Id[5, 5] = Id[5, 5] / 2.0;
                DoubleMatrix Cep_abstr2 = (4.0 * G * G * dgamma / norm_beta) * Id;
                DoubleMatrix Cep = D - Cep_abstr1 - Cep_abstr2;

                // Store
                gpsijk.SetValue(DataInGausspoint.currentStress, stress_tr - 2.0 * G * dgamma * n_vector);
                gpsijk.SetValue(DataInGausspoint.currentAlpha, alpha_n + Math.Sqrt(2.0 / 3.0) * dgamma);
                gpsijk.SetValue(DataInGausspoint.currentPlasticStrain, eps_pl_n + dgamma * n_vector);
                gpsijk.SetValue(DataInGausspoint.currentBackStress, backstress_n + 2.0 / 3.0 * K0 * dgamma * n_vector);
                gpsijk.SetValue(DataInGausspoint.currentDep, Cep);
                gpsijk.SetValue(DataInGausspoint.currentStrain, eps_total);

                //gps[i, j, k].currentStress = stress_tr - 2.0 * G * dgamma * n_vector;
                //gps[i, j, k].currentAlpha = alpha_n + Math.Sqrt(2.0 / 3.0) * dgamma;
                //gps[i, j, k].currentPlasticStrain = eps_pl_n + dgamma * n_vector;
                //gps[i, j, k].currentBackStress = backstress_n + 2.0 / 3.0 * K0 * dgamma * n_vector;
                //gps[i, j, k].currentDep = Cep;
                //gps[i, j, k].currentStrain = eps_total;
              }

              // Assembly
              DoubleMatrix Cijk = (DoubleMatrix)gpsijk.GetValue(DataInGausspoint.currentDep); // material matrix

              DoubleMatrix BTDB = NMathFunctions.TransposeProduct(B, MatrixFunctions.Product(Cijk, B));
              double detJ = Math.Abs(MatrixFunctions.Determinant(J));
              Ke += gpsijk.weight * detJbar * detJ * BTDB;
            }



            //int nm = 1;

            //if (nm == 0)
            //{
            //    double f_yield_tr;
            //    if (alpha_n==0)
            //    {
            //         f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * sigmaY;
            //    }
            //    else 
            //    {
            //       f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * (LinearInterpolation(alpha_n, Transformtable( table1), 'x'));
            //    }

            //   // vonMises yield function
            //    if (f_yield_tr < 0)
            //    // elastic
            //    {
            //        gps[i, j, k].currentStress = stress_tr.ToArray();
            //        gps[i, j, k].currentBackStress = backstress_n.ToArray();
            //        gps[i, j, k].currentAlpha = alpha_n;
            //        gps[i, j, k].currentPlasticStrain = eps_pl_n.ToArray();

            //        gps[i, j, k].Dep = D;
            //    }
            //    else
            //    {
            //        double dgamma = 0;
            //        double gdgama = -Math.Sqrt(2.0 / 3.0) * LinearInterpolation(alpha_n, Transformtable(table1), 'x') + norm_beta;
            //        double dgdgama;
            //        double alphak;
            //        alphak = alpha_n;
            //        double TOL = 10e-5;
            //        int maxIter = 10;
            //        int nIter = 0;

            //        //  double m = 5;
            //        do
            //        {

            //            nIter++;
            //            gdgama = -Math.Sqrt(2.0/3.0) * LinearInterpolation(alpha_n, Transformtable( table1), 'x') + norm_beta - 2.0 * G * dgamma;
            //            dgdgama = -2.0 * G * (1.0 + DplasticFun(alphak, Transformtable( table1)) / (3.0 * G));
            //            dgamma -= (gdgama / dgdgama);
            //            alphak = alpha_n + Math.Sqrt(2.0/3.0) * dgamma;
            //        }

            //        while ((Math.Abs(gdgama) > TOL) && (nIter < maxIter));
            //        DoubleVector stress = stress_tr - 2.0 * G * dgamma * n_vector;



            //        // Update consistent Elastoplastic material property tensor
            //        DoubleMatrix Cep_abstr1 = 2.0 * G * (1.0 / (1.0 + (DplasticFun(alpha_n, Transformtable(table1))) / (3.0 * G)) - 2.0 * G * dgamma / norm_beta) * ((DoubleMatrix)n_vector.OuterProduct(n_vector));


            //        DoubleMatrix Id = (DoubleMatrix)Idev.Clone();
            //        Id[3, 3] = Id[3, 3] / 2.0;
            //        Id[4, 4] = Id[4, 4] / 2.0;
            //        Id[5, 5] = Id[5, 5] / 2.0;
            //        DoubleMatrix Cep_abstr2 = (4.0 * G * G * dgamma / norm_beta) * Id;
            //        DoubleMatrix Cep = D - Cep_abstr1 - Cep_abstr2;

            //        double alpha = alpha_n + Math.Sqrt(2.0 / 3.0) * dgamma;
            //        DoubleVector eps_pl = eps_pl_n + dgamma * n_vector;

            //        // Store
            //        gps[i, j, k].currentAlpha = alpha;
            //        gps[i, j, k].currentPlasticStrain = eps_pl.ToArray();
            //        gps[i, j, k].currentStress = stress.ToArray();
            //        gps[i, j, k].Dep = Cep;
            //        //gps[i, j].Dep = D;
            //    }

            //    // Assembly
            //    DoubleMatrix Cijk = gps[i, j, k].Dep; // material matrix
            //    //DoubleVector stressS = gps[i, j].lastStress; // stress

            //    DoubleMatrix BTDB = (DoubleMatrix)B.Transpose() * Cijk * B;
            //    double detJ = Math.Abs(J.Determinant());
            //    Ke += gps[i, j, k].weight * detJbar * detJ * BTDB;
            //}

            //if (nm == 1)
            //{
            //    double f_yield_tr = norm_beta - Math.Sqrt(2.0 / 3.0) * (sigmaY + H0 * alpha_n); // vonMises yield function
            //    ////double f_yield_tr = norm_beta - (sigmaY + H0 * alpha_n);
            //    if (f_yield_tr < 0)
            //    // elastic
            //    {
            //        gps[i, j, k].currentStress = stress_tr;
            //        gps[i, j, k].currentAlpha = alpha_n;
            //        gps[i, j, k].currentPlasticStrain = eps_pl_n;
            //        gps[i, j, k].currentBackStress = backstress_n;
            //        gps[i, j, k].currentDep = D;
            //    }
            //    else
            //    {
            //        // plastic: Update stress by Return-Mapping algorithm
            //        double dgamma = f_yield_tr / (2.0 * G + 2.0 * K0 / 3.0 + 2.0 * H0 / 3.0);
            //        DoubleVector stress = stress_tr - 2.0 * G * dgamma * n_vector;



            //        // Update consistent Elastoplastic material property tensor
            //        DoubleMatrix Cep_abstr1 = 2.0 * G * (1.0 / (1.0 + (H0 + K0) / (3.0 * G)) - 2.0 * G * dgamma / norm_beta) * n_vector.OuterProduct(n_vector);

            //        //DoubleMatrix Cep_abstr1 = 4.0 * Gmod * Gmod * (1.0 / (2.0 * Gmod + 2.0 * H0 / 3.0 + 2.0 * K0 / 3.0) - dgamma / norm_beta) * (n_vecM * (DoubleMatrix)(n_vecM.Transpose()));
            //        DoubleMatrix Id = Idev.Clone();
            //        Id[3, 3] = Id[3, 3] / 2.0;
            //        Id[4, 4] = Id[4, 4] / 2.0;
            //        Id[5, 5] = Id[5, 5] / 2.0;
            //        DoubleMatrix Cep_abstr2 = (4.0 * G * G * dgamma / norm_beta) * Id;
            //        DoubleMatrix Cep = D - Cep_abstr1 - Cep_abstr2;

            //        double alpha = alpha_n + Math.Sqrt(2.0 / 3.0) * dgamma;
            //        DoubleVector eps_pl = eps_pl_n + dgamma * n_vector;
            //        DoubleVector qCurrent = backstress_n + 2.0 * K0 / 3.0 * dgamma * n_vector;


            //        // Store
            //        gps[i, j, k].currentBackStress = qCurrent;
            //        gps[i, j, k].currentAlpha = alpha;
            //        gps[i, j, k].currentPlasticStrain = eps_pl;
            //        gps[i, j, k].currentStress = stress;
            //        gps[i, j, k].currentDep = Cep;
            //    }

            //    // Assembly
            //    DoubleMatrix Cijk = gps[i, j, k].currentDep; // material matrix
            //    //DoubleVector stressS = gps[i, j, k].lastStress; // stress

            //    DoubleMatrix BTDB = (DoubleMatrix)B.Transpose() * Cijk * B;
            //    double detJ = Math.Abs(J.Determinant());
            //    Ke += gps[i, j, k].weight * detJbar * detJ * BTDB;
            //    // Notice:
            //    // In the iterative solver (i.e SolverLoadControlled or SolverArclengthControlled), when it is converged, the HistoryValuesLast have to be updated by: HistoryValuesLast = HistoryValuesCurrent


            //}
          }
        }
      }
    }


    public override void ComputeInternalForceElement(out DoubleVector fi)
    {
      PatchStructure3D patch = (PatchStructure3D)this.patch;
      int d = patch.GetCountField();
      TrivariateNURBSBasisFunction basis = (TrivariateNURBSBasisFunction)((NURBSVolume)(patch.GetGeometry())).Basis;
      int p = basis.GetDegree(0);
      int q = basis.GetDegree(1);
      int r = basis.GetDegree(2);
      int nen = (p + 1) * (q + 1) * (r + 1); // number of local basis functions
      var kvNoMulticiply1 = basis.GetKnotVector(0).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply2 = basis.GetKnotVector(1).GetKnotVectorNoMultiplicity();
      var kvNoMulticiply3 = basis.GetKnotVector(2).GetKnotVectorNoMultiplicity();
      int idx1 = patch.GetIPN(id, 0);
      int idx2 = patch.GetIPN(id, 1);
      int idx3 = patch.GetIPN(id, 2);
      fi = new DoubleVector(d * nen);
      DoubleMatrix B = new DoubleMatrix(6, d * nen);

      for (int i = 0; i < numberOfGaussPointOnEachDirection; i++)
      {
        for (int j = 0; j < numberOfGaussPointOnEachDirection; j++)
        {
          for (int k = 0; k < numberOfGaussPointOnEachDirection; k++)
          {

            double detJbar = 1.0 / 8.0
                    * (kvNoMulticiply1[idx1 + 1] - kvNoMulticiply1[idx1])
                    * (kvNoMulticiply2[idx2 + 1] - kvNoMulticiply2[idx2])
                    * (kvNoMulticiply3[idx3 + 1] - kvNoMulticiply3[idx3]);
            DoubleMatrix dNdxi = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdxi);
            DoubleMatrix J = JacobianAt(dNdxi);
            //DoubleMatrix invertJ = (DoubleMatrix)J.Inverse();
            DoubleMatrix dNdx = (DoubleMatrix)gps[i, j, k].GetValue(DataInGausspoint.dNdX);//dNdxi * invertJ;
            for (int kk = 0; kk < nen; kk++)
            {
              B[0, 3 * kk] = dNdx[kk, 0];
              B[1, 3 * kk + 1] = dNdx[kk, 1];
              B[2, 3 * kk + 2] = dNdx[kk, 2];
              B[3, 3 * kk + 1] = dNdx[kk, 2];
              B[3, 3 * kk + 2] = dNdx[kk, 1];
              B[4, 3 * kk] = dNdx[kk, 2];
              B[4, 3 * kk + 2] = dNdx[kk, 0];
              B[5, 3 * kk] = dNdx[kk, 1];
              B[5, 3 * kk + 1] = dNdx[kk, 0];
            }
            DoubleVector stressS = (DoubleVector)gps[i, j, k].GetValue(DataInGausspoint.currentStress);//gps[i, j, k].currentStress; // stress

            double detJ = Math.Abs(MatrixFunctions.Determinant(J));
            fi += gps[i, j, k].weight * detJbar * detJ * MatrixFunctions.Product(B.Transpose(), stressS);
            // Notice:
            // In the iterative solver (i.e SolverLoadControlled or SolverArclengthControlled), when it is converged, the HistoryValuesLast have to be updated by: HistoryValuesLast = HistoryValuesCurrent
          }
        }
      }
    }
    public double[,] Transformtable(double[,] table)
    {
      ///// ham nay dung cho bai toan multiliner////
      int legh = table.Length / 2;
      if (legh < 2)
      {
        throw new NotImplementedException("phai nhap nhieu hon hai diem");
      }
      //strain[1] = 0.0513;
      double Emodun = table[1, 1] / table[0, 1];
      double sigma_offset = table[1, 1];
      double[,] str_tressPlas = new double[2, legh];
      str_tressPlas[0, 0] = 0.0;
      str_tressPlas[1, 0] = sigma_offset;
      //  str_tressPlas[1, 0] = sigma_offset;

      for (int i = 1; i < legh; i++)
      {
        str_tressPlas[0, i] = table[0, i] - sigma_offset / Emodun;//// chu y ?????
        str_tressPlas[1, i] = table[1, i];
      }

      return str_tressPlas;
    }
    public double[,] Transformtable(double E, double[,] table, double sigma_offset)
    {
      ///// ham nay dung cho bai toan Ramber-Osgood va ham mu////
      int legh = table.Length / 2;
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

    //public double EqvlasticFun( double[,] table1)
    //{
    //    //double result = 0;

    //    double t1 = table1[0, 0];
    //    double t2 = table1[1, 0];
    //    double t3 = table1[2, 0];
    //    double t4 = table1[3, 0];
    //    double result = (Math.Sqrt(2) / 3) * Math.Sqrt((t1 - t2) * (t1 - t2) + (t1 - t3) * (t1 - t3) + (t2 - t3) * (t2 - t3) + (3 / 2) * 4 * t4 *t4);


    //    return result;
    //}


    public double sima_offset(double[,] table, double x_offset)
    {
      double ss1 = table[1, 1];
      double e1 = table[1, 0];

      double x_offstotal = x_offset + ss1;
      double E = ss1 / e1;
      double result = 0;

      int i = 1;
      while (table[0, i] < x_offstotal)
      {
        i++;
      }
      double x1 = table[0, i - 1];
      double y1 = table[1, i - 1];
      double x2 = table[0, i];
      double y2 = table[1, i];

      return result = ((y2 - y1) * (x_offstotal - x1)) / (x2 - x1) + y1;
    }

    //public double [,] Readfile (string ten)
    //{
    //    double[,] result = null;

    //    using (System.IO.StreamReader file = new System.IO.StreamReader(ten, true))
    //    {
    //        string text = file.ReadLine();
    //        int d = 0;
    //        while (text != null)
    //        {
    //            string[] temp = text.Split(new Char[] { ' ' });
    //            if (result == null)
    //            {
    //                result = new double[2, temp.Length];
    //            }
    //            for (int i = 0; i < temp.Length; i++)
    //            {
    //                result[d, i] = double.Parse(temp[i].Trim());
    //            }
    //            d++;

    //            text = file.ReadLine();
    //        }
    //        file.Close();
    //        file.Dispose();  

    //    }
    //    return result;
    //}






    public DoubleMatrix hp { get; set; }
    public DoubleMatrix data2 { get; set; }
    public DoubleMatrix e_s { get; set; }
    public void plmodata(double[,] table)
    {
      int m = 10;
      for (int i = 0; i < table.Length / 2; i++)
      {
        List<double> x = new List<double>();
        List<double> y = new List<double>();
        double n = (Math.Log(table[0, i] / table[0, i + 1]) / Math.Log(table[1, i] / table[1, i + 1]) - 1) / 2;

        double d = (table[0, i] / table[0, i + 1]) / (Math.Pow(table[1, i], (2 * n + 1)) + Math.Pow(table[1, i + 1], (2 * n + 1)));

        double delta = (table[1, i + 1] - table[1, i]) / m;

        y.Add(table[1, i]);
        double temp = table[1, i];
        while (true)
        {
          if (temp == table[1, i + 1]) break;
          temp += delta;
          y.Add(temp);
        }

        double[,] temp1 = new double[y.Count, 1];
        for (int ii = 0; ii < y.Count; ii++)
        {
          temp1[ii, 0] = y[ii];
          x.Add(d * Math.Pow(table[1, i], (2 * n + 1)));
        }
        data2 = new DoubleMatrix(temp1);
        temp1 = new double[2, y.Count];

        for (int ii = 0; ii < x.Count; ii++)
        {
          temp1[0, ii] = x[ii];
          temp1[1, ii] = y[ii];
        }

        e_s = new DoubleMatrix(temp1);
        temp1 = new double[y.Count - 1, 1];
        for (int ii = 0; ii < y.Count - 1; ii++)
        {
          double dx = table[0, ii + 1] - table[0, ii];
          double dy = table[1, ii + 1] - table[1, ii];
          temp1[ii, 0] = Math.Abs(dy / dx);
        }
        hp = new DoubleMatrix(temp1);
      }
    }
    //public double[,] Hp(DoubleMatrix hp, DoubleMatrix data2, DoubleMatrix e_s)
    //{
    //    for (int i = 0; i < hp.ColumnCount; i++)
    //    {
    //        for (int j = 0; j < data2.RowCount-1; j++)
    //        {

    //        }
    //    }
    //}
    public override DoubleVector StrainAt(double[] xi)
    {

      throw new NotImplementedException();
    }

    public override DoubleVector StressAt(double[] xi)
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
