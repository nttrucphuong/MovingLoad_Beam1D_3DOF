using DEMSoft.NURBS;
using CenterSpace.NMath.Core;

namespace DEMSoft.IGA
{
    public class PatchThermoelastic3D : AbstractPatch3D
    {
        public PatchThermoelastic3D(NURBSVolume volume)
                                   : base(volume, 4)
        {

        }


        /////////////////// Doi Luu
        public void ComputeHeatTransferConvectionPatch(ref DoubleVector rGlobal)
        {
            int nel = CountElements();//(n - p) * (m - q);//number of elements
                                        //if (!IsParallelProcesing)
                                        //{
            for (int i = 0; i < nel; i++)//Loop through elements
            {
                ElementThermoelastic3D elem = (ElementThermoelastic3D)listElement[i];
                if (elem.GetHeatTransferConvection() != null)
                {
                    elem.ComputeHeatTransferConvectionLoadVectorElement(ref rGlobal);
                }
            }
        }
        //////////////////////




        //public double[] GetDisplacementAt(double xi, double eta, double zeta)
        //{
        //    double[] disp = new double[3];
        //    var volume = GetVolume();
        //    var cps = volume.ControlPoints;
        //    var basis = (TrivariateNURBSBasisFunction)volume.Basis;
        //    int p = basis.GetDegree(0);
        //    int uSpan = basis.GetKnotVector(0).FindSpan(xi, p);
        //    int q = basis.GetDegree(1);
        //    int vSpan = basis.GetKnotVector(1).FindSpan(eta, q);
        //    int r = basis.GetDegree(2);
        //    int wSpan = basis.GetKnotVector(2).FindSpan(zeta, r);

        //    double[,,] Nuv = basis.GetValueTrivariateBasisFunctions(xi, eta, zeta);

        //    for (int i = 0; i <= p; i++)
        //        for (int j = 0; j <= q; j++)
        //            for (int k = 0; k <= r; k++)
        //            {

        //                disp[0] += Nuv[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.UX);
        //                disp[1] += Nuv[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.UY);
        //                disp[2] += Nuv[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.UZ);
        //            }
        //    return disp;
        //}

        //public double GetTempAt(double xi, double eta, double zeta)
        //{
        //    return GetApproximateAt(Result.TEMP, xi, eta, zeta);
        //}
        //public double GetTempOnGlobalAt(double x, double y, double z)
        //{
        //    return GetApproximateOnGlobalAt(Result.TEMP, x, y, z);
        //}
        //public double[] GetDisplacementOnGlobalAt(double x, double y, double z)
        //{
        //    double[] param = ((NURBSVolume)geometry[0]).Projection(x, y, z);
        //    return GetDisplacementAt(param[0], param[1], param[2]);
        //}


        //public double GetStressAt(double xi, double eta, double zeta, Result re, bool isFromGaussPoints)
        //{
        //    double disp = 0;
        //    if (isFromGaussPoints)
        //    {
        //        var vol = GetVolume();
        //        var cps = vol.ControlPoints;
        //        var basis = (TrivariateNURBSBasisFunction)vol.Basis;
        //        int p = basis.GetDegree(0);
        //        int uSpan = basis.GetKnotVector(0).FindSpan(xi, p);
        //        int q = basis.GetDegree(1);
        //        int vSpan = basis.GetKnotVector(1).FindSpan(eta, q);
        //        int r = basis.GetDegree(2);
        //        int wSpan = basis.GetKnotVector(2).FindSpan(zeta, r);
        //        double[,,] Nuvw = basis.GetValueTrivariateBasisFunctions(xi, eta, zeta);
        //        double sxx = 0, syy = 0, sxy = 0, szz = 0, syz = 0, sxz = 0;
        //        for (int i = 0; i <= p; i++)
        //            for (int j = 0; j <= q; j++)
        //                for (int k = 0; k <= r; k++)
        //                {
        //                    if (re != Result.SIGMAEQV)
        //                        disp += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(re);
        //                    else
        //                    {
        //                        sxx += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.SIGMAXX);
        //                        syy += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.SIGMAYY);
        //                        szz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.SIGMAZZ);
        //                        sxy += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.SIGMAXY);
        //                        syz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.SIGMAYZ);
        //                        sxz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.SIGMAXZ);
        //                    }
        //                }
        //        if (re == Result.SIGMAEQV)
        //            disp = Math.Sqrt(1.0 / 2.0 * (Math.Pow(sxx - syy, 2) + Math.Pow(syy - szz, 2) + Math.Pow(sxx - szz, 2) + 6.0 * (Math.Pow(sxy, 2) + Math.Pow(syz, 2) + Math.Pow(sxz, 2))));
        //    }
        //    //else
        //    //{
        //    //    int numElem = FindIndexOfElementAt(xi, eta, zeta);
        //    //    AbstractElementStructure3D elem = (AbstractElementStructure3D)listElement[numElem];
        //    //    DoubleVector stress = elem.StressAt(xi, eta, zeta);
        //    //    switch (re)
        //    //    {
        //    //        case Result.SIGMAXX:
        //    //            disp = stress[0];
        //    //            break;
        //    //        case Result.SIGMAYY:
        //    //            disp = stress[1];
        //    //            break;
        //    //        case Result.SIGMAZZ:
        //    //            disp = stress[2];
        //    //            break;
        //    //        case Result.SIGMAXY:
        //    //            disp = stress[3];
        //    //            break;
        //    //        case Result.SIGMAYZ:
        //    //            disp = stress[4];
        //    //            break;
        //    //        case Result.EPSILONXZ:
        //    //            disp = stress[5];
        //    //            break;
        //    //        case Result.SIGMAEQV:
        //    //            double sxx = stress[0], syy = stress[1], sxy = stress[2], szz = 0;

        //    //                DoubleVector strain = elem.StrainAt(xi, eta);
        //    //                double E = elem.Material.GetProperty(MaterialPropertyName.YoungModulus).GetValueProperty();
        //    //                double nu = elem.Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();

        //    //                szz = E / (1.0 + nu) * nu / (1.0 - 2.0 * nu) * (strain[0] + strain[1]);
        //    //            disp = Math.Sqrt(1.0 / 2.0 * (Math.Pow(sxx - syy, 2) + Math.Pow(syy - szz, 2) + Math.Pow(sxx - szz, 2) + 6.0 * Math.Pow(sxy, 2)));
        //    //            break;
        //    //    }
        //    //}
        //    return disp;
        //    //int numElem = FindIndexOfElementAt(xi, eta, zeta);
        //    //return ((AbstractStructure3DElement)(listElement[numElem])).StressAt(xi, eta, zeta);
        //}

        //public double GetStressOnGlobalAt(double x, double y, double z, Result re, bool isFromGaussPoints)
        //{
        //    double[] param = ((NURBSVolume)geometry[0]).Projection(x, y, z);
        //    return GetStressAt(param[0], param[1], param[2], re, isFromGaussPoints);
        //}

        //public double GetStrainAt(double xi, double eta, double zeta, Result re)
        //{
        //    double disp = 0;
        //    var vol = GetVolume();
        //    var cps = vol.ControlPoints;
        //    var basis = (TrivariateNURBSBasisFunction)vol.Basis;
        //    int p = basis.GetDegree(0);
        //    int uSpan = basis.GetKnotVector(0).FindSpan(xi, p);
        //    int q = basis.GetDegree(1);
        //    int vSpan = basis.GetKnotVector(1).FindSpan(eta, q);
        //    int r = basis.GetDegree(2);
        //    int wSpan = basis.GetKnotVector(2).FindSpan(zeta, r);
        //    double[,,] Nuvw = basis.GetValueTrivariateBasisFunctions(xi, eta, zeta);
        //    double sxx = 0, syy = 0, sxy = 0, szz = 0, syz = 0, sxz = 0;
        //    for (int i = 0; i <= p; i++)
        //        for (int j = 0; j <= q; j++)
        //            for (int k = 0; k <= r; k++)
        //            {
        //                if (re != Result.EPSILONEQV)
        //                    disp += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(re);
        //                else
        //                {
        //                    sxx += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.EPSILONXX);
        //                    syy += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.EPSILONYY);
        //                    szz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.EPSILONZZ);
        //                    sxy += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.EPSILONXY);
        //                    syz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.EPSILONYZ);
        //                    sxz += Nuvw[i, j, k] * cps[uSpan - p + i, vSpan - q + j, wSpan - r + k].GetResult(Result.EPSILONXZ);
        //                }
        //            }
        //    if (re == Result.EPSILONEQV)
        //    {
        //        double xx = 2.0 / 3.0 * sxx - 1.0 / 3.0 * syy - 1.0 / 3.0 * szz;
        //        double yy = -1.0 / 3.0 * sxx + 2.0 / 3.0 * syy - 1.0 / 3.0 * szz;
        //        double zz = -1.0 / 3.0 * sxx - 1.0 / 3.0 * syy + 2.0 / 3.0 * szz;

        //        disp = 2.0 / 3.0 * Math.Sqrt(3.0 / 2.0 * (Math.Pow(xx, 2) + Math.Pow(yy, 2) + Math.Pow(zz, 2)) + 3.0 / 4.0 * (Math.Pow(2 * sxy, 2) + Math.Pow(2 * syz, 2) + Math.Pow(2 * sxz, 2)));
        //    }
        //    return disp;
        //    //int numElem = FindIndexOfElementAt(xi, eta, zeta);
        //    //return ((AbstractStructure3DElement)(listElement[numElem])).StrainAt(xi, eta, zeta);
        //}



        ////public double[,,] CalculateResult(Result re, int resolution1, int resolution2, int resolution3)
        ////{
        ////    var volume = GetVolume();
        ////    var data = volume.GetParametricOnVolume(resolution1, resolution2, resolution3);
        ////    double[,,] val = new double[data[0].Length, data[1].Length, data[2].Length];
        ////    switch (re)
        ////    {
        ////        case Result.UX:
        ////        case Result.UY:
        ////        case Result.UZ:
        ////            if (!IsParallelProcesing)
        ////                for (int i = 0; i < data[0].Length; i++)
        ////                    for (int j = 0; j < data[1].Length; j++)
        ////                        for (int k = 0; k < data[2].Length; k++)
        ////                            val[i, j, k] = GetApproximateAt(re, data[0][i], data[1][j], data[2][k]);
        ////            else
        ////                Parallel.For(0, data[0].Length, i =>
        ////                {
        ////                    for (int j = 0; j < data[1].Length; j++)
        ////                        for (int k = 0; k < data[2].Length; k++)
        ////                            val[i, j, k] = GetApproximateAt(re, data[0][i], data[1][j], data[2][k]);
        ////                });
        ////            break;
        ////        case Result.USUM:
        ////            if (!IsParallelProcesing)
        ////                for (int k = 0; k < data[2].Length; k++)
        ////                    for (int j = 0; j < data[1].Length; j++)
        ////                        for (int i = 0; i < data[0].Length; i++)

        ////                        {
        ////                            double[] temp = new double[3];
        ////                            temp[0] = GetApproximateAt(Result.UX, data[0][i], data[1][j], data[2][k]);
        ////                            temp[1] = GetApproximateAt(Result.UY, data[0][i], data[1][j], data[2][k]);
        ////                            temp[2] = GetApproximateAt(Result.UZ, data[0][i], data[1][j], data[2][k]);
        ////                            val[i, j, k] = Math.Sqrt(temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2]);
        ////                        }
        ////            else
        ////                Parallel.For(0, data[0].Length, i =>
        ////                {
        ////                    for (int j = 0; j < data[1].Length; j++)
        ////                        for (int k = 0; k < data[2].Length; k++)
        ////                        {
        ////                            double[] temp = new double[3];
        ////                            temp[0] = GetApproximateAt(Result.UX, data[0][i], data[1][j], data[2][k]);
        ////                            temp[1] = GetApproximateAt(Result.UY, data[0][i], data[1][j], data[2][k]);
        ////                            temp[2] = GetApproximateAt(Result.UZ, data[0][i], data[1][j], data[2][k]);
        ////                            val[i, j, k] = Math.Sqrt(temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2]);
        ////                        }
        ////                });
        ////            break;
        ////        case Result.SIGMAXX:
        ////        case Result.SIGMAYY:
        ////        case Result.SIGMAZZ:
        ////        case Result.SIGMAXY:
        ////        case Result.SIGMAYZ:
        ////        case Result.SIGMAXZ:
        ////        case Result.SIGMAEQV:
        ////            if (!IsParallelProcesing)
        ////                for (int i = 0; i < data[0].Length; i++)
        ////                    for (int j = 0; j < data[1].Length; j++)
        ////                        for (int k = 0; k < data[2].Length; k++)
        ////                            val[i, j, k] = GetStressAt(data[0][i], data[1][j], data[2][k], re, true);
        ////            else
        ////                Parallel.For(0, data[0].Length, i =>
        ////                {
        ////                    for (int j = 0; j < data[1].Length; j++)
        ////                        for (int k = 0; k < data[2].Length; k++)
        ////                            val[i, j, k] = GetStressAt(data[0][i], data[1][j], data[2][k], re, true);
        ////                });
        ////            break;
        ////        case Result.EPSILONXX:
        ////        case Result.EPSILONYY:
        ////        case Result.EPSILONZZ:
        ////        case Result.EPSILONXY:
        ////        case Result.EPSILONYZ:
        ////        case Result.EPSILONXZ:
        ////        case Result.EPSILONEQV:
        ////            if (!IsParallelProcesing)
        ////                for (int i = 0; i < data[0].Length; i++)
        ////                    for (int j = 0; j < data[1].Length; j++)
        ////                        for (int k = 0; k < data[2].Length; k++)
        ////                            val[i, j, k] = GetStrainAt(data[0][i], data[1][j], data[2][k], re);
        ////            else
        ////                Parallel.For(0, data[0].Length, i =>
        ////                {
        ////                    for (int j = 0; j < data[1].Length; j++)
        ////                        for (int k = 0; k < data[2].Length; k++)
        ////                            val[i, j, k] = GetStrainAt(data[0][i], data[1][j], data[2][k], re);
        ////                });
        ////            break;
        ////        case Result.TEMP:
        ////            for (int i = 0; i < data[0].Length; i++)
        ////                for (int j = 0; j < data[1].Length; j++)
        ////                    for (int k = 0; k < data[2].Length; k++)
        ////                    {
        ////                        val[i, j, k] = GetTempAt(data[0][i], data[1][j], data[2][k]);
        ////                    }
        ////            break;
        ////    }
        ////    return val;
        ////}



        ///// <summary>
        ///// Get deformation surface to draw result with scale factor
        ///// </summary>
        ///// <param name="scale">scale factor</param>
        ///// <returns></returns>
        //public NURBSVolume GetDeformationVolume(double scale)
        //{
        //    var volumeCopy = (NURBSVolume)GetVolume().Clone();
        //    var cpsCopy = volumeCopy.ControlPoints;
        //    if (!IsParallelProcesing)
        //        for (int i = 0; i < cpsCopy.GetLength(0); i++)
        //            for (int j = 0; j < cpsCopy.GetLength(1); j++)
        //                for (int k = 0; k < cpsCopy.GetLength(2); k++)
        //                {
        //                    for (int kk = 0; kk < 3; kk++)
        //                        if (kk < GetCountField())
        //                            cpsCopy[i, j, k].SetCoordinate(kk, cpsCopy[i, j, k].GetCoordinate(kk) + scale * cpsCopy[i, j, k].GetULocal(kk));
        //                }
        //    else
        //        Parallel.For(0, cpsCopy.GetLength(0), i =>
        //        {
        //            for (int j = 0; j < cpsCopy.GetLength(1); j++)
        //                for (int k = 0; k < cpsCopy.GetLength(2); k++)
        //                {
        //                    for (int kk = 0; kk < 3; kk++)
        //                        if (kk < GetCountField())
        //                            cpsCopy[i, j, k].SetCoordinate(kk, cpsCopy[i, j, k].GetCoordinate(kk) + scale * cpsCopy[i, j, k].GetULocal(kk));
        //                }
        //        });
        //    return volumeCopy;
        //}
    }
}
