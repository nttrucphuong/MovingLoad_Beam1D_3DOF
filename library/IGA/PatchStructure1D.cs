using System;
using DEMSoft.NURBS;
using DEMSoft.EngineeringData;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
    public class PatchStructure1D : AbstractPatch1D
    {
        public Structure1DState StateStress
        { get; set; }

        public PatchStructure1D(NURBSCurve curve, Structure1DState state)
                    : base(curve)
        {
            StateStress = state;
        }
        public double StressAt(Result re, params double[] xi)
        {
            double disp = 0;
            int numElem = FindIndexOfElementAt(xi[0], xi[1]);
            AbstractElementStructure1D elem = (AbstractElementStructure1D)listElement[numElem];
            DoubleVector stress = elem.StressAt(xi[0], xi[1]);
            double nu = elem.Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
            switch (re)
            {
                case Result.SIGMAXX:
                    disp = stress[0];
                    break;
                case Result.SIGMAYY:
                    disp = stress[1];
                    break;
                case Result.SIGMAXY:
                    disp = stress[3];
                    break;
                case Result.SIGMAEQV:
                    double sxx = stress[0], syy = stress[1], szz = stress[2], sxy = stress[3];
                    if (StateStress == Structure1DState.Beam)
                    {
                        //double nu = elem.Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
                        szz = nu * (sxx + syy);
                    }
                    disp = Math.Sqrt(1.0 / 2.0 * (Math.Pow(sxx - syy, 2) + Math.Pow(syy - szz, 2) + Math.Pow(sxx - szz, 2) + 6.0 * Math.Pow(sxy, 2)));
                    break;
            }
            return disp;
        }
        public double StressOnGlobalAt(Result re, TypeCoordination coord, params double[] x)
        {
            double[] param = null;
            double stress = 0;
            switch (coord)
            {
                case TypeCoordination.Cartesian:
                    param = ((NURBSCurve)geometry[0]).Projection(x[0], x[1], 0);
                    stress = StressAt(re, param[0], param[1]);
                    break;
                case TypeCoordination.Cylindrical:
                    double X = x[0] * Math.Cos(x[1]);
                    double Y = x[0] * Math.Sin(x[1]);
                    param = ((NURBSCurve)geometry[0]).Projection(X, Y, 0);
                    double sXX = StressAt(Result.SIGMAXX, param[0], param[1]);
                    double sYY = StressAt(Result.SIGMAYY, param[0], param[1]);
                    double sZZ = StressAt(Result.SIGMAZZ, param[0], param[1]);
                    double sXY = StressAt(Result.SIGMAXY, param[0], param[1]);
                    DoubleMatrix stressXY = new DoubleMatrix(3, 3);
                    DoubleMatrix T = new DoubleMatrix(3, 3);
                    stressXY[0, 0] = sXX;
                    stressXY[1, 0] = stressXY[0, 1] = sXY;
                    stressXY[1, 1] = sYY;
                    stressXY[2, 2] = sZZ;
                    T[0, 0] = T[1, 1] = Math.Cos(x[1]);
                    T[0, 1] = -Math.Sin(x[1]);
                    T[1, 0] = Math.Sin(x[1]);
                    T[2, 2] = 1;
                    DoubleMatrix stressRT = NMathFunctions.TransposeProduct(T, NMathFunctions.Product(stressXY, T));
                    switch (re)
                    {
                        case Result.SIGMARR:
                            stress = stressRT[0, 0];
                            break;
                        case Result.SIGMATT:
                            stress = stressRT[1, 1];
                            break;
                        case Result.SIGMAZZ:
                            stress = stressRT[2, 2];
                            break;
                        case Result.SIGMART:
                            stress = stressRT[0, 1];
                            break;
                    }
                    break;
            }
            return stress;
        }
        public double StrainAt(Result re, params double[] xi)
        {
            double disp = 0;
            int numElem = FindIndexOfElementAt(xi[0], xi[1]);
            AbstractElementStructure1D elem = (AbstractElementStructure1D)listElement[numElem];
            DoubleVector strain = elem.StrainAt(xi[0], xi[1]);
            switch (re)
            {
                case Result.EPSILONXX:
                    disp = strain[0];
                    break;
                case Result.EPSILONYY:
                    disp = strain[1];
                    break;
                case Result.EPSILONZZ:
                    disp = strain[2];
                    if (StateStress == Structure1DState.Beam)
                    {
                        double nu = elem.Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
                        disp = -nu * (1.0 - nu) * (strain[0] + strain[1]);
                    }
                    break;
                case Result.EPSILONXY:
                    disp = strain[3];
                    break;
                case Result.EPSILONEQV:
                    double sxx = strain[0], syy = strain[1], szz = strain[2], sxy = strain[3];
                    if (StateStress == Structure1DState.Beam)
                    {
                        double nu = elem.Material.GetProperty(MaterialPropertyName.PoissonRatio).GetValueProperty();
                        szz = -nu * (1.0 - nu) * (sxx + syy);
                    }
                    double xx = 2.0 / 3.0 * sxx - 1.0 / 3.0 * syy - 1.0 / 3.0 * szz;
                    double yy = -1.0 / 3.0 * sxx + 2.0 / 3.0 * syy - 1.0 / 3.0 * szz;
                    double zz = -1.0 / 3.0 * sxx - 1.0 / 3.0 * syy + 2.0 / 3.0 * szz;
                    disp = 2.0 / 3.0 * Math.Sqrt(3.0 / 2.0 * (Math.Pow(xx, 2) + Math.Pow(yy, 2) + Math.Pow(zz, 2)) + 3.0 / 4.0 * (Math.Pow(2 * sxy, 2)));
                    break;
            }
            return disp;
        }

        public double StrainOnGlobalAt(Result re, TypeCoordination coord, params double[] x)
        {
            double[] param = null;
            double strain = 0;
            switch (coord)
            {
                case TypeCoordination.Cartesian:
                    param = GetCurve().Projection(x[0], x[1], 0);
                    strain = StrainAt(re, param[0], param[1]);
                    break;
                case TypeCoordination.Cylindrical:
                    double X = x[0] * Math.Cos(x[1]);
                    double Y = x[0] * Math.Sin(x[1]);
                    param = GetCurve().Projection(X, Y, 0);
                    double sXX = StrainAt(Result.EPSILONXX, param[0], param[1]);
                    double sYY = StrainAt(Result.EPSILONYY, param[0], param[1]);
                    double sZZ = StrainAt(Result.EPSILONZZ, param[0], param[1]);
                    double sXY = StrainAt(Result.EPSILONXY, param[0], param[1]);
                    DoubleMatrix strainXY = new DoubleMatrix(3, 3);
                    DoubleMatrix T = new DoubleMatrix(3, 3);
                    strainXY[0, 0] = sXX;
                    strainXY[1, 0] = strainXY[0, 1] = sXY;
                    strainXY[1, 1] = sYY;
                    strainXY[2, 2] = sZZ;
                    T[0, 0] = T[1, 1] = Math.Cos(x[1]);
                    T[0, 1] = -Math.Sin(x[1]);
                    T[1, 0] = Math.Sin(x[1]);
                    T[2, 2] = 1;
                    DoubleMatrix strainRT = NMathFunctions.TransposeProduct(T, NMathFunctions.Product(strainXY, T));
                    switch (re)
                    {
                        case Result.SIGMARR:
                            strain = strainRT[0, 0];
                            break;
                        case Result.SIGMATT:
                            strain = strainRT[1, 1];
                            break;
                        case Result.SIGMAZZ:
                            strain = strainRT[2, 2];
                            break;
                        case Result.SIGMART:
                            strain = strainRT[0, 1];
                            break;
                    }
                    break;
            }
            return strain;
        }
    }
}
