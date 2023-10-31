﻿using System;
using DEMSoft.NURBS;
using DEMSoft.Function;
using CenterSpace.NMath.Core;
using DEMSoft.Common;

namespace DEMSoft.IGA
{
    public enum TypeBeam
    { FSDT, HSDT }
    public class PatchStructureBeam : AbstractPatch1D
    {
        public double Getzz
        { get; set; }
        public TypeVFunction GetTypeVF
        { get; set; }
        /// <summary>
        /// Kinematics function in HSDT plate
        /// </summary>
        public FunctionRToR KinematicsFunction
        { get; set; }
        public TypeBeam TypeBeam
        { get; set; }

        public Structure1DState StateStress
        {
            get
            {
                throw new NotImplementedException();
            }

            set
            {
                throw new NotImplementedException();
            }
        }

        public PatchStructureBeam(NURBSCurve curve, double thickness, TypeBeam typeBeam)
                                        : base(curve)
        {
            countDimension = 1;
            this.TypeBeam = typeBeam;
            Thickness = thickness;
        }

        public void ComputeInternalForcePatch(out DoubleVector residualGlobal)
        {
            throw new NotImplementedException();
        }

        public double StressAt(Result re, params double[] xi)
        {
            double disp = 0;
            int numElem = FindIndexOfElementAt(xi[0], xi[1]);
            AbstractElementStructure1D elem = (AbstractElementStructure1D)listElement[numElem];
            DoubleVector stress = elem.StressAt(xi[0], xi[1]);
            switch (re)
            {
                case Result.SIGMAXX:
                    disp = stress[0];
                    break;
                case Result.SIGMAYY:
                    disp = stress[1];
                    break;
                case Result.SIGMAZZ:
                    disp = 0;
                    break;
                case Result.SIGMAXY:
                    disp = stress[2];
                    break;
                case Result.SIGMAXZ:
                    disp = stress[3];
                    break;
                case Result.SIGMAYZ:
                    disp = stress[4];
                    break;
                case Result.SIGMAEQV:
                    double sxx = stress[0], syy = stress[1], szz = 0, sxy = stress[2], sxz = stress[3], syz = stress[4];
                    disp = Math.Sqrt(1.0 / 2.0 * (Math.Pow(sxx - syy, 2) + Math.Pow(syy - szz, 2) + Math.Pow(sxx - szz, 2) + 6.0 * (Math.Pow(sxy, 2) + Math.Pow(sxz, 2) + Math.Pow(syz, 2))));
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
                    double sXZ = StressAt(Result.SIGMAXZ, param[0], param[1]);
                    double sYZ = StressAt(Result.SIGMAYZ, param[0], param[1]);
                    DoubleMatrix stressXY = new DoubleMatrix(3, 3);
                    DoubleMatrix T = new DoubleMatrix(3, 3);
                    stressXY[0, 0] = sXX;
                    stressXY[1, 0] = stressXY[0, 1] = sXY;
                    stressXY[1, 1] = sYY;
                    stressXY[2, 0] = stressXY[0, 2] = sXZ;
                    stressXY[2, 1] = stressXY[1, 2] = sYZ;
                    stressXY[2, 2] = sZZ;
                    T[0, 0] = T[1, 1] = Math.Cos(x[1]);
                    T[0, 1] = -Math.Sin(x[1]);
                    T[1, 0] = Math.Sin(x[1]);
                    T[2, 2] = 1;
                    DoubleMatrix stressRT = T.Transpose() * stressXY * T;
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
            AbstractElementStructureBeam elem = (AbstractElementStructureBeam)listElement[numElem];
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
                    disp = 0;
                    break;
                case Result.EPSILONXY:
                    disp = strain[2];
                    break;
                case Result.EPSILONXZ:
                    disp = strain[3];
                    break;
                case Result.EPSILONYZ:
                    disp = strain[4];
                    break;
                case Result.EPSILONEQV:
                    double sxx = strain[0], syy = strain[1], sxy = strain[2], sxz = strain[3], syz = strain[4], szz = 0;
                    double xx = 2.0 / 3.0 * sxx - 1.0 / 3.0 * syy - 1.0 / 3.0 * szz;
                    double yy = -1.0 / 3.0 * sxx + 2.0 / 3.0 * syy - 1.0 / 3.0 * szz;
                    double zz = -1.0 / 3.0 * sxx - 1.0 / 3.0 * syy + 2.0 / 3.0 * szz;
                    disp = 2.0 / 3.0 * Math.Sqrt(3.0 / 2.0 * (Math.Pow(xx, 2) + Math.Pow(yy, 2) + Math.Pow(zz, 2)) + 3.0 / 4.0 * (Math.Pow(2 * sxy, 2) + Math.Pow(2 * sxz, 2) + Math.Pow(2 * syz, 2)));
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
                    double sXZ = StrainAt(Result.EPSILONXZ, param[0], param[1]);
                    double sYZ = StrainAt(Result.EPSILONYZ, param[0], param[1]);
                    DoubleMatrix strainXY = new DoubleMatrix(3, 3);
                    DoubleMatrix T = new DoubleMatrix(3, 3);
                    strainXY[0, 0] = sXX;
                    strainXY[1, 0] = strainXY[0, 1] = sXY;
                    strainXY[1, 1] = sYY;
                    strainXY[2, 0] = strainXY[0, 2] = sXZ;
                    strainXY[2, 1] = strainXY[1, 2] = sYZ;
                    strainXY[2, 2] = sZZ;/////////////////////////////// chua hieu
                    T[0, 0] = T[1, 1] = Math.Cos(x[1]);
                    T[0, 1] = -Math.Sin(x[1]);
                    T[1, 0] = Math.Sin(x[1]);
                    T[2, 2] = 1;
                    DoubleMatrix strainRT = T.Transpose() * strainXY * T;
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
