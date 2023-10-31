using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DEMSoft.IGA
{
  public class Convergence
  {
    /// <summary>
    /// Criterion convergence
    /// </summary>
    public CriterionConvergence criterionConvergence;
    /// <summary>
    /// Type of norm
    /// </summary>
    public Norm norm;
    /// <summary>
    /// Maximum of iterator
    /// </summary>
    public int maximumIteration;
    /// <summary>
    /// Tolorance
    /// </summary>
    public double TOL;
    public bool IsSkipUnconvergenceStep;
    public Convergence(CriterionConvergence criterionConvergence, int maximumIteration, double TOL, Norm norm)
    {
      this.criterionConvergence = criterionConvergence;
      this.TOL = TOL;
      this.norm = norm;
      this.maximumIteration = maximumIteration;
    }
    public Convergence()
    {
      this.criterionConvergence = CriterionConvergence.ResidualCriterion;
      this.TOL = 10e-4;
      this.norm = Norm.L2Norm;
      this.maximumIteration = 100;
    }
  }
  /// <summary>
  /// Convergence Criteria
  /// </summary>
  /// <see cref="https://dianafea.com/manuals/d944/Analys/node395.html"/>
  public enum CriterionConvergence
  {
    /// <summary>
    /// Displacement criterion
    /// </summary>
    /// <remarks>The displacement norm is the Euclidian norm of the iterative displacement increment. To check convergence, the displacement norm is checked against the norm of the displacement increments in the first prediction of the increment.
    /// Displacement norm ratio = $\displaystyle {\frac{{ \sqrt{ \delta\mathbf{u}_{i}^{\mathrm{\scriptscriptstyle{lta\mathbf{u}_{0}^{\mathrm{\scriptscriptstyle{T}}}\; \Delta\mathbf{u}_{0} } }}}$	(31.17)
    /// From(31.17) it is clear that the ratio of the displacement norm after the first prediction(iteration 0) equals 1 by definition.To check convergence, always one additional iteration is necessary.
    /// </remarks>
    DisplacementCriterion,
    /// <summary>
    /// Force criterion
    /// </summary>
    /// <remarks>
    /// The force norm is the Euclidian norm of the out-of-balance force vector g . To check convergence, the force norm after the current iteration is checked against the norm of the initial unbalance g0
    /// Force norm ratio = $\displaystyle {\frac{{ \sqrt{ \mathbf{g}_{i}^{\mathrm{\scriptscriptstyle{T}}}\;{ \sqrt{ \mathbf{g}_{0}^{\mathrm{\scriptscriptstyle{T}}}\; \mathbf{g}_{0} } }}}$	(31.16)
    /// Because the reference force norm is known before the first prediction of displacements, the force norm ratio can be calculated directly after the first prediction, i = 1 in (31.16). This means that if the first prediction is correct(nearly linear behavior) the force norm can detect convergence right away and no unnecessary iterations have to be performed.
    /// </remarks>
    ForceCriterion,
    /// <summary>
    /// Energy criterion
    /// </summary>
    /// <remarks>
    /// A third way to check convergence is the energy norm. This norm is composed of internal forces and relative displacements as indicated in Figure 31.8 with $ \Delta$E0 and $ \delta$E1 . To determine convergence, the energy ratio is calculated as
    /// Energy norm ratio = $\displaystyle \left\vert\vphantom{ \frac{ \delta\mathbf{u}_{i}^{\mathrm{\scriptT}}}\: ( \mathbf{f}_{\mathrm{int},1} + \mathbf{f}_{\mathrm{int},0} ) } }\right.$$\displaystyle {\frac{{ \delta\mathbf{u}_{i}^{\mathrm{\scriptscriptstyle{T}}}\: tstyle{T}}}\: ( \mathbf{f}_{\mathrm{int},1} + \mathbf{f}_{\mathrm{int},0} ) }}}$$\displaystyle \left.\vphantom{ \frac{ \delta\mathbf{u}_{i}^{\mathrm{\scriptscri\: ( \mathbf{f}_{\mathrm{int},1} + \mathbf{f}_{\mathrm{int},0} ) } }\right\vert$	(31.18)
    /// Note that here the internal force is used and not the out-of-balance force.Use of the out-of-balance force would be improper, for a Line Search procedure could then minimize the norm, see Equation (31.15), before the solution really converges to equilibrium.As with the displacement norm, the energy norm also requires an additional iteration to detect convergence.
    /// The choice of the proper norm and its convergence criterion depends on the type of analysis.Using a lot of prescribed displacements generally makes the displacement norm less useful. On the other hand, a structure that can expand freely will hardly build up any internal forces and the force norm may be less useful.Always be sure that the reference norm(the denominator in the ratios) has a reasonable value i.e., not close to zero.
    /// Experience shows that the convergence criterion for softening type behavior should be more strict than the criterion that can be used in a hardening type analysis. If there is any doubt about the criterion to be used, it is advisable to perform the analysis with two distinct criteria and check the differences in results.If large differences occur, at least the less strict norm was to large.
    /// </remarks>
    EnergyCriterion,
    /// <summary>
    /// Residual criterion
    /// </summary>
    /// <remarks>
    /// The residual norm is also a Euclidian norm of the out-of-balance force vector g . Contraray to the force norm, the residual norm also takes the values in constrained degrees of freedom (supports and tyings) into account. To check convergence, DIANA compares the change in the residual norm during the current iteration with the change in the residual norm during the first prediction of displacements in the current step.
    /// Residual norm ratio = $\displaystyle {\frac{{ \left\vert \sqrt{ \mathbf{g}_{i}^{\mathrm{\scriptscriptsthbf{g}_{n}^{\mathrm{\scriptscriptstyle{T}}}\mathbf{g}_{n} } \, \right\vert }}}$	(31.19)
    /// Where gn denotes the out-of-balance force vector in the last iteration of the previous step.In the first step DIANA takes its value as zero.
    /// </remarks>
    ResidualCriterion,
    /// <summary>
    /// Displacement criterion Nam Ho-Kim (2.13)
    /// </summary>
    /// <see cref="Introduction to Nonlinear Finite Element Analysis"/>
    DisplacementCriterion2,
    /// <summary>
    /// Displacement criterion Nam Ho-Kim (2.13) sqrt
    /// </summary>
    /// <see cref="Introduction to Nonlinear Finite Element Analysis"/>
    DisplacementCriterion22,
    /// <summary>
    /// Residual criterion Nam Ho-Kim (2.12)
    /// </summary>
    /// <see cref="Introduction to Nonlinear Finite Element Analysis"/>
    ResidualCriterion2,
    /// <summary>
    /// Displacement criterion (6.3.63)
    /// </summary>
    /// <see cref="Ted Belytschko"/>
    DisplacementCriterion3,
    /// <summary>
    /// Displacement criterion (6.3.63) remove constrain value
    /// </summary>
    /// <see cref="Ted Belytschko"/>
    DisplacementCriterion4,
    /// <summary>
    /// Displacement criterion not normalize
    /// </summary>
    /// <see cref="Wu"/>
    DisplacementCriterion5,
    /// <summary>
    /// Energy criterion (6.3.64)
    /// </summary>
    /// <see cref="Ted Belytschko"/>
    EnergyCriterion3,
    /// <summary>
    /// Residual criterion (6.3.62)
    /// </summary>
    /// <see cref="Ted Belytschko"/>
    ResidualCriterion3,
  }

  public enum Norm
  {
    InfiniteNorm,
    L1Norm,
    L2Norm
  }
}
