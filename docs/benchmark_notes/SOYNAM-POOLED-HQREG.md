# SoyNAM pooled quantile Lasso

The empirical prediction comparator uses the published `hqreg` 1.4-1 implementation of Yi and Huang's semismooth Newton coordinate descent. Its registered path has 40 penalty ratios between 1 and 0.01, with `alpha=1`, standardisation computed within the corresponding training sample, `eps=1e-8`, and at most 100000 iterations. Prespecified adjustment variables have zero penalty factors. The native intercept is mapped to the eight exhaustive block indicators without changing predictions. The comparator supplies predictions and pinball losses.

The package's quantile routine uses an adaptive Huber continuation internally. Its source selects the tenth percentile of the preceding absolute residuals as the continuation width, monotonically reducing it along the path with a lower bound of 0.001 for the three quantiles used here. For width $g$, its loss is

\[
 L_{\tau,g}(r)=
 \begin{cases}r^2/(4g)+(\tau-1/2)r,&|r|\le g,\\
 |r|/2-g/4+(\tau-1/2)r,&|r|>g.
 \end{cases}
\]

The package applies this continuation to the quantile Lasso objective. The returned solution is evaluated using the original unsmoothed pinball loss on the held-out lines. We retain the native implementation and report its version.

`check_pooled.R` uses 200 actual SoyNAM lines. It checks deterministic repetition, exact intercept mapping, and agreement with an independent CVXR/Clarabel convex formulation at three quantiles and three penalties, including zero. The objective discrepancy must be below $2\times10^{-6}$ in working units. The zero-penalty calculation also records its unsmoothed loss difference from `quantreg::rq.fit.fnb`; the source-derived bound is $g/4$ plus the same numerical tolerance because $0\le\rho_\tau-L_{\tau,g}\le g/4$. These checks keep numerical optimisation error separate from the native continuation.

The three methods share the same complete-line training and test roles. Their calibration uses only internal training validation predictions, separately by year. No test outcome selects a penalty, model size or calibration value.

Sources are the [official function documentation](https://search.r-project.org/CRAN/refmans/hqreg/html/hqreg.html) and the CRAN package source `src/hqreg.c` and `src/data_process.c`. The method reference is Yi and Huang, *Semismooth Newton Coordinate Descent Algorithm for Elastic-Net Penalized Huber Loss Regression and Quantile Regression*, Journal of Computational and Graphical Statistics, DOI 10.1080/10618600.2016.1256816.
