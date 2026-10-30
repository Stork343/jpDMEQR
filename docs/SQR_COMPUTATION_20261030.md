# Certified computation for SQR

The implementation evaluates the registered SQR estimator with its full
gamma grid and nested subject folds. The coefficient objective, common
CLIME tolerance and minimum-magnitude symmetrisation retain their declared
definitions. Interval construction and all simulation seeds are unchanged.
The numerical specification is
`config/submission_comparator_numerics_20261030.yml`.

## CLIME solution certificates

The compiled driver calls the HiGHS C++ library. It solves restricted
linear programmes, adds violated rows and dual-priced columns, and
verifies the resulting direction against the complete Hessian. For a
direction omega and dual vector y, the checks are

\[
 \|H\omega-e_k\|_\infty\leq\gamma+10^{-7},\qquad
 \|H^Ty\|_\infty\leq1+10^{-7},
\]

\[
 \big|\|\omega\|_1-y_k+\gamma\|y\|_1\big|
 \leq10^{-7}\max(1,\|\omega\|_1).
\]

HiGHS uses primal and dual feasibility tolerances of 1e-8. A restricted
infeasibility conclusion requires a dual ray checked against the complete
Hessian. Its bound can be reused at smaller tolerances for that same
Hessian. Cached matrices are compared directly before reuse.

The first direction's parametric path retains its declared budget of
4p pivots. A single path supplies its attained tolerance range to the
remaining grid values. A numerical failure retains its coordinate,
final knot and pivot budget in the candidate record.

## Validation bounds

Let V be a validation Hessian and W the symmetrised CLIME matrix. The
fold loss is

\[
 L(W)=p^{-1}\|VW-I\|_F^2.
\]

After solving raw row k, every unknown entry in its symmetrised column
has absolute value bounded by the corresponding absolute raw coefficient.
Write that column's residual as r+Cz, where r contains the known entries
and the negative canonical vector. The unknown coefficients satisfy
\(|z_j|\leq a_j\). For any vector v,

\[
 \|r+Cz\|_2^2\geq
 2r^Tv-\|v\|_2^2-2\sum_j a_j|C_j^Tv|.
\]

Optimising the nonnegative scale of v gives the lower bound

\[
 \frac{[r^Tv-\sum_j a_j|C_j^Tv|]_+^2}{\|v\|_2^2}.
\]

The implementation evaluates this expression using the residual and a
least-squares projection as candidate directions. The inequality holds
for either direction, including when the projection matrix is rank
deficient. A rounding allowance is subtracted before evaluating the
bound. Fully resolved columns contribute their exact squared residuals.

Adding the column bounds gives a lower bound for the complete fold
loss. Adding already evaluated fold losses and dividing by five gives a
lower bound for the candidate's complete validation criterion. A candidate
is eliminated only when this bound exceeds the incumbent with the stated
rounding margin. Every eligible candidate retains all required rows and
folds. The final selected fit uses its complete symmetrised matrix.

## Reproduction

`research/submission_completion_20261027/install_compiled_clime.R` builds
the official HiGHS 1.15.1 source in the project's ignored dependency
directory. It uses the documented
[HiGHS C++ build interface](https://ergo-code.github.io/HiGHS/dev/interfaces/cpp/).
The driver is compiled with RcppEigen and uses one solver thread per
worker process.

`check_compiled_clime.R` compares full solutions with the independent R
HiGHS interface and checks the validation bounds against box-constrained
quadratic programmes solved by OSQP. `check_sqr_accelerated_procedure.R`
repeats complete saved datasets and compares the chosen tuning values,
coefficients and interval endpoints. `check_sqr_registered_case.R`
checks a registered setting using independent full-matrix solutions for
every selected validation fold and the final correction.

The execution scheduler allocates at most twelve comparator processes.
Completed queues release their allocation to the remaining queues. Each
replication retains its original dataset and tuning seeds, and completed
results remain available during a numerical implementation update.
