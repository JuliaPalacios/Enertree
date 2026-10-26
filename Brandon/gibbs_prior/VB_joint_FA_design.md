# Joint F and labeled A variational method

## Representation

The ordinary `(n-1) x (n-1)` F-matrix in `VB_master.R` determines an
**unlabeled ranked shape**. The binary `(n-1) x n` A-matrix in
`VB_master_A.R` has one named column per tip and one row per merger clade;
it determines a **labeled ranked topology**. The youngest merger is in the
bottom row, and the root is in the top row. Fixed coalescent times complete
the timed tree for likelihood evaluation. The mapping `A -> F` forgets tip
labels. Combining the two matrices adds no tree information; it provides a
way to score shape geometry and label placement separately.

The “A block” in `Variational Ranked Methods/fmatrix_blocks/check_AB_bijection.R`
is a different object: the square top-left block of an augmented labeled
F-matrix of size `(2n-2) x (2n-2)`. Its B block encodes leaf attachment ranks.
That script verifies one-to-one reconstruction of the augmented F from its
A and B blocks, or from its A block and leaf parent ranks, for `n=3,...,7`.
It does not claim that the ordinary unlabeled F is one-to-one with the
named binary A-matrix.

If an ordinary F has `c(F)` cherries, it has `m(F)=n!/2^c(F)` distinct
labelings. Exchanging the two members of a cherry leaves the labeled tree
unchanged; unique internal ranks prevent other subtree automorphisms.
The total number of labeled ranked histories is
`N_n = prod(k=2..n) choose(k,2) = n!(n-1)!/2^(n-1)`.
Enumeration gives 2,700 histories but only 16 ordinary F shapes at six tips.

## Joint family and exact gradient

Use the named binary A as the sampler state, with `F=F(A)`. Define

```
D(A;M) = ||F(A)-M||_F^2,       beta = exp(g),
s(A;M,g,Theta) = -beta D(A;M) + <Theta,A> - log m(F(A)),
q(A) = exp(s(A))/Z,            Z = sum over valid A exp(s(A)).
```

The multiplicity offset means that when `Theta=0`, summing q over the
labelings of each F gives exactly the F Gibbs family in `VB_master.R`:
`q(F) proportional to exp(-beta D(F;M))`. Conditional on F, labels are then
uniform over the `m(F)` distinct labeled histories. With nonzero Theta,
`q(A|F) proportional to exp(<Theta,A>)`, and the F marginal is its Gibbs
weight times the **average** label weight over that F fiber. Without the
offset, the baseline is uniform over labeled histories and the F marginal
contains a factor `m(F)`.

For fixed coalescent times, sequence log likelihood `ell(A)`, and a uniform
prior over labeled ranked histories, let
`h(A)=ell(A)-log N_n-s(A)`. The absolute ELBO is

```
L = E_q[h(A)] + log Z.
```

The score-function identity is
`d log q(A)/dphi = ds(A)/dphi - E_q[ds(A)/dphi]`.
Differentiating the ELBO gives `dL/dphi = Cov_q(h, ds/dphi)` because the
derivative of `-s` in h cancels the derivative of `+log Z`. Therefore the
full gradients are

```
grad_M     =  2 beta Cov_q(F, h),
grad_g     = -beta Cov_q(D, h),
grad_Theta =         Cov_q(A, h).
```

For the F-only case, set Theta to zero, average `ell(A)` over the distinct
labelings of each F, and use the first two gradients. Its joint-tree reward
includes `log m(F)` because `q(A)=q(F)/m(F)`. This connects the F approach
to the labeled A covariance gradient in `VB_master_A.R`.

## What needs Z

Metropolis-Hastings uses only `s(A')-s(A)`, so Z cancels. The covariance
gradient likewise has no numerical Z and can be estimated from MCMC draws
at the current parameters. Reporting an **absolute ELBO** still requires
the explicit `+log Z` term above.

For small n, enumerate all labeled histories for exact Z. For larger n,
`fa_logZ_reference` uses uniform labeled-history importance sampling:
`log Z = log N_n + log E_uniform[exp(s)]`. Its reference ESS must be checked;
it may be very low when q is concentrated. `fa_logZ_ais` provides an
annealed importance estimate, starting with exact uniform draws and using
intermediate scores `t*s`. For a trustworthy reported ELBO, use more
particles, a slower or adaptive temperature schedule, and replicated AIS
runs until the estimates and importance ESS are stable. Neither estimator
is used by the optimizer.

## Sampler and limitations

The sampler mixes two symmetric kernels: a uniformly chosen pair of label
columns is swapped, or an adjacent merger rank is selected and a valid
rank/NNI neighbor proposed. The phylodyn adjacent-rank/NNI proposal has
equal forward and reverse probabilities at the selected row. The acceptance
probability is `min(1, exp(s(A')-s(A)))`. Self proposals are retained. Every
valid named binary A represents one labeled ranked history; no extra
multiplicity factor is added to proposal probabilities.

`fa_fit` updates a continuous lower-triangular F center, `g`, and `Theta`
with bounded AdaGrad steps. The F center is clipped to `[0,n]` entrywise
and need not itself be a valid F-matrix. The A root row and mean of the
bottom row are fixed as score gauges. MCMC ESS and agreement among chains
must be checked separately from topology fit.

`exact_compare_VB_joint_FA.R` uses six taxa from the archived ten-tip
alignment and evaluates all 2,700 fixed-time sequence posterior states.
The repository has exact F-prior utilities, but no earlier six-tip exact
sequence-posterior comparison was found. The ten-tip example runs in memory.
The exact comparison prints its summary and writes a four-panel PNG by
default at `Brandon/gibbs_prior/joint_FA_exact_diagnostics.png`. A positional
argument chooses a different PNG path.

The plot's four panels show: (1) `KL(q || exact posterior)` for four model
distributions; (2) exact versus learned clade probabilities; (3) probabilities
of the ten labeled ranked histories with the highest exact posterior mass;
and (4) the expected number of differing A entries between two independent
draws. If `p_ij` is the inclusion probability of A entry `(i,j)`, the last
quantity is `sum_ij 2 p_ij (1-p_ij)`.

The four model bars are the uniform labeled-history prior, an F-only model
with uniform labels conditional on F, the same fitted F parameters paired
with the initial A weights, and the fully learned F+A model. The F-only and
initial-A bars are **ablations of the joint fit**, not separately optimized
variational fits. No BEAST trees are used in these plots.

From the repository root:

```
Rscript Brandon/gibbs_prior/exact_compare_VB_joint_FA.R
Rscript Brandon/gibbs_prior/exact_compare_VB_joint_FA.R /private/tmp/my_plot.png
Rscript Brandon/gibbs_prior/VB_joint_FA.R
```
