# Variational ranked methods incorporation plan

This note is for the next round of experiments where the center tree or matrix
M is treated as known. The existing Boltzmann idea can stay:

```text
q_beta(T | M) proportional to exp(-beta * d(T, M)^2 / scale)

log p(Y | M, beta)
  approx logsumexp_i(log p(Y | T_i) - beta * d(T_i, M)^2 / scale)
        - logsumexp_i(              - beta * d(T_i, M)^2 / scale)
```

The old Brandon code already implements this pattern with F-matrix L2 distance:

- `Brandon/gibbs_prior/gradient_ascent_M_beta.R` builds a cache of F-matrices,
  likelihoods, and squared distances.
- `Brandon/gibbs_prior/utils.R` defines `Logdistance_prob()`, `sampleF()`, and
  `log_likelihood_given_tree()`.
- `Brandon/gibbs_prior/gibbs_prior_joint_standardized_helpers.R` contains the
  MH experiments, F-matrix posterior comparisons, and RF diagnostics.

For the known-M phase, the cleanest abstraction is:

```text
support item: one candidate ranked tree/topology T_i
center:       known M, represented in the same space as T_i
energy:       D2_i = d(T_i, M)^2 / scale
likelihood:   ll_i = log p(Y | T_i)
weight:       exp(ll_i - beta * D2_i)
```

Then each idea below is mostly a different `support item` representation and
`distance_to_center()` function.

## Option 1: fully heterochronous F-matrix with fake branch extension

What changes:

- Replace isochronous F-matrices of size `(n - 1) x (n - 1)` with fully
  heterochronous F-matrices of size `(2n - 2) x (2n - 2)`.
- The diagonal is no longer fixed as `2, 3, ..., n`. It is a `+1/-1` event path:
  `+1` for a branching event and `-1` for a sampling/leaf event.
- A valid heterochronous F-matrix still has monotone rows, columns that decrease
  by at most one, and the same four-entry local bounds as the old F matrix, but
  sequential construction needs the extra previous-diagonal rule from the paper.

How to create the fake heterochronous representation:

- Convert every candidate tree and the true M into a fully ordered event tree.
- If the tree already has meaningful root-to-node branch-length depths, rank all
  internal nodes and leaves by those depths, breaking ties with tiny epsilons.
- If the data are isochronous and all leaves tie, choose a leaf-rank policy:
  deterministic epsilon ladder by tip label, random/multiple completions, or
  marginalization over several leaf orderings.
- Use the same policy for the support and for true M. Otherwise distances will
  partly measure an arbitrary leaf-rank convention.

What we need to implement:

- `as_hetero_fmat(tree, rank_policy)`: rooted `phylo` tree to heterochronous F.
- `validate_hetero_fmat(F)`: check the theorem constraints.
- `hetero_fmat_to_event_tree(F, event_times)`: reconstruct a tree with branch
  lengths for likelihood evaluation.
- `enumerate_or_sample_hetero_fmats(n)`: support generation from the paper's
  single-pass rules, or a wrapper around the authors' code.
- `hetero_f_distance(F, M)`: initially squared L2 on the lower triangle.
- Optional later: `nearby_hetero_fmat()` if we infer M continuously again.

Impact on the old code:

- `build_gradient_cache()` changes from `rEncod()` and `Fmat_from_myencod()` to
  heterochronous support generation.
- `precompute_tree_chain_distance_cache()` can be reused if support items are
  still matrices of equal dimension.
- `log_likelihood_given_tree()` needs a heterochronous matrix-to-tree converter
  instead of `mytree_from_F(F, coal_times)`.
- `compute_log_Z_est()` must not use `ZigZag(n)` for this support. For the
  cached marginal likelihood ratio, the support-size constant cancels; for MH
  log-Z terms, use cache-only normalization or the fully heterochronous support
  size.

Main risk:

- Fake leaf ordering can inject artificial signal. If the fake extension is
  branch-length-derived, this is a feature. If it is only a tie-breaker, we
  should either keep it deterministic and document it, or average over multiple
  completions.

Best first experiment:

- Implement heterochronous F validation/enumeration for `n = 3` and reproduce
  the four example matrices from the paper.
- Lift a small known true tree into heterochronous F space, compute
  `exp(-beta * ||F - M||^2)`, and verify the cache-based beta likelihood still
  behaves sensibly.

## Option 2: Robinson-Foulds distance

What changes:

- Keep the Boltzmann kernel, but set `d(T, M)` to rooted Robinson-Foulds
  distance between tree topologies.
- Work with `phylo` trees directly instead of F-matrices as the primary support
  item.
- In the current environment, RF is available through `ape/phangorn::RF.dist`
  and richer variants are available through `TreeDist`.

What we need to implement:

- `build_tree_distance_cache(trees, M_tree, distance = "rf")`.
- `rf_distance_to_center(tree, M_tree, rooted = TRUE, normalize = TRUE)`.
- A tree-support generator. Easiest short-term: convert old F-list support to
  `phylo` with fixed event times. Longer-term: generate labeled rooted trees
  directly.
- A labeling policy. Standard RF is a labeled split distance, so labels must be
  fixed, or we need to average/minimize over labelings if we want unlabeled
  ranked-shape behavior.

Impact on the old code:

- `tree_chain_squared_distances()` becomes a vector of cached RF distances.
- `log_likelihood_given_tree()` becomes simpler because the support item is
  already a `phylo` tree.
- The gradient update for M no longer applies directly. For known M this is no
  problem; for inferred M later, use discrete MH/tree proposals or a categorical
  distribution over support.

Main risk:

- Plain RF ignores rank/order and branch lengths. It may be a useful topology
  baseline, but it is not really a ranked-tree distance unless we augment it.

Best first experiment:

- Reuse the existing RF diagnostics code, but move RF into the Boltzmann energy.
- Compare beta-weighted posterior summaries against the F-L2 version using the
  same true M and likelihood cache.

## Option 3: BHV distance

What changes:

- Use a metric tree-space distance that accounts for topology and branch lengths.
- The support item must be a branch-length-bearing rooted tree with consistent
  labels.
- If rank is the target, branch lengths must encode the same event order or fake
  extension used for M.

What we need to implement:

- Choose a BHV engine. `distory` is not installed locally; `TreeDist` is
  installed but does not expose a BHV function in this environment.
- Build `bhv_distance_to_center(tree, M_tree)` once the engine is chosen.
- Decide how candidate branch lengths are assigned: fixed true event-time grid,
  tree-specific inferred lengths, or optimized/integrated branch lengths.
- Add distance scaling because BHV is in branch-length units and can dominate
  likelihood weights if beta is not calibrated.

Impact on the old code:

- Same cache structure as RF: precompute `D2_i = BHV(T_i, M)^2`.
- No simple continuous-M gradient unless we bring in specialized geometry.
- Log-Z is cache-based unless we have exact finite support and feasible pairwise
  distance computation.

Main risk:

- BHV is the most faithful metric-tree idea but the heaviest to engineer. The
  dependency and branch-length policy are the important decisions.

Best first experiment:

- Defer until RF is working, then install or vendor a BHV implementation and run
  the same true-M beta-only experiment.

## Recommendation

Start with RF as the fastest distance-swap baseline, and in parallel implement
the heterochronous F validator/enumerator. RF will tell us whether a
distance-kernel prior helps even without rank information. The heterochronous F
route is the more on-target ranked-method extension and fits the old Boltzmann
cache best once the converter and support generator exist. BHV should come
third because it requires a dependency decision and careful branch-length policy.

Likely code to copy or adapt from Brandon later:

- `Brandon/gibbs_prior/utils.R`
- `Brandon/gibbs_prior/gibbs_prior_joint_standardized_helpers.R`
- `Brandon/gibbs_prior/gradient_ascent_M_beta.R`
- `Brandon/Fmatrix Bernoulli/utils.R` for sequential F-matrix sampling patterns

