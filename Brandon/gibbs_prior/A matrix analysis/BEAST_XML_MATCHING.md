# How the matched BEAST XML was made

This note records how `outputs/20261025_221200/beast_matched_A.xml` was made for
the six-tip A-matrix experiment. The same procedure was used earlier for
`outputs/20261024_184601/beast_matched_A.xml` (ten tips). The purpose is to
compare BEAST tree samples with `VB_master_A.R` on the **same conditional tree
model**, rather than changing the sequence model or allowing BEAST to estimate
the node times that VI holds fixed.

## Source files

- Fitted VI result: `../logs/VB_master_A/20261025_221200/result.rds`.
- Exported alignment: `outputs/20261025_221200/sequences.fasta`, created by
  `export_final_A_samples.R` from that result. It has six taxa (`t1` through
  `t6`) and 1,000 DNA sites. The XML embeds these sequences, so BEAST does not
  read the FASTA separately when it runs.
- XML template: `outputs/20261024_184601/beast_matched_A.xml`. I copied its
  structure, replaced all sequence entries with the six-tip FASTA, replaced
  the starting Newick tree, and set the clock rate from the six-tip fit.
  The ten-tip sequences and starting tree were **not** retained.

The fit's `coal_times` are
`2.65484128591686, 4.38592953009912, 4.8485895762994,
4.94371830246364, 6.41912658023266`, from youngest to oldest event.

## Match to `VB_master_A.R`

| VI setting | BEAST XML setting |
| --- | --- |
| `fit$sequences` | The same taxon sequences embedded as `<sequence>` elements. |
| `phangorn::pml(..., bf = rep(0.25, 4), Q = c(1, 2, 1, 1, 2, 1), rate = 0.01)` | HKY with fixed `kappa = 2`, fixed A/C/G/T frequencies of 0.25, and a fixed strict-clock rate of 0.01. |
| No among-site rate variation or invariant-site mixture | `gammaCategoryCount = 1`, `proportionInvariant = 0`, and site `mutationRate = 1`. |
| Every A tree uses the same `fit$coal_times` | The starting tree has those internal heights. Only BEAST `Exchange` topology operators are enabled; no operator changes a node height, the clock, or substitution parameters. |
| Uniform prior over valid labeled ranked A matrices, conditional on those times | A constant-population coalescent tree prior with fixed `popSize = 1`. For contemporaneous tips and fixed event times, its density depends on the event times and lineage counts, not on which taxa merge. It is therefore constant across the allowed topologies. The choice of fixed population size changes the absolute prior density, but not relative topology probabilities. |

The two methods use different samplers: the A exporter samples from the fitted
variational distribution, while BEAST targets the conditional posterior from
the sequence likelihood and uniform topology prior. Their raw `posterior` or
ELBO values have different additive constants and should not be compared
directly. Compare tree probabilities, clades, RF distances, or cherries instead.

## Starting tree

I compared the saved `initial_tree` and `estimated_tree` under the **same**
phangorn likelihood used in VI. From `Brandon/gibbs_prior` as the R working
directory:

```r
fit <- readRDS("logs/VB_master_A/20261025_221200/result.rds")
score <- function(tree) as.numeric(phangorn::pml(
  tree, fit$sequences, bf = rep(0.25, 4),
  Q = c(1, 2, 1, 1, 2, 1), rate = fit$config$rate
)$logLik)
score(fit$initial_tree)    # about -2979.807
score(fit$estimated_tree)  # about -2975.059
```

I used `ape::write.tree(fit$estimated_tree)` as the XML `TreeParser` Newick
because it had the higher starting likelihood. This choice only initializes
the BEAST chain; it does not define the target posterior. For the earlier
ten-tip XML, the saved `initial_tree` had the higher likelihood, so that XML
starts there instead.

## Sampling and output

The XML sets `chainLength = 10000000` and logs every 1,000 steps. A complete
run writes 10,001 tree states including state 0 to
`beast_matched_A-tree.trees`, plus `beast_matched_A.log` and a state file in the
same experiment output folder. Keeping the XML and BEAST `.trees` file there
lets `compare_tree_samples.R` find the BEAST sample for this experiment.
The comparison script currently reads every tree in both files; it applies no
additional BEAST burn-in or thinning.

From the repository root, run in the six-tip output folder (the BEAST launcher
here does not handle a full input path containing spaces reliably):

```sh
cd "Brandon/gibbs_prior/A matrix analysis/outputs/20261025_221200"
"/Applications/BEAST 2.7.7/bin/beast" beast_matched_A.xml
```

The installed launcher directory says `2.7.7`, but its BEAST runtime reports
`2.7.8`; the XML requests `BEAST.base v2.7.8`.

## Validation performed

I ran `beast -validate beast_matched_A.xml`. BEAST parsed the file, reporting
six taxa, 1,000 sites, 112 site patterns, equal starting frequencies, and HKY.
I then ran a separate ten-step smoke-test copy in `/private/tmp` so the full
experiment outputs were untouched. BEAST reported starting sequence log
likelihood `-2975.0591`, matching phangorn's score for the saved
`estimated_tree`. The full 10-million-step run was **not** started during XML
creation. A ten-step run is not a convergence or mixing test.

## Repeating this for another experiment

1. Export that fit with `export_final_A_samples.R`, and use its new
   `outputs/<experiment-id>/sequences.fasta` and
   `../logs/VB_master_A/<experiment-id>/result.rds` together. Check the export
   summary's `source_result` field to confirm the pairing.
2. Read `fit$config$rate`, `fit$coal_times`, and the two saved candidate
   starting trees. Use the higher-likelihood candidate under the fixed VI
   phangorn model, or another valid fixed-times starting tree if warranted.
3. Copy a matched XML as a template. Replace **every** `<sequence>` entry,
   replace `TreeParser`'s `newick`, and set `clockRate` to `fit$config$rate`.
   Keep HKY `kappa = 2`, frequencies `0.25` each, site rate `1`, one gamma
   category, invariant proportion `0`, the fixed-population coalescent prior,
   and only topology `Exchange` operators unless the VI model has changed.
4. Save as `outputs/<experiment-id>/beast_matched_A.xml`; run
   `beast -validate beast_matched_A.xml` from that folder. Before a long run,
   check that BEAST's state-0 **sequence likelihood** agrees with phangorn's
   likelihood for the chosen starting tree. Do not compare BEAST's total
   posterior with the phangorn sequence likelihood.

If `VB_master_A.R` changes its likelihood, rate, time treatment, or prior,
update the XML model accordingly; copying this file alone would then no
longer give a matched comparison.

## Experiment `20261025_231748`

`outputs/20261025_231748/beast_matched_A.xml` was made by the same procedure,
using `../logs/VB_master_A/20261025_231748/result.rds` and the matching
`outputs/20261025_231748/sequences.fasta`. This fit has ten taxa, 1,000 sites,
rate `0.01`, and the same fixed event times as the earlier ten-tip experiment.
The exported FASTA was checked against the saved fit's sequences.

The saved `initial_tree` scored about `-3928.134` under the VI sequence
likelihood, versus about `-3999.425` for `estimated_tree`, so the XML starts
at `initial_tree`. BEAST validation reported ten taxa and 210 site patterns.
A separate ten-step smoke-test copy in `/private/tmp` gave BEAST state-0
sequence log likelihood `-3928.1344`, matching phangorn. The full chain was
not started while creating this XML.
