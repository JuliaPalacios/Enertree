# A matrix analysis exports

`export_final_A_samples.R` reads the final fit saved by `VB_master_A.R`. New
fits are saved under `Brandon/gibbs_prior/logs/VB_master_A/<experiment-id>/`.
Without an argument, the exporter selects the newest completed `result.rds`
from that location only. Pass an experiment ID, its folder, or its `result.rds`
path to export an older fit.

From the Enertree repository root, run:

```sh
Rscript "Brandon/gibbs_prior/A matrix analysis/export_final_A_samples.R"
```

In RStudio, with the working directory set to `Brandon/gibbs_prior`, use the
Console commands below. Clicking **Source** only defines the function; it
does not run the export. Running the entire script with **Run** also starts
the export, using the newest completed fit in `gibbs_prior/logs/VB_master_A`.

```r
source("A matrix analysis/export_final_A_samples.R")
export_final_A_samples("EXPERIMENT_ID")
```

When running the whole file with RStudio's **Run** button, set
`experiment_id <- "EXPERIMENT_ID"` near the top of the exporter. A command-line
argument overrides that setting. To select a fit and custom output folder from
the repository root, run:

```sh
Rscript "Brandon/gibbs_prior/A matrix analysis/export_final_A_samples.R" \
  "EXPERIMENT_ID" \
  "Brandon/gibbs_prior/A matrix analysis/outputs/EXPERIMENT_ID"
```

The default output folder is `A matrix analysis/outputs/<run-folder-name>`. It contains:

- `sequences.fasta`: the saved nucleotide alignment, with original taxon labels and site order.
- `sampled_trees.newick`: 1,000 trees, one per line, with the fit's fixed coalescent times as branch lengths.
- `sampled_trees.trees`: the same 1,000 trees in NEXUS format, suitable for reading with `ape::read.nexus()`.
- `sampled_A_matrices.rds`: the exact A matrices used to construct the exported trees.
- `sample_manifest.csv`: sample IDs, chain IDs, within-chain draw indices, and unnormalized log weights.
- `export_summary.rds`: source path, RNG seed, sampling schedule, acceptance rates, and distinct-tree count.
- `true_tree.newick`: the simulated true tree, when the saved fit contains one.

The script starts from the two saved final chain states, uses the fitted `theta`, and retains 500 draws from each chain. It uses the fit's saved burn-in and thinning settings, with RNG seed `20261024` for reproducibility. These are MCMC draws, so the file can include repeated trees. The original 400 final trees in `result.rds` are preserved separately and are not reused as the 1,000 exported draws.

Rows of `sample_manifest.csv`, lines of `sampled_trees.newick`, and trees in `sampled_trees.trees` are in the same sample order. NEXUS tree names may differ from the manifest's sample IDs, so use their order to pair records. The manifest's `log_unnormalized_q` is the A-matrix score `sum(theta * A)`, not a normalized log probability or sequence log likelihood.

## Compare with BEAST

See `BEAST_XML_MATCHING.md` for how the six-tip and ten-tip BEAST XML files
were built to match the fixed-times A-matrix model.

Place the BEAST `.trees` file in the same experiment output folder as
`sampled_trees.newick` and `sequences.fasta`, then run:

```sh
Rscript "Brandon/gibbs_prior/A matrix analysis/compare_tree_samples.R" 20261024_175352
```

You can use any output folder name or pass its full path instead of the example
ID. If there is only one exported folder, you can omit the ID. With multiple
folders, specify one; the script lists available IDs if none is selected. The script
looks for one additional `.trees` file inside that folder, regardless of its name; it
never selects BEAST trees from another experiment. If needed, give the BEAST
file and true-tree file explicitly as arguments 2 and 3. The true tree is
discovered from the experiment folder or a source folder with an exactly
matching FASTA alignment. For older exports without `true_tree.newick`, the
comparison script can recover it from the matching saved fit's `result.rds`.

In RStudio, with `Brandon/gibbs_prior` as the working directory, you can run
the entire `compare_tree_samples.R` file with **Run**. To choose among multiple
experiments this way, edit `experiment_id <- NULL` near the top of that file to,
for example, `experiment_id <- "20261024_184601"`. Or use the Console:

```r
source("A matrix analysis/compare_tree_samples.R")
compare_tree_samples("20261024_184601")
```

The script writes three PNG plots, `tree_summary.csv`, `ess_summary.csv`, and
`clade_ess.csv` to
`A matrix analysis/plots/<experiment-folder-name>/`. The plots compare the
rooted RF distance from the true tree and the number of cherries, using
proportions so BEAST and A-matrix samples of different sizes are comparable.
It includes all trees in each file without applying additional burn-in or
thinning.

`ess_summary.csv` reports `coda::effectiveSize` for the rooted RF distance,
cherry count, and the least efficient varying clade indicator. It has one row
per A-matrix chain and one row for the BEAST chain. `clade_ess.csv` gives the
frequency and ESS for each clade indicator in each chain; a constant indicator
has `NA` ESS because its autocorrelation cannot be estimated. `ess_per_draw.png`
shows these ESS values divided by the number of retained trees, allowing the
500-draw A-matrix chains to be compared with the larger BEAST sample. These
ESS values measure autocorrelation of tree summaries, not agreement with the
true tree or BEAST distribution.
