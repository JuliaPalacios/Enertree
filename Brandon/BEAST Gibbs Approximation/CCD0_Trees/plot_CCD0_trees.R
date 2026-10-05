library(ape)

trees <- list(
  CCD0Beast = read.nexus("CCD0Beast.tree"),
  CCD0model = read.nexus("CCD0model.tree"),
  true_tree = read.tree("true_tree.newick")
)

pdf("CCD0_tree_plots.pdf", width = 15, height = 5)
par(mfrow = c(1, length(trees)))

for (tree_name in names(trees)) {
  plot.phylo(trees[[tree_name]], main = tree_name)
}

dev.off()
