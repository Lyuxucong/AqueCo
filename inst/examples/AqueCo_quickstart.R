# AqueCo quick start — bundled KO/OE/WT example
# This file is documentation only. No package source changes.
# Run in an R session with aqueco installed.

library(aqueco)
library(ggplot2)

example_path <- function(name) {
  path <- system.file("extdata", name, package = "aqueco")
  if (!nzchar(path)) stop("Missing example file: ", name)
  path
}

# Input OTU format: taxa in rows, samples in columns.
otu <- read.csv(example_path("otutab.csv"), row.names = 1, check.names = FALSE)
taxonomy <- read.csv(example_path("taxonomy.csv"), row.names = 1, check.names = FALSE)
metadata <- read.csv(example_path("metadata.csv"), row.names = 1, check.names = FALSE)
env <- read.csv(example_path("env.csv"), row.names = 1, check.names = FALSE)

# The tree is optional for the basic modules below.
tree <- if (requireNamespace("ape", quietly = TRUE)) {
  ape::read.tree(example_path("otus.tree"))
} else NULL

stopifnot(all(colnames(otu) %in% rownames(metadata)))
stopifnot(all(colnames(otu) %in% rownames(env)))

dat <- aq_import(otu = otu, taxonomy = taxonomy,
                 metadata = metadata, env = env, tree = tree)
aq_check(dat)
print(table(dat$metadata$Group))

# Taxonomic composition
tax <- aq_taxa(dat, level = "Phylum")
p_tax <- aq_comp(tax, type = "bar", summary = "mean", group = "Group")
print(p_tax)

# Alpha diversity
alpha <- aq_alpha(dat, group = "Group")
p_alpha <- aq_alpha_plot(alpha, index = "Shannon", group = "Group")
print(p_alpha)

# Beta diversity
beta <- aq_beta(dat, method = "bray", group = "Group")
p_pcoa <- aq_beta_plot(beta, method = "PCoA", color = "Group")
print(p_pcoa)

# Neutral community model
ncm <- aq_ncm(dat, group = "Group")
print(ncm$summary)
p_ncm <- aq_ncm_plot(ncm, facet = TRUE)
print(p_ncm)

# Export
dir.create("AqueCo_results", showWarnings = FALSE)
ggsave("AqueCo_results/alpha_shannon.pdf", p_alpha, width = 7, height = 5)
ggsave("AqueCo_results/beta_pcoa.pdf", p_pcoa, width = 7, height = 5)
ggsave("AqueCo_results/ncm.pdf", p_ncm, width = 10, height = 5)
aq_ncm_export(ncm, file = "AqueCo_results/NCM")
saveRDS(dat, "AqueCo_results/example_data.rds")
