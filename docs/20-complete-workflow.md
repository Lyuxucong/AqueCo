# Complete Example Workflow

[← Table of contents](README.md)

This walkthrough demonstrates data loading, diversity, composition, a neutral model, and figure export. Run each section in R before relying on the outputs.

```r
library(aqueco)
library(ggplot2)

p <- function(f) system.file("extdata", f, package = "aqueco")
otu <- read.csv(p("otutab.csv"), row.names = 1, check.names = FALSE)
taxonomy <- read.csv(p("taxonomy.csv"), row.names = 1, check.names = FALSE)
metadata <- read.csv(p("metadata.csv"), row.names = 1, check.names = FALSE)
env <- read.csv(p("env.csv"), row.names = 1, check.names = FALSE)
tree <- ape::read.tree(p("otus.tree"))

dat <- aq_import(otu = otu, taxonomy = taxonomy,
                 metadata = metadata, env = env, tree = tree)
aq_check(dat)

# Alpha diversity
alpha <- aq_alpha(dat, group = "Group")
p_alpha <- aq_alpha_plot(alpha, index = "Shannon", group = "Group")
print(p_alpha)

# Beta diversity
beta <- aq_beta(dat, method = "bray",
                ordination = "PCoA", group = "Group")
p_beta <- aq_beta_plot(beta, color = "Group")
print(p_beta)

# Taxonomic composition
p_comp <- aq_comp(dat, type = "bar", level = "Phylum",
                  summary = "mean", group = "Group")
print(p_comp)

# Neutral community model
ncm <- aq_ncm(dat, group = "Group")
p_ncm <- aq_ncm_plot(ncm, facet = TRUE)
print(p_ncm)

# Export
dir.create("results", showWarnings = FALSE)
ggsave("results/alpha.png", p_alpha, width = 7, height = 5, dpi = 300)
ggsave("results/beta.png", p_beta, width = 7, height = 5, dpi = 300)
ggsave("results/composition.png", p_comp, width = 9, height = 5, dpi = 300)
ggsave("results/ncm.png", p_ncm, width = 10, height = 5, dpi = 300)

saveRDS(list(alpha = alpha, beta = beta, ncm = ncm),
        "results/analysis_objects.rds")
capture.output(sessionInfo(), file = "results/sessionInfo.txt")
```

**Validation note:** This is an API-grounded example, not a claim that the code has been executed in the current environment. Confirm the object types and return values against your installed package.

## Suggested output directory

```text
results/
├── alpha.png
├── beta.png
├── composition.png
├── ncm.png
├── analysis_objects.rds
└── sessionInfo.txt
```
