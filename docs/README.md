# AqueCo User Manual

**A comprehensive English-language guide to AqueCo 0.1.0**

**Author and maintainer:** Xucong Lyu · lvxvcong@outlook.com

This manual provides practical instructions, function signatures, example code, methodological notes,  for the AqueCo R package. It is designed to be read directly on GitHub, with no separate website, rendering engine, or special build step.

## Contents

1. [Introduction](01-introduction.md)
2. [Installation and Data Import](02-install-data.md)
3. [Taxonomy and Community Composition](03-taxonomy.md)
4. [Alpha Diversity](04-alpha.md)
5. [Beta Diversity and Ordination](05-beta.md)
6. [Beta Diversity Partitioning](06-beta-decomp.md)
7. [Shared and Unique Taxa](07-sets.md)
8. [Environmental Variables](08-environment.md)
9. [Environmental Correlations](09-env-correlation.md)
10. [Environment–Microbe Associations](10-env-microbe.md)
11. [Mantel Analysis](11-mantel.md)
12. [RDA and CCA](12-ordination.md)
13. [Variation Partitioning Analysis](13-vpa.md)
14. [Neutral Community Model](14-ncm.md)
15. [Differential Feature Analysis](15-lefse.md)
16. [Co-occurrence Networks](16-network.md)
17. [Null Models and Community Assembly](17-nullmodel.md)
18. [Distance–Decay Relationships](18-ddr.md)
19. [Plot Styling and Export](19-visualization.md)
20. [Complete Example Workflow](20-complete-workflow.md)
21. [Troubleshooting and FAQ](21-troubleshooting.md)
22. [Statistical Reporting and Reproducibility](22-reproducibility.md)
23. [Function Reference](23-function-reference.md)

## Quick start

```r
library(aqueco)
p <- function(f) system.file("extdata", f, package = "aqueco")
otu <- read.csv(p("otutab.csv"), row.names = 1, check.names = FALSE)
taxonomy <- read.csv(p("taxonomy.csv"), row.names = 1, check.names = FALSE)
metadata <- read.csv(p("metadata.csv"), row.names = 1, check.names = FALSE)
env <- read.csv(p("env.csv"), row.names = 1, check.names = FALSE)
tree <- ape::read.tree(p("otus.tree"))
dat <- aq_import(otu, taxonomy, metadata, env, tree)
aq_check(dat)
```

## Accuracy and versioning

The function signatures were taken from the uploaded AqueCo source package. The tutorial was not executed end-to-end in an R runtime. Please verify code and results before publishing it as a validated analytical reference.

## How to use on GitHub

Place this `docs/` folder in the root of your AqueCo repository. Add the following line to the main repository `README.md`:

```markdown
📘 [Complete AqueCo User Manual (English)](docs/README.md)
```

GitHub renders the Markdown pages and links automatically. No website configuration is required.
