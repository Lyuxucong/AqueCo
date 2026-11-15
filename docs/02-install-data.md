# Installation, Data Preparation, and Import

[← Table of contents](README.md)

## Installation

The package declares R >= 4.1.0. Install from a local verified source directory:

```r
install.packages("remotes")
remotes::install_local("D:/AqueCo/aqueco", dependencies = TRUE)
library(aqueco)
packageVersion("aqueco")
```

After the GitHub repository is available, replace the placeholder username:

```r
remotes::install_github("YOUR_USERNAME/AqueCo", dependencies = TRUE)
```

Some workflows require packages listed in `Suggests` in `DESCRIPTION`; install these only as required for the analyses you intend to run.

## File layout and identifiers

| Input | File | Expected orientation |
|---|---|---|
| Abundance | `otutab.csv` | Taxa × samples |
| Taxonomy | `taxonomy.csv` | Taxa × taxonomic ranks |
| Metadata | `metadata.csv` | Samples × sample attributes |
| Environment | `env.csv` | Samples × environmental variables |
| Tree | `otus.tree` | Newick phylogenetic tree |

Sample IDs must match across abundance, metadata, and environment tables. Taxon IDs must match abundance and taxonomy; phylogenetic tips should match taxa for phylogenetic analyses.

## Read the bundled data

```r
library(aqueco)

p <- function(f) system.file("extdata", f, package = "aqueco")
otu <- read.csv(p("otutab.csv"), row.names = 1, check.names = FALSE)
taxonomy <- read.csv(p("taxonomy.csv"), row.names = 1, check.names = FALSE)
metadata <- read.csv(p("metadata.csv"), row.names = 1, check.names = FALSE)
env <- read.csv(p("env.csv"), row.names = 1, check.names = FALSE)
tree <- ape::read.tree(p("otus.tree"))

dat <- aq_import(
  otu = otu, taxonomy = taxonomy, metadata = metadata,
  env = env, tree = tree
)
aq_check(dat)
```

## Verify data alignment

```r
dim(otu)
dim(metadata)
table(metadata$Group)
stopifnot(setequal(colnames(otu), rownames(metadata)))
stopifnot(setequal(rownames(otu), rownames(taxonomy)))
stopifnot(setequal(colnames(otu), rownames(env)))
```

Check for duplicated identifiers, nonnumeric abundance values, negative abundances, zero-depth samples, missing covariates, and unexpected factor levels.

## Understanding the bundled example

The fixed dataset contains 18 samples, grouped into KO, OE, and WT (six samples per treatment). These data are suitable for learning how to call functions, but some complex models require greater replication or additional experimental controls.

## Reproducibility

```r
saveRDS(dat, "aqueco_data.rds")
capture.output(sessionInfo(), file = "sessionInfo.txt")
```

**Important:** The bundled CSV files are located in `inst/extdata` and should be accessed through `system.file()`. They are not automatically loaded through `data("metadata")`.
