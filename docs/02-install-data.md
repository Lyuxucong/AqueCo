# Installation, Data Preparation, and Import

[← Table of contents](README.md)

## Installation
Install AqueCo directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("Lyuxucong/AqueCo", dependencies = TRUE)
```

Load the package and check the version:

```r
library(aqueco)
packageVersion("aqueco")
```

> If the repository is private, installation requires GitHub access and authentication.

## Documentation

See the [Complete AqueCo User Manual](docs/README.md) for detailed instructions and examples.

## Author and Maintainer

**Xucong Lyu**  
Email: lvxvcong@outlook.com  
Repository: https://github.com/Lyuxucong/AqueCo



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
