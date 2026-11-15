# AqueCo

**AqueCo** is an R package for reproducible analysis and visualization of 16S/ITS amplicon microbiome data.

> **Status:** Development / beta preparation. The currently inspected package source declares version `0.1.0`. A future `v0.2.0-beta.1` release tag should only be created after final verification. This is not a CRAN release.

[中文快速入门](docs/QUICKSTART_CN.md) · [完整分析示例](inst/examples/AqueCo_quickstart.R)

## Features

- Data import and validation; taxonomy composition
- Alpha diversity and beta diversity, including NMDS and PCoA
- Beta diversity decomposition; Venn/UpSet visualizations
- Environmental summaries and correlations; Mantel tests
- RDA/CCA constrained ordination and variation partitioning (VPA)
- Sloan neutral community model (NCM)
- LEfSe differential abundance analysis (optional `microeco`)
- Co-occurrence networks and network topology (optional `igraph`)
- Phylogenetic null models and community assembly (optional `ape` / `picante`)
- Distance–decay relationships (DDR)

**Interpretation note:** These analyses are tools, not automatic evidence of causation, species interaction, or ecological mechanism. Validate settings and statistical assumptions for each study.

## Requirements and installation

- R >= 4.1.0
- `ggplot2`, `vegan`, and the package's declared dependencies
- Some features require additional packages listed in `DESCRIPTION` under `Suggests`.

### Install from GitHub

```r
install.packages("remotes")
remotes::install_github(
  "Lyuxucong/AqueCo",
  dependencies = TRUE
)
library(aqueco)
packageVersion("aqueco")
```

If the GitHub repository is private, you must have permission and authenticate with GitHub.

### Install from a Local Source Directory

```r
remotes::install_local(
  "D:/AqueCo/aqueco",
  dependencies = TRUE
)
```

Replace the path with the actual directory containing the AqueCo `DESCRIPTION` file.

### Update AqueCo

```r
remotes::install_github(
  "Lyuxucong/AqueCo",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE
)
```

### Verify the Installation

```r
library(aqueco)
packageVersion("aqueco")
help(package = "aqueco")
```

## Bundled example data

AqueCo includes five files in `inst/extdata/`:

| File | Purpose |
|---|---|
| `otutab.csv` | ASV/OTU counts (taxa × samples) |
| `taxonomy.csv` | Taxonomic annotation |
| `metadata.csv` | Sample IDs, groups, sampling metadata, coordinates |
| `env.csv` | Environmental variables |
| `otus.tree` | Phylogenetic tree in Newick format |

The bundled dataset contains **1,363 ASVs across 18 samples**, with **KO, OE, WT** groups (6 samples per group). These files are accessed via `system.file()`; they are **not** registered datasets for `data("metadata")`.

```r
library(aqueco)

example_path <- function(file) {
  p <- system.file("extdata", file, package = "aqueco")
  if (!nzchar(p)) stop("Bundled example file not found: ", file)
  p
}

otu <- read.csv(example_path("otutab.csv"), row.names = 1,
                check.names = FALSE)
taxonomy <- read.csv(example_path("taxonomy.csv"), row.names = 1,
                     check.names = FALSE)
metadata <- read.csv(example_path("metadata.csv"), row.names = 1,
                     check.names = FALSE)
env <- read.csv(example_path("env.csv"), row.names = 1,
                check.names = FALSE)

tree <- NULL
if (requireNamespace("ape", quietly = TRUE)) {
  tree <- ape::read.tree(example_path("otus.tree"))
}

dat <- aq_import(otu = otu, taxonomy = taxonomy,
                 metadata = metadata, env = env, tree = tree)
aq_check(dat)
table(dat$metadata$Group)
```

**Input convention:** `aq_import()` expects the external OTU table to have **taxa in rows and samples in columns**. It transposes internally. Metadata and environment tables use sample IDs as row names.

## Minimal analysis

```r
alpha <- aq_alpha(dat, group = "Group")
print(alpha)
p_alpha <- aq_alpha_plot(alpha, index = "Shannon", group = "Group")
print(p_alpha)

beta <- aq_beta(dat, method = "bray", group = "Group")
p_beta <- aq_beta_plot(beta, method = "PCoA", color = "Group")
print(p_beta)

ncm <- aq_ncm(dat, group = "Group")
print(ncm$summary)
print(aq_ncm_plot(ncm, facet = TRUE))
```

For a complete import and analysis walkthrough, see [中文快速入门](docs/QUICKSTART_CN.md) and [the executable quick-start script](inst/examples/AqueCo_quickstart.R).

## Reproducibility and verification

Run checks from the package source directory:

```r
devtools::load_all()
devtools::test()
devtools::check(error_on = "never")
```

A previous local test report recorded **182 passing tests and 0 failures** for a stabilization build; this does **not** establish that the future GitHub release has been tested or that all statistical results have been independently validated. The author reports that the documentation warnings have since disappeared; the current archive has not been rechecked in an R runtime here.

For randomized analyses, set and report the seed, permutation count, filtering thresholds, sample sizes, and package versions.

## Known limitations

- LEfSe, networks, phylogenetic null models, and some visualizations require optional dependencies.
- Six samples per treatment in the example data are insufficient for robust inference in many network and permutation-based analyses.
- Network correlations should not be interpreted as confirmed biological interactions.
- Null-model implementations and RCbray randomization may differ from other R packages; independent benchmarking is recommended.
- This repository is intended as a beta research tool pending full validation.

## License, citation, and maintenance

The package declares `MIT + file LICENSE`, with Xucong Lyu as the copyright holder. Verify third-party data permissions before making the repository public.

**Author and maintainer:** Xucong Lyu (lvxvcong@outlook.com).

For citation metadata, see `CITATION.cff`. Please cite the repository and exact version/tag used in your analysis.
