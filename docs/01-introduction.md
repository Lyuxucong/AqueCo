# Introduction to AqueCo

[← Table of contents](README.md)

AqueCo is an R package for analysis and visualization of amplicon-based microbial community data. It brings together data import and validation, taxonomic composition, alpha and beta diversity, environmental analyses, constrained ordination, variation partitioning, neutral community modeling, differential feature analysis, association networks, null-model workflows, and distance–decay analysis.

## Intended audience

This manual is written for microbial ecologists, graduate students, and researchers who need reproducible, documented workflows.

## Analysis roadmap

1. **Prepare the inputs.** Align sample IDs, taxon IDs, environmental measurements, and optional phylogenetic information.
2. **Explore the data.** Summarize composition, sample coverage, and diversity.
3. **Define a question.** Select the analysis based on the ecological hypothesis.
4. **Check assumptions.** Evaluate replication, model suitability, and data transformations.
5. **Analyze and visualize.** Save fitted objects and high-resolution figures.
6. **Report transparently.** Record software versions, parameters, random seeds, and limitations.

## Core object

`aq_import()` creates an AqueCo data object. Internally, abundance data are organized as samples × taxa. Read raw abundance tables with taxa as rows and samples as columns, then use `aq_import()` to align and construct the object.

## Tutorial conventions

- R code uses `aqueco 0.1.0`.
- The bundled example has 18 samples across KO, OE, and WT treatments.
- Placeholder figures are **not actual analytical outputs**.
- Code snippets were checked against source help signatures but have not been executed as a complete tutorial.
- Functions that need additional samples or study-specific variable definitions are presented as templates, not validated results.

## Getting help

```r
library(aqueco)
help(package = "aqueco")
?aq_import
?aq_ncm
```
