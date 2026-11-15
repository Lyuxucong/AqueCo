# Troubleshooting and Frequently Asked Questions

[← Table of contents](README.md)

## Why is the example dataset not found using `data()`?

The example is stored as CSV files under `inst/extdata`. Use `system.file("extdata", "metadata.csv", package = "aqueco")` and `read.csv()`.

## Why does a grouping column appear to be missing?

Check `names(dat$metadata)` and confirm sample IDs are row names. The examples use a column named `Group`.

## Why does a network analysis fail or return an empty network?

The bundled example has six samples per group. `aq_network()` defaults to `min_samples = 10`, so this example may not meet its requirements. Strict multiple-testing correction and correlation thresholds can also produce no edges. Increase replication rather than relaxing thresholds solely to create a plot.

## Why is the neutral model fit poor?

The model may not adequately describe the occurrence–abundance relationship, or the data may violate assumptions. Fit quality is not a direct measure of the percentage of stochastic assembly.

## What does an R package WARNING mean?

Read the exact message. Documentation, dependency, installation, and test warnings have different implications. Do not assume a green build status means all warnings have been resolved.

## What if an optional package is unavailable?

Consult the package's `DESCRIPTION` and the function error message. Install only the dependencies needed for the intended module.

## How should I report an issue?

Provide the AqueCo version, R version, `sessionInfo()`, a minimal reproducible example, expected behavior, actual output, and the full error message. Remove any private sample metadata first.
