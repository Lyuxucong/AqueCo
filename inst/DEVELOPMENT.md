# Development and validation checklist

1. Replace the placeholder maintainer email in DESCRIPTION.
2. In R, run:
   ```r
   install.packages(c("devtools", "roxygen2", "testthat", "vegan", "ggplot2"))
   devtools::document()
   devtools::test()
   devtools::check()
   ```
3. Install optional packages for optional functions:
   `ggVennDiagram`, `ComplexUpset`, `betapart`, `ggtern`, `ggExtra`.
4. Review the statistical design for your experiment before interpreting PERMANOVA.
   Check multivariate dispersion with `vegan::betadisper()`; PERMANOVA can be
   sensitive to unequal dispersion.
5. NMDS stress is stored in `aq_beta_result$stress`. PCoA axis percentages use
   positive eigenvalues only. PCoA with non-Euclidean distances may require
   correction; this implementation requests `cmdscale(add=TRUE)`.
6. `marginal=TRUE` currently uses ggExtra's marginal boxplots. Axis-wise
   tests are returned in `$axis_tests`; the marginal plot is not a fully
   annotated paired panel with significance brackets in this initial release.
7. Beta diversity decomposition uses betapart pairwise Sørensen/Jaccard
   partitioning. The ternary display is a compositional display of turnover
   and nestedness shares plus a residual component; interpret it as a
   visualization of relative components, not a new independent beta index.
8. No automatic rarefaction or normalization is performed. Record and justify
   transformations in the analysis plan.
9. Run `R CMD check --as-cran` before any CRAN submission. Bioconductor has
   separate package review and release requirements.

## Palette compatibility

Base R palette names can vary by R version. The built-in palette definitions use names supported by `grDevices::hcl.colors()`; inspect your installation with `grDevices::hcl.pals()` if adding new palettes.
