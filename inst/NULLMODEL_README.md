# AqueCo Null Model module (development build)

This module is an original AqueCo implementation inspired by the workflow of microeco's `trans_nullmodel`; it is not a verbatim port.

## Requirements
Install `ape` and `picante` for phylogenetic betaNRI / betaNTI; `vegan` and `ggplot2` for RCbray and visualization.

## Notes
- Input is an `aq_data` object with sample × taxon abundance internally.
- `filter_thres` applies to the mean sample-wise relative abundance (strict >).
- betaNRI = -SES(betaMPD), betaNTI = SES(betaMNTD).
- RCbray uses a richness- and total-read-preserving frequency-weighted multinomial null; results need not match microeco's RCbray algorithm.
- For RCbray, use nonnegative integer raw counts.
- Pairwise Wilcoxon p-values are exploratory because sample pairs are nonindependent.
- This build has not been run in an R environment here; run testthat and compare against microeco before publication.
