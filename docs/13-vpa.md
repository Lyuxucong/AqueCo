# Variation Partitioning Analysis

[← Table of contents](README.md) · [Data import](02-install-data.md)

## Overview

Partition explained community variation among predictor sets.

**Research question.** What ecological pattern does this method describe, and what assumptions are necessary for interpretation? Identify the response variable, sample unit, grouping factor, and analysis scale before proceeding.

## Data requirements

Start with the `dat` object created in [Data Import](02-install-data.md). Confirm the availability of relevant sample metadata, taxonomic assignments, environmental variables, and (where required) a phylogenetic tree. Use `aq_check(dat)` before each analysis.

## Main function: `aq_vpa()`

The following signature was extracted from the current package help file:

```r
aq_vpa(
  data,
  groups,
  include_space = FALSE,
  latitude = "latitude",
  longitude = "longitude",
  spatial_name = "Space",
  env_transform = c("scale", "none", "log1p_scale", "log10_scale"),
  community_transform = "hellinger",
  min_prevalence = 0,
  min_abundance = 0,
  permutations = 999,
  spatial_autocor = "positive"
)
```

### Parameter reference

| Argument | Description |
|---|---|
| `data` | An aq_data object. |
| `groups` | A named list of character vectors. Each element names columns in data$env, for example list(Physical=c("Temperature","Depth"), Chemical=c("TN","TP","pH")). |
| `include_space` | Include spatial variation from metadata coordinates. |
| `latitude` | Name of the latitude column in data$metadata. |
| `longitude` | Name of the longitude column in data$metadata. |
| `spatial_name` | Name used for the spatial explanatory set. |
| `env_transform` | Transformation for environmental variables: "scale" (center and scale; default), "none", "log1p_scale", or "log10_scale". |
| `community_transform` | Community transformation; default Hellinger. |
| `min_prevalence` | Minimum taxon prevalence before VPA. |
| `min_abundance` | Minimum abundance defining presence. |
| `permutations` | Number of permutations used by vegan::varpart when estimating CCA fractions; for the default Hellinger-RDA VPA this is not used to test individual fractions. |
| `spatial_autocor` | dbMEM eigenfunctions to retain: positive or non-null. |

### Return value

An aq_vpa object containing the vegan varpart result and input tables.

## Worked example

The example uses the fixed KO/OE/WT dataset from Chapter 2. Where the bundled data do not meet a method's assumptions, code is intentionally commented out rather than claiming a successful fit.

```r
library(aqueco)
# Define predictor sets from names(dat$env) before running:
# result <- aq_vpa(dat, groups = predictor_sets)
```

Inspect the result before further processing:

```r
# str(result)
# names(result)
# summary(result)
```

## Plotting

The current plotting/companion function signature is:

```r
aq_vpa_plot(
  data,
  cutoff = 0,
  digits = 2,
  alpha = 0.35,
  colors = NULL,
  cex = 1,
  ...
)
```

```r
# Run the analysis on a suitable dataset first.
# p <- aq_vpa_plot(result)
```

## Interpretation and reporting

Shared fractions are not causal effects; adjusted fractions may be negative.

For a manuscript, report the sample size, filtering criteria, input abundance scale, grouping variables, relevant model settings, test statistics, P-value adjustment, and package version. Do not infer causal mechanisms solely from a visualization.

### Example Methods wording

> We analyzed the relevant microbial community characteristics using AqueCo (version 0.1.0). The analysis settings, input scale, filtering thresholds, and statistical tests were specified according to the study design. [Replace this sentence with the actual method, parameters, and sample numbers.]

### Example Results wording

> The analysis revealed [insert verified result]. The corresponding figure displays [insert axes and groups]. Statistical support was assessed using [insert actual test and correction]. [Do not fill these placeholders without analyzing the study data.]

## Exporting results

When the plotting function returns a `ggplot` object, export it with:

```r
# ggplot2::ggsave("figure.png", plot = p,
#                 width = 8, height = 6, dpi = 300)
```

Some modules expose dedicated `aq_*_export()` functions. Consult their help pages before assuming output filenames or object members.

## Troubleshooting

- **Missing grouping variable:** confirm the column name in `dat$metadata`.
- **Missing optional dependency:** inspect the error and install the required package.
- **Empty output:** inspect sample counts, abundance thresholds, missing values, and assumptions.
- **Unexpected significance:** check multiple testing, permutations, and pseudoreplication.

## Function documentation

```r
?aq_vpa
?aq_vpa_plot
```

**Further reading:** Refer to the methods cited in the AqueCo package documentation and the primary statistical literature for this method.
