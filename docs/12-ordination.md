# RDA and CCA

[← Table of contents](README.md) · [Data import](02-install-data.md)

## Overview

Perform constrained ordination with environmental predictors.

**Research question.** What ecological pattern does this method describe, and what assumptions are necessary for interpretation? Identify the response variable, sample unit, grouping factor, and analysis scale before proceeding.

## Data requirements

Start with the `dat` object created in [Data Import](02-install-data.md). Confirm the availability of relevant sample metadata, taxonomic assignments, environmental variables, and (where required) a phylogenetic tree. Use `aq_check(dat)` before each analysis.

## Main function: `aq_rda_cca()`

The following signature was extracted from the current package help file:

```r
aq_rda_cca(
  data,
  method = c("RDA", "CCA"),
  env_variables = NULL,
  group = NULL,
  level = NULL,
  abundance = c("count", "relative"),
  transform = NULL,
  select = "forward",
  permutations = 999,
  p_cutoff = 0.05,
  prevalence = 0,
  min_abundance = 0,
  scaling = 2,
  mode = "publication",
  tests = c("global", "term", "axis")
)
```

### Parameter reference

| Argument | Description |
|---|---|
| `data` | An aq_data object containing otu, env, and optionally metadata. |
| `method` | Constrained ordination method: RDA or CCA. |
| `env_variables` | Optional environmental variables to consider. NULL selects all numeric columns in data$env. |
| `group` | Optional metadata column used to color sample points. |
| `level` | Optional taxonomic rank for aggregating taxa; NULL uses OTU/ASV. |
| `abundance` | Aggregated abundance scale when level is supplied: count or relative. |
| `transform` | Community transformation. RDA defaults to hellinger; CCA uses none. |
| `select` | Forward-selection method. Currently forward. |
| `permutations` | Number of permutations for selection and tests. |
| `p_cutoff` | Significance threshold used by forward selection. |
| `prevalence` | Minimum fraction of samples in which a taxon is present. |
| `min_abundance` | Minimum abundance defining presence. |
| `scaling` | Vegan ordination score scaling, 1 or 2. |
| `mode` | Computation mode: `publication` or `fast`. |
| `tests` | Permutation tests to request: `global`, `term`, and/or `axis`. |

### Return value

An aq_ordination object containing the selected model, scores and tests.

## Worked example

The example uses the fixed KO/OE/WT dataset from Chapter 2. Where the bundled data do not meet a method's assumptions, code is intentionally commented out rather than claiming a successful fit.

```r
library(aqueco)
result <- aq_rda_cca(dat, method = "RDA", group = "Group")
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
aq_rda_cca_plot(
  data,
  show_environment = TRUE,
  show_species = FALSE,
  show_ellipse = FALSE,
  ellipse_level = 0.95,
  point_size = 3,
  point_alpha = 0.85,
  label_samples = FALSE,
  palette = "nature",
  colors = NULL,
  arrow_multiplier = 1,
  theme = "publication"
)
```

```r
p <- aq_rda_cca_plot(result)
print(p)
```

## Interpretation and reporting

Assess collinearity and the ratio of predictors to samples.

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
?aq_rda_cca
?aq_rda_cca_plot
```

**Further reading:** Refer to the methods cited in the AqueCo package documentation and the primary statistical literature for this method.
