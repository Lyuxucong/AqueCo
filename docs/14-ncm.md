# Neutral Community Model

[← Table of contents](README.md) · [Data import](02-install-data.md)

## Overview

Fit an occurrence–abundance neutral community model.

**Research question.** What ecological pattern does this method describe, and what assumptions are necessary for interpretation? Identify the response variable, sample unit, grouping factor, and analysis scale before proceeding.

## Data requirements

Start with the `dat` object created in [Data Import](02-install-data.md). Confirm the availability of relevant sample metadata, taxonomic assignments, environmental variables, and (where required) a phylogenetic tree. Use `aq_check(dat)` before each analysis.

## Main function: `aq_ncm()`

The following signature was extracted from the current package help file:

```r
aq_ncm(
  data,
  group = NULL,
  groups = NULL,
  level = NULL,
  min_prevalence = 0,
  min_abundance = 0,
  N = NULL
)
```

### Parameter reference

| Argument | Description |
|---|---|
| `data` | An aq_data object. |
| `group` | Optional metadata column defining independent communities. |
| `groups` | Optional subset of levels to analyze. |
| `level` | Optional taxonomic rank to aggregate with aq_taxa. |
| `min_prevalence` | Minimum within-community prevalence for retaining taxa. |
| `min_abundance` | Minimum count defining taxon presence. |
| `N` | Optional community size used in the Sloan model. Defaults to mean library size within each analysis set. |

### Return value

An aq_ncm object containing summary, taxon-level results and fits.

## Worked example

The example uses the fixed KO/OE/WT dataset from Chapter 2. Where the bundled data do not meet a method's assumptions, code is intentionally commented out rather than claiming a successful fit.

```r
library(aqueco)
result <- aq_ncm(dat, group = "Group")
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
aq_ncm_plot(
  data,
  facet = TRUE,
  show_ci = TRUE,
  show_class = TRUE,
  point_size = 2,
  point_alpha = 0.65,
  palette = c(Above = "#D95F5F", Neutral = "#8A8A8A", Below = "#4C91C6"),
  theme = "publication"
)
```

```r
p <- aq_ncm_plot(result)
print(p)
```

## Interpretation and reporting

R-squared is a model-fit statistic, not the percentage of neutral assembly.

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
?aq_ncm
?aq_ncm_plot
```

**Further reading:** Refer to the methods cited in the AqueCo package documentation and the primary statistical literature for this method.
