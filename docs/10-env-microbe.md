# Environment–Microbe Associations

[← Table of contents](README.md) · [Data import](02-install-data.md)

## Overview

Associate taxon abundances with environmental variables.

**Research question.** What ecological pattern does this method describe, and what assumptions are necessary for interpretation? Identify the response variable, sample unit, grouping factor, and analysis scale before proceeding.

## Data requirements

Start with the `dat` object created in [Data Import](02-install-data.md). Confirm the availability of relevant sample metadata, taxonomic assignments, environmental variables, and (where required) a phylogenetic tree. Use `aq_check(dat)` before each analysis.

## Main function: `aq_env_microbe_cor()`

The following signature was extracted from the current package help file:

```r
aq_env_microbe_cor(
  data,
  env_variables = NULL,
  level = "Genus",
  abundance = c("relative", "count"),
  method = c("spearman", "pearson", "kendall"),
  group = NULL,
  prevalence = 0.1,
  min_abundance = 0,
  top = NULL,
  p_adjust = "BH",
  use = c("pairwise.complete.obs", "complete.obs")
)
```

### Parameter reference

| Argument | Description |
|---|---|
| `data` | An aq_data object containing env and taxonomy. |
| `env_variables` | Optional environmental columns; NULL selects all numeric. |
| `level` | Taxonomic rank used to aggregate OTUs/ASVs. |
| `abundance` | Taxon abundance scale: relative or count. |
| `method` | Pearson, Spearman, or Kendall correlation. |
| `group` | Optional metadata column; correlations are calculated separately within each group. |
| `prevalence` | Minimum fraction of samples in the analysis subset in which a taxon must be present (0 to 1). |
| `min_abundance` | Minimum abundance used to define presence. |
| `top` | Optional maximum number of taxa, ranked by mean abundance. |
| `p_adjust` | P-value adjustment method across all tested taxon-environment pairs within each analysis group. |
| `use` | Missing-data handling: pairwise.complete.obs or complete.obs. |

### Return value

An aq_env_microbe_cor object with pairwise results and abundance data.

## Worked example

The example uses the fixed KO/OE/WT dataset from Chapter 2. Where the bundled data do not meet a method's assumptions, code is intentionally commented out rather than claiming a successful fit.

```r
library(aqueco)
result <- aq_env_microbe_cor(dat, group = "Group")
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
aq_env_microbe_cor_plot(
  data,
  show_value = TRUE,
  show_sig = TRUE,
  only_significant = FALSE,
  p_cutoff = 0.05,
  label_size = 3,
  low = "#3C5488",
  mid = "#F7F7F7",
  high = "#E64B35",
  theme = "publication"
)
```

```r
p <- aq_env_microbe_cor_plot(result)
print(p)
```

## Interpretation and reporting

Compositionality and multiple testing can produce spurious associations.

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
?aq_env_microbe_cor
?aq_env_microbe_cor_plot
```

**Further reading:** Refer to the methods cited in the AqueCo package documentation and the primary statistical literature for this method.
