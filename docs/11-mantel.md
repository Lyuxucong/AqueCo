# Mantel Analysis

[← Table of contents](README.md) · [Data import](02-install-data.md)

## Overview

Assess relationships between community and environmental distance matrices.

**Research question.** What ecological pattern does this method describe, and what assumptions are necessary for interpretation? Identify the response variable, sample unit, grouping factor, and analysis scale before proceeding.

## Data requirements

Start with the `dat` object created in [Data Import](02-install-data.md). Confirm the availability of relevant sample metadata, taxonomic assignments, environmental variables, and (where required) a phylogenetic tree. Use `aq_check(dat)` before each analysis.

## Main function: `aq_mantel()`

The following signature was extracted from the current package help file:

```r
aq_mantel(
  data,
  group,
  env_variables = NULL,
  level = NULL,
  abundance = c("relative", "count"),
  community_method = "bray",
  env_method = "euclidean",
  permutations = 999,
  p_adjust = "BH",
  prevalence = 0.1,
  min_abundance = 0,
  cor_method = "pearson"
)
```

### Parameter reference

| Argument | Description |
|---|---|
| `data` | An aq_data object containing otu, metadata, and env. |
| `group` | Metadata column defining treatment groups. |
| `env_variables` | Environmental variables to test; NULL selects all numeric. |
| `level` | Taxonomic rank for aggregating microbial abundances; NULL uses OTUs. |
| `abundance` | Relative abundance or counts when level is not NULL. |
| `community_method` | Community dissimilarity, e.g. bray, jaccard, euclidean. |
| `env_method` | Distance for each environmental variable: euclidean or manhattan. |
| `permutations` | Number of permutations. |
| `p_adjust` | P-value adjustment method within each group. |
| `prevalence` | Minimum taxon prevalence within each group. |
| `min_abundance` | Minimum abundance defining presence. |
| `cor_method` | Correlation method for the environmental matrix. |

### Return value

An aq_mantel object with Mantel results and environmental correlations.

## Worked example

The example uses the fixed KO/OE/WT dataset from Chapter 2. Where the bundled data do not meet a method's assumptions, code is intentionally commented out rather than claiming a successful fit.

```r
library(aqueco)
result <- aq_mantel(dat, group = "Group")
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
aq_mantel_plot(
  data,
  cor_show_sig = TRUE,
  cor_p_cutoff = 0.05,
  mantel_p_breaks = c(0.01, 0.05),
  mantel_p_labels = c("< 0.01", "0.01 – 0.05", ">= 0.05"),
  mantel_r_breaks = c(0.2, 0.3),
  show_mantel_label = FALSE,
  show_group_labels = TRUE,
  group_label_size = 3.8,
  low = "#3C80B8",
  mid = "#F7F7F7",
  high = "#F04438",
  theme = "publication",
  x_label_angle = 90,
  x_label_size = 10
)
```

```r
p <- aq_mantel_plot(result)
print(p)
```

## Interpretation and reporting

Distance-matrix observations are not independent. Report permutation settings.

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
?aq_mantel
?aq_mantel_plot
```

**Further reading:** Refer to the methods cited in the AqueCo package documentation and the primary statistical literature for this method.
