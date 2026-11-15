# Beta Diversity Partitioning

[← Table of contents](README.md) · [Data import](02-install-data.md)

## Overview

Partition pairwise dissimilarity into turnover and nestedness components.

**Research question.** What ecological pattern does this method describe, and what assumptions are necessary for interpretation? Identify the response variable, sample unit, grouping factor, and analysis scale before proceeding.

## Data requirements

Start with the `dat` object created in [Data Import](02-install-data.md). Confirm the availability of relevant sample metadata, taxonomic assignments, environmental variables, and (where required) a phylogenetic tree. Use `aq_check(dat)` before each analysis.

## Main function: `aq_beta_decomp()`

The following signature was extracted from the current package help file:

```r
aq_beta_decomp(data, group, method = c("sorensen", "jaccard"), prevalence = 0)
```

### Parameter reference

| Argument | Description |
|---|---|
| `data` | An `aq_data` object (samples in rows, taxa in columns). |
| `group` | Metadata column defining groups. Required. |
| `method` | Partitioning family: `sorensen` or `jaccard`. |
| `prevalence` | Optional within-group prevalence threshold between 0 and 1 (inclusive). A taxon is retained in a group when its prevalence is >= this value. The default 0 retains any taxon present in at least one sample. |

### Return value

An `aq_beta_decomp` object containing pairwise components, group means, and settings.

## Worked example

The example uses the fixed KO/OE/WT dataset from Chapter 2. Where the bundled data do not meet a method's assumptions, code is intentionally commented out rather than claiming a successful fit.

```r
library(aqueco)
result <- aq_beta_decomp(dat, group = "Group", method = "sorensen")
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
aq_beta_decomp_plot(
  data,
  point_size = 1.8,
  point_alpha = 0.3,
  show_mean = TRUE,
  mean_size = 4,
  palette = "nature",
  colors = NULL,
  theme = "publication"
)
```

```r
p <- aq_beta_decomp_plot(result)
print(p)
```

## Interpretation and reporting

Prevalence filters and presence/absence conversion influence partitioning.

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
?aq_beta_decomp
?aq_beta_decomp_plot
```

**Further reading:** Refer to the methods cited in the AqueCo package documentation and the primary statistical literature for this method.
