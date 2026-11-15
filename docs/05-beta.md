# Beta Diversity and Ordination

[← Table of contents](README.md) · [Data import](02-install-data.md)

## Overview

Evaluate between-sample community dissimilarity and ordination.

**Research question.** What ecological pattern does this method describe, and what assumptions are necessary for interpretation? Identify the response variable, sample unit, grouping factor, and analysis scale before proceeding.

## Data requirements

Start with the `dat` object created in [Data Import](02-install-data.md). Confirm the availability of relevant sample metadata, taxonomic assignments, environmental variables, and (where required) a phylogenetic tree. Use `aq_check(dat)` before each analysis.

## Main function: `aq_beta()`

The following signature was extracted from the current package help file:

```r
aq_beta(
  data,
  method = "bray",
  ordination = c("NMDS", "PCoA"),
  group = NULL,
  permutations = 999,
  transform = c("none", "hellinger", "log", "relative"),
  k = 2,
  trymax = 100,
  binary = FALSE,
  formula = NULL
)
```

### Parameter reference

| Argument | Description |
|---|---|
| `data` | aq_data. |
| `method` | Distance method passed to vegan::vegdist. |
| `ordination` | NMDS or PCoA. |
| `group` | Optional metadata variable for PERMANOVA. |
| `permutations` | Number of permutations. |
| `transform` | none, hellinger, log, or relative. |
| `k` | Number of ordination dimensions. |
| `trymax` | NMDS maximum random starts. |
| `binary` | Presence/absence distance option. |
| `formula` | Optional PERMANOVA formula, e.g. `~ Treatment + Time`. |

### Return value

Inspect the returned object with `str(result)` and consult the installed help page for the documented return structure.

## Worked example

The example uses the fixed KO/OE/WT dataset from Chapter 2. Where the bundled data do not meet a method's assumptions, code is intentionally commented out rather than claiming a successful fit.

```r
library(aqueco)
result <- aq_beta(dat, method = "bray", ordination = "PCoA", group = "Group")
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
aq_beta_plot(
  data,
  method = NULL,
  color = NULL,
  shape = NULL,
  size = 3,
  point_alpha = 0.85,
  palette = "nature",
  colors = NULL,
  centroid = FALSE,
  ellipse = FALSE,
  ellipse_level = 0.95,
  label = FALSE,
  marginal = FALSE,
  test = c("kruskal", "wilcox", "none"),
  theme = "publication"
)
```

```r
p <- aq_beta_plot(result)
print(p)
```

## Interpretation and reporting

Ordination displays dissimilarities; apparent group separation alone does not establish statistical significance.

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
?aq_beta
?aq_beta_plot
```

**Further reading:** Refer to the methods cited in the AqueCo package documentation and the primary statistical literature for this method.
