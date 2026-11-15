# Taxonomy and Community Composition

[← Table of contents](README.md) · [Data import](02-install-data.md)

## Overview

Summarize taxonomic abundance and create composition plots.

**Research question.** What ecological pattern does this method describe, and what assumptions are necessary for interpretation? Identify the response variable, sample unit, grouping factor, and analysis scale before proceeding.

## Data requirements

Start with the `dat` object created in [Data Import](02-install-data.md). Confirm the availability of relevant sample metadata, taxonomic assignments, environmental variables, and (where required) a phylogenetic tree. Use `aq_check(dat)` before each analysis.

## Main function: `aq_taxa()`

The following signature was extracted from the current package help file:

```r
aq_taxa(
  data,
  level = "Genus",
  abundance = c("count", "relative"),
  unclassified = c("keep", "drop")
)
```

### Parameter reference

| Argument | Description |
|---|---|
| `data` | An `aq_data` object. |
| `level` | Taxonomic level. |
| `abundance` | `count` or `relative`. |
| `unclassified` | How to handle missing labels: keep or drop. |

### Return value

An `aq_taxa` object.

## Worked example

The example uses the fixed KO/OE/WT dataset from Chapter 2. Where the bundled data do not meet a method's assumptions, code is intentionally commented out rather than claiming a successful fit.

```r
library(aqueco)
taxa <- aq_taxa(dat, level = "Genus", abundance = "relative")
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
aq_comp(
  data,
  type = c("bar", "heatmap", "bubble", "venn"),
  level = "Genus",
  summary = c("sample", "mean"),
  group = NULL,
  facet = NULL,
  abundance = "relative",
  top = NULL,
  other = FALSE,
  palette = "nature",
  colors = NULL,
  position = c("stack", "fill"),
  theme = "publication",
  facet_scales = c("free_x", "fixed", "free", "free_y"),
  x_angle = 90,
  x_size = 8,
  x_labels = TRUE,
  heatmap_scale = c("none", "row_z"),
  heatmap_colors = c("#3C5488", "#F7F7F7", "#E64B35"),
  cluster_samples = FALSE,
  cluster_taxa = FALSE,
  bubble_color = c("taxon", "group", "abundance"),
  bubble_size_range = c(1.5, 9),
  alpha = 0.85,
  show_values = FALSE,
  prevalence = 0,
  min_abundance = 0
)
```

```r
p <- aq_comp(dat, type = "bar", level = "Genus", group = "Group")
print(p)
```

## Interpretation and reporting

Taxonomic composition is relative to the data processing choices. Relative abundance is not absolute cell abundance.

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
?aq_taxa
?aq_comp
```

**Further reading:** Refer to the methods cited in the AqueCo package documentation and the primary statistical literature for this method.
