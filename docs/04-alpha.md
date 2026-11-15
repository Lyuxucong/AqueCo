# Alpha Diversity

[← Table of contents](README.md) · [Data import](02-install-data.md)

## Overview

Quantify within-sample diversity and compare treatment groups.

**Research question.** What ecological pattern does this method describe, and what assumptions are necessary for interpretation? Identify the response variable, sample unit, grouping factor, and analysis scale before proceeding.

## Data requirements

Start with the `dat` object created in [Data Import](02-install-data.md). Confirm the availability of relevant sample metadata, taxonomic assignments, environmental variables, and (where required) a phylogenetic tree. Use `aq_check(dat)` before each analysis.

## Main function: `aq_alpha()`

The following signature was extracted from the current package help file:

```r
aq_alpha(data, index = "all", group = NULL)
```

### Parameter reference

| Argument | Description |
|---|---|
| `data` | aq_data. |
| `index` | Indices or "all". |
| `group` | Metadata columns to retain. |

### Return value

Inspect the returned object with `str(result)` and consult the installed help page for the documented return structure.

## Worked example

The example uses the fixed KO/OE/WT dataset from Chapter 2. Where the bundled data do not meet a method's assumptions, code is intentionally commented out rather than claiming a successful fit.

```r
library(aqueco)
result <- aq_alpha(dat, group = "Group")
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
aq_alpha_plot(
  data,
  index = "Shannon",
  group,
  facet = NULL,
  palette = "nature",
  colors = NULL,
  fill = TRUE,
  alpha = 0.7,
  point = TRUE,
  point_size = 1.8,
  test = c("kruskal", "wilcox", "anova", "none"),
  p_adjust = "BH",
  show_p = TRUE,
  theme = "publication"
)
```

```r
p <- aq_alpha_plot(result)
print(p)
```

## Interpretation and reporting

Consider sequencing depth, unequal sample sizes, and multiple comparisons.

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
?aq_alpha
?aq_alpha_plot
```

**Further reading:** Refer to the methods cited in the AqueCo package documentation and the primary statistical literature for this method.
