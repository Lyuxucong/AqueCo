# Environmental Variables

[← Table of contents](README.md) · [Data import](02-install-data.md)

## Overview

Summarize and compare environmental measurements.

**Research question.** What ecological pattern does this method describe, and what assumptions are necessary for interpretation? Identify the response variable, sample unit, grouping factor, and analysis scale before proceeding.

## Data requirements

Start with the `dat` object created in [Data Import](02-install-data.md). Confirm the availability of relevant sample metadata, taxonomic assignments, environmental variables, and (where required) a phylogenetic tree. Use `aq_check(dat)` before each analysis.

## Main function: `aq_env()`

The following signature was extracted from the current package help file:

```r
aq_env(
  data,
  group,
  variables = NULL,
  test = c("kruskal", "anova", "none"),
  posthoc = c("both", "dunn", "wilcox", "none"),
  p_adjust = "BH"
)
```

### Parameter reference

| Argument | Description |
|---|---|
| `data` | An `aq_data` object containing `env` and `metadata`. |
| `group` | One metadata column defining the groups. |
| `variables` | Optional character vector of environmental variables. `NULL` selects all numeric columns in `data$env`. |
| `test` | Overall test: `kruskal`, `anova`, or `none`. |
| `posthoc` | Pairwise test: `dunn`, `wilcox`, `both`, or `none`. |
| `p_adjust` | P-value adjustment method, e.g. `BH`, `holm`, or `bonferroni`. |

### Return value

An `aq_env` object containing long-format data, overall tests, pairwise tests, and analysis settings.

## Worked example

The example uses the fixed KO/OE/WT dataset from Chapter 2. Where the bundled data do not meet a method's assumptions, code is intentionally commented out rather than claiming a successful fit.

```r
library(aqueco)
result <- aq_env(dat, group = "Group")
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
aq_env_plot(
  data,
  type = c("boxplot", "bar"),
  variables = NULL,
  facet = TRUE,
  facet_scales = "free_y",
  bar_error = c("SE", "SD", "CI"),
  point = TRUE,
  point_size = 1.8,
  alpha = 0.75,
  palette = "nature",
  colors = NULL,
  show_p = TRUE,
  show_sig = TRUE,
  pairwise_label = FALSE,
  pairwise_method = c("dunn", "wilcox"),
  pairwise_label_type = c("p", "signif"),
  theme = "publication"
)
```

```r
p <- aq_env_plot(result)
print(p)
```

## Interpretation and reporting

Check units, missing values, and distributional assumptions.

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
?aq_env
?aq_env_plot
```

**Further reading:** Refer to the methods cited in the AqueCo package documentation and the primary statistical literature for this method.
