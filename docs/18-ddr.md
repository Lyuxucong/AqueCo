# Distance–Decay Relationships

[← Table of contents](README.md) · [Data import](02-install-data.md)

## Overview

Analyze community dissimilarity against geographic or environmental distance.

**Research question.** What ecological pattern does this method describe, and what assumptions are necessary for interpretation? Identify the response variable, sample unit, grouping factor, and analysis scale before proceeding.

## Data requirements

Start with the `dat` object created in [Data Import](02-install-data.md). Confirm the availability of relevant sample metadata, taxonomic assignments, environmental variables, and (where required) a phylogenetic tree. Use `aq_check(dat)` before each analysis.

## Main function: `aq_ddr()`

The following signature was extracted from the current package help file:

```r
aq_ddr(data, distance=c("geographic","environmental"), group=NULL, groups=NULL,
  scope=c("overall","group","both"), lon="Longitude", lat="Latitude",
  env_variables=NULL, standardize=TRUE, permutations=999, seed=123)
```

### Parameter reference

| Argument | Description |
|---|---|
| `data` | aq_data object. |
| `distance` | geographic, environmental, or both. |
| `group` | Metadata group column. |
| `groups` | Optional selected group names. |
| `scope` | overall, group, or both. |
| `lon` | Longitude column name in metadata. |
| `lat` | Latitude column name in metadata. |
| `env_variables` | Numeric environmental columns; by default all numeric except coordinates. |
| `standardize` | Standardize environmental variables before Euclidean distance. |
| `permutations` | Mantel permutation count. |
| `seed` | Random seed. |

### Return value

Inspect the returned object with `str(result)` and consult the installed help page for the documented return structure.

## Worked example

The example uses the fixed KO/OE/WT dataset from Chapter 2. Where the bundled data do not meet a method's assumptions, code is intentionally commented out rather than claiming a successful fit.

```r
library(aqueco)
result <- aq_ddr(dat, distance = "environmental", group = "Group", scope = "both")
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
aq_ddr_plot(data, distance = NULL, groups = NULL, compare = TRUE,
  show_ci = TRUE, show_points = TRUE)
```

```r
p <- aq_ddr_plot(result)
print(p)
```

## Interpretation and reporting

Linear fit and permutation-based Mantel tests answer different questions.

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
?aq_ddr
?aq_ddr_plot
```

**Further reading:** Refer to the methods cited in the AqueCo package documentation and the primary statistical literature for this method.
