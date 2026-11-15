# Co-occurrence Networks

[← Table of contents](README.md) · [Data import](02-install-data.md)

## Overview

Construct and describe association networks.

**Research question.** What ecological pattern does this method describe, and what assumptions are necessary for interpretation? Identify the response variable, sample unit, grouping factor, and analysis scale before proceeding.

## Data requirements

Start with the `dat` object created in [Data Import](02-install-data.md). Confirm the availability of relevant sample metadata, taxonomic assignments, environmental variables, and (where required) a phylogenetic tree. Use `aq_check(dat)` before each analysis.

## Main function: `aq_network()`

The following signature was extracted from the current package help file:

```r
aq_network(data, group="Group", groups=NULL, method="spearman",
  min_abundance=0.0005, min_prevalence=0.3, cor_cutoff=0.6, p_cutoff=0.05,
  p_adjust="BH", min_samples=10)
```

### Parameter reference

| Argument | Description |
|---|---|
| `group` | Input data or grouping parameter. |
| `data` | See module documentation in inst/NETWORK_README.md. |
| `groups` | See module documentation in inst/NETWORK_README.md. |
| `method` | See module documentation in inst/NETWORK_README.md. |
| `min_abundance` | See module documentation in inst/NETWORK_README.md. |
| `min_prevalence` | See module documentation in inst/NETWORK_README.md. |
| `cor_cutoff` | See module documentation in inst/NETWORK_README.md. |
| `p_cutoff` | See module documentation in inst/NETWORK_README.md. |
| `p_adjust` | See module documentation in inst/NETWORK_README.md. |
| `min_samples` | See module documentation in inst/NETWORK_README.md. |

### Return value

Inspect the returned object with `str(result)` and consult the installed help page for the documented return structure.

## Worked example

The example uses the fixed KO/OE/WT dataset from Chapter 2. Where the bundled data do not meet a method's assumptions, code is intentionally commented out rather than claiming a successful fit.

```r
library(aqueco)
# The example has six samples per group; default min_samples = 10.
# result <- aq_network(dat, group = "Group")
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
aq_net_plot(data, layout="fr", node_color="Phylum", node_size="Degree", top_n=10,
  edge_color="sign", show_legend=TRUE, remove_isolates=FALSE, seed=123)
```

```r
# p <- aq_net_plot(result)  # After constructing a suitable network.
```

## Interpretation and reporting

Correlation-based edges do not establish direct ecological interactions.

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
?aq_network
?aq_net_plot
```

**Further reading:** Refer to the methods cited in the AqueCo package documentation and the primary statistical literature for this method.
