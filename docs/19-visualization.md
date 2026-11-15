# Plot Styling, Colors, and Publication Export

[← Table of contents](README.md)

## Available functions

AqueCo provides `aq_pal()`, `aq_pal_list()`, `aq_pal_show()`, `theme_aq()`, `aq_theme_set()`, and `aq_save()`.

```r
library(aqueco)
?aq_pal
?aq_pal_list
?aq_pal_show
?theme_aq
?aq_save
```

## Figure design

Keep the same treatment colors across all plots. Specify whether intervals represent standard errors, standard deviations, confidence intervals, or distributions. Use legible axis labels and export at the journal's required resolution.

## Export a plot

```r
# After constructing a ggplot object named p:
# aq_save(p, filename = "figure.tiff", width = 8, height = 6, dpi = 600)
# ggplot2::ggsave("figure.png", plot = p,
#                 width = 8, height = 6, dpi = 300)
```

## Recommended figure caption checklist

A figure caption should identify the data type, number of biological replicates, transformation, test, multiple-testing correction, confidence interval definition, and any thresholds.
