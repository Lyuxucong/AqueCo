# Declare ggplot2 aesthetic symbols used with non-standard evaluation.
utils::globalVariables(c(
  ".data", "Abundance", "Axis1", "Axis2", "Group", "SampleID",
  "SharedPart", "Taxon", "Unique1Part", "Unique2Part",
  "Group", "Similarity_pct", "Turnover_pct", "Nestedness_pct",
  "Value", "p_adjusted"
))

utils::globalVariables(c("lda", "feature", "enrich_group"))

# Data-frame columns referenced in ggplot2 aesthetics and dplyr NSE.
utils::globalVariables(c(
  "percentage", "process", ".aq_marginal_group", "Variable1", "Variable2",
  "Correlation", "label", "Environment", "GroupValue", "Mean", "LabelX",
  "y_pos", "Group1", "from", "to", "p_category", "r_category",
  "Classification", "CI_lower", "CI_upper", "x", "y", "xend", "yend",
  "category", "Pi", "Zi", "Role", "value", "Sample1", "Sample2",
  "Variable"
))
