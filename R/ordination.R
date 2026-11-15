#' Run RDA or CCA with forward selection of environmental variables
#'
#' Samples are rows and taxa are columns in the AqueCo object. RDA uses
#' Hellinger-transformed community data by default. CCA uses the original
#' non-negative abundance matrix because CCA has its own chi-square geometry.
#' @param data An aq_data object containing otu, env, and optionally metadata.
#' @param method Constrained ordination method: RDA or CCA.
#' @param env_variables Optional environmental variables to consider. NULL selects
#'   all numeric columns in data$env.
#' @param group Optional metadata column used to color sample points.
#' @param level Optional taxonomic rank for aggregating taxa; NULL uses OTU/ASV.
#' @param abundance Aggregated abundance scale when level is supplied: count or relative.
#' @param transform Community transformation. RDA defaults to hellinger; CCA uses none.
#' @param select Forward-selection method. Currently forward.
#' @param permutations Number of permutations for selection and tests.
#' @param p_cutoff Significance threshold used by forward selection.
#' @param prevalence Minimum fraction of samples in which a taxon is present.
#' @param min_abundance Minimum abundance defining presence.
#' @param scaling Vegan ordination score scaling, 1 or 2.
#' @param mode Computation mode: `publication` or `fast`.
#' @param tests Permutation tests to request: `global`, `term`, and/or `axis`.
#' @return An aq_ordination object containing the selected model, scores and tests.
#' @export
aq_rda_cca <- function(data, method = c("RDA", "CCA"),
                       env_variables = NULL, group = NULL,
                       level = NULL, abundance = c("count", "relative"),
                       transform = NULL, select = "forward",
                       permutations = 999, p_cutoff = 0.05,
                       prevalence = 0, min_abundance = 0,
                       scaling = 2, mode = "publication", tests = c("global", "term", "axis")) {
  if (!inherits(data, "aq_data")) stop("`data` must be an aq_data object.", call. = FALSE)
  if (is.null(data$env)) stop("Environmental data are required in data$env.", call. = FALSE)
  method <- match.arg(method)
  abundance <- match.arg(abundance)
  select <- match.arg(select, c("forward", "none"))
  mode <- match.arg(mode, c("publication", "fast"))
  tests <- match.arg(tests, c("global","term","axis"), several.ok=TRUE)
  if (mode=="fast" && identical(permutations,999)) permutations <- 99L
  if (!is.numeric(permutations) || length(permutations) != 1L ||
      permutations < 1 || permutations != as.integer(permutations))
    stop("`permutations` must be a positive integer.")
  if (!is.numeric(p_cutoff) || length(p_cutoff) != 1L ||
      p_cutoff <= 0 || p_cutoff > 1) stop("`p_cutoff` must be between 0 and 1.")
  if (!is.numeric(prevalence) || prevalence < 0 || prevalence > 1)
    stop("`prevalence` must be between 0 and 1.")
  if (!scaling %in% c(1, 2)) stop("`scaling` must be 1 or 2.")

  env <- as.data.frame(data$env, check.names = FALSE)
  samples <- rownames(data$otu)
  if (is.null(rownames(env)) || !all(samples %in% rownames(env)))
    stop("data$env row names must contain all OTU sample IDs.")
  env <- env[samples, , drop = FALSE]
  if (is.null(env_variables))
    env_variables <- names(env)[vapply(env, is.numeric, logical(1))]
  if (!is.character(env_variables) || !length(env_variables))
    stop("Select at least one numeric environmental variable.")
  unknown <- setdiff(env_variables, names(env))
  if (length(unknown)) stop("Unknown environmental variables: ",
                            paste(unknown, collapse = ", "))
  bad <- env_variables[!vapply(env[env_variables], is.numeric, logical(1))]
  if (length(bad)) stop("Environmental variables must be numeric: ",
                        paste(bad, collapse = ", "))
  env <- env[, unique(env_variables), drop = FALSE]

  if (is.null(level)) {
    comm <- as.matrix(data$otu)
    level_used <- "OTU/ASV"
  } else {
    tx <- aq_taxa(data, level = level, abundance = abundance, unclassified = "keep")
    comm <- as.matrix(tx$abundance)
    level_used <- tx$parameters$level
  }
  comm <- comm[samples, , drop = FALSE]
  prev <- colMeans(comm > min_abundance, na.rm = TRUE)
  comm <- comm[, is.finite(prev) & prev >= prevalence, drop = FALSE]
  comm <- comm[, colSums(comm, na.rm = TRUE) > 0, drop = FALSE]
  if (!ncol(comm)) stop("No taxa remain after prevalence filtering.")
  if (any(!is.finite(comm)) || any(comm < 0))
    stop("Community abundance contains invalid values.")

  transform_used <- if (is.null(transform)) {
    if (method == "RDA") "hellinger" else "none"
  } else match.arg(transform, c("none", "hellinger", "log", "total"))
  if (method == "CCA" && transform_used == "hellinger")
    stop("Hellinger transformation is intended for RDA. For CCA use transform = 'none'.")
  if (method == "RDA" && transform_used == "hellinger")
    comm_use <- vegan::decostand(comm, method = "hellinger")
  else if (transform_used == "total")
    comm_use <- vegan::decostand(comm, method = "total")
  else if (transform_used == "log")
    comm_use <- log1p(comm)
  else comm_use <- comm

  # Keep only samples complete for the selected environmental variables.
  ok <- stats::complete.cases(env)
  env <- env[ok, , drop = FALSE]
  comm_use <- comm_use[rownames(env), , drop = FALSE]
  if (nrow(env) < 4L) stop("At least four complete samples are required.")
  varying <- vapply(env, function(x) length(unique(x)) > 1L, logical(1))
  env <- env[, varying, drop = FALSE]
  if (!ncol(env)) stop("Selected environmental variables have no variation.")
  if (nrow(env) <= ncol(env) + 1L)
    warning("There are few samples relative to the number of candidate environmental variables.")

  md <- NULL
  if (!is.null(group)) {
    if (is.null(data$metadata) || length(group) != 1L ||
        !group %in% names(data$metadata))
      stop("`group` must name a column in data$metadata.")
    md <- data$metadata[rownames(env), , drop = FALSE]
    grp <- droplevels(factor(md[[group]]))
    if (anyNA(grp)) stop("Grouping variable contains missing values.")
  } else {
    grp <- factor(rep("All", nrow(env)), levels = "All")
  }

  Y <- comm_use
  full_formula <- stats::reformulate(names(env), response = "Y")
  null_formula <- Y ~ 1
  if (method == "RDA") {
    null_model <- vegan::rda(null_formula, data = env)
    full_model <- vegan::rda(full_formula, data = env)
    selected <- if (select=="none") full_model else tryCatch(
      vegan::ordiR2step(null_model, scope = stats::formula(full_model),
                        direction = "forward", R2scope = TRUE,
                        Pin = p_cutoff,
                        permutations = as.integer(permutations), trace = FALSE),
      error = function(e) stop("RDA forward selection failed: ", conditionMessage(e))
    )
  } else {
    null_model <- vegan::cca(null_formula, data = env)
    full_model <- vegan::cca(full_formula, data = env)
    selected <- if (select=="none") full_model else tryCatch(
      vegan::ordistep(null_model, scope = stats::formula(full_model),
                      direction = "forward", Pin = p_cutoff,
                      permutations = as.integer(permutations), trace = FALSE),
      error = function(e) stop("CCA forward selection failed: ", conditionMessage(e))
    )
  }
  selected_vars <- attr(stats::terms(selected), "term.labels")
  if (!length(selected_vars)) {
    model <- null_model
    if (select=="forward")
      warning("Forward selection retained no environmental variables at the specified significance threshold.")
  } else model <- selected

  global_test <- if ("global" %in% tests) tryCatch(
    vegan::anova.cca(model, permutations = as.integer(permutations)),
    error = function(e) NULL
  ) else NULL
  term_test <- if (length(selected_vars) && "term" %in% tests) tryCatch(
    vegan::anova.cca(model, by = "term", permutations = as.integer(permutations)),
    error = function(e) NULL
  ) else NULL
  axis_test <- if (length(selected_vars) && "axis" %in% tests) tryCatch(
    vegan::anova.cca(model, by = "axis", permutations = as.integer(permutations)),
    error = function(e) NULL
  ) else NULL

  eig <- if (length(selected_vars)) vegan::eigenvals(model, model = "constrained") else numeric()
  total_eig <- sum(eig, na.rm = TRUE)
  axis_df <- data.frame(
    Axis = paste0(if (method == "RDA") "RDA" else "CCA", seq_along(eig)),
    Eigenvalue = as.numeric(eig),
    Explained_constrained_pct = if (total_eig > 0) as.numeric(eig) / total_eig * 100 else NA_real_,
    Explained_total_pct = if (!is.null(model$tot.chi) && model$tot.chi > 0)
      as.numeric(eig) / model$tot.chi * 100 else NA_real_
  )
  if (nrow(axis_df) > 0) axis_df <- utils::head(axis_df, 2)

  if (method == "RDA") {
    rsq <- vegan::RsquareAdj(model)
    r2 <- unname(rsq$r.squared)
    adj_r2 <- unname(rsq$adj.r.squared)
  } else {
    r2 <- if (!is.null(model$tot.chi) && model$tot.chi > 0)
      sum(eig) / model$tot.chi else NA_real_
    adj_r2 <- NA_real_
  }
  p_global <- if (!is.null(global_test) && nrow(global_test))
    as.numeric(global_test[1, "Pr(>F)"]) else NA_real_

  get_scores <- function(display) {
    if (!length(selected_vars)) return(data.frame())
    sc <- tryCatch(vegan::scores(model, display = display,
                                 choices = 1:2, scaling = scaling),
                   error = function(e) NULL)
    if (is.null(sc) || is.null(dim(sc))) return(data.frame())
    as.data.frame(sc)
  }
  site_scores <- get_scores("sites")
  if (nrow(site_scores)) {
    site_scores$SampleID <- rownames(site_scores)
    site_scores$Group <- as.character(grp[match(site_scores$SampleID, rownames(env))])
  }
  species_scores <- get_scores("species")
  if (nrow(species_scores)) species_scores$Taxon <- rownames(species_scores)
  bp_scores <- get_scores("bp")
  if (nrow(bp_scores)) bp_scores$Variable <- rownames(bp_scores)

  structure(list(
    model = model, null_model = null_model, full_model = full_model,
    site_scores = site_scores, species_scores = species_scores,
    biplot_scores = bp_scores, selected_variables = selected_vars,
    candidate_variables = names(env), global_test = global_test,
    term_test = term_test, axis_test = axis_test, axis_explained = axis_df,
    r_squared = r2, adjusted_r_squared = adj_r2, p_value = p_global,
    env_data = env, community = comm_use, group = grp,
    parameters = list(method = method, group = group, level = level_used,
                      abundance = abundance, transform = transform_used,
                      permutations = as.integer(permutations),
                      p_cutoff = p_cutoff, prevalence = prevalence,
                      scaling = scaling, select = select, mode = mode, tests = tests)
  ), class = "aq_ordination")
}

#' Plot AqueCo RDA or CCA results
#'
#' @param data An aq_ordination result.
#' @param show_environment Show selected environmental vectors.
#' @param show_species Show species/taxon labels.
#' @param show_ellipse Draw group confidence ellipses when possible.
#' @param ellipse_level Confidence level for group ellipses.
#' @param point_size Sample point size.
#' @param point_alpha Sample point transparency.
#' @param label_samples Label sample points.
#' @param palette AqueCo discrete palette.
#' @param colors Optional named group colors.
#' @param arrow_multiplier Scaling multiplier for environmental arrows.
#' @param theme AqueCo theme.
#' @return An aq_plot object.
#' @export
aq_rda_cca_plot <- function(data, show_environment = TRUE,
                            show_species = FALSE, show_ellipse = FALSE,
                            ellipse_level = 0.95, point_size = 3,
                            point_alpha = 0.85, label_samples = FALSE,
                            palette = "nature", colors = NULL,
                            arrow_multiplier = 1, theme = "publication") {
  if (!inherits(data, "aq_ordination")) stop("`data` must be an aq_ordination object.")
  sc <- data$site_scores
  if (!nrow(sc) || ncol(sc) < 2L)
    stop("The selected model has fewer than two constrained axes; no 2D ordination plot is available.")
  ax <- data$axis_explained
  method <- data$parameters$method
  xlab <- if (nrow(ax) >= 1) sprintf("%s1 (%.1f%%)", method, ax$Explained_total_pct[1]) else paste0(method, "1")
  ylab <- if (nrow(ax) >= 2) sprintf("%s2 (%.1f%%)", method, ax$Explained_total_pct[2]) else paste0(method, "2")
  if (!"Group" %in% names(sc)) sc$Group <- "All"
  sc$Group <- factor(sc$Group)
  if (is.null(colors)) colors <- stats::setNames(aq_pal(palette, nlevels(sc$Group)),
                                                  levels(sc$Group))
  axis1 <- names(sc)[1]
  axis2 <- names(sc)[2]
  p <- ggplot2::ggplot(sc, ggplot2::aes(x = .data[[axis1]], y = .data[[axis2]],
                                         color = Group)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey80") +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_color_manual(values = colors, name = if (is.null(data$parameters$group)) "Group" else data$parameters$group) +
    ggplot2::labs(x = xlab, y = ylab,
                  title = paste0(method, " ordination"),
                  subtitle = paste0("R\u00b2 = ", formatC(data$r_squared, digits = 3, format = "f"),
                                    if (is.finite(data$adjusted_r_squared))
                                      paste0("; adjusted R\u00b2 = ", formatC(data$adjusted_r_squared, digits = 3, format = "f"))
                                    else "",
                                    "; P = ", ifelse(is.na(data$p_value), "NA",
                                    format.pval(data$p_value, digits = 3)))) +
    theme_aq(theme)

  if (isTRUE(show_ellipse) && nlevels(sc$Group) > 1) {
    p <- p + ggplot2::stat_ellipse(
      ggplot2::aes(group = Group), type = "t", level = ellipse_level,
      linewidth = 0.6, show.legend = FALSE
    )
  }
  bp <- data$biplot_scores
  if (isTRUE(show_environment) && nrow(bp)) {
    bp$Axis1 <- bp[[1]] * arrow_multiplier
    bp$Axis2 <- bp[[2]] * arrow_multiplier
    p <- p +
      ggplot2::geom_segment(
        data = bp,
        ggplot2::aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
        inherit.aes = FALSE,
        arrow = grid::arrow(length = grid::unit(0.18, "cm")),
        color = "#B24A33", linewidth = 0.6
      ) +
      ggplot2::geom_text(
        data = bp,
        ggplot2::aes(x = Axis1, y = Axis2, label = Variable),
        inherit.aes = FALSE, color = "#B24A33", size = 3.5,
        vjust = -0.5
      )
  }
  sp <- data$species_scores
  if (isTRUE(show_species) && nrow(sp)) {
    sp$Axis1 <- sp[[1]] * arrow_multiplier
    sp$Axis2 <- sp[[2]] * arrow_multiplier
    p <- p + ggplot2::geom_text(
      data = sp,
      ggplot2::aes(x = Axis1, y = Axis2, label = Taxon),
      inherit.aes = FALSE, color = "grey30", size = 2.8
    )
  }
  if (isTRUE(label_samples)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = SampleID), size = 2.8,
      nudge_y = 0.02, show.legend = FALSE
    )
  }
  structure(list(plot = p, results = data,
                 parameters = data$parameters), class = "aq_plot")
}

#' @exportS3Method print aq_ordination
print.aq_ordination <- function(x, ...) {
  cat("AqueCo ", x$parameters$method, " analysis\n",
      "Transformation: ", x$parameters$transform, "\n",
      "Selected variables: ",
      if (length(x$selected_variables)) paste(x$selected_variables, collapse = ", ") else "None",
      "\nR-squared: ", signif(x$r_squared, 4),
      "\nAdjusted R-squared: ", ifelse(is.finite(x$adjusted_r_squared),
                                       signif(x$adjusted_r_squared, 4), "NA"),
      "\nGlobal permutation P: ", signif(x$p_value, 4), "\n", sep = "")
  invisible(x)
}
