#' Variation partitioning analysis (VPA) of microbial communities
#'
#' Explanatory variables are supplied as named groups (for example Physical
#' and Chemical). If spatial analysis is requested, longitude/latitude are
#' converted to geographic distances and positive dbMEM eigenfunctions are
#' added as a separate spatial explanatory table.
#'
#' @param data An aq_data object.
#' @param groups A named list of character vectors. Each element names columns
#'   in data$env, for example list(Physical=c("Temperature","Depth"),
#'   Chemical=c("TN","TP","pH")).
#' @param include_space Include spatial variation from metadata coordinates.
#' @param latitude Name of the latitude column in data$metadata.
#' @param longitude Name of the longitude column in data$metadata.
#' @param spatial_name Name used for the spatial explanatory set.
#' @param env_transform Transformation for environmental variables: "scale"
#'   (center and scale; default), "none", "log1p_scale", or "log10_scale".
#' @param community_transform Community transformation; default Hellinger.
#' @param min_prevalence Minimum taxon prevalence before VPA.
#' @param min_abundance Minimum abundance defining presence.
#' @param permutations Number of permutations used by vegan::varpart when
#'   estimating CCA fractions; for the default Hellinger-RDA VPA this is not
#'   used to test individual fractions.
#' @param spatial_autocor dbMEM eigenfunctions to retain: positive or non-null.
#' @return An aq_vpa object containing the vegan varpart result and input tables.
#' @export
aq_vpa <- function(data, groups, include_space = FALSE,
                   latitude = "latitude", longitude = "longitude",
                   spatial_name = "Space",
                   env_transform = c("scale", "none", "log1p_scale", "log10_scale"),
                   community_transform = "hellinger",
                   min_prevalence = 0, min_abundance = 0,
                   permutations = 999,
                   spatial_autocor = "positive") {
  if (!inherits(data, "aq_data")) stop("`data` must be an aq_data object.")
  if (is.null(data$env)) stop("Environmental data are required in data$env.")
  if (!is.list(groups) || length(groups) < 1L ||
      is.null(names(groups)) || any(!nzchar(names(groups))))
    stop("`groups` must be a named list of environmental variable groups.")
  if (length(groups) + as.integer(isTRUE(include_space)) < 2L ||
      length(groups) + as.integer(isTRUE(include_space)) > 4L)
    stop("VPA requires 2 to 4 explanatory sets (environmental groups plus optional space).")
  if (anyDuplicated(names(groups))) stop("Group names must be unique.")
  env_transform <- match.arg(env_transform)
  if (!identical(community_transform, "hellinger"))
    stop("For this VPA module, community_transform must be 'hellinger'.")
  if (!is.numeric(min_prevalence) || min_prevalence < 0 || min_prevalence > 1)
    stop("`min_prevalence` must be between 0 and 1.")
  if (!is.numeric(permutations) || length(permutations) != 1L || permutations < 1)
    stop("`permutations` must be a positive integer.")

  samples <- rownames(data$otu)
  env_all <- as.data.frame(data$env, check.names = FALSE)
  if (is.null(rownames(env_all)) || !all(samples %in% rownames(env_all)))
    stop("data$env must contain all sample IDs as row names.")
  env_all <- env_all[samples, , drop = FALSE]

  selected_names <- unique(unlist(groups, use.names = FALSE))
  if (!length(selected_names)) stop("At least one environmental variable is required.")
  missing_vars <- setdiff(selected_names, names(env_all))
  if (length(missing_vars))
    stop("Environmental variables not found in data$env: ",
         paste(missing_vars, collapse = ", "))
  non_numeric <- selected_names[!vapply(env_all[selected_names], is.numeric, logical(1))]
  if (length(non_numeric))
    stop("VPA environmental variables must be numeric: ",
         paste(non_numeric, collapse = ", "))
  if (anyDuplicated(selected_names))
    stop("An environmental variable cannot belong to more than one group.")

  md <- data$metadata
  if (isTRUE(include_space)) {
    if (is.null(md)) stop("Spatial VPA requires data$metadata.")
    if (!all(c(latitude, longitude) %in% names(md)))
      stop("Spatial columns not found in metadata: ",
           paste(setdiff(c(latitude, longitude), names(md)), collapse = ", "))
    md <- md[samples, , drop = FALSE]
    if (!is.numeric(md[[latitude]]) || !is.numeric(md[[longitude]]))
      stop("Latitude and longitude columns must be numeric decimal degrees.")
    lat <- md[[latitude]]
    lon <- md[[longitude]]
    if (any(!is.finite(lat) | lat < -90 | lat > 90) ||
        any(!is.finite(lon) | lon < -180 | lon > 180))
      stop("Invalid geographic coordinates. Latitude must be [-90, 90] and longitude [-180, 180].")
  }

  # Align all inputs and drop samples with missing values in any selected predictor.
  complete <- stats::complete.cases(env_all[, selected_names, drop = FALSE])
  if (isTRUE(include_space))
    complete <- complete & stats::complete.cases(md[, c(latitude, longitude), drop = FALSE])
  keep_samples <- samples[complete]
  if (length(keep_samples) < 5L)
    stop("At least five complete samples are required for VPA.")
  env <- env_all[keep_samples, , drop = FALSE]
  otu <- as.matrix(data$otu[keep_samples, , drop = FALSE])
  prev <- colMeans(otu > min_abundance)
  otu <- otu[, is.finite(prev) & prev >= min_prevalence, drop = FALSE]
  otu <- otu[, colSums(otu) > 0, drop = FALSE]
  if (ncol(otu) < 2L) stop("At least two taxa must remain after prevalence filtering.")
  Y <- vegan::decostand(otu, method = "hellinger")

  transform_group <- function(x, label) {
    x <- as.data.frame(x, check.names = FALSE)
    if (env_transform == "log1p_scale") {
      if (any(as.matrix(x) < 0)) stop("log1p_scale requires non-negative environmental variables in group ", label)
      x[] <- lapply(x, log1p)
    } else if (env_transform == "log10_scale") {
      if (any(as.matrix(x) <= 0)) stop("log10_scale requires strictly positive environmental variables in group ", label)
      x[] <- lapply(x, log10)
    }
    varying <- vapply(x, function(v) is.finite(stats::sd(v)) && stats::sd(v) > 0, logical(1))
    x <- x[, varying, drop = FALSE]
    if (!ncol(x)) stop("No variable with non-zero variance remains in group ", label)
    if (env_transform != "none") {
      x[] <- lapply(x, function(v) as.numeric(scale(v)))
    }
    x
  }

  X <- lapply(names(groups), function(nm) {
    transform_group(env[, groups[[nm]], drop = FALSE], nm)
  })
  names(X) <- names(groups)
  if (any(vapply(X, ncol, integer(1)) < 1L))
    stop("Each explanatory group must contain at least one usable variable.")

  spatial_info <- NULL
  if (isTRUE(include_space)) {
    if (!requireNamespace("adespatial", quietly = TRUE)) {
      stop("Package `adespatial` is required for spatial dbMEM. Install with install.packages('adespatial').",
           call. = FALSE)
    }
    if (!requireNamespace("geosphere", quietly = TRUE)) {
      stop("Package `geosphere` is required to calculate geographic distances. Install with install.packages('geosphere').",
           call. = FALSE)
    }
    lat <- md[keep_samples, latitude]
    lon <- md[keep_samples, longitude]
    # Use great-circle (Haversine) distances in kilometres, not raw degrees.
    coords <- cbind(longitude = lon, latitude = lat)
    rownames(coords) <- keep_samples
    if (nrow(unique(as.data.frame(coords))) < 3L)
      stop("At least three distinct sampling locations are required for spatial dbMEM.")
    D_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
    D <- stats::as.dist(D_km)
    spatial <- adespatial::dbmem(
      D, MEM.autocor = match.arg(spatial_autocor, c("positive", "non-null")),
      silent = TRUE
    )
    spatial <- as.matrix(spatial)
    if (is.null(dim(spatial)) || ncol(spatial) < 1L)
      stop("No usable spatial dbMEM eigenfunctions could be generated.")
    rownames(spatial) <- keep_samples
    colnames(spatial) <- paste0("dbMEM", seq_len(ncol(spatial)))
    X[[spatial_name]] <- as.data.frame(spatial, check.names = FALSE)
    spatial_info <- list(
      latitude = latitude, longitude = longitude,
      coordinate_method = "great-circle Haversine distance (km)",
      coordinates = coords, distance_km = D_km, dbmem = spatial,
      autocorrelation = spatial_autocor
    )
  }

  if (length(X) < 2L || length(X) > 4L)
    stop("After adding spatial predictors, VPA must contain 2 to 4 explanatory sets.")
  if (nrow(Y) <= sum(vapply(X, ncol, integer(1))) + 1L)
    warning("The number of explanatory variables is high relative to sample size; adjusted fractions may be unstable.")

  # vegan::varpart expects two to four explanatory tables after Y.
  # Pass each table explicitly as a numeric matrix. Y is already Hellinger
  # transformed above, so do not request a second transformation in varpart.
  X_matrices <- lapply(X, function(z) {
    z <- as.matrix(z)
    storage.mode(z) <- "double"
    z
  })
  # Call varpart directly by the number of sets. Using do.call() here makes
  # varpart's match.call()/deparse bookkeeping inspect a nested list call,
  # which can trigger "number of items to replace..." warnings.
  if (length(X_matrices) == 2L) {
    vp <- vegan::varpart(Y, X_matrices[[1L]], X_matrices[[2L]])
  } else if (length(X_matrices) == 3L) {
    vp <- vegan::varpart(Y, X_matrices[[1L]], X_matrices[[2L]],
                         X_matrices[[3L]])
  } else {
    vp <- vegan::varpart(Y, X_matrices[[1L]], X_matrices[[2L]],
                         X_matrices[[3L]], X_matrices[[4L]])
  }
  # Preserve the names used by AqueCo rather than varpart's deparsed call.
  vp$tables <- names(X_matrices)
  ind <- as.data.frame(vp$part$indfract, check.names = FALSE)
  fract <- as.data.frame(vp$part$fract, check.names = FALSE)
  # vegan versions differ: some expose Adj.R.squared, others Adj.R.square.
  # Keep the original field and provide AqueCo's stable Adj.R.square alias.
  normalize_fraction_names <- function(x) {
    nm <- names(x)
    adj <- intersect(c("Adj.R.square", "Adj.R.squared"), nm)
    if (length(adj) && !"Adj.R.square" %in% nm)
      names(x)[match(adj[1L], names(x))] <- "Adj.R.square"
    raw <- intersect(c("R.square", "R.squared"), names(x))
    if (length(raw) && !"R.square" %in% names(x))
      names(x)[match(raw[1L], names(x))] <- "R.square"
    x
  }
  ind <- normalize_fraction_names(ind)
  fract <- normalize_fraction_names(fract)
  structure(list(
    varpart = vp,
    fractions = ind,
    table_fractions = fract,
    individual_fractions = ind,
    overall_fractions = fract,
    explanatory = X,
    community = Y,
    samples = keep_samples,
    spatial = spatial_info,
    group_names = names(X),
    parameters = list(
      env_groups = groups, include_space = include_space,
      spatial_name = if (include_space) spatial_name else NULL,
      env_transform = env_transform,
      community_transform = "hellinger",
      min_prevalence = min_prevalence,
      min_abundance = min_abundance,
      permutations = as.integer(permutations)
    )
  ), class = "aq_vpa")
}

#' Plot VPA as a Venn diagram
#' @param data An aq_vpa result.
#' @param cutoff Minimum adjusted R-squared fraction to display.
#' @param digits Number of digits in fraction labels.
#' @param alpha Circle fill transparency.
#' @param colors Optional fill colors, one per explanatory group.
#' @param cex Label size.
#' @param ... Additional arguments passed to vegan's plot.varpart.
#' @return The Venn diagram is drawn; the aq_vpa object is returned invisibly.
#' @export
aq_vpa_plot <- function(data, cutoff = 0, digits = 2,
                        alpha = 0.35, colors = NULL, cex = 1,
                        ...) {
  if (!inherits(data, "aq_vpa")) stop("`data` must be an aq_vpa object.")
  if (is.null(colors)) colors <- aq_pal("nature", length(data$group_names))
  if (length(colors) != length(data$group_names))
    stop("`colors` must have one color per explanatory group.")
  graphics::plot(
    data$varpart,
    Xnames = data$group_names,
    bg = colors,
    alpha = alpha,
    cutoff = cutoff,
    digits = digits,
    cex = cex,
    ...
  )
  invisible(data)
}

#' Print VPA summary
#' @param x An aq_vpa object.
#' @param ... Ignored.
#' @export
print.aq_vpa <- function(x, ...) {
  cat("AqueCo variation partitioning analysis\n")
  cat("Samples:", length(x$samples), "\n")
  cat("Community transformation: Hellinger\n")
  cat("Explanatory sets:", paste(x$group_names, collapse = ", "), "\n\n")
  print(x$varpart)
  invisible(x)
}
