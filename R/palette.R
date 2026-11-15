.aq_palettes <- list(
  nature = c("#386641", "#6A994E", "#A7C957", "#F2E8CF", "#BC6C25", "#8C4A3C", "#577590", "#7B6D8D"),
  npg = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"),
  lancet = c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6"),
  set2 = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),
  viridis = grDevices::hcl.colors(8, "viridis"),
  inferno = grDevices::hcl.colors(8, "Inferno")
)

#' Get an AqueCo color palette
#' @param name Palette name.
#' @param n Number of colors.
#' @param type Discrete or continuous.
#' @export
aq_pal <- function(name = "nature", n = 8, type = c("discrete", "continuous")) {
  type <- match.arg(type)
  if (!name %in% names(.aq_palettes)) stop("Unknown palette: ", name)
  base <- .aq_palettes[[name]]
  if (length(n) != 1 || n < 1) stop("n must be a positive number.")
  if (type == "continuous") grDevices::colorRampPalette(base)(n) else {
    if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
  }
}

#' List available AqueCo palettes
#' @export
aq_pal_list <- function() names(.aq_palettes)

#' Display palette colors
#' @param name Palette name.
#' @param n Number of colors.
#' @export
aq_pal_show <- function(name = "nature", n = 8) {
  cols <- aq_pal(name, n)
  graphics::barplot(rep(1, length(cols)), col = cols, border = NA, axes = FALSE,
                    names.arg = seq_along(cols), ylab = "", main = name)
  invisible(cols)
}
