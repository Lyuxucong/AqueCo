#' Save an AqueCo plot
#' @param plot ggplot or AqueCo plot result.
#' @param filename Output filename.
#' @param width Width in inches.
#' @param height Height in inches.
#' @param dpi Resolution.
#' @export
aq_save <- function(plot, filename, width=8, height=6, dpi=600) {
  p <- if (inherits(plot,"aq_plot")) plot$plot else if (inherits(plot,"aq_beta")) NULL else plot
  if (is.null(p)) stop("Please supply a ggplot or an AqueCo plot result containing $plot.")
  ggplot2::ggsave(filename,p,width=width,height=height,dpi=dpi)
  invisible(filename)
}

#' @exportS3Method print aq_plot
print.aq_plot <- function(x, ...) {
  if (!is.null(x$plot)) print(x$plot)
  invisible(x)
}
