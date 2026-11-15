.aq_alpha_one <- function(x) {
  x <- as.numeric(x); x <- x[!is.na(x)]
  if (any(x < 0) || any(!is.finite(x))) stop("Abundances must be finite and non-negative.")
  x <- x[x > 0]; S <- length(x); N <- sum(x)
  ids <- c("Observed","Chao1","ACE","Fisher","Shannon","Simpson","InvSimpson",
           "GiniSimpson","Pielou","Heip","Evar","Coverage","BergerParker")
  out <- stats::setNames(rep(NA_real_, length(ids)), ids)
  if (!S || !N) return(out)
  p <- x/N; H <- -sum(p*log(p)); D <- sum(p^2)
  out["Observed"] <- S
  est <- vegan::estimateR(x)
  # vegan may return a vector or matrix depending on input shape/version.
  if (is.matrix(est) || is.data.frame(est)) {
    out["Chao1"] <- unname(est["S.chao1", 1])
    out["ACE"] <- unname(est["S.ACE", 1])
  } else {
    out["Chao1"] <- unname(est["S.chao1"])
    out["ACE"] <- unname(est["S.ACE"])
  }
  out["Fisher"] <- tryCatch(vegan::fisher.alpha(x), error=function(e) NA_real_)
  out["Shannon"] <- H; out["Simpson"] <- D; out["InvSimpson"] <- 1/D
  out["GiniSimpson"] <- 1-D
  out["Pielou"] <- if (S > 1) H/log(S) else NA_real_
  out["Heip"] <- if (S > 1) (exp(H)-1)/(S-1) else NA_real_
  out["Evar"] <- if (S > 1) 1-(2/pi)*atan(mean((log(p)-mean(log(p)))^2)) else NA_real_
  out["Coverage"] <- 1-sum(x==1)/N
  out["BergerParker"] <- max(p)
  out
}

#' Calculate alpha diversity indices
#' @param data aq_data.
#' @param index Indices or "all".
#' @param group Metadata columns to retain.
#' @export
aq_alpha <- function(data, index = "all", group = NULL) {
  if (!inherits(data, "aq_data")) stop("data must be aq_data.")
  valid <- c("Observed","Chao1","ACE","Fisher","Shannon","Simpson","InvSimpson",
             "GiniSimpson","Pielou","Heip","Evar","Coverage","BergerParker")
  if (identical(index, "all")) index <- valid
  if (!all(index %in% valid)) stop("Unknown alpha index.")
  if (anyNA(data$otu) || any(data$otu < 0) || any(!is.finite(data$otu)))
    stop("OTU must contain finite non-negative values without NA.")
  if (any(abs(data$otu-round(data$otu)) > 1e-8))
    warning("Some abundance values are non-integers; richness estimators are designed for count data.")
  m <- t(vapply(seq_len(nrow(data$otu)), function(i) .aq_alpha_one(data$otu[i,]), numeric(length(valid))))
  colnames(m) <- valid; rownames(m) <- rownames(data$otu)
  d <- data.frame(SampleID=rownames(m), m[,index,drop=FALSE], check.names=FALSE)
  if (!is.null(data$metadata)) d <- cbind(d, data$metadata[match(d$SampleID, rownames(data$metadata)),,drop=FALSE])
  if (!is.null(group) && !all(group %in% names(d))) stop("Unknown metadata grouping variable.")
  long <- do.call(rbind, lapply(index, function(i) data.frame(
    SampleID=d$SampleID, Index=i, Value=d[[i]],
    if (!is.null(data$metadata)) data$metadata[match(d$SampleID, rownames(data$metadata)),,drop=FALSE],
    check.names=FALSE
  )))
  structure(list(data=d, long=long, parameters=list(index=index, group=group)), class="aq_alpha")
}

#' Plot alpha diversity boxplots
#' @param data aq_alpha object.
#' @param index One or more indices.
#' @param group Metadata grouping variable.
#' @param facet Facet variables, including `Index`.
#' @param palette Palette.
#' @param colors Named custom colors.
#' @param fill Fill boxes.
#' @param point Show sample points.
#' @param test none, wilcox, kruskal, or anova.
#' @param p_adjust P-value adjustment.
#' @param show_p Add overall test labels.
#' @param alpha Box fill transparency.
#' @param point_size Sample point size.
#' @param theme AqueCo theme style.
#' @export
aq_alpha_plot <- function(data, index="Shannon", group, facet=NULL, palette="nature",
                          colors=NULL, fill=TRUE, alpha=0.7, point=TRUE,
                          point_size=1.8, test=c("kruskal","wilcox","anova","none"),
                          p_adjust="BH", show_p=TRUE, theme="publication") {
  if (!inherits(data,"aq_alpha")) stop("data must be aq_alpha.")
  test <- match.arg(test)
  d <- data$long[data$long$Index %in% index,,drop=FALSE]
  if (!group %in% names(d)) stop("Unknown group variable.")
  if (length(index)>1 && is.null(facet)) stop("Use facet='Index' for multiple indices.")
  if (!is.null(facet) && !all(facet %in% names(d))) stop("Unknown facet variable.")
  d[[group]] <- factor(d[[group]])
  lev <- levels(d[[group]])
  cols <- if (is.null(colors)) aq_pal(palette,length(lev)) else colors
  if (is.null(names(cols))) names(cols) <- lev
  p <- ggplot2::ggplot(d,ggplot2::aes(x=.data[[group]],y=Value,color=.data[[group]]))
  if (fill) p <- p + ggplot2::geom_boxplot(ggplot2::aes(fill=.data[[group]]),alpha=alpha,outlier.shape=NA)
  else p <- p + ggplot2::geom_boxplot(fill="white",outlier.shape=NA)
  if (point) p <- p + ggplot2::geom_jitter(width=.15,size=point_size,alpha=.65)
  p <- p + ggplot2::scale_color_manual(values=cols) + theme_aq(theme)
  if (fill) p <- p + ggplot2::scale_fill_manual(values=cols)
  if (!is.null(facet)) p <- p + ggplot2::facet_wrap(stats::reformulate(facet),scales="free_y")
  stats_out <- NULL
  if (test != "none") {
    # Build the formula from the grouping-column name. This ensures that
    # Value and the grouping variable are both evaluated inside data = z.
    test_formula <- stats::reformulate(group, response = "Value")

    split_key <- if (is.null(facet)) {
      rep("all", nrow(d))
    } else {
      interaction(d[facet], drop = TRUE, lex.order = TRUE)
    }
    groups <- split(d, split_key, drop = TRUE)

    stat_list <- lapply(groups, function(z) {
      z <- z[!is.na(z$Value) & !is.na(z[[group]]), , drop = FALSE]
      z[[group]] <- droplevels(factor(z[[group]]))

      if (length(unique(z[[group]])) < 2L) return(NULL)

      pv <- if (test == "wilcox" && nlevels(z[[group]]) == 2L) {
        stats::wilcox.test(test_formula, data = z)$p.value
      } else if (test == "anova") {
        stats::oneway.test(test_formula, data = z, var.equal = TRUE)$p.value
      } else {
        # Kruskal-Wallis is also the fallback for wilcox with >2 groups.
        stats::kruskal.test(test_formula, data = z)$p.value
      }

      row <- data.frame(
        Index = if ("Index" %in% names(z)) as.character(z$Index[1]) else NA_character_,
        p_value = as.numeric(pv),
        stringsAsFactors = FALSE
      )

      # Preserve facet columns so annotations are drawn in the right panels.
      if (!is.null(facet)) {
        for (f in facet) row[[f]] <- z[[f]][1]
      }
      row
    })

    stat_list <- Filter(Negate(is.null), stat_list)
    if (length(stat_list)) {
      stats_out <- do.call(rbind, stat_list)
      rownames(stats_out) <- NULL
      stats_out$p_adjusted <- stats::p.adjust(
        stats_out$p_value,
        method = p_adjust
      )

      if (show_p) {
        p <- p + ggplot2::geom_text(
          data = stats_out,
          mapping = ggplot2::aes(
            x = Inf,
            y = Inf,
            label = paste0("p=", signif(p_adjusted, 3))
          ),
          inherit.aes = FALSE,
          hjust = 1.1,
          vjust = 1.5
        )
      }
    }
  }
  structure(list(plot=p,statistics=stats_out,data=d),class="aq_plot")
}

#' Export alpha diversity table
#' @param x aq_alpha object.
#' @param file CSV path.
#' @export
aq_alpha_export <- function(x,file="alpha_diversity.csv") {
  if (!inherits(x,"aq_alpha")) stop("x must be aq_alpha.")
  utils::write.csv(x$data,file,row.names=FALSE,na="")
  invisible(file)
}

#' @exportS3Method print aq_alpha
print.aq_alpha <- function(x,...) {
  cat("AqueCo alpha diversity:",nrow(x$data),"samples;",
      length(x$parameters$index),"indices\n")
  print(utils::head(x$data)); invisible(x)
}


#' @exportS3Method print aq_plot
print.aq_plot <- function(x, ...) {
  if (is.null(x$plot)) {
    stop("Invalid aq_plot object: missing ggplot component.")
  }
  print(x$plot)
  invisible(x)
}
