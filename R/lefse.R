#' Run LEfSe using microeco
#'
#' Uses microeco::microtable and microeco::trans_diff (method = "lefse").
#' OTU abundances are stored as samples x taxa in aq_data; conversion to
#' microeco's taxa x samples convention is automatic.
#' @param data An aq_data object with taxonomy and metadata.
#' @param group Metadata column containing groups.
#' @param min_total_abundance Minimum total ASV count (inclusive).
#' @param alpha Kruskal-Wallis significance cutoff.
#' @param p_adjust_method microeco p-value adjustment method.
#' @param ... Additional arguments to microeco::trans_diff$new().
#' @return aq_lefse object with res_diff, model and dataset.
#' @export
aq_lefse <- function(data, group, min_total_abundance = 10,
                     alpha = 0.05, p_adjust_method = "none", ...) {
  if (!inherits(data,"aq_data")) stop("data must be an aq_data object.")
  if (!requireNamespace("microeco",quietly=TRUE))
    stop("Install microeco: install.packages('microeco').")
  if (is.null(data$taxonomy) || is.null(data$metadata))
    stop("LEfSe requires taxonomy and metadata.")
  if (length(group)!=1L || !is.character(group) ||
      !group %in% names(data$metadata))
    stop("group must name a metadata column.")
  if (length(min_total_abundance)!=1L || !is.numeric(min_total_abundance) ||
      !is.finite(min_total_abundance) || min_total_abundance<0)
    stop("min_total_abundance must be non-negative.")
  otu <- as.matrix(data$otu)
  otu <- otu[,colSums(otu)>=min_total_abundance,drop=FALSE]
  if (!ncol(otu)) stop("No ASVs remain after filtering.")
  if (any(rowSums(otu)==0)) stop("Filtering left samples with zero abundance.")
  if (any(abs(otu-round(otu))>1e-8))
    warning("LEfSe input is not integer counts; verify normalization and filtering.")
  taxa <- as.data.frame(data$taxonomy,check.names=FALSE)
  if (!all(colnames(otu) %in% rownames(taxa)))
    stop("Taxonomy lacks some retained ASV IDs.")
  taxa <- taxa[colnames(otu),,drop=FALSE]
  taxa[] <- lapply(taxa,function(x) {
    x <- as.character(x)
    x <- sub("^[A-Za-z]__", "",x)
    x[is.na(x) | trimws(x)==""] <- "unclassified"
    x
  })
  # microeco's own taxonomy tidier also resolves hierarchical formatting.
  taxa <- microeco::tidy_taxonomy(taxa)
  meta <- as.data.frame(data$metadata[rownames(otu),,drop=FALSE])
  if (anyNA(meta[[group]]) || length(unique(meta[[group]]))<2L)
    stop("Group column must contain at least two non-missing groups.")
  mt <- microeco::microtable$new(
    otu_table = as.data.frame(t(otu), check.names = FALSE),
    tax_table = taxa,
    sample_table = meta,
    auto_tidy = FALSE
  )
  # LEfSe method is supplied explicitly to avoid conflicting user dots.
  extra <- list(...)
  reserved <- c("dataset","method","group","alpha","p_adjust_method")
  if (any(names(extra) %in% reserved))
    stop("Do not pass reserved arguments in ...: ",
         paste(intersect(names(extra),reserved),collapse=", "))
  model <- do.call(microeco::trans_diff$new,
    c(list(dataset=mt,method="lefse",group=group,alpha=alpha,
           p_adjust_method=p_adjust_method),extra))
  structure(list(res_diff=model$res_diff,model=model,dataset=mt,
    group=group,parameters=list(min_total_abundance=min_total_abundance,
      alpha=alpha,p_adjust_method=p_adjust_method)),
    class="aq_lefse")
}

#' Plot LEfSe results with explicit feature selection
#'
#' LDA bars are constructed from res_diff rather than delegating feature
#' selection to microeco. The cladogram still uses microeco's plot method.
#' @param data aq_lefse result.
#' @param type lda or cladogram.
#' @param threshold Minimum absolute LDA score.
#' @param group_order Optional group order.
#' @param colors Optional group colors, ideally named by group.
#' @param top_n Maximum number of features across all groups; NULL for no cap.
#' @param top_n_per_group Maximum features per group; NULL for no cap.
#' @param tax_level Optional taxonomy rank (e.g. Genus), based on Taxa prefixes.
#' @param short_labels Whether to shorten hierarchical taxonomy labels.
#' @param ... For lda, supports width and use_number (alias of top_n).
#'   For cladogram, passed to microeco.
#' @return ggplot object.
#' @export
aq_lefse_plot <- function(data,type=c("lda","cladogram"),
                          threshold=2,group_order=NULL,colors=NULL,
                          top_n=NULL,top_n_per_group=NULL,
                          tax_level=NULL,short_labels=TRUE,...) {
  if (!inherits(data,"aq_lefse")) stop("data must be aq_lefse.")
  type <- match.arg(type)
  extra <- list(...)
  if (type=="cladogram") {
    if (!is.null(top_n) || !is.null(top_n_per_group) || !is.null(tax_level))
      warning("top_n, top_n_per_group and tax_level only apply to LDA bars.")
    if (!is.null(group_order)) extra$group_order <- group_order
    if (!is.null(colors)) extra$color <- colors
    return(do.call(data$model$plot_diff_cladogram,extra))
  }
  if (!is.null(extra$use_number)) {
    if (is.null(top_n)) top_n <- extra$use_number
    extra$use_number <- NULL
  }
  width <- if (is.null(extra$width)) 0.75 else extra$width
  extra$width <- NULL
  if (length(extra)) warning("Unused LDA plotting arguments: ",
                             paste(names(extra),collapse=", "))
  check_n <- function(x,nm) {
    if (is.null(x)) return(NULL)
    if (!is.numeric(x) || length(x)!=1L || !is.finite(x) ||
        x<1 || x!=as.integer(x)) stop(nm," must be a positive integer.")
    as.integer(x)
  }
  top_n <- check_n(top_n,"top_n")
  top_n_per_group <- check_n(top_n_per_group,"top_n_per_group")
  df <- as.data.frame(data$res_diff)
  required <- c("Taxa","Group","LDA")
  if (!all(required %in% names(df)))
    stop("LEfSe res_diff must contain Taxa, Group, and LDA. Found: ",
         paste(names(df),collapse=", "))
  df$Taxa <- as.character(df$Taxa)
  df$Group <- as.character(df$Group)
  df$LDA <- suppressWarnings(as.numeric(as.character(df$LDA)))
  df <- df[!is.na(df$Taxa) & nzchar(df$Taxa) &
             !is.na(df$Group) & is.finite(df$LDA) &
             abs(df$LDA)>=threshold,,drop=FALSE]
  if (!is.null(tax_level)) {
    rank_prefix <- c(Kingdom="k",Phylum="p",Class="c",Order="o",
                     Family="f",Genus="g",Species="s")
    if (length(tax_level)!=1L || !tax_level %in% names(rank_prefix))
      stop("tax_level must be one of: ",paste(names(rank_prefix),collapse=", "))
    # Taxa may be a single rank token or a lineage separated by | or ;.
    # Inspect the terminal rank, not an ancestor's prefix.
    terminal <- trimws(sub("^.*[|;]","",df$Taxa))
    df <- df[grepl(paste0("^",rank_prefix[[tax_level]],"__"),terminal),
             ,drop=FALSE]
  }
  if (!nrow(df)) stop("No LEfSe features satisfy the plotting filters.")
  if (is.null(group_order)) group_order <- unique(df$Group)
  group_order <- unique(c(group_order,setdiff(unique(df$Group),group_order)))
  df$Group <- factor(df$Group,levels=group_order)
  df <- df[order(-abs(df$LDA),df$Taxa),,drop=FALSE]
  if (!is.null(top_n_per_group)) {
    df <- do.call(rbind,lapply(split(df,df$Group,drop=TRUE),
                              function(x) utils::head(x,top_n_per_group)))
  }
  if (!is.null(top_n)) df <- utils::head(df[order(-abs(df$LDA),df$Taxa),,drop=FALSE],top_n)
  if (!nrow(df)) stop("No LEfSe features remain.")
  if (isTRUE(short_labels)) {
    label <- trimws(sub("^.*[|;]","",df$Taxa))
    label <- sub("^[a-z]__", "",label)
    label[!nzchar(label)] <- df$Taxa[!nzchar(label)]
  } else label <- df$Taxa
  # Distinct feature identifiers prevent duplicate labels from collapsing bars.
  df$.feature <- paste0("feature_",seq_len(nrow(df)))
  df$.label <- make.unique(label,sep="_")
  # Display each group as a contiguous block, strongest first.
  df <- df[order(df$Group,-abs(df$LDA)),,drop=FALSE]
  df$.feature <- factor(df$.feature,levels=rev(as.character(df$.feature)))
  breaks <- as.character(df$.feature)
  labels <- df$.label
  names(labels) <- breaks
  p <- ggplot2::ggplot(df,ggplot2::aes(x=.data$LDA,y=.data$.feature,
                                       fill=.data$Group))+
    ggplot2::geom_col(width=width)+
    ggplot2::scale_y_discrete(breaks=breaks,labels=labels)+
    ggplot2::labs(x="LDA score",y=NULL,fill="Group")+
    ggplot2::theme_bw(base_size=11)+
    ggplot2::theme(panel.grid.major.y=ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_text(size=9))
  if (!is.null(colors)) {
    if (is.null(names(colors))) {
      if (length(colors)<length(group_order))
        stop("colors must have at least one value per group.")
      names(colors) <- group_order[seq_along(group_order)]
    }
    missing <- setdiff(as.character(unique(df$Group)),names(colors))
    if (length(missing)) stop("Missing colors for groups: ",
                              paste(missing,collapse=", "))
    p <- p+ggplot2::scale_fill_manual(values=colors,drop=TRUE)
  }
  p
}

#' Export microeco LEfSe results
#' @param data aq_lefse result.
#' @param file CSV output path.
#' @return Invisibly returns the CSV output path.
#' @details Exports the differential-feature table stored in `data$res_diff`.
#' @examples
#' \dontrun{
#' aq_lefse_export(lefse_result, file = "LEfSe_results.csv")
#' }
#' @export
aq_lefse_export <- function(data,file="LEfSe_results.csv") {
  if (!inherits(data,"aq_lefse")) stop("data must be aq_lefse.")
  utils::write.csv(data$res_diff,file,row.names=FALSE)
  invisible(file)
}

#' @export
print.aq_lefse <- function(x,...) {
  cat("AqueCo LEfSe (microeco)\n",
      "Group:",x$group,"\n",
      "Differential features:",nrow(x$res_diff),"\n")
  invisible(x)
}
