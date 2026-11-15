#' Distance-decay analysis of Bray-Curtis dissimilarity
#'
#' Computes pairwise Bray-Curtis dissimilarity against geographic (haversine km)
#' or standardized environmental Euclidean distance. Regression R2 describes
#' the displayed linear fit; Mantel permutation P accounts for dependent pairs.
#' @param data aq_data object (samples x taxa).
#' @param distance geographic, environmental, or both.
#' @param group Metadata grouping column; NULL for overall only.
#' @param groups Optional subset of group names; NULL selects all.
#' @param scope overall, group, or both.
#' @param lon Longitude column name in metadata.
#' @param lat Latitude column name in metadata.
#' @param env_variables Numeric environmental columns to include; NULL all
#'   numeric env columns excluding coordinate columns.
#' @param standardize Standardize environmental variables to unit variance.
#' @param permutations Mantel permutations (at least 1).
#' @param seed Random seed for reproducible Mantel tests.
#' @return aq_ddr with pairs, summary, and settings.
#' @export
aq_ddr <- function(data, distance=c("geographic","environmental"),
                   group=NULL, groups=NULL, scope=c("overall","group","both"),
                   lon="Longitude",lat="Latitude",env_variables=NULL,
                   standardize=TRUE,permutations=999,seed=123) {
  if (!inherits(data,"aq_data")) stop("data must be aq_data.")
  if (!requireNamespace("vegan",quietly=TRUE)) stop("Install vegan.")
  distance <- match.arg(distance,c("geographic","environmental"),several.ok=TRUE)
  scope <- match.arg(scope)
  if (is.null(group) && scope!="overall")
    stop("Provide group for scope='group' or 'both'.")
  if (length(permutations)!=1L || !is.finite(permutations) ||
      permutations<1 || permutations!=as.integer(permutations))
    stop("permutations must be a positive integer.")
  otu <- as.matrix(data$otu)
  if (any(rowSums(otu)<=0)) stop("Samples with zero abundance are not supported.")
  samples <- rownames(otu)
  md <- data$metadata
  if (!is.null(group)) {
    if (is.null(md) || !group %in% names(md)) stop("Grouping column not found.")
    grp <- as.character(md[samples,group])
    if (anyNA(grp)) stop("Missing group values.")
    if (!is.null(groups)) {
      if (!all(groups %in% grp)) stop("Unknown groups: ",
                                     paste(setdiff(groups,grp),collapse=", "))
      keep <- grp %in% groups
      otu <- otu[keep,,drop=FALSE]; samples <- rownames(otu); grp <- grp[keep]
    }
  } else grp <- rep("Overall",length(samples))
  if (nrow(otu)<3L) stop("At least three samples required.")
  # Drop taxa absent in the selected samples.
  otu <- otu[,colSums(otu)>0,drop=FALSE]
  if (!ncol(otu)) stop("No taxa available.")
  bc <- as.matrix(vegan::vegdist(otu,method="bray"))
  all_dist <- list()
  if ("geographic" %in% distance) {
    if (is.null(md) || !all(c(lon,lat) %in% names(md)))
      stop("Longitude/latitude columns must be in metadata: ",lon,", ",lat)
    coords <- md[samples,c(lon,lat),drop=FALSE]
    xy <- suppressWarnings(sapply(coords,as.numeric))
    xy <- matrix(xy,ncol=2,dimnames=list(samples,c(lon,lat)))
    if (any(!is.finite(xy)) || any(abs(xy[,1])>180) || any(abs(xy[,2])>90))
      stop("Invalid or missing longitude/latitude values.")
    rad <- pi/180
    lo <- xy[,1]*rad; la <- xy[,2]*rad
    dlon <- outer(lo,lo,"-"); dlat <- outer(la,la,"-")
    h <- sin(dlat/2)^2 + outer(cos(la),cos(la))*sin(dlon/2)^2
    geo <- matrix(6371.0088*2*asin(sqrt(pmin(1,pmax(0,h)))),
                  nrow=length(samples),ncol=length(samples),
                  dimnames=list(samples,samples))
    diag(geo) <- 0
    all_dist$Geographic <- geo
  }
  if ("environmental" %in% distance) {
    if (is.null(data$env)) stop("Environmental data are required.")
    env <- data$env[samples,,drop=FALSE]
    if (is.null(env_variables)) {
      env_variables <- names(env)[vapply(env,is.numeric,logical(1))]
      env_variables <- setdiff(env_variables,c(lon,lat))
    }
    if (!length(env_variables) || !all(env_variables %in% names(env)))
      stop("No valid environmental variables selected.")
    env <- env[,env_variables,drop=FALSE]
    if (!all(vapply(env,is.numeric,logical(1))) ||
        any(!is.finite(as.matrix(env))))
      stop("Environmental variables must be finite numeric values.")
    varying <- vapply(env,stats::sd,numeric(1))>0
    if (!all(varying))
      warning("Dropping constant environmental variables: ",
              paste(names(env)[!varying],collapse=", "))
    env <- env[,varying,drop=FALSE]
    if (!ncol(env)) stop("No varying environmental variables.")
    mat <- as.matrix(env)
    if (standardize) mat <- scale(mat)
    all_dist$Environmental <- as.matrix(stats::dist(mat,method="euclidean"))
  }
  groups_to_run <- if (scope=="overall") list(Overall=samples) else
    if (scope=="group") split(samples,grp) else
      c(list(Overall=samples),split(samples,grp))
  results <- list(); summaries <- list()
  if (!is.null(seed)) set.seed(seed)
  for (nm in names(groups_to_run)) {
    ids <- groups_to_run[[nm]]
    if (length(ids)<3L) {
      warning("Skipping ",nm,": fewer than three samples.")
      next
    }
    for (type in names(all_dist)) {
      dx <- all_dist[[type]][ids,ids,drop=FALSE]
      dy <- bc[ids,ids,drop=FALSE]
      ix <- lower.tri(dx)
      x <- dx[ix]; y <- dy[ix]
      fit <- if (length(unique(x))>1L) stats::lm(y~x) else NULL
      if (is.null(fit)) {
        warning("Skipping ",nm," ",type,": all predictor distances equal.")
        next
      }
      r2 <- summary(fit)$r.squared
      slope <- unname(stats::coef(fit)[2])
      # Mantel p is a matrix-permutation p-value, not an lm p-value.
      mt <- vegan::mantel(stats::as.dist(dx),stats::as.dist(dy),
                          method="pearson",permutations=as.integer(permutations))
      key <- paste(nm,type,sep="::")
      results[[key]] <- data.frame(Group=nm,Distance=type,
        Sample1=rownames(dx)[row(dx)[ix]],Sample2=colnames(dx)[col(dx)[ix]],
        Predictor=x,BC=y)
      summaries[[key]] <- data.frame(Group=nm,Distance=type,
        Samples=length(ids),Pairs=length(x),R2=r2,Slope=slope,
        Mantel_r=unname(mt$statistic),Mantel_p=mt$signif,
        Permutations=as.integer(permutations))
    }
  }
  if (!length(summaries)) stop("No valid group-distance combinations.")
  structure(list(pairs=do.call(rbind,results),
    summary=do.call(rbind,summaries),
    parameters=list(distance=distance,group=group,groups=groups,scope=scope,
      lon=lon,lat=lat,env_variables=env_variables,
      standardize=standardize,permutations=permutations,seed=seed)),
    class="aq_ddr")
}

#' Plot distance-decay relationships
#' @param data aq_ddr result.
#' @param distance geographic or environmental; NULL uses all.
#' @param groups Optional subset of groups to display.
#' @param compare Overlay curves with shared axes; FALSE facets by group.
#' @param show_ci Show 95 percent lm confidence ribbons.
#' @param show_points Display pairwise points.
#' @return A ggplot object.
#' @details Confidence ribbons are descriptive linear-model intervals, not
#'   confidence intervals for the Mantel permutation test.
#' @examples
#' \dontrun{
#' aq_ddr_plot(ddr_result)
#' }
#' @export
aq_ddr_plot <- function(data,distance=NULL,groups=NULL,compare=TRUE,
                        show_ci=TRUE,show_points=TRUE) {
  if (!inherits(data,"aq_ddr")) stop("data must be aq_ddr.")
  df <- data$pairs; sm <- data$summary
  if (!is.null(distance)) {
    distance <- match.arg(distance,c("geographic","environmental"),several.ok=TRUE)
    wanted <- ifelse(distance=="geographic","Geographic","Environmental")
    df <- df[df$Distance %in% wanted,,drop=FALSE]
    sm <- sm[sm$Distance %in% wanted,,drop=FALSE]
  }
  if (!is.null(groups)) {
    df <- df[df$Group %in% groups,,drop=FALSE]
    sm <- sm[sm$Group %in% groups,,drop=FALSE]
  }
  if (!nrow(df)) stop("No distance-decay data for selected filters.")
  sm$Label <- sprintf("%s: R\u00b2 = %.3f, Mantel p = %s",sm$Group,sm$R2,
    ifelse(sm$Mantel_p<0.001,"< 0.001",sprintf("%.3f",sm$Mantel_p)))
  df$Group <- factor(df$Group,levels=unique(sm$Group))
  sm$Group <- factor(sm$Group,levels=levels(df$Group))
  p <- ggplot2::ggplot(df,ggplot2::aes(x=.data$Predictor,y=.data$BC,
                                       color=.data$Group))
  if (show_points) p <- p+ggplot2::geom_point(alpha=.22,size=1)
  if (show_ci) p <- p+ggplot2::geom_smooth(method="lm",formula=y~x,
    se=TRUE,alpha=.12,linewidth=.85)
  else p <- p+ggplot2::geom_smooth(method="lm",formula=y~x,
    se=FALSE,linewidth=.85)
  p <- p+ggplot2::labs(x="Distance",y="Bray-Curtis dissimilarity",
    color="Group",caption="R\u00b2: linear fit; p: Mantel matrix permutation test")+
    ggplot2::theme_bw(base_size=12)+
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank())
  if (compare) {
    sm <- sm[order(sm$Distance,sm$Group),,drop=FALSE]
    sm$.rank <- stats::ave(seq_len(nrow(sm)),sm$Distance,FUN=seq_along)
    sm$.y <- 1.02 - (sm$.rank-1)*0.065
    p <- p+ggplot2::facet_wrap(~Distance,scales="free_x")+
      ggplot2::geom_text(data=sm,
        mapping=ggplot2::aes(x=-Inf,y=.data$.y,label=.data$Label,
                              color=.data$Group),
        inherit.aes=FALSE,hjust=-0.02,vjust=1,size=3.1,
        show.legend=FALSE)+
      ggplot2::coord_cartesian(ylim=c(0,1.08))
  } else {
    p <- p+ggplot2::facet_grid(Group~Distance,scales="free_x")+
      ggplot2::geom_text(data=sm,
        mapping=ggplot2::aes(x=-Inf,y=Inf,label=.data$Label),
        inherit.aes=FALSE,hjust=-0.02,vjust=1.2,size=3.1)
  }
  p
}

#' Export distance-decay pairwise data and statistics
#' @param data aq_ddr result.
#' @param prefix Output filename prefix, optionally including a directory path.
#' @return Invisibly returns the two CSV output paths.
#' @details Writes `<prefix>_pairs.csv` and `<prefix>_summary.csv`.
#'   The destination directory must already exist.
#' @examples
#' \dontrun{
#' aq_ddr_export(ddr_result, prefix = "results/DDR")
#' }
#' @export
aq_ddr_export <- function(data,prefix="DDR") {
  if (!inherits(data,"aq_ddr")) stop("data must be aq_ddr.")
  utils::write.csv(data$pairs,paste0(prefix,"_pairs.csv"),row.names=FALSE)
  utils::write.csv(data$summary,paste0(prefix,"_summary.csv"),row.names=FALSE)
  invisible(c(paste0(prefix,"_pairs.csv"),paste0(prefix,"_summary.csv")))
}
