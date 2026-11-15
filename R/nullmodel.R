
.aq_null_matrix <- function(x) {
  if (!inherits(x, "aq_null")) stop("Expected aq_null object.")
  x$otu
}
.aq_null_seed <- function(seed) {
  if (!is.null(seed)) set.seed(seed)
}
.aq_null_dist <- function(x, ids) {
  x[ids, ids, drop=FALSE]
}
#' Prepare AqueCo null-model data
#' @param data aq_data object.
#' @param filter_thres Minimum mean relative abundance, strictly greater.
#' @param group Metadata column containing treatment labels.
#' @param seed Random seed.
#' @param tree Optional phylo tree overriding data$tree.
#' @return aq_null object.
#' @export
aq_null <- function(data, filter_thres=0.0005, group="Group", seed=123, tree=NULL, cache=FALSE) {
  if (!inherits(data,"aq_data")) stop("data must be aq_data.")
  otu <- data$otu
  if (any(rowSums(otu)<=0)) stop("Samples with zero total abundance.")
  if (!is.numeric(filter_thres) || length(filter_thres)!=1L ||
      !is.finite(filter_thres) || filter_thres<0) stop("Invalid filter_thres.")
  rel <- sweep(otu,1,rowSums(otu),"/")
  keep <- colMeans(rel)>filter_thres
  otu <- otu[,keep,drop=FALSE]
  if (ncol(otu)<2L) stop("Fewer than two taxa after filtering.")
  if (any(rowSums(otu)==0)) stop("Filtering produced empty samples.")
  tr <- if (is.null(tree)) data$tree else tree
  if (!is.null(tr)) {
    if (!requireNamespace("ape",quietly=TRUE)) stop("Install ape for phylogenetic analyses.")
    if (!inherits(tr,"phylo")) stop("tree must be a phylo object.")
    if (anyDuplicated(tr$tip.label)) stop("Duplicate tree tip labels.")
    shared <- intersect(colnames(otu),tr$tip.label)
    if (length(shared)<2) stop("Fewer than two shared taxa with tree.")
    otu <- otu[,shared,drop=FALSE]
    tr <- ape::keep.tip(tr,shared)
    if (is.null(tr$edge.length) || any(!is.finite(tr$edge.length)) ||
        any(tr$edge.length<0)) stop("Tree requires finite nonnegative edge lengths.")
    if (any(rowSums(otu)==0)) stop("Tree matching produced empty samples.")
  }
  meta <- data$metadata
  if (!is.null(group) && (is.null(meta) || !group %in% names(meta)))
    stop("Grouping column missing in metadata.")
  phylo_dist <- if (isTRUE(cache) && !is.null(tr)) {
    d <- ape::cophenetic.phylo(tr)
    d[colnames(otu),colnames(otu),drop=FALSE]
  } else NULL
  structure(list(otu=otu, tree=tr, phylo_dist=phylo_dist,
    metadata=meta, group=group,
    seed=seed, filter_thres=filter_thres,
    input_taxa=ncol(data$otu), retained_taxa=ncol(otu)),class="aq_null")
}
.aq_null_map <- function(runs, workers, seed, fun) {
  if (length(workers)!=1L || !is.finite(workers) || workers<1 ||
      workers!=as.integer(workers)) stop("workers must be a positive integer.")
  if (!is.null(seed)) {
    set.seed(seed)
    seeds <- sample.int(.Machine$integer.max,runs,replace=TRUE)
  } else seeds <- rep(NA_integer_,runs)
  task <- function(k) {
    if (!is.na(seeds[k])) set.seed(seeds[k])
    fun(k)
  }
  if (workers>1L && .Platform$OS.type!="windows") {
    parallel::mclapply(seq_len(runs),task,mc.cores=workers,mc.set.seed=FALSE)
  } else {
    if (workers>1L) warning("Parallel workers currently require Unix/macOS; using sequential execution on Windows.")
    lapply(seq_len(runs),task)
  }
}
.aq_phylo_distance <- function(null, runs, weighted, model, method, seed, workers=1) {
  if (is.null(null$tree)) stop("A phylogenetic tree is required.")
  if (!identical(model,"taxa.labels")) stop("Only taxa.labels is supported.")
  if (runs<2 || runs!=as.integer(runs)) stop("runs must be integer >= 2.")
  if (length(workers)!=1L || !is.finite(workers) ||
      workers<1 || workers!=as.integer(workers)) stop("Invalid workers.")
  if (!requireNamespace("picante",quietly=TRUE)) stop("Install picante.")
  if (!requireNamespace("ape",quietly=TRUE)) stop("Install ape.")
  .aq_null_seed(seed)
  comm <- null$otu
  d <- if (!is.null(null$phylo_dist)) null$phylo_dist else
    ape::cophenetic.phylo(null$tree)
  d <- d[colnames(comm),colnames(comm),drop=FALSE]
  fn <- if (method=="betampd") picante::comdist else picante::comdistnt
  calculate <- function(dist) as.matrix(fn(comm,dist,abundance.weighted=weighted))
  observed <- calculate(d)
  n <- nrow(observed)
  sum_null <- matrix(0,n,n,dimnames=dimnames(observed))
  sumsq <- sum_null
  # Use chunked reductions to avoid retaining every full null matrix.
  chunk_size <- if (workers>1L) max(1L,ceiling(runs/(workers*4L))) else min(50L,runs)
  if (!is.null(seed)) {
    set.seed(seed)
    seeds <- sample.int(.Machine$integer.max,runs,replace=TRUE)
  } else seeds <- rep(NA_integer_,runs)
  calc_chunk <- function(indices) {
    local_sum <- matrix(0,n,n)
    local_sq <- matrix(0,n,n)
    for (k in indices) {
      if (!is.na(seeds[k])) set.seed(seeds[k])
      shuffled <- sample(rownames(d))
      perm <- d
      dimnames(perm) <- list(shuffled,shuffled)
      perm <- perm[colnames(comm),colnames(comm),drop=FALSE]
      m <- calculate(perm)
      local_sum <- local_sum+m
      local_sq <- local_sq+m*m
    }
    list(sum=local_sum,sq=local_sq)
  }
  chunks <- split(seq_len(runs),ceiling(seq_len(runs)/chunk_size))
  if (workers>1L && .Platform$OS.type!="windows") {
    parts <- parallel::mclapply(chunks,calc_chunk,mc.cores=workers,
                               mc.set.seed=FALSE)
  } else {
    if (workers>1L) warning("Parallel workers currently require Unix/macOS; using sequential execution on Windows.")
    parts <- lapply(chunks,calc_chunk)
  }
  for (part in parts) {
    sum_null <- sum_null+part$sum
    sumsq <- sumsq+part$sq
  }
  mu <- sum_null/runs
  sdv <- sqrt(pmax(0,(sumsq-runs*mu*mu)/(runs-1)))
  ses <- (observed-mu)/sdv
  ses[sdv<1e-12] <- NA_real_
  if (method=="betampd") ses <- -ses # betaNRI = -SES(betaMPD)
  diag(ses) <- NA_real_
  structure(list(matrix=ses, observed=observed, null_mean=mu,
    null_sd=sdv, metadata=null$metadata, group=null$group,
    measure=if(method=="betampd") "betaNRI" else "betaNTI",
    runs=runs, seed=seed, null_model=model),class="aq_null_result")
}
#' Calculate betaNRI (negative SES of betaMPD)
#' @param data aq_null object.
#' @param runs Number of permutations.
#' @param abundance.weighted Whether to weight by abundance.
#' @param null.model Only taxa.labels.
#' @param seed Random seed; defaults to the seed in data.
#' @export
aq_beta_nri <- function(data,runs=500,abundance.weighted=TRUE,
                        null.model="taxa.labels",seed=data$seed,workers=1) {
  .aq_phylo_distance(data,runs,abundance.weighted,null.model,"betampd",seed,workers)
}
#' Calculate betaNTI (SES of betaMNTD)
#' @inheritParams aq_beta_nri
#' @export
aq_beta_nti <- function(data,runs=500,abundance.weighted=TRUE,
                        null.model="taxa.labels",seed=data$seed,workers=1) {
  .aq_phylo_distance(data,runs,abundance.weighted,null.model,"betamntd",seed,workers)
}
#' Calculate RCbray with a fixed-richness and abundance-informed null model
#' @param data aq_null object.
#' @param runs Number of randomizations.
#' @param seed Random seed.
#' @return aq_null_result; RCbray in [-1,1].
#' @export
aq_rcbray <- function(data,runs=1000,seed=data$seed,workers=1) {
  comm <- .aq_null_matrix(data)
  if (runs<2 || runs!=as.integer(runs)) stop("runs must be integer >= 2.")
  if (length(workers)!=1L || !is.finite(workers) ||
      workers<1 || workers!=as.integer(workers)) stop("Invalid workers.")
  if (!requireNamespace("vegan",quietly=TRUE)) stop("Install vegan.")
  .aq_null_seed(seed)
  observed <- as.matrix(vegan::vegdist(comm,method="bray"))
  n <- nrow(comm); count_less <- matrix(0,n,n)
  count_equal <- matrix(0,n,n)
  freq <- colSums(comm>0); prob <- freq/sum(freq)
  pool <- colSums(comm); pool <- pool/sum(pool)
  # Retains each sample's richness and read total. Taxa selection is
  # weighted by occurrence frequency; reads follow pooled abundance.
  one_iter <- function(iter) {
    sim <- matrix(0,n,ncol(comm))
    for (i in seq_len(n)) {
      richness <- sum(comm[i,]>0)
      ids <- sample.int(ncol(comm),richness,replace=FALSE,prob=prob)
      total <- round(sum(comm[i,]))
      if (total<richness) stop("Noninteger or insufficient counts; RCbray requires raw counts.")
      sim[i,ids] <- 1
      if (total>richness)
        sim[i,ids] <- sim[i,ids]+as.vector(rmultinom(1,total-richness,prob=pool[ids]))
    }
    d <- as.matrix(vegan::vegdist(sim,method="bray"))
    list(less=(d<observed),equal=(abs(d-observed)<1e-12))
  }
  # Chunk work to reduce parallel scheduling overhead.
  if (!is.null(seed)) {
    set.seed(seed)
    seeds <- sample.int(.Machine$integer.max,runs,replace=TRUE)
  } else seeds <- rep(NA_integer_,runs)
  chunk_size <- if(workers>1L) max(1L,ceiling(runs/(workers*4L))) else min(50L,runs)
  chunks <- split(seq_len(runs),ceiling(seq_len(runs)/chunk_size))
  run_chunk <- function(indices) {
    less <- matrix(0,n,n); equal <- matrix(0,n,n)
    for (iter in indices) {
      if (!is.na(seeds[iter])) set.seed(seeds[iter])
      z <- one_iter(iter)
      less <- less+z$less; equal <- equal+z$equal
    }
    list(less=less,equal=equal)
  }
  if (workers>1L && .Platform$OS.type!="windows") {
    parts <- parallel::mclapply(chunks,run_chunk,mc.cores=workers,
                               mc.set.seed=FALSE)
  } else {
    if (workers>1L) warning("Parallel workers currently require Unix/macOS; using sequential execution on Windows.")
    parts <- lapply(chunks,run_chunk)
  }
  for (part in parts) {
    count_less <- count_less+part$less
    count_equal <- count_equal+part$equal
  }
  rc <- 2*(count_less+0.5*count_equal)/runs-1
  dimnames(rc) <- dimnames(observed)
  diag(rc) <- NA_real_
  structure(list(matrix=rc,observed=observed,metadata=data$metadata,
    group=data$group,measure="RCbray",runs=runs,seed=seed,
    null_model="frequency-weighted richness-preserving multinomial"),
    class="aq_null_result")
}
#' Infer ecological assembly processes from betaNTI and RCbray
#' @param betanti aq_beta_nti result.
#' @param rcbray aq_rcbray result.
#' @param group Metadata group column.
#' @param comparison within, between, or all.
#' @export
aq_assembly <- function(betanti,rcbray,group="Group",comparison="within") {
  comparison <- match.arg(comparison,c("within","between","all"))
  if (!identical(rownames(betanti$matrix),rownames(rcbray$matrix)))
    stop("Sample order mismatch.")
  ids <- rownames(betanti$matrix)
  meta <- betanti$metadata
  if (is.null(meta) || !group %in% names(meta)) stop("Group metadata missing.")
  ix <- which(upper.tri(betanti$matrix),arr.ind=TRUE)
  a <- betanti$matrix[ix]; b <- rcbray$matrix[ix]
  g1 <- as.character(meta[ids[ix[,1]],group])
  g2 <- as.character(meta[ids[ix[,2]],group])
  within <- g1==g2
  take <- switch(comparison,within=within,between=!within,all=rep(TRUE,length(within)))
  process <- ifelse(a>2,"Heterogeneous Selection",
    ifelse(a< -2,"Homogeneous Selection",
    ifelse(b>0.95,"Dispersal Limitation",
    ifelse(b< -0.95,"Homogenizing Dispersal","Undominated"))))
  pair_group <- ifelse(within,g1,paste(pmin(g1,g2),pmax(g1,g2),sep=" vs "))
  pairs <- data.frame(Sample1=ids[ix[,1]],Sample2=ids[ix[,2]],
    Group=pair_group,betaNTI=a,RCbray=b,process=process,
    stringsAsFactors=FALSE)
  pairs <- pairs[take & is.finite(a) & is.finite(b) & !is.na(g1) & !is.na(g2),,drop=FALSE]
  levels <- c("Heterogeneous Selection","Homogeneous Selection",
    "Dispersal Limitation","Homogenizing Dispersal","Undominated")
  if (nrow(pairs)) {
    tab <- as.data.frame(table(Group=pairs$Group,
      process=factor(pairs$process,levels=levels)),stringsAsFactors=FALSE)
    names(tab)[3] <- "Count"
    tab$percentage <- stats::ave(tab$Count,tab$Group,FUN=function(x) if(sum(x)) 100*x/sum(x) else 0)
    tab <- tab[tab$Count>0 | tab$Group %in% unique(pairs$Group),,drop=FALSE]
  } else tab <- data.frame(Group=character(),process=character(),
                            Count=integer(),percentage=numeric())
  structure(list(summary=tab,pairs=pairs,comparison=comparison),class="aq_assembly")
}
#' Summarize within- and between-group null-model distances
#' @param data aq_null_result.
#' @param group Metadata group column.
#' @param comparison within, between, all.
#' @export
aq_null_distance <- function(data,group="Group",comparison="all") {
  comparison <- match.arg(comparison,c("within","between","all"))
  ids <- rownames(data$matrix)
  meta <- data$metadata
  if (is.null(meta) || !group %in% names(meta)) stop("Group column missing.")
  ix <- which(upper.tri(data$matrix),arr.ind=TRUE)
  a <- as.character(meta[ids[ix[,1]],group])
  b <- as.character(meta[ids[ix[,2]],group])
  same <- a==b
  keep <- switch(comparison,within=same,between=!same,all=rep(TRUE,length(same)))
  out <- data.frame(Sample1=ids[ix[,1]],Sample2=ids[ix[,2]],
    Group=ifelse(same,a,paste(pmin(a,b),pmax(a,b),sep=" vs ")),
    value=data$matrix[ix],stringsAsFactors=FALSE)
  out[keep & is.finite(out$value) & !is.na(a) & !is.na(b),,drop=FALSE]
}
#' Pairwise Wilcoxon comparison (exploratory; pairs are not independent)
#' @param data aq_null_distance output or aq_null_result.
#' @param group Group column if data is aq_null_result.
#' @export
aq_null_test <- function(data,group="Group") {
  if (inherits(data,"aq_null_result")) data <- aq_null_distance(data,group,"within")
  if (!all(c("Group","value") %in% names(data))) stop("Expected Group and value.")
  groups <- unique(data$Group)
  if (length(groups)<2) stop("At least two groups required.")
  comb <- combn(groups,2,simplify=FALSE)
  out <- do.call(rbind,lapply(comb,function(x) {
    p <- tryCatch(stats::wilcox.test(data$value[data$Group==x[1]],
      data$value[data$Group==x[2]],exact=FALSE)$p.value,error=function(e) NA_real_)
    data.frame(Group1=x[1],Group2=x[2],p=p)
  }))
  out$p_adjust <- stats::p.adjust(out$p,method="BH")
  warning("Sample pairs share samples; Wilcoxon p-values are exploratory, not independent-pair inference.")
  out
}
#' Plot null-model group comparisons or heatmap
#' @param data aq_null_result.
#' @param type boxplot or heatmap.
#' @param group Group column.
#' @param add mean or none.
#' @param threshold Optional horizontal thresholds.
#' @param test Optional exploratory wilcox; use aq_null_test separately.
#' @param cluster Whether to cluster heatmap.
#' @export
aq_null_plot <- function(data,type="boxplot",group="Group",add="mean",
                         threshold=NULL,test=NULL,cluster=FALSE) {
  type <- match.arg(type,c("boxplot","heatmap"))
  if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2.")
  if (type=="boxplot") {
    d <- aq_null_distance(data,group,"within")
    p <- ggplot2::ggplot(d,ggplot2::aes(x=Group,y=value,fill=Group))+
      ggplot2::geom_boxplot(outlier.shape=NA,alpha=.6)+
      ggplot2::geom_jitter(width=.12,height=0,alpha=.55,size=1.6)+
      ggplot2::theme_bw()+ggplot2::labs(y=data$measure,x=NULL)
    if (identical(add,"mean")) p <- p+ggplot2::stat_summary(fun=mean,geom="point",shape=4,size=3)
    if (length(threshold)) for (v in threshold)
      p <- p+ggplot2::geom_hline(yintercept=v,linetype=2)
    if (!is.null(test) && test!="wilcox") stop("Only exploratory wilcox is supported.")
    return(p)
  }
  m <- data$matrix
  diag(m) <- 0
  if (cluster) {
    ord <- stats::hclust(stats::as.dist(replace(m,!is.finite(m),0)))$labels
    ord <- rownames(m)[stats::hclust(stats::as.dist(replace(m,!is.finite(m),0)))$order]
    m <- m[ord,ord,drop=FALSE]
  }
  d <- as.data.frame(as.table(m),stringsAsFactors=FALSE)
  names(d) <- c("Sample1","Sample2","value")
  d$Sample1 <- factor(d$Sample1,levels=rownames(m))
  d$Sample2 <- factor(d$Sample2,levels=rev(colnames(m)))
  ggplot2::ggplot(d,ggplot2::aes(Sample1,Sample2,fill=value))+
    ggplot2::geom_tile()+ggplot2::scale_fill_gradient2(low="#4477AA",
      mid="white",high="#C75B50",midpoint=0,na.value="grey90")+
    ggplot2::coord_equal()+ggplot2::theme_bw()+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1))+
    ggplot2::labs(x=NULL,y=NULL,fill=data$measure)
}
#' Plot ecological assembly composition
#' @param data aq_assembly result.
#' @param type stacked.
#' @param palette publication.
#' @export
aq_assembly_plot <- function(data,type="stacked",palette="publication") {
  if (!inherits(data,"aq_assembly")) stop("Expected aq_assembly.")
  if (type!="stacked") stop("Only stacked plots supported.")
  p <- ggplot2::ggplot(data$summary,
    ggplot2::aes(x=Group,y=percentage,fill=process))+
    ggplot2::geom_col(width=.7)+ggplot2::theme_bw()+
    ggplot2::labs(x=NULL,y="Ecological Processes (%)",fill="Process")
  if (palette=="publication") p <- p+ggplot2::scale_fill_manual(
    values=c("Heterogeneous Selection"="#D66A57",
      "Homogeneous Selection"="#E9B56C",
      "Dispersal Limitation"="#5B8DB8",
      "Homogenizing Dispersal"="#65A58A","Undominated"="#A8A5B5"),
    drop=FALSE)
  p
}
#' Export null model and assembly results
#' @param data aq_null_result or aq_assembly.
#' @param file File prefix.
#' @export
aq_null_export <- function(data,file="AqueCo_Null") {
  if (inherits(data,"aq_null_result")) {
    utils::write.csv(data$matrix,paste0(file,"_",data$measure,"_matrix.csv"))
    utils::write.csv(aq_null_distance(data),paste0(file,"_",data$measure,"_pairs.csv"),row.names=FALSE)
  } else if (inherits(data,"aq_assembly")) {
    utils::write.csv(data$summary,paste0(file,"_assembly_summary.csv"),row.names=FALSE)
    utils::write.csv(data$pairs,paste0(file,"_assembly_pairs.csv"),row.names=FALSE)
  } else stop("Unsupported result type.")
  invisible(TRUE)
}
