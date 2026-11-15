
# AqueCo network module: correlation networks and graph diagnostics.
.aq_net_require <- function() {
  if (!requireNamespace("igraph", quietly=TRUE))
    stop("Package 'igraph' is required. Install with install.packages('igraph').")
}
.aq_net_one <- function(otu, taxonomy, method, cor_cutoff, p_cutoff, p_adjust,
                        min_abundance, min_prevalence, min_samples, group_name) {
  if (nrow(otu)<min_samples)
    warning(sprintf("Network %s has only %s samples; inference may be unstable.",
                    group_name,nrow(otu)))
  if (nrow(otu)<4) stop("At least four samples are required.")
  rel <- sweep(otu,1,rowSums(otu),"/")
  keep <- colMeans(rel)>min_abundance &
    colMeans(otu>0)>=min_prevalence &
    apply(otu,2,stats::sd)>0
  rel <- rel[,keep,drop=FALSE]
  if (ncol(rel)<2) stop("Fewer than two taxa remain after filtering.")
  taxa <- colnames(rel)
  pairs <- utils::combn(taxa,2)
  rho <- p <- rep(NA_real_,ncol(pairs))
  for (k in seq_len(ncol(pairs))) {
    a <- rel[,pairs[1,k]]; b <- rel[,pairs[2,k]]
    ct <- suppressWarnings(stats::cor.test(a,b,method=method,
                exact=FALSE))
    rho[k] <- unname(ct$estimate)
    p[k] <- ct$p.value
  }
  padj <- stats::p.adjust(p,method=p_adjust)
  edges_all <- data.frame(from=pairs[1,],to=pairs[2,],
    weight=rho,p=p,p_adjust=padj,stringsAsFactors=FALSE)
  good <- is.finite(rho) & is.finite(padj) &
    abs(rho)>=cor_cutoff & padj<p_cutoff
  edges <- edges_all[good,,drop=FALSE]
  edges$sign <- ifelse(edges$weight>0,"Positive","Negative")
  nodes <- data.frame(name=taxa,stringsAsFactors=FALSE)
  if (!is.null(taxonomy)) {
    tx <- taxonomy[match(taxa,rownames(taxonomy)),,drop=FALSE]
    nodes <- cbind(nodes,tx)
  }
  .aq_net_require()
  graph <- igraph::graph_from_data_frame(edges,directed=FALSE,vertices=nodes)
  structure(list(graph=graph,nodes=nodes,edges=edges,
    edges_all=edges_all,group=group_name,method=method,
    samples=rownames(otu),n_samples=nrow(otu),
    retained_taxa=length(taxa),
    parameters=list(cor_cutoff=cor_cutoff,p_cutoff=p_cutoff,
      p_adjust=p_adjust,min_abundance=min_abundance,
      min_prevalence=min_prevalence)),class="aq_network")
}
#' Build correlation-based microbial co-occurrence networks
#' @param data aq_data object with sample-by-taxon abundance matrix.
#' @param group Metadata grouping column, or NULL for pooled network.
#' @param groups Optional group labels.
#' @param method spearman or pearson.
#' @param min_abundance Minimum mean relative abundance (strictly greater).
#' @param min_prevalence Minimum fraction of samples with nonzero counts.
#' @param cor_cutoff Minimum absolute correlation.
#' @param p_cutoff Adjusted p-value cutoff.
#' @param p_adjust Multiple testing correction method.
#' @param min_samples Threshold for small-sample warning.
#' @export
aq_network <- function(data,group="Group",groups=NULL,method="spearman",
                       min_abundance=0.0005,min_prevalence=0.3,
                       cor_cutoff=0.6,p_cutoff=0.05,p_adjust="BH",
                       min_samples=10) {
  if (!inherits(data,"aq_data")) stop("Expected aq_data.")
  method <- match.arg(method,c("spearman","pearson"))
  if (any(!is.finite(data$otu)) || any(data$otu<0) ||
      any(rowSums(data$otu)<=0)) stop("Invalid OTU counts.")
  if (any(!is.finite(c(min_abundance,min_prevalence,cor_cutoff,p_cutoff))) ||
      min_abundance<0 || min_prevalence<0 || min_prevalence>1 ||
      cor_cutoff<0 || cor_cutoff>1 || p_cutoff<=0 || p_cutoff>1)
    stop("Invalid network thresholds.")
  if (length(min_samples)!=1 || !is.finite(min_samples) || min_samples<4) stop("Invalid min_samples.")
  if (!p_adjust %in% stats::p.adjust.methods) stop("Invalid p_adjust.")
  if (is.null(group)) {
    return(.aq_net_one(data$otu,data$taxonomy,method,cor_cutoff,p_cutoff,
      p_adjust,min_abundance,min_prevalence,min_samples,"All"))
  }
  if (is.null(data$metadata) || !group %in% names(data$metadata))
    stop("Grouping column not found.")
  g <- as.character(data$metadata[rownames(data$otu),group])
  if (anyNA(g)) stop("Missing group labels.")
  if (is.null(groups)) groups <- unique(g)
  if (!all(groups %in% g)) stop("Unknown group label.")
  nets <- lapply(groups,function(z) .aq_net_one(data$otu[g==z,,drop=FALSE],
    data$taxonomy,method,cor_cutoff,p_cutoff,p_adjust,
    min_abundance,min_prevalence,min_samples,z))
  names(nets) <- groups
  structure(nets,class=c("aq_network_set","list"))
}
.aq_net_each <- function(x,fn) {
  if (inherits(x,"aq_network_set")) return(lapply(x,fn))
  if (inherits(x,"aq_network")) return(fn(x))
  stop("Expected aq_network or aq_network_set.")
}
#' Summarize network topology
#' @param data aq_network or aq_network_set.
#' @export
aq_net_topology <- function(data) {
  .aq_net_require()
  f <- function(x) {
    g <- x$graph; n <- igraph::vcount(g); e <- igraph::ecount(g)
    components <- igraph::components(g)
    largest <- igraph::induced_subgraph(g,which(components$membership==
                  which.max(components$csize)))
    data.frame(Group=x$group,Samples=x$n_samples,Nodes=n,Edges=e,
      Positive_edges=sum(x$edges$weight>0),
      Negative_edges=sum(x$edges$weight<0),
      Average_degree=if(n) mean(igraph::degree(g)) else NA_real_,
      Density=if(n>1) igraph::edge_density(g) else NA_real_,
      Clustering=if(n>2 && e) igraph::transitivity(g,type="global") else NA_real_,
      Components=components$no,
      Mean_path_largest=if(igraph::vcount(largest)>1 && igraph::ecount(largest)>0)
        mean(igraph::distances(largest, weights=NA)[upper.tri(igraph::distances(largest, weights=NA))])
        else NA_real_,stringsAsFactors=FALSE)
  }
  if (inherits(data,"aq_network_set"))
    return(do.call(rbind,lapply(data,f)))
  f(data)
}
#' Compute node-level centrality measures
#' @param data aq_network or aq_network_set.
#' @export
aq_net_node <- function(data) {
  .aq_net_require()
  f <- function(x) {
    g <- x$graph
    data.frame(Group=x$group,name=igraph::V(g)$name,
      Degree=igraph::degree(g),
      Strength=igraph::strength(g,weights=abs(igraph::E(g)$weight)),
      Betweenness=igraph::betweenness(g,directed=FALSE,weights=NA),
      Closeness=igraph::closeness(g,normalized=TRUE,weights=NA),
      stringsAsFactors=FALSE)
  }
  if (inherits(data,"aq_network_set"))
    return(do.call(rbind,lapply(data,f)))
  f(data)
}
#' Detect modules using Louvain on absolute edge weights
#' @param data aq_network or aq_network_set.
#' @param seed Seed for stochastic community detection.
#' @export
aq_net_module <- function(data,seed=123) {
  .aq_net_require()
  f <- function(x) {
    g <- x$graph
    if (!is.null(seed)) set.seed(seed)
    if (!igraph::ecount(g)) {
      membership <- seq_len(igraph::vcount(g))
      modularity <- NA_real_
    } else {
      fit <- igraph::cluster_louvain(g,weights=abs(igraph::E(g)$weight))
      membership <- igraph::membership(fit)
      modularity <- igraph::modularity(fit)
    }
    data.frame(Group=x$group,name=igraph::V(g)$name,
      Module=as.integer(membership),Modularity=modularity)
  }
  if (inherits(data,"aq_network_set"))
    return(do.call(rbind,lapply(data,f)))
  f(data)
}
#' Compute within-module Zi and participation Pi roles
#' @param data aq_network or aq_network_set.
#' @param modules Optional module table from aq_net_module.
#' @export
aq_net_zp <- function(data,modules=NULL) {
  .aq_net_require()
  f <- function(x) {
    g <- x$graph
    mod <- if (is.null(modules)) aq_net_module(x) else
      modules[modules$Group==x$group,,drop=FALSE]
    m <- mod$Module[match(igraph::V(g)$name,mod$name)]
    if (anyNA(m)) stop("Module membership incomplete.")
    adj <- as.matrix(igraph::as_adjacency_matrix(g,sparse=FALSE))
    ki <- rowSums(adj)
    kin <- vapply(seq_along(m),function(i) sum(adj[i,m==m[i]]),numeric(1))
    zi <- numeric(length(m))
    for (z in unique(m)) {
      ids <- which(m==z)
      s <- stats::sd(kin[ids])
      zi[ids] <- if(length(ids)>1 && is.finite(s) && s>0)
        (kin[ids]-mean(kin[ids]))/s else 0
    }
    pi <- vapply(seq_along(m),function(i) {
      if (ki[i]==0) return(0)
      1-sum(vapply(unique(m),function(z)
        (sum(adj[i,m==z])/ki[i])^2,numeric(1)))
    },numeric(1))
    role <- ifelse(zi>2.5 & pi>0.62,"Network hub",
            ifelse(zi>2.5,"Module hub",
            ifelse(pi>0.62,"Connector","Peripheral")))
    data.frame(Group=x$group,name=igraph::V(g)$name,
      Module=m,Zi=zi,Pi=pi,Role=role)
  }
  if (inherits(data,"aq_network_set"))
    return(do.call(rbind,lapply(data,f)))
  f(data)
}
#' Identify candidate hub/connector taxa from Zi-Pi
#' @param data aq_network or aq_network_set.
#' @export
aq_net_keystone <- function(data) {
  z <- aq_net_zp(data)
  z[z$Role!="Peripheral",,drop=FALSE]
}
#' Compare group network topology
#' @param data aq_network_set.
#' @export
aq_net_compare <- function(data) {
  if (!inherits(data,"aq_network_set")) stop("Expected grouped networks.")
  aq_net_topology(data)
}
#' Plot microbial networks using igraph layout and ggplot2
#' @param data Single aq_network object.
#' @param layout fr or circle.
#' @export
aq_net_plot <- function(data,layout="fr",node_color="Phylum",
                        node_size="Degree",top_n=10,
                        edge_color="sign",show_legend=TRUE,
                        remove_isolates=FALSE,seed=123) {
  .aq_net_require()
  if (!inherits(data,"aq_network")) stop("Select one group: net[['KO']].")
  layout <- match.arg(layout,c("fr","circle"))
  if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2.")
  g <- data$graph
  if (remove_isolates) {
    g <- igraph::induced_subgraph(g,igraph::V(g)[igraph::degree(g)>0])
  }
  if (igraph::vcount(g)==0) {
    return(ggplot2::ggplot()+ggplot2::theme_void()+
      ggplot2::labs(title=paste0("Network: ",data$group),
                    subtitle="No connected nodes at selected thresholds"))
  }
  if (!is.null(seed)) set.seed(seed)
  xy <- if (layout=="fr") igraph::layout_with_fr(g,
    weights=if (igraph::ecount(g)) abs(igraph::E(g)$weight) else NULL)
    else igraph::layout_in_circle(g)
  nd <- data.frame(name=igraph::V(g)$name,x=xy[,1],y=xy[,2],
                   Degree=as.numeric(igraph::degree(g)),
                   stringsAsFactors=FALSE)
  if (!node_size %in% c("Degree","Strength")) stop("node_size must be Degree or Strength.")
  nd$Strength <- as.numeric(igraph::strength(g,
    weights=if(igraph::ecount(g)) abs(igraph::E(g)$weight) else numeric(0)))
  if (node_color %in% c("Module","Role")) {
    z <- aq_net_zp(data)
    nd$category <- as.character(z[[node_color]][match(nd$name,z$name)])
  } else if (node_color=="none") {
    nd$category <- "Taxa"
  } else {
    tax <- data$nodes
    candidates <- names(tax)
    exact <- match(node_color,candidates)
    if (is.na(exact)) exact <- match(tolower(node_color),tolower(candidates))
    if (is.na(exact)) {
      if (identical(node_color,"Phylum") && length(setdiff(candidates,"name"))==0L) {
        nd$category <- "Taxa"
      } else {
        stop("Taxonomic level not found in network nodes: ",node_color,
             ". Available: ",paste(candidates,collapse=", "))
      }
    } else {
      nd$category <- as.character(tax[[exact]][match(nd$name,tax$name)])
    }
  }
  nd$category[is.na(nd$category) | trimws(nd$category)==""] <- "Unclassified"
  if (length(top_n)!=1 || !is.numeric(top_n) || is.na(top_n) ||
      top_n<1 || top_n!=as.integer(top_n)) stop("top_n must be a positive integer.")
  counts <- sort(table(nd$category),decreasing=TRUE)
  # Always retain Unclassified as an explicit category.
  leading <- utils::head(setdiff(names(counts),"Unclassified"),top_n)
  nd$category[!nd$category %in% c(leading,"Unclassified")] <- "Others"
  nd$category <- factor(nd$category,
    levels=c(leading,intersect("Unclassified",unique(nd$category)),
             intersect("Others",unique(nd$category))))
  ed <- igraph::as_data_frame(g,what="edges")
  if (nrow(ed)) {
    ed$x <- nd$x[match(ed$from,nd$name)]
    ed$y <- nd$y[match(ed$from,nd$name)]
    ed$xend <- nd$x[match(ed$to,nd$name)]
    ed$yend <- nd$y[match(ed$to,nd$name)]
    ed$sign <- ifelse(ed$weight>0,"Positive","Negative")
  }
  palette <- grDevices::hcl.colors(length(levels(nd$category)),palette="Dark 3")
  names(palette) <- levels(nd$category)
  if ("Others" %in% names(palette)) palette["Others"] <- "#B5B5B5"
  if ("Unclassified" %in% names(palette)) palette["Unclassified"] <- "#D4C6B0"
  p <- ggplot2::ggplot()
  if (nrow(ed)) {
    if (edge_color=="sign") {
      p <- p+ggplot2::geom_segment(data=ed,
        ggplot2::aes(x=x,y=y,xend=xend,yend=yend,color=sign),
        alpha=.35,linewidth=.35)+
        ggplot2::scale_color_manual(name="Correlation",
          values=c(Positive="#D66A57",Negative="#4B82B2"))
    } else if (edge_color=="none") {
      p <- p+ggplot2::geom_segment(data=ed,
        ggplot2::aes(x=x,y=y,xend=xend,yend=yend),
        color="#B7B7B7",alpha=.35,linewidth=.35)
    } else stop("edge_color must be sign or none.")
  }
  p <- p+ggplot2::geom_point(data=nd,
      ggplot2::aes(x=x,y=y,fill=category,size=.data[[node_size]]),
      shape=21,color="white",stroke=.18,alpha=.95)+
    ggplot2::scale_fill_manual(values=palette,name=node_color,drop=TRUE)+
    ggplot2::scale_size_continuous(range=c(1.5,8),name=node_size)+
    ggplot2::coord_equal()+ggplot2::theme_void(base_size=12)+
    ggplot2::labs(title=paste0("Network: ",data$group))+
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=15),
                   legend.position=if(show_legend) "right" else "none")
  p
}

#' Plot topology metrics across groups
#' @param data Output from aq_net_compare or aq_network_set.
#' @param metric Numeric topology metric column.
#' @export
aq_net_compare_plot <- function(data,metric="Average_degree") {
  if (inherits(data,"aq_network_set")) data <- aq_net_compare(data)
  if (!metric %in% names(data) || !is.numeric(data[[metric]]))
    stop("Invalid numeric metric.")
  ggplot2::ggplot(data,ggplot2::aes(x=Group,y=rlang::.data[[metric]],fill=Group))+
    ggplot2::geom_col(width=.65)+ggplot2::theme_bw()+
    ggplot2::labs(y=metric,x=NULL)
}
#' Plot Zi-Pi node roles
#' @param data Output of aq_net_zp or aq_network.
#' @export
aq_net_zp_plot <- function(data) {
  if (inherits(data,"aq_network") || inherits(data,"aq_network_set"))
    data <- aq_net_zp(data)
  ggplot2::ggplot(data,ggplot2::aes(x=Pi,y=Zi,color=Role))+
    ggplot2::geom_point(alpha=.8)+
    ggplot2::geom_hline(yintercept=2.5,linetype=2)+
    ggplot2::geom_vline(xintercept=.62,linetype=2)+
    ggplot2::facet_wrap(~Group)+ggplot2::theme_bw()
}
#' Export network tables and GraphML
#' @param data aq_network or aq_network_set.
#' @param file Output prefix.
#' @export
aq_net_export <- function(data,file="AqueCo_Network") {
  .aq_net_require()
  nets <- if (inherits(data,"aq_network_set")) data else list(data)
  for (x in nets) {
    prefix <- paste0(file,"_",gsub("[^A-Za-z0-9_-]","_",x$group))
    utils::write.csv(x$nodes,paste0(prefix,"_nodes.csv"),row.names=FALSE)
    utils::write.csv(x$edges,paste0(prefix,"_edges.csv"),row.names=FALSE)
    utils::write.csv(aq_net_node(x),paste0(prefix,"_centrality.csv"),row.names=FALSE)
    utils::write.csv(aq_net_zp(x),paste0(prefix,"_ZiPi.csv"),row.names=FALSE)
    igraph::write_graph(x$graph,paste0(prefix,".graphml"),format="graphml")
  }
  utils::write.csv(aq_net_topology(data),
                   paste0(file,"_topology.csv"),row.names=FALSE)
  invisible(TRUE)
}
