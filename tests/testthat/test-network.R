
test_that("network construction and metrics support grouped data", {
  skip_if_not_installed("igraph")
  set.seed(2026)
  otu <- matrix(sample(1:40, 12*6, replace=TRUE),12,6,
     dimnames=list(paste0("S",1:12),paste0("ASV",1:6)))
  md <- data.frame(Group=rep(c("KO","WT"),each=6),
                   row.names=rownames(otu))
  dat <- aq_data(otu,metadata=md)
  net <- suppressWarnings(aq_network(dat,groups=c("KO","WT"),
     min_abundance=0,min_prevalence=0,cor_cutoff=.3,p_cutoff=1))
  expect_s3_class(net,"aq_network_set")
  expect_true(all(c("KO","WT") %in% names(net)))
  top <- aq_net_topology(net)
  expect_equal(nrow(top),2L)
  expect_true(all(top$Nodes>=2))
  expect_equal(nrow(aq_net_zp(net)),sum(top$Nodes))
  expect_s3_class(aq_net_plot(net[["KO"]]),"ggplot")
  expect_s3_class(aq_net_compare_plot(net),"ggplot")
})
test_that("network rejects invalid thresholds", {
  skip_if_not_installed("igraph")
  otu <- matrix(1:24,6,4,dimnames=list(paste0("S",1:6),
    paste0("T",1:4)))
  dat <- aq_data(otu,metadata=data.frame(Group=rep("KO",6),
    row.names=rownames(otu)))
  expect_error(aq_network(dat,cor_cutoff=2),"thresholds")
})

test_that("negative correlation weights do not break topology shortest paths", {
  skip_if_not_installed("igraph")
  g <- igraph::make_ring(4)
  igraph::E(g)$weight <- c(-0.8, 0.7, -0.6, 0.9)
  x <- structure(list(graph=g, group="KO", n_samples=6,
    edges=data.frame(weight=igraph::E(g)$weight)), class="aq_network")
  result <- aq_net_topology(x)
  expect_equal(result$Mean_path_largest, 4/3)
})

test_that("network plotting supports taxonomic fill and isolate filtering", {
  skip_if_not_installed("igraph")
  g <- igraph::make_empty_graph(4,directed=FALSE)
  g <- igraph::add_edges(g,c(1,2,2,3))
  igraph::V(g)$name <- paste0("ASV",1:4)
  igraph::E(g)$weight <- c(0.7,-0.8)
  nd <- data.frame(name=paste0("ASV",1:4),
    Phylum=c("Firmicutes","Proteobacteria","Firmicutes",NA))
  edges <- data.frame(from=c("ASV1","ASV2"),to=c("ASV2","ASV3"),
    weight=c(.7,-.8),sign=c("Positive","Negative"))
  x <- structure(list(graph=g,nodes=nd,edges=edges,group="KO"),
    class="aq_network")
  expect_s3_class(aq_net_plot(x,node_color="Phylum"),"ggplot")
  expect_s3_class(aq_net_plot(x,remove_isolates=TRUE),"ggplot")
  expect_s3_class(aq_net_plot(x,node_color="Module"),"ggplot")
  expect_error(aq_net_plot(x,node_color="Nonexistent"),"not found")
})

test_that("default plot supports networks without taxonomy", {
  skip_if_not_installed("igraph")
  g <- igraph::make_ring(4)
  igraph::V(g)$name <- paste0("ASV",1:4)
  igraph::E(g)$weight <- c(.7,-.6,.8,-.7)
  x <- structure(list(graph=g,nodes=data.frame(name=paste0("ASV",1:4)),
    edges=data.frame(from=character(),to=character()),group="KO"),
    class="aq_network")
  expect_s3_class(aq_net_plot(x),"ggplot")
  expect_error(aq_net_plot(x,node_color="Class"),"not found")
})
