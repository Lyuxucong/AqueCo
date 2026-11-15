# AqueCo network module — development release

## Implemented
`aq_network`, `aq_net_topology`, `aq_net_node`, `aq_net_module`,
`aq_net_zp`, `aq_net_keystone`, `aq_net_compare`, `aq_net_plot`,
`aq_net_compare_plot`, `aq_net_zp_plot`, `aq_net_export`.

Correlation networks use sample-wise relative abundances, Spearman/Pearson
tests and BH adjustment across all candidate taxon pairs within each group.
Edges satisfy `abs(r) >= cor_cutoff` and `adjusted_p < p_cutoff`.
Disconnected and isolated nodes are retained.

Community detection uses Louvain with absolute edge weights; Zi-Pi is
calculated on the unweighted graph. Network plotting is via ggplot2.
GraphML is exported for Cytoscape.

## Scientific limitations
The included dataset has 6 samples per treatment (KO/OE/WT).
Within-group correlations are underpowered, and multiple-testing correction
can yield zero edges. No claim of biological interaction should be inferred.
Topological group comparisons are descriptive, not inferential.
SPIEC-EASI, SparCC, differential edges, bootstrapped stability and
network robustness are not implemented in this development build.
Run tests and R CMD check in an R environment before publication.
