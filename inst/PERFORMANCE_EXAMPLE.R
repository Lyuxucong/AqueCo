
library(aqueco)
# dat is an aq_data object with otu, env, metadata, tree.
rda_quick <- aq_rda_cca(dat,method="RDA",select="none",
  mode="fast",tests="global",group="Group")
rda_full <- aq_rda_cca(dat,method="RDA",select="forward",
  mode="publication",permutations=999,group="Group")
null <- aq_null(dat,group="Group",filter_thres=.0005,cache=TRUE,seed=123)
print(dim(null$otu))
# On Windows use workers=1; on Unix/macOS workers=2 or 4.
nti <- aq_beta_nti(null,runs=99,workers=1)
nri <- aq_beta_nri(null,runs=99,workers=1)
rc <- aq_rcbray(null,runs=99,workers=1)
