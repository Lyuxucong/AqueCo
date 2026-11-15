library(aqueco)
library(ggplot2)
# dat <- aq_import(otu=otu,taxonomy=taxonomy,metadata=metadata)
# Change Group to the actual metadata column name.
res <- aq_lefse(dat,group="Group",min_total_abundance=10,
                alpha=0.05,p_adjust_method="none")
head(res$res_diff)
aq_lefse_export(res,"LEfSe_results.csv")
p_bar <- aq_lefse_plot(res,type="lda",threshold=2,width=0.7,
  group_order=c("IW","CW","TW"),
  colors=c("#ffc168","#7f181b","#537b35"))
ggsave("LEfSe_bar.png",p_bar,width=9,height=9,dpi=300,bg="white")
p_clado <- aq_lefse_plot(res,type="cladogram",
  use_taxa_num=150,use_feature_num=50,clade_label_level=5,
  group_order=c("IW","CW","TW"),
  colors=c("#ffc168","#7f181b","#537b35"),
  sep="|",branch_size=0.5,alpha=0.20,clade_label_size=3)
ggsave("LEfSe_cladogram.tiff",p_clado,width=18,height=12,dpi=300,bg="white")
