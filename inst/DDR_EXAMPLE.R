
library(aqueco)
library(ggplot2)
# dat <- aq_import(otu=otu,taxonomy=taxonomy,metadata=metadata,env=env)
# Longitude/latitude must be in metadata, in decimal degrees (WGS84).
# If coordinate column names differ, set lon and lat explicitly.
ddr <- aq_ddr(dat,group="Group",scope="both",
              distance=c("geographic","environmental"),
              lon="Longitude",lat="Latitude",
              permutations=999,seed=123)
print(ddr$summary)
p_all <- aq_ddr_plot(ddr)
p_iw <- aq_ddr_plot(ddr,groups="IW",distance="geographic")
p_compare <- aq_ddr_plot(ddr,groups=c("IW","CW","TW"),
                         distance="environmental",compare=TRUE)
ggsave("DDR_comparison.tiff",p_compare,width=10,height=6,
       dpi=600,compression="lzw")
aq_ddr_export(ddr,"DDR")
