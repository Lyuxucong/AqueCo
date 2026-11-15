# =============================================================================
# AqueCo | 微生物生态分析完整示例
# =============================================================================
#
# 文件来源：用户提供的 AqueCo_full_workflow_example.R
# 整理原则：保留原有函数、参数、分析顺序和数据路径；
#           仅优化章节、空行、注释及可读性，不擅自改变统计设置。
#
# 使用方式：
#   1. 根据本机实际情况调整 setwd() 和五个数据文件路径。
#   2. 在 RStudio 中按章节运行，建议先检查数据导入结果。
#   3. devtools::document() / test() / check() 属于开发检查，
#      不必在每次常规分析时执行。
#   4. 部分模块（如 LEfSe、网络、零模型）计算量较大。
#
# 重要提醒：
#   - 原脚本同时出现 KO/OE/WT 与 IW/CW/TW 分组设置；
#     LEfSe 的 group_order 和颜色应与实际 metadata$Group 一致。
#   - 原脚本使用 facet = "Type"、variables = "pH" 等字段；
#     请确认相应列确实存在于上传的数据中。
#   - 本文件为排版与注释整理稿，未在 R 环境中逐项执行。
# =============================================================================

# ============================================================================
# 00 | 环境设置与包检查
# ============================================================================
# 开发阶段使用 load_all() 加载源码；若仅使用已安装包，可直接 library(aqueco)。
# 注意：document()/test()/check() 可能耗时，并可能改变 man/ 和 NAMESPACE。

setwd("D:/aqueco2/aqueco")

# 核心依赖加载
library(devtools)
library(ggplot2)
library(vegan)
library(ape)

devtools::load_all()
devtools::document()
devtools::test()
devtools::check()
library(aqueco)
#设置风格
aq_theme_set("publication")


# ============================================================================
# 01 | 读取数据与构建 AqueCo 对象
# ============================================================================
# 数据路径为原脚本中的本地路径，请按实际存储位置调整。
# 重点检查 OTU 行列方向、样本 ID 是否与 metadata/env 一致。

metadata <- read.csv("D:/aqueco2/data2/metadata.csv",row.names = 1)

env <- read.csv("D:/aqueco2/data2/env.csv",row.names = 1)

otu <- read.csv("D:/aqueco2/data2/otutab.csv",row.names = 1)

taxonomy <- read.csv("D:/aqueco2/data2/taxonomy.csv",row.names = 1)

tree <- read.tree("D:/aqueco2/data2/tree.nwk", tree.names = NULL)

# 查看模拟数据
#dim(otu)
#head(otu[, 1:6])
#head(taxonomy)
#head(metadata)
#head(env)

#
# 2. 导入 AqueCo 数据对象
#

dat <- aq_import(
  otu = otu,
  taxonomy = taxonomy,
  metadata = metadata,
  env = env,
  tree = tree
)

# 检查数据
aq_check(dat)

#

# ============================================================================
# 02 | 配色方案
# ============================================================================
# 颜色向量必须与真实分组水平对应；若 Group 不是 factor，levels() 可能返回 NULL。

# 3. AqueCo 配色系统
#

# 查看内置配色方案
#aq_pal_list()

# 提取 3 个分组需要的颜色
group_colors <- aq_pal("nature", n = 3)
names(group_colors) <- levels(metadata$Group)
group_colors

# 展示配色
#aq_pal_show("nature", n = 8)

# 也可以使用：
# aq_pal_show("npg", n = 8)
# aq_pal_show("lancet", n = 8)
# aq_pal_show("set2", n = 8)
# aq_pal_show("viridis", n = 8)
# aq_pal_show("inferno", n = 8)

#
# 4. 分类水平汇总
#

# ============================================================================
# 03 | 物种组成与可视化
# ============================================================================
# 组成图分别展示样本、分组均值、热图、气泡图及韦恩图。
# facet='Type' 只有在元数据存在 Type 列时才可使用。

# 门水平汇总
#tax_phylum <- aq_taxa(
#  data = dat,
#  level = "phylum",
#  abundance = "relative",
#  unclassified = "keep" )

#tax_phylum

# 属水平汇总
#tax_genus <- aq_taxa(
#  data = dat,
#  level = "Genus",
#  abundance = "relative",
#  unclassified = "keep")

#tax_genus

# 查看汇总矩阵
#dim(tax_phylum$relative)
#head(tax_phylum$relative[, 1:min(6, ncol(tax_phylum$relative))])

# 导出属水平相对丰度表
#write.csv(
#  tax_genus$relative,
#  file = "AqueCo_genus_relative_abundance.csv",
#  quote = FALSE
#)

# 5. 物种组成可视化

# 5.1 每个样本的门水平组成
p_comp_sample <- aq_comp(
  data = dat,
  type = "bar",
  level = "Phylum",
  summary = "sample",
  abundance = "relative",
  top = 6,
  other = TRUE,
 facet = "Group",
  palette = "nature",
  position = "fill",
  theme = "publication"
)

print(p_comp_sample)

#aq_save(
#  p_comp_sample,
#  filename = "AqueCo_composition_phylum_sample.pdf",
#  width = 11,
#  height = 6
#)

# 5.2 按处理组汇总属水平组成，并按时间分面
p_comp_group <- aq_comp(
  data = dat,
  type = "bar",
  level = "Genus",
  summary = "mean",
  group = "Group",
  facet = "Type",
  abundance = "relative",
  top = 10,
  other = TRUE,
  palette = "npg",
  position = "fill",
  theme = "publication"
)

print(p_comp_group)

#aq_save(
#  p_comp_group,
#  filename = "AqueCo_composition_genus_group.pdf",
#  width = 10,
#  height = 6
#)

#热图
p2 <- aq_comp(
  data = dat,
  type = "heatmap",
  level = "Genus",
  summary = "sample",
  abundance = "relative",
  top = 15,
  facet = "Group",
  palette = "nature",
  heatmap_scale = "none",
  cluster_samples = FALSE,
  cluster_taxa = FALSE,
  theme = "publication"
)

print(p2)

#气泡图
p3 <- aq_comp(
  data = dat,
  type = "bubble",
  level = "Genus",
  summary = "sample",
  abundance = "relative",
  top = 15,
  facet = "Group",
  palette = "nature",
  bubble_color = "taxon",
  bubble_size_range = c(1.5, 9),
  alpha = 0.85,
  theme = "publication"
)

print(p3)

#韦恩图
p_venn <- aq_comp(
  data = dat,
  type = "venn",
  level = "Genus",
  group = "Group",
  prevalence = 0.1,
  min_abundance = 0,
  palette = "nature"
)

print(p_venn$plot)


# ============================================================================
# 04 | LEfSe 差异分析
# ============================================================================
# LEfSe 的 group_order 与颜色使用 IW/CW/TW；确认与实际分组一致。
# 差异结果受丰度过滤、显著性阈值及多重检验设置影响。

lefse <- aq_lefse(
  data = dat,
  group = "Group",
  min_total_abundance = 10,
  alpha = 0.05,
  p_adjust_method = "none"
)

# 查看差异结果
#head(lefse$res_diff)

# 导出结果
#aq_lefse_export(lefse, "LEfSe_results.csv")

p_genus <- aq_lefse_plot(
  lefse,
  type = "lda",
  threshold = 2,
  tax_level = "Genus",
# top_n_per_group = 15,
  top_n_per_group = 10,
  short_labels = TRUE,
  group_order = c("IW", "CW", "TW"),
  colors = c(
    IW = "#ffc168",
    CW = "#7f181b",
    TW = "#537b35"
  )
)

print(p_genus)

p_clado <- aq_lefse_plot(
  lefse,
  type = "cladogram",
  use_taxa_num = 150,
  use_feature_num = 50,
  clade_label_level = 5,
  group_order = c("IW", "CW", "TW"),
  colors = c("#ffc168", "#7f181b", "#537b35"),
  sep = "|",
  branch_size = 0.5,
  alpha = 0.20,
  clade_label_size = 3
)

print(p_clado)

ggsave(
  "LEfSe_cladogram.tiff",
  p_clado,
  width = 18,
  height = 12,
  dpi = 300,
  bg = "white"
)


# ============================================================================
# 05 | Alpha 多样性
# ============================================================================
# Alpha 多样性包括丰富度与均匀度指标；组间检验以原脚本设置为准。

# 7.1 一次计算全部 Alpha 指数
alpha_res <- aq_alpha(
  data = dat,
  index = "all",
  group = c("Group")
)

alpha_res

# 查看全部结果
head(alpha_res$data)

# 导出所有 Alpha 指数
#aq_alpha_export(
#  alpha_res,
#  file = "AqueCo_alpha_all_indices.csv")

# 7.2 Shannon 指数箱线图
p_alpha_shannon <- aq_alpha_plot(
  data = alpha_res,
  index = "Shannon",
  group = "Group",
  palette = "nature",
  colors = group_colors,
  fill = TRUE,
  alpha = 0.65,
  point = TRUE,
  point_size = 2,
  test = "kruskal",
  p_adjust = "BH",
  show_p = TRUE,
  theme = "publication"
)

print(p_alpha_shannon)

#aq_save(
#  p_alpha_shannon,
#  filename = "AqueCo_alpha_shannon.pdf",
#  width = 7,
#  height = 6)
#

# 查看统计检验结果
p_alpha_shannon$statistics

# 7.3 多个 Alpha 指数分面箱线图
p_alpha_multi <- aq_alpha_plot(
  data = alpha_res,
  index = c("Observed", "Chao1", "Shannon", "Simpson", "Pielou"),
  group = "Group",
  facet = "Index",
  palette = "npg",
  fill = TRUE,
  alpha = 0.65,
  point = TRUE,
  point_size = 1.6,
  test = "kruskal",
  p_adjust = "BH",
  show_p = TRUE,
  theme = "publication"
)

print(p_alpha_multi)

#aq_save(
#  p_alpha_multi,
#  filename = "AqueCo_alpha_multiple_indices.pdf",
#  width = 11,
#  height = 8
#)


# ============================================================================
# 06 | Beta 多样性：NMDS 与 PCoA
# ============================================================================
# NMDS stress 用于评价排序拟合；PCoA 坐标轴解释率需要结合距离定义理解。
# PERMANOVA 的显著性与组内离散度不同，必要时单独检验离散度。

beta_nmds <- aq_beta(
  data = dat,
  method = "bray",
  ordination = "NMDS",
  group = "Group",
  permutations = 999,
  transform = "relative",
  k = 2,
  trymax = 100,
  binary = FALSE
)

# 查看 NMDS 结果和 stress
beta_nmds
beta_nmds$stress

# NMDS 排序图：分组颜色、形状、中心点和置信椭圆
p_nmds <- aq_beta_plot(
  data = beta_nmds,
  method = "NMDS",
  color = "Group",
  shape = "Group",
  size = 3.2,
  point_alpha = 0.85,
  palette = "nature",
  colors = group_colors,
  centroid = TRUE,
  ellipse = TRUE,
  ellipse_level = 0.95,
  label = FALSE,
  marginal = T,   #改为T，可添加边缘图
  test = "kruskal",
  theme = "publication"
)

print(p_nmds)

# 查看 NMDS 两轴的组间检验
p_nmds$axis_tests

#aq_save(
#  p_nmds,
#  filename = "AqueCo_NMDS_bray.pdf",
 # width = 8,
 # height = 7
#)

# 查看 PERMANOVA
beta_nmds$permanova

# 9. Beta 多样性：PCoA

beta_pcoa <- aq_beta(
  data = dat,
  method = "bray",
  ordination = "PCoA",
  group = "Group",
  permutations = 999,
  transform = "relative",
  k = 2,
  binary = FALSE
)

# PCoA 坐标轴解释率
beta_pcoa$variance[1:2]

# PERMANOVA 结果
beta_pcoa$permanova

# 9.1 标准 PCoA 图
p_pcoa <- aq_beta_plot(
  data = beta_pcoa,
  method = "PCoA",
  color = "Group",
  shape = "Group",
  size = 3.2,
  point_alpha = 0.85,
  palette = "npg",
  colors = group_colors,
  centroid = TRUE,
  ellipse = TRUE,
  ellipse_level = 0.95,
  label = FALSE,
  marginal = T,
  test = "kruskal",
  theme = "publication"
)

print(p_pcoa)

# 查看轴向显著性检验
p_pcoa$axis_tests

#aq_save(
#  p_pcoa,
#  filename = "AqueCo_PCoA_bray.pdf",
#  width = 8,
#  height = 7
#)


# ============================================================================
# 07 | Beta 多样性分解
# ============================================================================
# Sørensen 分解的组成部分需按包内定义解释，注意相对比例与绝对量区别。

beta_decomp <- aq_beta_decomp(
  data = dat,
  group = "Group",
  method = "sorensen",
  prevalence = 0
)

p_decomp <- aq_beta_decomp_plot(
  data = beta_decomp,
  point_size = 5,
  point_alpha = 0.30,
  show_mean = TRUE,
  mean_size = 4,
  palette = "nature",
  theme = "publication"
)

print(p_decomp)

  # 所有组的两两比较结果
  head(beta_decomp$data)

  # 各组的平均组成
beta_decomp$means

# 导出 CSV
#  write.csv(
#    beta_decomp$data,
#    "Beta_decomposition_pairwise.csv",
#    row.names = FALSE
#  )

#  write.csv(
#    beta_decomp$means,
#    "Beta_decomposition_group_means.csv",
#    row.names = FALSE
#  )


# ============================================================================
# 08 | 环境因子差异与相关性
# ============================================================================
# 环境因子名称必须与 env.csv 的列名完全一致（大小写敏感）。
# 原脚本单变量示例为 pH；若不存在请换成真实列名。

#计算所有环境因子
env_res <- aq_env(
  data = dat,
  group = "Group",
  variables = NULL,#可以选定特定因子 variables = c("pH", "TN", "TP")
  test = "kruskal",
  posthoc = "both",
  p_adjust = "BH"
)

#env_res

#env_res$statistics
#env_res$pairwise

#环境因子箱线图
#一次性展示所有环境因子

p_env <- aq_env_plot(
  data = env_res,
  type = "boxplot",
  facet = TRUE,
  facet_scales = "free_y",
  point = TRUE,
  point_size = 1.8,
  palette = "nature",
  show_p = TRUE,
  show_sig = TRUE,
  theme = "publication"
)

print(p_env)

#单独做个别
p_do <- aq_env_plot(
  data = env_res,
  type = "boxplot",
  variables = "pH",
  facet = FALSE,
  palette = "nature",
  show_p = TRUE,
  theme = "publication"
)

print(p_do)

#柱状图
p_env_bar <- aq_env_plot(
  data = env_res,
  type = "bar",
  facet = TRUE,
  bar_error = "SE",
  point = TRUE,
  point_size = 1.8,
  palette = "nature",
  show_p = TRUE,
  theme = "publication"
)

print(p_env_bar)

#事后检验
# 同时计算
#posthoc = "both"

# 只使用 Dunn
#posthoc = "dunn"

# 只使用 Wilcoxon
#posthoc = "wilcox"

#结果导出
#aq_env_export(
#  env_res,
#  file = "Environmental_Analysis"
#)

#环境因子自相关

env_cor <- aq_env_cor(
  data = dat,
  variables = NULL,
  method = "spearman",
  p_adjust = "BH"
)

env_cor

#绘制相关性热图
p_cor <- aq_env_cor_plot(
  data = env_cor,
  show_value = TRUE,
  show_sig = TRUE,
  label_size = 3,
  theme = "publication"
)

print(p_cor)

#分组相关
env_cor_group <- aq_env_cor(
  data = dat,
  variables = NULL,
  method = "spearman",
  group = "Group",
  p_adjust = "BH"
)

p_cor_group <- aq_env_cor_plot(
  data = env_cor_group,
  show_value = TRUE,
  show_sig = TRUE,
  theme = "publication"
)

print(p_cor_group)

#结果导出
#write.csv(
#  env_cor$results,
#  "Environmental_Correlation.csv",
#  row.names = FALSE
#)

#

#计算相关性
env_microbe_cor <- aq_env_microbe_cor(
  data = dat,
  env_variables = NULL,
  level = "Genus",
  abundance = "relative",
  method = "spearman",
  group = NULL,
  prevalence = 0.10,
  min_abundance = 0,
  top = 10,
  p_adjust = "BH"
)

env_microbe_cor

#绘制相关性气泡图

p_env_microbe <- aq_env_microbe_cor_plot(
  data = env_microbe_cor,
  show_value = F,
  show_sig = T,
  only_significant = F,
  p_cutoff = 0.05,
  label_size = 3,
  theme = "publication"
)

print(p_env_microbe)

#按处理组分别分析
env_microbe_cor_group <- aq_env_microbe_cor(
  data = dat,
  level = "Genus",
  abundance = "relative",
  method = "spearman",
  group = "Group",
  prevalence = 0.10,
  p_adjust = "BH"
)

env_microbe_cor_group$results
env_microbe_cor_group$abundance

p_group <- aq_env_microbe_cor_plot(
  data = env_microbe_cor_group,
  show_value = TRUE,
  show_sig = TRUE,
  theme = "publication"
)

print(p_group)


# ============================================================================
# 09 | Mantel 检验
# ============================================================================
# Mantel 置换检验的结果受距离定义、置换次数及样本量影响。

mantel_res <- aq_mantel(
  data = dat,
  group = "Group",
  level = "Genus",
  abundance = "relative",
  community_method = "bray",
  permutations = 999,
  p_adjust = "BH",
  prevalence = 0.10
)

p_mantel <- aq_mantel_plot(
  data = mantel_res,
  cor_show_sig = TRUE,
  cor_p_cutoff = 0.05,
  mantel_p_breaks = c(0.01, 0.05),
  mantel_r_breaks = c(0.2, 0.3),

  show_group_labels = TRUE,
  group_label_size = 3.8,

  theme = "publication"
)

print(p_mantel)


# ============================================================================
# 10 | RDA / CCA 排序
# ============================================================================
# RDA / CCA 的环境变量、转换方式、筛选策略需与研究设计一致。

#RDA

rda_res <- aq_rda_cca(
  data = dat,
  method = "RDA",
  # env_variables = c("pH", "TN", "TP", "NH4", "NO3"  ),
  mode = "fast",
  group = "Group",

  transform = "hellinger",

  select = "forward",
  permutations = 99,
  p_cutoff = 0.05,

  prevalence = 0.05
)

# 最终保留的环境变量
#rda_res$selected_variables

# 候选环境变量
#rda_res$candidate_variables

# 整体模型检验
#rda_res$global_test

# 各变量显著性检验
#rda_res$term_test

# 各排序轴解释率
#rda_res$axis_explained

# R2
#rda_res$r_squared

# 调整R2
#rda_res$adjusted_r_squared

# P值
#rda_res$p_value

#绘制 RDA 图
p_rda <- aq_rda_cca_plot(
  data = rda_res,

  show_environment = TRUE,
  show_species = FALSE,

  show_ellipse = TRUE,
  ellipse_level = 0.95,

  point_size = 3,
  point_alpha = 0.85,

  palette = "nature",

  theme = "publication"
)

print(p_rda)

#CCA
cca_res <- aq_rda_cca(
  data = dat,
  method = "CCA",

#  env_variables = c("pH", "TN", "TP", "NH4", "NO3"),

  group = "Group",

  transform = "none",

  select = "forward",
  permutations = 999,
  p_cutoff = 0.05
)

p_cca <- aq_rda_cca_plot(
  data = cca_res,
  show_environment = TRUE,
  show_ellipse = TRUE,
  palette = "nature",
  theme = "publication"
)

print(p_cca)


# ============================================================================
# 11 | VPA 方差分解
# ============================================================================
# VPA 分组变量应有明确生态学含义；注意共线性与负解释分量。

vpa_res <- aq_vpa(
  data = dat,

  groups = list(

    Physical = c(
      "Altitude",
      "Temperature",
      "Precipitation",
      "TOC"
    ),
    Chemical = c(
      "NH4",
      "NO3",
      "pH",
      "Conductivity",
      "TN"
    )
  ),

  # 纳入空间变量
  include_space = TRUE,

  latitude = "latitude",
  longitude = "longitude",
  # 环境变量标准化
  env_transform = "scale",
  # 物种数据 Hellinger 转换
  community_transform = "hellinger",
  min_prevalence = 0.05
)

# 查看完整分析结果
#vpa_res

# 查看各个解释变量组
#vpa_res$group_names

# 查看 VPA 具体分解结果
#vpa_res$fractions

# 查看各解释变量组的总体解释结果
#vpa_res$table_fractions

# 查看使用的空间变量
#vpa_res$spatial$dbmem

#绘制 Venn 图
aq_vpa_plot(
  data = vpa_res,
  cutoff = 0,
  digits = 2,
  alpha = 60,
  colors = c(
    "#6BAED6",
    "#74C69D",
    "#F4A261"
  )
)


# ============================================================================
# 12 | 群落构建过程：NCM 与零模型
# ============================================================================
# NCM 的 Nm 与 R² 不代表中性过程贡献率；组间比较需考虑测序深度。
# 零模型通常计算耗时，建议先用较小 runs 检查流程，再做正式分析。

#整体分析

ncm_all <- aq_ncm(
  data = dat
)

# 查看整体模型结果
ncm_all$summary

# 查看物种水平结果
head(ncm_all$taxa)

#按Group分组
ncm_group <- aq_ncm(
  data = dat,
  group = "Group"
)

# 查看三个组的模型结果
ncm_group$summary

ncm_group$summary[, c(
  "Group",
  "R2",
  "Nm",
  "m",
  "N",
  "Samples",
  "Taxa"
)]

#NCM 分面绘图
p_ncm <- aq_ncm_plot(
  data = ncm_group,

  facet = TRUE,

  show_ci = TRUE,
  show_class = TRUE,

  point_size = 2,
  point_alpha = 0.65,

  theme = "publication"
)

print(p_ncm)

#分别绘制

p_list <- aq_ncm_plot(
  data = ncm_group,
  facet = FALSE,
  show_ci = TRUE,
  show_class = TRUE,
  point_size = 2,
  point_alpha = 0.65
)

# 分别显示
print(p_list$KO)
print(p_list$OE)
print(p_list$WT)

#初始化零模型
null <- aq_null(
  data = dat,
  filter_thres = 0.0005,
  group = "Group",
  seed = 123
)

#betaNRI 分析
nri <- aq_beta_nri(
  data = null,
  runs = 500,
  abundance.weighted = TRUE,
  workers = 1,
  null.model = "taxa.labels"
)

nri$matrix

#betaNTI 分析
nti <- aq_beta_nti(
  data = null,
  runs = 500,
  abundance.weighted = TRUE,
  null.model = "taxa.labels"
)

nti$matrix

#RCbray 分析
rc <- aq_rcbray(
  data = null,
  runs = 1000,
  seed = 123
)

rc$matrix

#五类群落构建过程
assembly <- aq_assembly(
  betanti = nti,
  rcbray = rc,
  group = "Group",
  comparison = "within"
)

assembly$summary
assembly$pairs

#绘图
# betaNTI 组间比较
p1 <- aq_null_plot(
  data = nti,
  type = "boxplot",
  group = "Group",
  add = "mean",
  threshold = c(-2, 2),
  test = "wilcox"
)

# 五类生态过程
p2 <- aq_assembly_plot(
  data = assembly,
  type = "stacked",
  palette = "publication"
)

# betaNTI 热图
p3 <- aq_null_plot(
  data = nti,
  type = "heatmap",
  cluster = TRUE
)

print(p1)
print(p2)
print(p3)


# ============================================================================
# 13 | 微生物共现网络
# ============================================================================
# 网络阈值与过滤规则会影响节点、边和模块结构。
# 不同组的网络规模不同，拓扑指标比较需谨慎解释。

##构建 KO、OE、WT 共现网络
net <- aq_network(
  data = dat,
  group = "Group",
  groups = c("IW", "CW", "TW"),
  method = "spearman",
  min_abundance = 0.0005,
  min_prevalence = 0.3,
  cor_cutoff = 0.6,
  p_cutoff = 0.05,
  p_adjust = "BH"
)

# 查看网络对象
names(net)

# 查看各组保留的边数
sapply(net, function(x) nrow(x$edges))

##计算网络拓扑指标
topology <- aq_net_topology(net)

print(topology)

# 导出
#write.csv(
#  topology,
#  "Network_Topology.csv",
#  row.names = FALSE
#)

##节点中心性分析
node_metrics <- aq_net_node(net)

head(node_metrics)

# 按 Degree 从高到低排序
node_metrics <- node_metrics[
  order(-node_metrics$Degree),
]

# 查看连接度最高的节点
head(node_metrics, 20)

#write.csv(
#  node_metrics,
#  "Network_Node_Centrality.csv",
#  row.names = FALSE
#)

##网络模块分析
modules <- aq_net_module(
  data = net,
  seed = 123
)

head(modules)

# 查看各组模块数量
aggregate(
  Module ~ Group,
  data = modules,
  FUN = function(x) length(unique(x))
)

#write.csv(
#  modules,
#  "Network_Modules.csv",
#  row.names = FALSE
#)

##Zi–Pi 关键物种分析

zp <- aq_net_zp(net)

head(zp)

# 各类节点数量
table(zp$Group, zp$Role)

# 提取候选关键物种
keystone <- aq_net_keystone(net)

print(keystone)

#write.csv(
#  zp,
#  "Network_ZiPi.csv",
#  row.names = FALSE
#)

#write.csv(
#  keystone,
#  "Network_Keystone.csv",
#  row.names = FALSE
#)

##绘制三组共现网络图

p_ko <- aq_net_plot(
  net$IW,
  layout = "fr",
  node_color = "Phylum",
  node_size = "Degree",
  top_n = 10,
  edge_color = "sign",
  remove_isolates = TRUE
)

print(p_ko)
print(p_oe)
print(p_wt)

##绘制网络拓扑比较图
comparison <- aq_net_compare(net)

# 节点数
p_nodes <- aq_net_compare_plot(
  comparison,
  metric = "Nodes"
)

# 边数
p_edges <- aq_net_compare_plot(
  comparison,
  metric = "Edges"
)

# 平均度
p_degree <- aq_net_compare_plot(
  comparison,
  metric = "Average_degree"
)

# 网络密度
p_density <- aq_net_compare_plot(
  comparison,
  metric = "Density"
)

print(p_nodes)
print(p_edges)
print(p_degree)
print(p_density)

#绘制 Zi–Pi 图

p_zp <- aq_net_zp_plot(zp)

print(p_zp)

ggsave(
  "Network_Figures/Network_ZiPi.tiff",
  p_zp,
  width = 10,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

#导出所有结果

#aq_net_export(
#  data = net,
#  file = "AqueCo_Network"
#)


# ============================================================================
# 14 | 距离衰减关系 DDR
# ============================================================================
# 地理距离计算依赖正确的经纬度列与单位。
# DDR 斜率、显著性及拟合优度应结合距离范围和置换方法报告。

ddr <- aq_ddr(
  data = dat,
  group = "Group",
  groups = c("IW", "CW", "TW"),
  scope = "both",
  distance = c("geographic", "environmental"),
  lon = "longitude",
  lat = "latitude",
# env_variables = c("NH4", "NO3", "pH","Conductivity", "TN" ),
  standardize = TRUE,
  permutations = 999,
  seed = 123
)

ddr$summary

# 三组地理距离衰减
p_geo <- aq_ddr_plot(
  ddr,
  distance = "geographic",
  groups = c("IW", "CW", "TW"),
  compare = TRUE
)

print(p_geo)

# 三组环境距离衰减
p_env <- aq_ddr_plot(
  ddr,
  distance = "environmental",
  groups = c("IW", "CW", "TW"),
  compare = TRUE
)

print(p_env)

p_iw <- aq_ddr_plot(
  ddr,
  distance = "geographic",
  groups = "IW"
)

print(p_iw)
