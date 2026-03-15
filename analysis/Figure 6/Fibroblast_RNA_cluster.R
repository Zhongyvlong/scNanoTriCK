Fibroblast single-Cell RNA-seq Analysis

# ==============================================================================
# Pipeline: Fibroblast Single-Cell RNA-seq Analysis
# Description: QC, Integration (Harmony), Clustering, and Cell Type Identification
# ==============================================================================

library(Seurat)
library(Harmony)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(clusterProfiler)
library(org.Mm.eg.db)

# --- 1. 数据加载与预处理 ---
# 建议：在 GitHub 上将路径改为相对路径或变量
data_dir <- "./data/counts/"
samples <- c("Hu25001", "Hu25004", "Hu25011", "Hu25012")
metadata_map <- list(
  "Hu25001" = "rep2_p1", "Hu25004" = "rep1_p1", 
  "Hu25011" = "rep1_p2", "Hu25012" = "rep2_p2"
)

# 批量读取并创建 Seurat 对象
obj_list <- lapply(samples, function(s) {
  counts <- read.table(paste0(data_dir, s, "_RNA_count.txt"))
  obj <- CreateSeuratObject(counts = counts, project = s, min.cells = 10, min.features = 200)
  obj$sample <- metadata_map[[s]]
  return(obj)
})

# 合并数据
fibroblast <- merge(obj_list[[1]], y = obj_list[2:4], project = "fibroblast")

# --- 2. 质量控制 (QC) ---
fibroblast[["percent.mt"]] <- PercentageFeatureSet(fibroblast, pattern = "^mt-")

# 根据测序阶段标记 Stage (S3=E14.5, S5=P1, S6=Adult, S7=MI, S8=Aging)
fibroblast$stage <- case_when(
  grepl("_S3$", colnames(fibroblast)) ~ "E14.5",
  grepl("_S5$", colnames(fibroblast)) ~ "P1",
  grepl("_S6$", colnames(fibroblast)) ~ "Adult",
  grepl("_S8$", colnames(fibroblast)) ~ "Aging",
  TRUE ~ "Other"
)

# 过滤低质量细胞
fibroblast.filter <- subset(fibroblast, 
                            subset = nCount_RNA > 1500 & 
                              nCount_RNA < 10000 & 
                              percent.mt < 10 & 
                              stage != "Other")

# --- 3. 标准化与降维 ---
fibroblast.filter <- JoinLayers(fibroblast.filter)
fibroblast.filter <- NormalizeData(fibroblast.filter) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
  RunPCA()

# 使用 Harmony 整合批次 (按 sample 整合)
fibroblast.filter <- RunHarmony(fibroblast.filter, group.by.vars = "sample", dims.use = 1:15)

# --- 4. DoubletFinder 去双细胞 ---
# 注意：实际使用时需根据不同样本分别运行，此处为简化示例
sweep.res <- paramSweep(fibroblast.filter, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_val <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

nExp_poi <- round(0.075 * nrow(fibroblast.filter@meta.data)) # 假设 7.5% 双细胞率
fibroblast.filter <- doubletFinder(fibroblast.filter, PCs = 1:10, pN = 0.25, pK = pK_val, nExp = nExp_poi)

# 过滤掉双细胞 (需根据具体列名修改)
df_col <- grep("DF.classifications", colnames(fibroblast.filter@meta.data), value = TRUE)
fibroblast.filter <- subset(fibroblast.filter, cells = colnames(fibroblast.filter)[fibroblast.filter@meta.data[[df_col]] == "Singlet"])

# --- 5. 聚类与可视化 ---
fibroblast.filter <- FindNeighbors(fibroblast.filter, dims = 1:15, reduction = "harmony") %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(dims = 1:15, reduction = "harmony")

# --- 6. 标记基因与注释 ---
# 常用 Marker 列表
markers_list <- list(
  Fibroblast = c("Pdgfra", "Postn", "Col1a1", "Dcn"),
  Endothelial = c("Pecam1", "Vwf", "Cdh5"),
  Valve_derived = c("Hapln1", "Lef1", "Sall3"),
  Active_Fibro = c("Meox1", "Col8a1", "Cilp")
)

# 绘制标记基因点图
DotPlot(fibroblast.filter, features = unlist(markers_list), assay = "RNA") + 
  coord_flip() + 
  scale_color_gradient2(low = "#5891BF", high = "#AA3538", mid = "white")

# --- 7. 保存结果 ---
saveRDS(fibroblast.filter, file = "output/fibroblast_final_processed.rds")