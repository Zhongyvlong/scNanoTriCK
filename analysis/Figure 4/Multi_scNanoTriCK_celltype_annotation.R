



# ==============================================================================
# Pipeline: Single-Cell Reference Projection & Annotation Transfer (E5.0 Embryo)
# Logic: Project new Dual_RNA (Batch 4) onto E5.0 Reference (Batch 1-3)
# ==============================================================================

library(Seurat)
library(Harmony)
library(SingleR)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)

# --- 1. 数据加载与环境准备 ---
# e5.0.old: 包含 rna_v8 标签的参考集
# Dual_RNA: 待加入的查询集 (Batch 4)
ref <- readRDS("E5.0_rna_cluster_v8.rds")
query <- readRDS("Dual_RNA.rds")

# --- 2. 统一特征过滤 (删除线粒体、Gm 基因等) ---
# 建议在合并前对 Query 进行初步 QC
query[["percent.mt"]] <- PercentageFeatureSet(query, pattern = "^mt-")

# 合并并标记 Batch
e5.0 <- merge(ref, y = query, add.cell.ids = c("ref", "query"))
e5.0$batch[is.na(e5.0$batch)] <- "batch4" # 假设原参考集已有 batch1-3

# 过滤干扰基因 (Gm-, mt-, Rp-)
genes_to_keep <- grep("^(Gm|mt-|Rp[sl])", rownames(e5.0), value = TRUE, invert = TRUE)
e5.0 <- subset(e5.0, features = genes_to_keep)
e5.0 <- JoinLayers(e5.0)

# --- 3. 引用投影 (Reference PCA Projection) ---

e5.0 <- NormalizeData(e5.0, normalization.method = "LogNormalize", scale.factor = 10000)
VariableFeatures(e5.0) <- VariableFeatures(ref) 
pca_loadings <- Loadings(ref[["pca"]])
common_genes <- intersect(rownames(pca_loadings), rownames(e5.0))
pca_loadings <- pca_loadings[common_genes, ]
scaled_data <- GetAssayData(e5.0, slot = "data")[common_genes, ] 
new_pca_embeddings <- t(scaled_data) %*% pca_loadings
e5.0[["pca_ref"]] <- CreateDimReducObject(
  embeddings = as.matrix(new_pca_embeddings),
  loadings = pca_loadings,
  key = "PC_",
  assay = "RNA"
)

# --- 4. Harmony 批次校正 ---
e5.0 <- RunHarmony(e5.0, group.by.vars = "batch", reduction.use = "pca_ref", 
                   dims.use = 1:30, max_iter = 5)

# --- 5. UMAP 投影可视化 ---
# 将 Query 投影到 Ref 的 UMAP 模型上
ref <- RunUMAP(ref, dims = 1:30, return.model = TRUE)
query <- MapQuery(anchorset = anchors, reference = ref, query = query,
                  refdata = list(celltype = "rna_v8"), 
                  reference.reduction = "pca", reduction.model = "umap")

# --- 7. 结果保存与导出 ---
DimPlot(e5.0, group.by = "SingleR_label", label = TRUE) + ggtitle("SingleR Annotation on Merged Object")
saveRDS(e5.0, "E5.0_Integrated_with_Batch4.rds")