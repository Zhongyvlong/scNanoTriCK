########## fibroblast DNA clustering ------------------------------------
##### merged peak in signac
subpair<-read.table("Analysis/250910_fibro/metadata/fibro_pairedcell_dna200.txt",header = T)
subpair$dna_archr<-paste(subpair$Nano,subpair$CB,sep = "#")
rownames(subpair)<-subpair$dna_archr

peak<-import("Analysis/250910_fibro/frip/merged_allpeak.bed",format = "BED")
count_peak<-matrix(nrow = 264143)
for(i in c("HuNACHIP25001","HuNACHIP25002","HuNACHIP25004","HuNACHIP25009")){
  file=list.files(path = paste0("/data/csy/zyl_analysis/demultiplex/fibroblast_highthroughput/",i,"/demultiplex/fragment"), pattern = "*_signac_fragfromflank.bed.gz$", full.names = T)
  fragments <- CreateFragmentObject(file)
  mat<-FeatureMatrix(fragments = fragments,features = peak )
  count_peak<-cbind(count_peak,mat[,intersect(colnames(mat),subpair$dna_archr)])
}
count_peak<-count_peak[,2:ncol(count_peak)]

fibroOb<-readRDS("fibroOb_peakcount_filter_dnadb.RDS")
fibroOb <- fibroOb[rowSums(fibroOb)!=0,]
fibroOb <- RunTFIDF(fibroOb)
fibroOb <- FindTopFeatures(fibroOb, min.cutoff = 'q5')
fibroOb <- RunSVD(object = fibroOb)
DepthCor(fibroOb,n=20)
set.seed(1234)
fibroOb <- RunHarmony(fibroOb,group.by.vars=c("replicate"),early_stop=F,sigma=5,theta=10,lambda=10,reduction.use = "lsi",assay.use="peaks",project.dim=F) #,sigma=0.3,lambda=1,early_stop=F
DepthCor(fibroOb,reduction = "harmony",n=20)
fibroOb <- RunUMAP(object = fibroOb,reduction = 'harmony',dims = c(3:15),min.dist = 0.5,spread = 0.8,return.model = T) #,spread = 0.5
fibroOb <- FindNeighbors(object = fibroOb,reduction = 'harmony',annoy.metric = "cosine",dims = c(3:15))
fibroOb <- FindClusters(object = fibroOb,algorithm = 3,resolution = 1.2,verbose = FALSE)
p1<-DimPlot(object = fibroOb, label = TRUE,pt.size=0.1)+labs(title = c("3:15 | 1.2"))
p2<-DimPlot(object = fibroOb, group.by = "stage",pt.size=0.1)
p3<-FeaturePlot(object = fibroOb, features = "nCount_peaks",pt.size=0.3)
p4<-FeaturePlot(object = fibroOb, features = "frip",pt.size=0.3)
(p1+p2)/(p3|p4)
fibroOb<-readRDS("fibroOb_cluster_v3.RDS")
VlnPlot(fibroOb,features = c("nCount_peaks","nFeature_peaks","dup_reads","frip"),log = T,ncol = 2)
dna_meta<-as.data.frame(fibroOb@meta.data)

########## clustree
library(clustree)
resol<-FindClusters(object = fibroOb,algorithm = 3,resolution = seq(0.4,1.5,0.1))
clutree<-clustree(resol)
print(clutree)+scale_color_manual(values = brewer.pal(12,"Paired"))

source("Analysis/TSSGeneActivity.R")
library(data.table)
DefaultAssay(fibroOb) <- "peaks"
gene.activities <- TSSGeneActivity(fibroOb, extend.upstream = 30000, extend.downstream = 30000)
gene.activities <- gene.activities[rowSums(gene.activities)!=0,]
fibroOb[['geneact']] <- CreateAssayObject(counts = gene.activities)
fibroOb <- NormalizeData(object = fibroOb,assay = 'geneact',normalization.method = 'LogNormalize',scale.factor = median(fibroOb$nCount_geneact) )

markers <- FindAllMarkers( fibroOb,logfc.threshold = 0.25,min.pct = 0.1,only.pos =T,assay ="geneact")%>% dplyr::filter(p_val<0.05)%>%arrange(-avg_log2FC)
markers<-markers[-grep("^Gm",markers$gene),]
table(markers$cluster)
DotPlot(fibroOb, features = unique((markers%>%group_by(cluster)%>%slice_head(n = 6) %>%ungroup())$gene), cluster.idents =F,assay = "geneact",dot.scale = 10) + coord_flip() + scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+theme(axis.text.y=element_text(size=10))

genelist<-c("Col9a1","Col9a3","Sall3","Hapln1","Lef1","Spint2","Ociad2", # valve
            "Ccnb1","Pbk","Top2a","Ube2c", # cycle
            "Cdh5","Flt1","Tie1","Kdr","Rasip1","Myct1", # endo
            "Tbx18","Wt1","Aldh1a2", # epi
            #"Gata4","Tagln2","Socs3","Nr4a1",
            "Ly6c1","Pi16","Ly6a","Gfpt2","Cd248","Meox1","Cilp","Ckb","Mt2", # F-SH, F-Act
            "Apoe","Wif1","Dkk3","Cytl1", # F-Wnt
            "Postn","Pdgfa","Ddr2","Cthrc1","Fn1","Acta2") # F-MYO
DotPlot(fibroOb, features = genelist, cluster.idents =F,assay = "geneact",dot.scale = 8) + coord_flip() + scale_y_discrete(limits=c("0","6","8","5","3","7","2","4","1","10","9")) + scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+theme(axis.text.y=element_text(size=10))

FeaturePlot(fibroOb,features = c("Gata4","Tagln2","Socs3","Nr4a1"),max.cutoff = "q95",min.cutoff = "q30",order = T,ncol=2,pt.size = 0.1)


########## fibroblast WNN clustering ------------------------------------
fibroOb<-readRDS("fibroOb_cluster_DNA.RDS")
fibro_rna<-readRDS("fibro_rna_R1.rds")
DimPlot(object = fibro_rna, label = TRUE,pt.size=0.1,group.by = "seurat_clusters",reduction = "umap")+labs(title = c("fibroblast RNA"))
DimPlot(object = fibroOb, label = TRUE,pt.size=0.1,group.by = "seurat_clusters",reduction = "umap")+labs(title = c("fibroblast DNA"))
fibro_rna<-fibro_rna[,fibroOb$rna_id]
fibro_rna$dna_archr<-fibroOb@meta.data[match(colnames(fibro_rna),fibroOb$rna_id),"dna_archr"]

fibro_merged<-fibroOb
fibro_merged$dna.clu<-fibro_merged$seurat_clusters
fibro_merged$rna.clu<-fibro_rna@meta.data[match(fibro_merged$rna_id,colnames(fibro_rna)),"seurat_clusters"]
fibro_merged@assays$RNA<-CreateAssay5Object(counts = fibro_rna@assays$RNA$counts%>%`[`(, fibro_merged$rna_id)%>%`colnames<-`(fibro_merged$dna_archr), data = fibro_rna@assays$RNA$data%>%`[`(, fibro_merged$rna_id)%>%`colnames<-`(fibro_merged$dna_archr), key="rna_" )
fibro_merged@assays$RNA$scale.data<-fibro_rna@assays$RNA$scale.data %>% `[`(, fibro_merged$rna_id) %>% `colnames<-`(fibro_merged$dna_archr)
fibro_merged@assays$RNA@meta.data<-fibro_rna@assays$RNA@meta.data
fibro_merged@reductions$rnapca<-CreateDimReducObject(embeddings = fibro_rna@reductions$pca@cell.embeddings%>%`[`(fibro_merged$rna_id, )%>%`rownames<-`(fibro_merged$dna_archr),assay = "RNA",key = "rnapca_")
fibro_merged@reductions$rnaumap<-CreateDimReducObject(embeddings = fibro_rna@reductions$umap@cell.embeddings%>%`[`(fibro_merged$rna_id, )%>%`rownames<-`(fibro_merged$dna_archr),assay = "RNA",key = "rnaumap_")

fibro_merged <- FindMultiModalNeighbors(fibro_merged, reduction.list = list("rnapca","harmony"), dims.list = list(1:10, 3:15))
fibro_merged <- RunUMAP(fibro_merged, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",spread = 0.5)
fibro_merged <- FindClusters(fibro_merged, graph.name = "wsnn",resolution = 1.0, algorithm = 3, verbose = FALSE)
fibro_merged$wnn.clu<-paste0("C",fibro_merged$seurat_clusters)
fibro_merged$wnn.clu<-dplyr::recode(fibro_merged$wnn.clu,"C9"="C1","C10"="C2","C0"="C3","C4"="C4","C12"="C4","C7"="C5","C3"="C7","C8"="C8","C1"="C9","C2"="C10","C11"="C11","C6"="C12","C5"="C13")
fibro_merged$wnn.clu<-factor(fibro_merged$wnn.clu,levels = paste0("C",1:13))
DimPlot(fibro_merged, reduction = "wnn.umap", group.by = "wnn.clu", label = TRUE, repel = TRUE,cols = wnncol) + ggtitle("WNN cluster")
table(fibro_merged$wnn.clu)

p1 <- DimPlot(fibro_merged, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("WNN cluster")
p2 <- DimPlot(fibro_merged, reduction = "wnn.umap", group.by = "stage", label = TRUE, repel = TRUE) + ggtitle("WNN stage")
p3 <- DimPlot(fibro_merged, reduction = "umap", group.by = "dna.clu", label = TRUE, repel = TRUE) + ggtitle("DNA")
p4 <- DimPlot(fibro_merged, reduction = "wnn.umap", group.by = "rna.clu", label = TRUE, repel = TRUE) + ggtitle("RNA inWNN")
p5 <- DimPlot(fibro_merged, reduction = "rnaumap", group.by = "rna.clu", label = TRUE, repel = TRUE) + ggtitle("RNA")
p6 <- DimPlot(fibro_merged, reduction = "wnn.umap", group.by = "dna.clu", label = TRUE, repel = TRUE) + ggtitle("DNA inWNN")
(p1+p2+p3)/(p4+p5+p6) & theme(plot.title = element_text(hjust = 0.5))

gene_showing_paper <- c(
  "Wt1", "Aldh1a2", "Nfatc1", 
  "Gata4", "Tbx5", "Hapln1",   "Lef1", "Col9a1", "Cdh5", "Kdr", "Adgrl4",  "Mki67", "Top2a", "H2afz", "Cenpa", "Hmgb2",
  "Pi16", "Ly6a", "Gfpt2", "Cxcl14", "Hsd11b1", "Lpl", "Dpep1", 
  "G0s2", "Apoe", "Ccl19", "Fgl2", "Inmt", "Wif1", "Prg4", "Dkk3", "Meox1", 
  "Eln", "Cilp", "Col8a1","Cthrc1", "Spp1", "Acta2", "Tpm2"
)
p1<-DotPlot(fibro_merged, features = gene_showing_paper, cluster.idents =F,assay = "RNA",dot.scale = 8) + coord_flip() + scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+theme(axis.text.y=element_text(size=10))
p2<-DotPlot(fibro_merged, features = gene_showing_paper, cluster.idents =F,assay = "geneact",dot.scale = 8) + coord_flip() + scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+theme(axis.text.y=element_text(size=10))
p1+p2

#####sankey or allulium
library(sankeyD3)
library(networkD3)
library(webshot)
fibro_merged<-readRDS("fibro_mergedOb_cluster_v1.RDS")
DimPlot(fibro_merged, reduction = "umap", group.by = "dna.clu", label = TRUE, repel = TRUE) + ggtitle("DNA cluster")
table(fibro_merged$dna.clu)

cross_tab <- as.data.frame(fibro_merged@meta.data[,c("dna.clu","wnn.clu","rna.clu")])%>%
  dplyr::mutate(dna.clu=paste0("dna_",dna.clu),wnn.clu=paste0("wnn_",wnn.clu),rna.clu=paste0("rna_C",rna.clu))%>%
  dplyr::mutate(dna.clu=factor(dna.clu,levels=paste0("dna_C",1:11)), wnn.clu=factor(wnn.clu,levels=paste0("wnn_C",1:13)), rna.clu=factor(rna.clu,levels=paste0("rna_C",1:13)) )

##### three group
links<-rbind.data.frame(cross_tab%>%group_by(dna.clu,wnn.clu)%>%dplyr::summarise(weight=n())%>%setNames(c("source","target","weight")), 
                           cross_tab%>%group_by(wnn.clu,rna.clu)%>%dplyr::summarise(weight=n())%>%setNames(c("source","target","weight")) )

nodes <- data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique())
nodes$name<-factor(nodes$name,levels = c(paste0("dna_C",1:11), paste0("wnn_C",1:13), paste0("rna_C",1:13) ))
nodes<-nodes[order(nodes$name),, drop = FALSE]
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1
nodes$color<-c(dnacol,wnncol,rnacol)
links$linkcol<-c(links[1:105,]$target,links[106:181,]$source)
colour='d3.scaleOrdinal().domain(["dna_C1","dna_C2","dna_C3","dna_C4","dna_C5","dna_C6","dna_C7","dna_C8","dna_C9","dna_C10","dna_C11","wnn_C1","wnn_C2","wnn_C3","wnn_C4","wnn_C5","wnn_C7","wnn_C8","wnn_C9","wnn_C10","wnn_C11","wnn_C12","wnn_C13",
"rna_C1","rna_C2","rna_C3","rna_C4","rna_C5","rna_C6","rna_C7","rna_C8","rna_C9","rna_C10","rna_C11","rna_C12","rna_C13"]).range(["#F5CFE4","#FCED82","#D2EBC8","#7DBFA7","#8FA4AE","#9B5B33","#EE934E","#3C77AF","#F5D2A8","#B383B9","#D1352B","#F5CFE4","#FCED82",
"#D2EBC8","#BBDD78","#7DBFA7","#9B5B33","#AECDE1","#EE934E","#3C77AF","#F5D2A8","#B383B9","#D1352B","#F5CFE4","#FCED82","#D2EBC8","#BBDD78","#7DBFA7","#8FA4AE","#9B5B33","#AECDE1","#EE934E","#3C77AF","#F5D2A8","#B383B9","#D1352B"])'
p<-sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget",
              
              Value = "weight", NodeID = "name",nodeWidth =10,units = 'TWh',iterations = 0,#order node
              
              height=600,width=600,colourScale=colour,fontSize = 8,LinkGroup = "linkcol")

dnacol=c("#F5CFE4","#FCED82","#D2EBC8","#7DBFA7","#8FA4AE","#9B5B33","#EE934E","#3C77AF","#F5D2A8","#B383B9","#D1352B")
wnncol=c("#F5CFE4","#FCED82","#D2EBC8","#BBDD78","#7DBFA7","#9B5B33","#AECDE1","#EE934E","#3C77AF","#F5D2A8","#B383B9","#D1352B")
rnacol=c("#F5CFE4","#FCED82","#D2EBC8","#BBDD78","#7DBFA7","#8FA4AE","#9B5B33","#AECDE1","#EE934E","#3C77AF","#F5D2A8","#B383B9","#D1352B")

library(webshot)
setwd("Photo/")
saveNetwork(p,"./Photo/WNNcluster_sankey.html")
webshot("./Photo/WNNcluster_sankey.html" , "./Photo/WNNcluster_sankey.pdf")
webshot2::webshot("./Photo/WNNcluster_sankey.html" , "./Photo/WNNcluster_sankey.pdf")



