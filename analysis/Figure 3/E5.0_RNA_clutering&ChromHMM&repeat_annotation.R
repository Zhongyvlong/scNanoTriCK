library(Seurat) #5.0.1
library(Signac) #1.12.0
library(SeuratWrappers)
library(harmony)
library(cisTopic) #0.3.0
library(dplyr)
library(stringr)
library(ggplot2)
library(clustree)
library(GenomicRanges)
library(RColorBrewer)
set.seed(1234)
setwd("./Nanopore/")
#
########## E5.0 RNA clustering ------------------------------------
##### paired RNA and DNA
e5.0_rna<-readRDS("./Nanopore/E5.0_RNA_finally_cells.rds")
e5.0_rna<-JoinLayers(e5.0_rna)
cell_meta<-read.table("./Nanopore/DNA_RNA_paired_info.txt",header = T)
cell_meta<-subset(cell_meta,nFeature_RNA>800&nFeature_RNA<8000&nCount_RNA<40000&nCount_RNA>800&percent.mt<4&dedup_read<70000&dedup_read>1000&frip>0.2)
e5.0_rna<-e5.0_rna[,cell_meta$RNA_cell]
e5.0_meta<-as.data.frame(e5.0_rna@meta.data)
e5.0_rna@meta.data<-e5.0_meta
e5.0_rna[["percent.rb"]]<-PercentageFeatureSet(e5.0_rna,pattern = "^Rp[sl]")
e5.0_rna<-e5.0_rna[-grep("^Gm[0-9]",rownames(e5.0_rna)),]

##### Normalize
e5.0_rna <- NormalizeData(e5.0_rna, normalization.method = "LogNormalize", scale.factor = median(e5.0_rna$nCount_RNA))
e5.0_rna <- FindVariableFeatures(e5.0_rna, selection.method = "vst", nfeatures = 2000)
LabelPoints(plot = VariableFeaturePlot(e5.0_rna), points = head(VariableFeatures(e5.0_rna), 10), repel = TRUE)
e5.0_rna <- ScaleData(e5.0_rna, features = rownames(e5.0_rna),assay = "RNA",vars.to.regress = c("nCount_RNA","orig.ident")) 
e5.0_rna <- RunPCA(e5.0_rna, features = VariableFeatures(e5.0_rna),assay = "RNA")
ElbowPlot(e5.0_rna,n=50)
e5.0_rna <- RunUMAP(e5.0_rna, dims = 1:20,assay = "RNA",reduction = "pca",return.model = T)

DefaultAssay(e5.0_rna)<-"RNA"
set.seed(1234)
e5.0_rna <- RunHarmony(e5.0_rna,group.by.vars=c("batch"),reduction.use = "pca",dims.use=1:30,early_stop=F,max_iter=1)
ElbowPlot(e5.0_rna,reduction = "harmony",n=50)
e5.0_rna <- RunUMAP(e5.0_rna, dims = 1:10,min.dist = 0.1,reduction = "harmony",return.model = T)
e5.0_rna <- FindNeighbors(e5.0_rna, dims = 1:10,annoy.metric = "cosine",reduction = "harmony")
e5.0_rna <- FindClusters(e5.0_rna,algorithm = 1, resolution = 0.3)
p1<-DimPlot(e5.0_rna,label=T,pt.size = 0.5)+labs(title="1:10,0.3")
p2<-DimPlot(e5.0_rna,label=T,group.by = "batch",pt.size = 0.5)
p3<-DimPlot(e5.0_rna,label=T,group.by = "orig.ident",pt.size = 0.5)
p4<-FeaturePlot(e5.0_rna,features = "nCount_RNA",pt.size = 0.5)
(p1+p2)/(p3+p4)
e5.0_rna <- RunALRA(e5.0_rna,assay = "RNA")
genes_PaperShow <- c ( "Pou3f1","Dnmt3a","Lef1","Tdgf1","Pdgfra","Sox7","Lama1","Gata4","Col4a1","Nodal","Ttr","Hnf4a","Apob","Lhx1","Esrrb", "Cdx2","Eomes", "Krt8","Hand1","Gata3","Plac1","Cdcp1","Pecam1","Prl3d1","Cd72", "Cxcl2", "S100a9")
DotPlot(e5.0_rna,features = genes_PaperShow, cluster.idents =F,assay = "alra") + coord_flip() + scale_color_gradientn(colours = c("navy", "white", "firebrick"))+theme(axis.text.y=element_text(size=12))
DimPlot(e5.0_rna,label=T,group.by = "Histone",pt.size = 0.5)


##########chromHMM QC ------------------------------------
##### State Number
state_num<-read.table("state_num_incluster.txt",sep = "\t",header = T)
df<-reshape2::melt(as.matrix(state_num))
colnames(df)<-c("Cluster","State","Number")
df$Cluster<-as.character(df$Cluster)
df$State<-as.character(df$State)
df<-subset(df,State%in%c("E1","E2","E4","E5","E7"))
p<-ggplot(df,aes(color=Cluster))+geom_segment(aes(x=State, y=0, yend=log10(Number)),position = position_dodge(width = 0.7))+
  geom_point(aes(x=State, y=log10(Number)),position = position_dodge(width = 0.7),size=5)+
  scale_color_manual(values = c("#6FD3F2","#FCB359","#C49A6C","#B2798F"))+
  scale_y_continuous(limits = c(0,5))+
  theme_bw()+#theme(axis.text.x = element_text(angle=45,hjust = 1))+
  labs(title = "Log10 Number of Cluster States",x="State",y="Log10(State Number)")
pdf("StateNum_inCluster.pdf",height = 4,width = 8) 
print(p)
dev.off()

##### State Percentage
df<-reshape2::melt(as.matrix(state_num))
colnames(df)<-c("Cluster","State","Number")
df<-subset(df,State%in%c("E1","E2","E4","E5","E7"))

df<-df%>%group_by(Cluster)%>%mutate(Percentage=100*(Number/sum(Number)))
p<-ggplot(df,aes(color=State))+geom_segment(aes(x=Cluster, y=0, yend=Percentage),position = position_dodge(width = 0.6))+
  geom_point(aes(x=Cluster, y=Percentage),position = position_dodge(width = 0.6),size=5)+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD"))+
  scale_y_continuous(limits = c(0,60))+
  theme_bw()+#theme(axis.text.x = element_text(angle=45,hjust = 1))+
  labs(title = "Percentage of State in Cluster",x="State",y="Percentage(%)")
pdf("StatePct_inCluster.pdf",height = 4,width = 8) # 25.09.02
print(p)
dev.off()

df<-df%>%group_by(State)%>%mutate(Percentage=100*(Number/sum(Number)))
p<-ggplot(df,aes(color=Cluster))+geom_segment(aes(x=State, y=0, yend=Percentage),position = position_dodge(width = 0.8))+
  geom_point(aes(x=State, y=Percentage),position = position_dodge(width = 0.8),size=5)+
  scale_color_manual(values = c("#6FD3F2","#FCB359","#C49A6C","#B2798F"))+
  scale_y_continuous(limits = c(0,62))+
  theme_bw()+#theme(axis.text.x = element_text(angle=45,hjust = 1))+
  labs(title = "Percentage of Cluster in State",x="State",y="Percentage(%)")
pdf("ClusterPct_inState.pdf",height = 3,width = 8) # 25.09.02
print(p)
dev.off()

##### State genomic Annotation
BiocManager::install("ChIPseeker")
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_downstream = T)
options(ChIPseeker.ignore_promoter_subcategory = T)
peak<-makeGRangesFromDataFrame(C1Epi_state[C1Epi_state$state=="E7", 1:3],keep.extra.columns=T,ignore.strand=T,seqinfo=NULL)
peak_anno<-as.data.frame(annotatePeak(peak, tssRegion = c(-3000, 3000), TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,annoDb = "org.Mm.eg.db"))

result<-data.frame()
for(clu in c("C1Epi","C2PaE","C3VE","C4ExE")){
  state_seq<-read.table(paste0("Analysis/250829_chromHMM/celltypecon_hmm_2kb_7/celltype_segment/",clu,"_7_newsegments.txt"),header = F)
  colnames(state_seq)<-c("chr","start","end","state","line")
  for(st in c("E1","E2","E4","E5","E7")){
    peak<-state_seq[state_seq$state==st, 1:3]
    peakgs<-makeGRangesFromDataFrame(peak,keep.extra.columns=T,ignore.strand=T,seqinfo=NULL)
    peak_anno<-annotatePeak(peakgs, tssRegion = c(-3000, 3000), TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,annoDb = "org.Mm.eg.db")
    annostat<-data.frame(Feature=c("Promoter","Intron","Exon","UTR","Intergenic"), 
                         Number=c(sum(grepl("Promoter",peak_anno@anno$annotation)), sum(grepl("Intron",peak_anno@anno$annotation)), sum(grepl("Exon",peak_anno@anno$annotation)), sum(grepl("UTR",peak_anno@anno$annotation)), sum(grepl("Intergenic|Downstream",peak_anno@anno$annotation))), 
                         Cluster=clu, State=st )
    result<-rbind(result,annostat)
  }
}
result<-result%>%group_by(Cluster,State)%>%mutate(Frequency=round(100*(Number/sum(Number)),digits = 2) )
result<-read.table("State_genomeAnno_num_inCluster.txt",header = T)
result$Group<-paste(result$Cluster,result$State,sep = "_")
result$Feature<-factor(result$Feature,levels = c("Promoter","Intron","Exon","UTR","Intergenic"))
sum_data <- result %>% group_by(Group) %>% summarise(Total = sum(Number)) %>% mutate(num=label_number(accuracy = 0.1, scale = 1/1000, suffix = "K")(Total))
p<-ggplot()+geom_bar(data=result,aes(x=Group,y=Number,fill = Feature),stat = "identity",alpha=0.9)+
  scale_fill_manual(values = c("#399335","#256ea2","#603c87","#ea841e","#b4403e"))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1),plot.title = element_text(hjust = 0.5))+
  #geom_text(data = sum_data,aes(x=Group,y=Total,label=Total),vjust=-0.5,size=4)+
  geom_text(data = sum_data,aes(x=Group,y=Total,label=num),vjust=-0.5,size=4)+
  scale_y_continuous(limits = c(0,42000))+
  labs(title = "Peak annotation in State",x="Cluster","State Number")
pdf("State_genomeAnno_num_inCluster.pdf",height = 4,width = 8) 
print(p)
dev.off()

pct_data <- result %>% group_by(Group) %>% mutate(Label = sprintf("%.0f%%", Frequency),Position = 100-(cumsum(Frequency)-(Frequency/2)) ) %>%ungroup()
p<-ggplot()+geom_bar(data=result,aes(x=Group,y=Frequency,fill = Feature),stat = "identity",alpha=0.9)+
  scale_fill_manual(values = c("#399335","#256ea2","#603c87","#ea841e","#b4403e"))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1),plot.title = element_text(hjust = 0.5))+
  geom_text(data = pct_data,aes(x=Group,y=Position,label=Label),size=2.8,fontface="bold",color="white")+
  labs(title = "Peak annotation in State",x="Cluster","State Frequency")
pdf("State_genomeAnno_pct_inCluster.pdf",height = 4,width = 8) 
print(p)
dev.off()

##### State repeat Annotation
state_repeat<-read.table("C1Epi_overRepeat.bed",header = F,sep = "\t")
colnames(state_repeat)<-c("chr","start","end","state","line","genoName","genoStart","genoEnd","genoLeft","strand","repName","repClass","repFamily","overBase")
state_repeat<-subset(state_repeat,state%in%c("E1","E2","E4","E5","E7"))
state_repeat<-C1Epi_repeat%>%mutate(stateseq=paste0(chr,":",start,"-",end), repeatSign=case_when(grepl("DNA",repClass) ~ "DNA", grepl("LINE",repClass) ~ "LINE", grepl("SINE",repClass) ~ "SINE", grepl("LTR",repClass) ~ "LTR", grepl("Low_complexity",repClass) ~ "Low_complexity", repClass=="." ~ "NonRepeats", TRUE ~ "Other" )) 
Cstate_repgroup<-state_repeat%>%group_by(stateseq,state)%>%arrange(desc(overBase))%>%mutate(newSign=if_else(sum(overBase)<666,"NonRepeats","Repeats"))%>%ungroup() 
state_repgroup<-state_repgroup%>%group_by(stateseq,state,repeatSign,newSign)%>%arrange(desc(overBase))%>%summarise(overBase_g=sum(overBase))%>%group_by(stateseq,state)%>%dplyr::slice(which.max(overBase_g))%>%mutate(repeatSign=if_else(newSign=="NonRepeats","NonRepeats",repeatSign))%>%mutate(Cluster="C1Epi")%>%ungroup() 

state_rep<-rbind(C1Epi_repgroup%>%group_by(Cluster,state,repeatSign)%>%summarise(Number=n()), C2PaE_repgroup%>%group_by(Cluster,state,repeatSign)%>%summarise(Number=n()), C3VE_repgroup%>%group_by(Cluster,state,repeatSign)%>%summarise(Number=n()), C4ExE_repgroup%>%group_by(Cluster,state,repeatSign)%>%summarise(Number=n()))
state_rep$Group<-paste(state_rep$Cluster,state_rep$state,sep = "_")
state_rep<-state_rep%>%group_by(Group)%>%mutate(Frequency=100*(Number/sum(Number)))
state_rep$repeatSign<-factor(state_rep$repeatSign,levels = c("DNA","LINE","SINE","LTR","Low_complexity","Other","NonRepeats"))

p<-ggplot()+geom_bar(data=state_rep,aes(x=Group,y=Frequency,fill = repeatSign),stat = "identity",alpha=0.9,width = 0.5)+
  scale_fill_manual(values = c("#a0c9e5","#E8C7D7","#FCF6B5","#99CCCC","#e867a7","#FBB463","#BCBEC0"))+ # "#a0c9e5","#ffc94a","#ff914e","#95CC5E","#e867a7","#947ff2","#9d5c39"
  theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1),plot.title = element_text(hjust = 0.5))+
  labs(title = "Repeat annotation in Cluster",x="Cluster","Repeat Frequency")
pdf("State_repeatAnno_pct_inCluster_allonethird.pdf",height = 4,width = 8) 
print(p)
dev.off()

