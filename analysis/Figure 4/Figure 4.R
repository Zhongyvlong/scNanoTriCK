########## Bivalent ggalluvial ------------------------------------
library(ggalluvial)
celltype_segment<-read.table("celltype_newsegments.txt",header = F)
colnames(celltype_segment)<-c("chr","start","end","state","line","cluster")
celltype_segment[1:5,]
celltype_segment$stateseq<-paste0(celltype_segment$chr,":",celltype_segment$start,"-",celltype_segment$end)
celltype_segment$group<-paste(celltype_segment$cluster,celltype_segment$state,sep = "_")
table(cellt)
celltype_state<-cbind(celltype_segment[celltype_segment$cluster=="C1Epi",c(7,4)], celltype_segment[celltype_segment$cluster=="C2PaE",c(7,4)], celltype_segment[celltype_segment$cluster=="C3VE",c(7,4)], celltype_segment[celltype_segment$cluster=="C4ExE",c(7,4)] )[,c(1,2,4,6,8)]
colnames(celltype_state)<-c("stateseq","C1Epi","C2PaE","C3VE","C4ExE")
celltype_state<-as.data.frame(apply(celltype_state,2,function(x){recode(x,"E1"="Others","E2"="Others","E3"="Others","E6"="Others")}))
celltype_state$Ncount<-apply(celltype_state,1,function(x){sum(str_count(x,pattern = "Others"))})
celltype_state<-subset(celltype_state,Ncount<4)

c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")
cols = c(C1Epi_E4="#2CA02C", C1Epi_E5="#D62728", C1Epi_E7="#9467BD", C1Epi_Others="lightgrey",C2PaE_E4="#2CA02C", C2PaE_E5="#D62728", C2PaE_E7="#9467BD", C2PaE_Others="lightgrey",C3VE_E4="#2CA02C", C3VE_E5="#D62728", C3VE_E7="#9467BD", C3VE_Others="lightgrey",C4ExE_E4="#2CA02C", C4ExE_E5="#D62728", C4ExE_E7="#9467BD", C4ExE_Others="lightgrey")

##### cluster biv state
state_comp<-celltype_state[,2:3]
state_comp$Ncount<-apply(state_comp,1,function(x){sum(str_count(x,pattern = "Others"))})
state_comp<-subset(state_comp,Ncount<2)
state_comp<-state_comp%>%group_by(C1Epi,C2PaE)%>%dplyr::summarise(number=n())
state_comp$C1Epi<-paste0("C1Epi_",state_comp$C1Epi)
state_comp$C2PaE<-paste0("C2PaE_",state_comp$C2PaE)
state_comp<-subset(state_comp,(C1Epi!="C1Epi_Others"&C2PaE!="C2PaE_Others")&!(C1Epi=="C1Epi_E4"&C2PaE=="C2PaE_E4")&!(C1Epi=="C1Epi_E7"&C2PaE=="C2PaE_E7") )

pdf("C1Epi_ToC2PaE_alluvium.pdf",height = 4.5,width = 6.5)
ggplot(state_comp, aes(y = number, axis1 = C1Epi,axis2 = C2PaE)) +
  geom_alluvium(aes(fill = C1Epi),alpha=0.4)+geom_stratum(aes(fill = after_stat(stratum)),width = 0.2)+geom_text(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_fill_manual(values = cols)+scale_x_discrete(limits = c("C1Epi", "C2PaE"), expand = c(0.05, 0.05))+theme_bw()
dev.off()

##### cluster biv related state 
state_comp<-celltype_state[,2:3]
state_comp$Ncount<-apply(state_comp,1,function(x){sum(str_count(x,pattern = "Others"))})
state_comp<-subset(state_comp,Ncount<2)
state_comp<-state_comp%>%group_by(C1Epi,C2PaE)%>%dplyr::summarise(number=n())
state_comp$C1Epi<-paste0("C1Epi_",state_comp$C1Epi)
state_comp$C2PaE<-paste0("C2PaE_",state_comp$C2PaE)
pdf("C1Epi_ToC2PaE_bivRelated_alluvium.pdf",height = 4.5,width = 6.5)
ggplot(state_comp, aes(y = number, axis1 = C1Epi,axis2 = C2PaE)) +
  geom_alluvium(aes(fill = C1Epi),alpha=0.4)+geom_stratum(aes(fill = after_stat(stratum)),width = 0.2)+geom_text(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_fill_manual(values = cols)+scale_x_discrete(limits = c("C1Epi", "C2PaE"), expand = c(0.05, 0.05))+theme_bw()
dev.off()


########## EK4_To_VBi regions ------------------------------------
##### identify EK4_To_VBi regions
celltype_state<-read.table("celltype_state.txt",header = T)
subbiv<-subset(celltype_state,C1Epi=="E7"&C3VE=="E5")
subbiv[,c("chr","start","end")]<-str_split_fixed(subbiv$stateseq,"[:-]",n=3)
subbiv[,6:7]<-apply(subbiv[,6:7],2,function(x){as.numeric(x)})

##### EK4_To_VBi signal in Epi and VE
biv_peak<-makeGRangesFromDataFrame(read.table("EK4_To_VBi.bed",header = F)%>%setNames(c("chr","start","end")), keep.extra.columns=T,ignore.strand=T,seqinfo=NULL)
fragments <- CreateFragmentObject("E5.0_dna_signac_fragfromflank.bed.gz")
mat<-FeatureMatrix(fragments = fragments,features = biv_peak )
mat<-mat[,intersect(e5.0_rna$cell,colnames(mat))]
mat<-mat[rowSums(mat)!=0,]
mat[1:5,1:5]
mat<-NormalizeData(mat)
mat_df<-reshape2::melt(as.matrix(mat))%>%setNames(c("bivRegion","cell","signal"))%>%left_join(e5.0_rna@meta.data[,c("cell","celltype","protein")],by="cell")
biv_score<-mat_df%>%group_by(bivRegion,celltype,protein)%>%dplyr::summarise(score=mean(signal))
biv_score_df<-reshape2::dcast(biv_score,bivRegion+celltype~protein,value.var = "score")
ggplot(subset(biv_score_df,celltype=="C1Epi"), aes(x=H3K27me3,y=H3K4me3))+geom_point()+theme_bw()+labs(title = "EK4_To_VBi state in Epi")+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed")+
  scale_x_continuous(limits = c(0,0.6))+scale_y_continuous(limits = c(0,0.6))

#####Epi/VE subset
mat<-as.data.frame(t(as.matrix(NormalizeData(readRDS("EK4_To_VBi_bivmat.RDS")))))%>%rownames_to_column(var = "cell")
e5.0_table<-left_join(e5.0_rna@meta.data[,c("cell","protein","celltype")], mat, by="cell" )
e5.0_table[is.na(e5.0_table)]<-0
subEpi<-subset(e5.0_table,celltype=="C1Epi")
subEpi<-data.frame(K4mean=colMeans(subEpi[subEpi$protein=="H3K4me3",4:4741]), K27mean=colMeans(subEpi[subEpi$protein=="H3K27me3",4:4741]), K9mean=colMeans(subEpi[subEpi$protein=="H3K9me3",4:4741]) ) %>% dplyr::mutate(K27rate=(K4mean+1)/(K27mean+1), K9rate=(K4mean+1)/(K9mean+1) )
subVE<-subset(e5.0_table,celltype=="C3VE")
subVE<-data.frame(K4mean=colMeans(subVE[subVE$protein=="H3K4me3",4:4741]), K27mean=colMeans(subVE[subVE$protein=="H3K27me3",4:4741]), K9mean=colMeans(subVE[subVE$protein=="H3K9me3",4:4741]) ) %>% dplyr::mutate(K27rate=(K27mean+1)/(K9mean+1), K4rate=(K4mean+1)/(K9mean+1) )

##### no rate
biv_score<-mat_df%>%group_by(bivRegion,celltype,protein)%>%dplyr::summarise(score=mean(signal))
biv_score_df<-reshape2::dcast(biv_score,bivRegion+celltype~protein,value.var = "score")
biv_score_df<-subset(biv_score_df,bivRegion%in%intersect(rownames(subset(subEpi,K4mean>0)), rownames(subset(subVE,K4mean>0&K27mean>0)) ))
#saveRDS(biv_score_df,"e5.0_bivscoredf_rate0.RDS")
##### rate > 1
biv_score<-mat_df%>%group_by(bivRegion,celltype,protein)%>%dplyr::summarise(score=mean(signal))
biv_score_df<-reshape2::dcast(biv_score,bivRegion+celltype~protein,value.var = "score")
biv_score_df<-subset(biv_score_df,bivRegion%in%intersect(rownames(subset(subEpi,K4mean>0&K27rate>1&K9rate>1)), rownames(subset(subVE,K4mean>0&K27mean>0&K27rate>1&K4rate>1)) ))
#saveRDS(biv_score_df,"e5.0_bivscoredf_rate1.RDS")
##### rate > 1.1
biv_score<-mat_df%>%group_by(bivRegion,celltype,protein)%>%dplyr::summarise(score=mean(signal))
biv_score_df<-reshape2::dcast(biv_score,bivRegion+celltype~protein,value.var = "score")
biv_score_df<-subset(biv_score_df,bivRegion%in%intersect(rownames(subset(subEpi,K4mean>0&K27rate>1.1&K9rate>1.1)), rownames(subset(subVE,K4mean>0&K27mean>0&K27rate>1.1&K4rate>1.1)) ))
#saveRDS(biv_score_df,"e5.0_bivscoredf_rate11.RDS")

biv_score<-mat_df%>%group_by(bivRegion,celltype,protein)%>%dplyr::summarise(score=mean(signal))
biv_score_df<-reshape2::dcast(biv_score,bivRegion+celltype~protein,value.var = "score")
biv_score_df<-subset(biv_score_df,bivRegion%in%intersect(rownames(subset(subEpi,K4mean>0.05&K27rate>1&K9rate>1)), rownames(subset(subVE,K4mean>0&K27mean>0&K27rate>1&K4rate>1)) ))
p1<-ggplot(subset(biv_score_df,celltype=="C1Epi"), aes(x=H3K27me3,y=H3K4me3))+geom_point(color="grey")+theme_bw()+labs(title = "E7ToE5 state in Epi")+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed")+
  scale_x_continuous(limits = c(0,0.6))+scale_y_continuous(limits = c(0,0.6))+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=4)+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")
p2<-ggplot(subset(biv_score_df,celltype=="C3VE"), aes(x=H3K27me3,y=H3K4me3))+geom_point(color="grey")+theme_bw()+labs(title = "E7ToE5 state in VE")+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed")+
  scale_x_continuous(limits = c(0,0.65))+scale_y_continuous(limits = c(0,0.65))+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=4)+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")
pdf("e5.0_E7ToE5_inEpiVE_rate1_dot.pdf",height = 4.5,width = 10)
p1+p2
dev.off()

##### GO
library(readxl)
golist<-read_excel("e5.0_bivscoredf_GO.xlsx",sheet = "e5.0_bivscoredf_rate11_geneGO")[c(1,2,7,9,19,24,28,29,31,32),]
golist<-golist[order(golist$Count),]%>%rowwise()%>%dplyr::mutate(GeneRatio=eval(parse(text=GeneRatio)))
golist$Description<-factor(golist$Description,levels = golist$Description)
pdf("e5.0_E7ToE5_inEpiVE_rate11_geneGO.pdf",height = 5.5,width = 6.5)
ggplot(golist,aes(x=GeneRatio,y=Description,color=-log10(p.adjust),size =Count ))+geom_point()+theme_bw()+#scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )+
  labs(title = "EK4_To_VBi peaks annotation in GO",x="GeneRatio",y="GO BP Term")+ scale_size(range=c(3, 8),breaks = pretty(golist[["Count"]], n=4))+
  scale_color_gradientn(colours = c("#327eba","#e06663"),limits=c(1.35,1.75),oob = scales::squish )
dev.off()

##### dual chip
E5.0_dualrna<-readRDS("Dual_RNA_assigned.rds")
E5.0_dualmeta<-read.table("Dual_DNA_RNA_paired_info_filtered.txt",header = T)
E5.0_dualmeta$dna_id<-str_split_fixed(E5.0_dualmeta$DNA_cell,"_Du",n=2)[,1]
E5.0_dualmeta$dual_type<-E5.0_dualrna@meta.data[match(E5.0_dualmeta$RNA_cell,colnames(E5.0_dualrna)),"seurat_clusters"]
E5.0_dualmeta$dna_cell<-substr(E5.0_dualmeta$dna_id,0,nchar(E5.0_dualmeta$dna_id)-2)
fragments_dual <- CreateFragmentObject("E5.0_dual_signac_fragfromflank.bed.gz")
mat_dual<-FeatureMatrix(fragments = fragments_dual,features = biv_peak )
mat_dual<-mat_dual[,intersect(E5.0_dualmeta$dna_id,colnames(mat_dual))]
mat_dual<-NormalizeData(mat_dual)
mat_dual_df<-reshape2::melt(as.matrix(mat_dual))%>%setNames(c("bivRegion","cell","signal"))%>%left_join(E5.0_dualmeta[,c("dna_id","dual_type","protein","dedup_read")],by=c("cell"="dna_id") )
mat_dual_df$dna_cell<-substr(mat_dual_df$cell,0,nchar(mat_dual_df$cell)-2)

biv_score_dual<-mat_dual_df%>%group_by(dna_cell,dual_type,protein)%>%dplyr::summarise(score=mean(signal))
biv_score_dual_df<-reshape2::dcast(biv_score_dual,dna_cell+dual_type~protein,value.var = "score")
ggplot(subset(biv_score_dual_df,dual_type=="0"), aes(x=H3K27me3,y=H3K4me3))+geom_point()+theme_bw()+labs(title = "E7ToE5 state in dual VE cells")+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed")+
  scale_x_continuous(limits = c(0,0.8))+scale_y_continuous(limits = c(0,0.8))

##### dual filter
mat_dual<-as.data.frame(t(as.matrix(NormalizeData(mat_dual))))%>%rownames_to_column(var = "dna_id")
dual_table<-left_join(E5.0_dualmeta[,c("dna_id","dual_type","dna_cell","protein","dedup_read")],mat_dual,by=c("dna_id"))
dual_table[is.na(dual_table)]<-0
rownames(dual_table)<-dual_table$dna_id
subregion<-data.frame(K4ncount=colSums(subset(dual_table,protein=="H3K4me3"&dual_type==0)[,6:4743]>0, na.rm = TRUE), K27ncount=colSums(subset(dual_table,protein=="H3K27me3"&dual_type==0)[,6:4743]>0, na.rm = TRUE) )
subcell<-data.frame(K4ncount=rowSums(subset(dual_table,protein=="H3K4me3"&dual_type==0)[,6:4743]>0, na.rm = TRUE), K27ncount=rowSums(subset(dual_table,protein=="H3K27me3"&dual_type==0)[,6:4743]>0, na.rm = TRUE) ) %>% dplyr::mutate(dna_cell=substr(rownames(.),0,nchar(rownames(.))-2))
dual_biv<-reshape2::melt(dual_table,id.vars=c("dna_id","dual_type","dna_cell","protein","dedup_read"),value.name = "signal",variable.name="bivRegion")%>%dplyr::filter(bivRegion%in% rownames(subset(subregion,K4ncount>=3&K27ncount>=3)) )
dual_biv<-dual_biv%>%group_by(dna_cell,dual_type,protein)%>%dplyr::summarise(score=mean(signal))
dual_biv<-reshape2::dcast(dual_biv,dna_cell+dual_type~protein,value.var = "score")
ggplot(subset(dual_biv,dual_type=="0"), aes(x=H3K27me3,y=H3K4me3))+geom_point()+theme_bw()+labs(title = "EK4_To_VBi region in dual VE cells")+geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed")+
  #scale_x_continuous(limits = c(0,0.15))+scale_y_continuous(limits = c(0,0.15))+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=4)+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")

pdf("dual_VE_E7ToE5_bivpeak_dot.pdf",height = 4.5,width = 5)
ggplot(subset(dual_biv,dual_type=="0")%>%dplyr::filter(rowSums(.[,3:4]>0)>1), aes(x=H3K27me3,y=H3K4me3))+geom_point(size=3)+theme_bw()+labs(title = "E7ToE5 state in dual VE cells")+geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed")+
  #scale_x_continuous(limits = c(0,0.15))+scale_y_continuous(limits = c(0,0.15))+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=4)+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")
dev.off()









