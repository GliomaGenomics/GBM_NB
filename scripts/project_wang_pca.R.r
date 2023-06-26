library(ggplot2)
library(ggrepel)
library(patchwork)
library(gridExtra)
library(cowplot)

pdf("wang/wang-diaz_C1.pdf", width=6, height=6)


pc <- read.table("wang/GBM_C1_PCA_loadings.txt", header=TRUE, row.names=1,check.names=FALSE)
# pc <- read.table("seurat/seurat_pca_weights.txt", header=TRUE, row.names=1,check.names=FALSE)
data <- read.table("seurat/seurat_wc_norm.txt", header=TRUE, row.names=1,check.names=FALSE)
data <- read.table("seurat/seurat_wc_scaled.txt", header=TRUE, row.names=1,check.names=FALSE)
samples<-colnames(data)
pc_genes<-data.frame(row.names(pc))
datafilt<-merge(data,pc_genes,by.x='row.names',by.y=1,all.y=TRUE)
row.names(datafilt)<-datafilt$Row.names
data<-subset(datafilt, select = -Row.names)
data[is.na(data)] <- 0
new=t(data)
pc<-as.matrix(pc)
pcx<- new%*%pc

intercept=-3.19

meta<-read.table('data/nano-biopsy_whole-cell/star_gfp_min0.4_uniq/bam_stats/metrics_combined_wc_0.4.txt',sep='\t',header=1)
row.names(meta)<-meta$sample
meta<-meta[!is.na(meta$neftel),]
meta<-meta[meta$category =="M059K_T" |meta$category =="M059K_U",]
samples<-meta$sample
pcxf<-pcx[samples,] #needed for correct order as well as filtering
comb1<-as.data.frame(cbind(pcxf[,1],meta[samples,]$category))
comb1<-split(comb1,comb1$V2)
eqvar_pval1<-t.test(as.numeric(noquote(comb1[[1]]$V1)),as.numeric(noquote(comb1[[2]]$V1)),var.equal = T)$p.value
wilcox_pval1<-wilcox.test(as.numeric(noquote(comb1[[1]]$V1)),as.numeric(noquote(comb1[[2]]$V1)))$p.value
comb2<-as.data.frame(cbind(pcxf[,2],meta[samples,]$category))
comb2<-split(comb2,comb2$V2)
eqvar_pval2<-t.test(as.numeric(noquote(comb2[[1]]$V1)),as.numeric(noquote(comb2[[2]]$V1)),var.equal = T)$p.value
wilcox_pval2<-wilcox.test(as.numeric(noquote(comb2[[1]]$V1)),as.numeric(noquote(comb2[[2]]$V1)))$p.value
ggplot(as.data.frame(pcxf), aes(x=pcxf[,1], y=pcxf[,2])) +
  theme_classic() +
  geom_point(aes(colour = meta[samples,]$category), stroke=0,size = 4, alpha=1/2) +
  theme(legend.title=element_blank()) +
  labs(x = "PC1", y= "PC2")+
  geom_vline(xintercept=intercept,linetype=2)
ggplot(as.data.frame(pcxf), aes(x=meta[samples,]$category, y=pcxf[,1])) +  
  geom_boxplot()+
  ggtitle(label = paste("P vlaue for t-test or wilcox:",eqvar_pval1,", ",wilcox_pval1,sep=""))
ggplot(as.data.frame(pcxf), aes(x=meta[samples,]$category, y=pcxf[,2])) +  
  geom_boxplot()+
  ggtitle(label = paste("P vlaue for t-test or wilcox:",eqvar_pval2,", ",wilcox_pval2,sep=""))
















#########nb

pc <- read.table("wang/GBM_C1_PCA_loadings.txt", header=TRUE, row.names=1,check.names=FALSE)
# pc <- read.table("seurat/seurat_pca_weights.txt", header=TRUE, row.names=1,check.names=FALSE)
# pc <- read.table("seurat/seurat_pca_weights_a0.25.txt", header=TRUE, row.names=1,check.names=FALSE)
# pc <- read.table("seurat/seurat_pca_weights_u0.25.txt", header=TRUE, row.names=1,check.names=FALSE)
# pc <- read.table("seurat/seurat_pca_weights_scaled_regressed.txt", header=TRUE, row.names=1,check.names=FALSE)
# pc <- read.table("seurat/seurat_pca_weights_sct_scaled.txt", header=TRUE, row.names=1,check.names=FALSE)
# data <- read.table("seurat/seurat_nb_scaled_regressed.txt", header=TRUE, row.names=1,check.names=FALSE)
# data <- read.table("seurat/seurat_nb_sct_norm.txt", header=TRUE, row.names=1,check.names=FALSE)
# data <- read.table("seurat/seurat_nb_sct_scaled.txt", header=TRUE, row.names=1,check.names=FALSE)

data <- read.table("seurat/seurat_nb_norm.txt", header=TRUE, row.names=1,check.names=FALSE)
data <- read.table("seurat/seurat_nb_scaled.txt", header=TRUE, row.names=1,check.names=FALSE)

samples<-colnames(data)
pc_genes<-data.frame(row.names(pc))
datafilt<-merge(data,pc_genes,by.x='row.names',by.y=1,all.y=TRUE)
row.names(datafilt)<-datafilt$Row.names
data<-subset(datafilt, select = -Row.names)
data[is.na(data)] <- 0
new=t(data)
pc<-as.matrix(pc)
pcx<- new%*%pc

intercept=3.51


meta<-read.table('data/nano-biopsy/star_gfp_min0.1_uniq/bam_stats/metrics_combined_nb_0.1.txt',sep='\t',header=1)
row.names(meta)<-meta$sample
meta<-meta[!is.na(meta$neftel),]
meta<-meta[meta$best_replicate =="Y",]
meta<-meta[meta$initial_treatment =="T" & !is.na(meta$initial_treatment),]
meta<-meta[meta$category=="I",]
meta$initial_survived[meta$initial_survived=='Y']<-'initial_treated_survived'
meta$initial_survived[meta$initial_survived=='N']<-'initial_treated_died'
samples<-meta$sample
pcxf<-pcx[samples,] #needed for correct order as well as filtering
comb1<-as.data.frame(cbind(pcxf[,1],meta[samples,]$initial_survived))
comb1<-split(comb1,comb1$V2)
eqvar_pval1<-t.test(as.numeric(noquote(comb1[[1]]$V1)),as.numeric(noquote(comb1[[2]]$V1)),var.equal = T)$p.value
wilcox_pval1<-wilcox.test(as.numeric(noquote(comb1[[1]]$V1)),as.numeric(noquote(comb1[[2]]$V1)))$p.value
comb2<-as.data.frame(cbind(pcxf[,2],meta[samples,]$initial_survived))
comb2<-split(comb2,comb2$V2)
eqvar_pval2<-t.test(as.numeric(noquote(comb2[[1]]$V1)),as.numeric(noquote(comb2[[2]]$V1)),var.equal = T)$p.value
wilcox_pval2<-wilcox.test(as.numeric(noquote(comb2[[1]]$V1)),as.numeric(noquote(comb2[[2]]$V1)))$p.value
ggplot(as.data.frame(pcxf), aes(x=pcxf[,1], y=pcxf[,2])) +
  theme_classic() +
  geom_point(aes(colour = meta[samples,]$initial_survived), stroke=0,size = 4, alpha=1/2) +
  theme(legend.title=element_blank()) +
  labs(x = "PC1", y= "PC2")+
  geom_vline(xintercept=intercept,linetype=2)
ggplot(as.data.frame(pcxf), aes(x=meta[samples,]$initial_survived, y=pcxf[,1])) +  
  geom_boxplot()+
  ggtitle(label = paste("P vlaue for t-test or wilcox:",eqvar_pval1,", ",wilcox_pval1,sep=""))
ggplot(as.data.frame(pcxf), aes(x=meta[samples,]$initial_survived, y=pcxf[,2])) +  
  geom_boxplot()+
  ggtitle(label = paste("P vlaue for t-test or wilcox:",eqvar_pval2,", ",wilcox_pval2,sep=""))




meta<-read.table('data/nano-biopsy/star_gfp_min0.1_uniq/bam_stats/metrics_combined_nb_0.1.txt',sep='\t',header=1)
row.names(meta)<-meta$sample
meta<-meta[!is.na(meta$neftel),]
meta<-meta[meta$best_replicate =="Y",]
meta<-meta[meta$initial_treatment =="U" & !is.na(meta$initial_treatment),]
meta<-meta[meta$category=="I",]
meta$initial_survived[meta$initial_survived=='Y']<-'initial_untreated_survived'
meta$initial_survived[meta$initial_survived=='N']<-'initial_untreated_died'
samples<-meta$sample
pcxf<-pcx[samples,] #needed for correct order as well as filtering
comb1<-as.data.frame(cbind(pcxf[,1],meta[samples,]$initial_survived))
comb1<-split(comb1,comb1$V2)
eqvar_pval1<-t.test(as.numeric(noquote(comb1[[1]]$V1)),as.numeric(noquote(comb1[[2]]$V1)),var.equal = T)$p.value
wilcox_pval1<-wilcox.test(as.numeric(noquote(comb1[[1]]$V1)),as.numeric(noquote(comb1[[2]]$V1)))$p.value
comb2<-as.data.frame(cbind(pcxf[,2],meta[samples,]$initial_survived))
comb2<-split(comb2,comb2$V2)
eqvar_pval2<-t.test(as.numeric(noquote(comb2[[1]]$V1)),as.numeric(noquote(comb2[[2]]$V1)),var.equal = T)$p.value
wilcox_pval2<-wilcox.test(as.numeric(noquote(comb2[[1]]$V1)),as.numeric(noquote(comb2[[2]]$V1)))$p.value
ggplot(as.data.frame(pcxf), aes(x=pcxf[,1], y=pcxf[,2])) +
  theme_classic() +
  geom_point(aes(colour = meta[samples,]$initial_survived), stroke=0,size = 4, alpha=1/2) +
  theme(legend.title=element_blank()) +
  labs(x = "PC1", y= "PC2")+
  geom_vline(xintercept=intercept,linetype=2)
ggplot(as.data.frame(pcxf), aes(x=meta[samples,]$initial_survived, y=pcxf[,1])) +  
  geom_boxplot()+
  ggtitle(label = paste("P vlaue for t-test or wilcox:",eqvar_pval1,", ",wilcox_pval1,sep=""))
ggplot(as.data.frame(pcxf), aes(x=meta[samples,]$initial_survived, y=pcxf[,2])) +  
  geom_boxplot()+
  ggtitle(label = paste("P vlaue for t-test or wilcox:",eqvar_pval2,", ",wilcox_pval2,sep=""))


meta<-read.table('data/nano-biopsy/star_gfp_min0.1_uniq/bam_stats/metrics_combined_nb_0.1.txt',sep='\t',header=1)
row.names(meta)<-meta$sample
meta<-meta[!is.na(meta$neftel),]
meta<-meta[meta$best_replicate =="Y",]
meta<-meta[meta$MEDIAN_CV_COVERAGE <3,]
meta<-meta[meta$category =="L_T" | meta$category =="L_U",]
meta$category[meta$category=='L_T']<-'longitudinal_treated'
meta$category[meta$category=='L_U']<-'longitudinal_untreated'
samples<-meta$sample
pcxf<-pcx[samples,] #needed for correct order as well as filtering
comb1<-as.data.frame(cbind(pcxf[,1],meta[samples,]$category))
comb1<-split(comb1,comb1$V2)
eqvar_pval1<-t.test(as.numeric(noquote(comb1[[1]]$V1)),as.numeric(noquote(comb1[[2]]$V1)),var.equal = T)$p.value
wilcox_pval1<-wilcox.test(as.numeric(noquote(comb1[[1]]$V1)),as.numeric(noquote(comb1[[2]]$V1)))$p.value
comb2<-as.data.frame(cbind(pcxf[,2],meta[samples,]$category))
comb2<-split(comb2,comb2$V2)
eqvar_pval2<-t.test(as.numeric(noquote(comb2[[1]]$V1)),as.numeric(noquote(comb2[[2]]$V1)),var.equal = T)$p.value
wilcox_pval2<-wilcox.test(as.numeric(noquote(comb2[[1]]$V1)),as.numeric(noquote(comb2[[2]]$V1)))$p.value
ggplot(as.data.frame(pcxf), aes(x=pcxf[,1], y=pcxf[,2])) +
  theme_classic() +
  geom_point(aes(colour = meta[samples,]$category), stroke=0,size = 4, alpha=1/2) +
  theme(legend.title=element_blank()) +
  labs(x = "PC1", y= "PC2")+
  geom_vline(xintercept=intercept,linetype=2)
ggplot(as.data.frame(pcxf), aes(x=meta[samples,]$category, y=pcxf[,1])) +  
  geom_boxplot()+
  ggtitle(label = paste("P vlaue for t-test or wilcox:",eqvar_pval1,", ",wilcox_pval1,sep=""))
ggplot(as.data.frame(pcxf), aes(x=meta[samples,]$category, y=pcxf[,2])) +  
  geom_boxplot()+
  ggtitle(label = paste("P vlaue for t-test or wilcox:",eqvar_pval2,", ",wilcox_pval2,sep=""))


meta<-read.table('data/nano-biopsy/star_gfp_min0.1_uniq/bam_stats/metrics_combined_nb_0.1.txt',sep='\t',header=1)
row.names(meta)<-meta$sample
meta<-meta[!is.na(meta$neftel),]
meta<-meta[meta$best_replicate =="Y",]
meta<-meta[!is.na(meta$neftel) & !is.na(meta$initial_treatment),]
meta<-meta[meta$initial_treatment =="T",]
meta<-meta[meta$category =="L_T" | (meta$category =="I" & meta$initial_survived=='Y'),]
meta$category[meta$category=='L_T']<-'longitudinal_treated'
meta$category[meta$category=='I']<-'initial_treated_survived'
samples<-meta$sample
pcxf<-pcx[samples,] #needed for correct order as well as filtering
comb1<-as.data.frame(cbind(pcxf[,1],meta[samples,]$category))
comb1<-split(comb1,comb1$V2)
eqvar_pval1<-t.test(as.numeric(noquote(comb1[[1]]$V1)),as.numeric(noquote(comb1[[2]]$V1)),var.equal = T)$p.value
wilcox_pval1<-wilcox.test(as.numeric(noquote(comb1[[1]]$V1)),as.numeric(noquote(comb1[[2]]$V1)))$p.value
comb2<-as.data.frame(cbind(pcxf[,2],meta[samples,]$category))
comb2<-split(comb2,comb2$V2)
eqvar_pval2<-t.test(as.numeric(noquote(comb2[[1]]$V1)),as.numeric(noquote(comb2[[2]]$V1)),var.equal = T)$p.value
wilcox_pval2<-wilcox.test(as.numeric(noquote(comb2[[1]]$V1)),as.numeric(noquote(comb2[[2]]$V1)))$p.value
ggplot(as.data.frame(pcxf), aes(x=pcxf[,1], y=pcxf[,2])) +
  theme_classic() +
  geom_point(aes(colour = meta[samples,]$category), stroke=0,size = 4, alpha=1/2) +
  theme(legend.title=element_blank()) +
  labs(x = "PC1", y= "PC2")+
  geom_vline(xintercept=intercept,linetype=2)
ggplot(as.data.frame(pcxf), aes(x=meta[samples,]$category, y=pcxf[,1])) +  
  geom_boxplot()+
  ggtitle(label = paste("P vlaue for t-test or wilcox:",eqvar_pval1,", ",wilcox_pval1,sep=""))
ggplot(as.data.frame(pcxf), aes(x=meta[samples,]$category, y=pcxf[,2])) +  
  geom_boxplot()+
  ggtitle(label = paste("P vlaue for t-test or wilcox:",eqvar_pval2,", ",wilcox_pval2,sep=""))

meta<-read.table('data/nano-biopsy/star_gfp_min0.1_uniq/bam_stats/metrics_combined_nb_0.1.txt',sep='\t',header=1)
row.names(meta)<-meta$sample
meta<-meta[!is.na(meta$neftel),]
meta<-meta[meta$best_replicate =="Y",]
meta<-meta[!is.na(meta$neftel) & !is.na(meta$initial_treatment),]
meta<-meta[meta$initial_treatment =="U",]
meta<-meta[meta$category =="L_U" | (meta$category =="I" & meta$initial_survived=='Y'),]
meta$category[meta$category=='L_U']<-'longitudinal_untreated'
meta$category[meta$category=='I']<-'initial_untreated_survived'
samples<-meta$sample
pcxf<-pcx[samples,] #needed for correct order as well as filtering
comb1<-as.data.frame(cbind(pcxf[,1],meta[samples,]$category))
comb1<-split(comb1,comb1$V2)
eqvar_pval1<-t.test(as.numeric(noquote(comb1[[1]]$V1)),as.numeric(noquote(comb1[[2]]$V1)),var.equal = T)$p.value
wilcox_pval1<-wilcox.test(as.numeric(noquote(comb1[[1]]$V1)),as.numeric(noquote(comb1[[2]]$V1)))$p.value
comb2<-as.data.frame(cbind(pcxf[,2],meta[samples,]$category))
comb2<-split(comb2,comb2$V2)
eqvar_pval2<-t.test(as.numeric(noquote(comb2[[1]]$V1)),as.numeric(noquote(comb2[[2]]$V1)),var.equal = T)$p.value
wilcox_pval2<-wilcox.test(as.numeric(noquote(comb2[[1]]$V1)),as.numeric(noquote(comb2[[2]]$V1)))$p.value
ggplot(as.data.frame(pcxf), aes(x=pcxf[,1], y=pcxf[,2])) +
  theme_classic() +
  geom_point(aes(colour = meta[samples,]$category), stroke=0,size = 4, alpha=1/2) +
  theme(legend.title=element_blank()) +
  labs(x = "PC1", y= "PC2")+
  geom_vline(xintercept=intercept,linetype=2)
ggplot(as.data.frame(pcxf), aes(x=meta[samples,]$category, y=pcxf[,1])) +  
  geom_boxplot()+
  ggtitle(label = paste("P vlaue for t-test or wilcox:",eqvar_pval1,", ",wilcox_pval1,sep=""))
ggplot(as.data.frame(pcxf), aes(x=meta[samples,]$category, y=pcxf[,2])) +  
  geom_boxplot()+
  ggtitle(label = paste("P vlaue for t-test or wilcox:",eqvar_pval2,", ",wilcox_pval2,sep=""))



meta<-read.table('data/nano-biopsy/star_gfp_min0.1_uniq/bam_stats/metrics_combined_nb_0.1.txt',sep='\t',header=1)
row.names(meta)<-meta$sample
meta<-meta[!is.na(meta$neftel),]
meta<-meta[!is.na(meta$neftel) & !is.na(meta$initial_treatment),]
meta<-meta[meta$initial_treatment =="T",]
meta<-meta[meta$best_replicate =="Y",]
meta<-meta[meta$category =="L_T" | (meta$category =="I" & meta$initial_survived=='Y'),]
meta<-meta[meta$category =="L_T" | (meta$category =="I" & meta$initial_survived=='Y'),]
meta$category[meta$category=='L_T']<-'longitudinal_treated'
meta$category[meta$category=='I']<-'initial_treated_survived'
samples<-meta$sample
pcxf<-pcx[samples,] #needed for correct order as well as filtering
pairs<-read.table('meta/biopsy_pairs.txt',header=TRUE)
vals<- data.frame(matrix(ncol = 4, nrow = 0))
colnames(vals)<-c('i_1','i_2','r_1','r_2')
lins<- data.frame(matrix(ncol = 4, nrow = 0))
colnames(lins)<-c('xx','yy','long','direction')
for(i in 1:nrow(pairs)) {
    initial <- pairs[i,2]
    recurrent <- pairs[i,1]
        if (initial %in% rownames(meta) & recurrent %in% rownames(meta)){
        i_x<-pcxf[,1][initial][[1]]
        i_y<-pcxf[,2][initial][[1]]
        r_x<-pcxf[,1][recurrent][[1]]
        r_y<-pcxf[,2][recurrent][[1]]
        tre<-ifelse(i_x<r_x,'mesenchymal','proneural')
        lins[nrow(lins) + 1,]<-list(i_x, i_y, recurrent,tre)
        lins[nrow(lins) + 1,]<-list(r_x, r_y, recurrent,tre)
        vals[nrow(vals) + 1,]<-list(i_x, i_y, r_x,r_y)
}
}
pval1<-t.test(vals$i_1,vals$r_1,paired = T)$p.value
pval2<-t.test(vals$i_2,vals$r_2,paired = T)$p.value

l<-ggplot(lins, aes(x = xx, y = yy,group=long,color=direction)) +
  geom_path(arrow = arrow(length = unit(0.35, "cm"),angle=20,type="closed"))+
  labs(color = "Direction")+ 
  scale_colour_manual(values=c("royalblue2", "red3"),breaks=c("mesenchymal","proneural"))+
  ggtitle(label = paste("P vlaues for paired t-test:",pval1,pval2,sep=""))

p<-ggplot(lins, aes(x = xx, y = yy,group=long,color=direction)) +
  geom_point(data=as.data.frame(pcxf),aes(x=pcxf[,1], y=pcxf[,2],color = meta[samples,]$category), stroke=0,size = 4, alpha=3/4,inherit.aes = FALSE)+
  labs(color = "Biopsy")+ 
  scale_color_manual(values=c("green4", "purple2"),breaks=c('longitudinal_treated','initial_treated_survived'))
lp<-ggplot(lins, aes(x = xx, y = yy,group=long,color=direction)) +
  geom_path(arrow = arrow(length = unit(0.35, "cm"),angle=20,type="closed"))+
  theme_classic() +
  geom_point(data=as.data.frame(pcxf),aes(x=pcxf[,1], y=pcxf[,2],color = meta[samples,]$category), stroke=0,size = 4, alpha=3/4,inherit.aes = FALSE) +
  labs(x = "PC1", y= "PC2")+
  geom_vline(xintercept=intercept,linetype=2)+
  theme(legend.position = "none")+ 
  scale_color_manual(values=c("green4", "purple2","royalblue2", "red3"),breaks=c('longitudinal_treated','initial_treated_survived',"mesenchymal","proneural"))
grid.arrange(l,p,nrow = 2)

meta<-read.table('data/nano-biopsy/star_gfp_min0.1_uniq/bam_stats/metrics_combined_nb_0.1.txt',sep='\t',header=1)
row.names(meta)<-meta$sample
meta<-meta[!is.na(meta$neftel),]
meta<-meta[!is.na(meta$neftel) & !is.na(meta$initial_treatment),]
meta<-meta[meta$initial_treatment =="U",]
meta<-meta[meta$best_replicate =="Y",]
meta<-meta[meta$category =="L_U" | (meta$category =="I" & meta$initial_survived=='Y'),]
meta<-meta[meta$category =="L_U" | (meta$category =="I" & meta$initial_survived=='Y'),]
meta$category[meta$category=='L_U']<-'longitudinal_untreated'
meta$category[meta$category=='I']<-'initial_untreated_survived'
samples<-meta$sample
pcxf<-pcx[samples,] #needed for correct order as well as filtering
pairs<-read.table('meta/biopsy_pairs.txt',header=TRUE)
vals<- data.frame(matrix(ncol = 4, nrow = 0))
colnames(vals)<-c('i_1','i_2','r_1','r_2')
lins<- data.frame(matrix(ncol = 4, nrow = 0))
colnames(lins)<-c('xx','yy','long','direction')
for(i in 1:nrow(pairs)) {
    initial <- pairs[i,2]
    recurrent <- pairs[i,1]
        if (initial %in% rownames(meta) & recurrent %in% rownames(meta)){
        i_x<-pcxf[,1][initial][[1]]
        i_y<-pcxf[,2][initial][[1]]
        r_x<-pcxf[,1][recurrent][[1]]
        r_y<-pcxf[,2][recurrent][[1]]
        tre<-ifelse(i_x<r_x,'mesenchymal','proneural')
        lins[nrow(lins) + 1,]<-list(i_x, i_y, recurrent,tre)
        lins[nrow(lins) + 1,]<-list(r_x, r_y, recurrent,tre)
        vals[nrow(vals) + 1,]<-list(i_x, i_y, r_x,r_y)
}
}

pval1<-t.test(vals$i_1,vals$r_1,paired = T)$p.value
pval2<-t.test(vals$i_2,vals$r_2,paired = T)$p.value

l<-ggplot(lins, aes(x = xx, y = yy,group=long,color=direction)) +
  geom_path(arrow = arrow(length = unit(0.35, "cm"),angle=20,type="closed"))+
  labs(color = "Direction")+ 
  scale_colour_manual(values=c("royalblue2", "red3"),breaks=c("mesenchymal","proneural"))+
  ggtitle(label = paste("P vlaues for paired t-test:",pval1,pval2,sep=""))

p<-ggplot(lins, aes(x = xx, y = yy,group=long,color=direction)) +
  geom_point(data=as.data.frame(pcxf),aes(x=pcxf[,1], y=pcxf[,2],color = meta[samples,]$category), stroke=0,size = 4, alpha=3/4,inherit.aes = FALSE)+
  labs(color = "Biopsy")+ 
  scale_color_manual(values=c("green4", "purple2"),breaks=c('longitudinal_untreated','initial_untreated_survived'))
lp<-ggplot(lins, aes(x = xx, y = yy,group=long,color=direction)) +
  geom_path(arrow = arrow(length = unit(0.35, "cm"),angle=20,type="closed"))+
  theme_classic() +
  geom_point(data=as.data.frame(pcxf),aes(x=pcxf[,1], y=pcxf[,2],color = meta[samples,]$category), stroke=0,size = 4, alpha=3/4,inherit.aes = FALSE) +
  labs(x = "PC1", y= "PC2")+
  geom_vline(xintercept=intercept,linetype=2)+
  theme(legend.position = "none")+ 
  scale_color_manual(values=c("green4", "purple2","royalblue2", "red3"),breaks=c('longitudinal_untreated','initial_untreated_survived',"mesenchymal","proneural"))
grid.arrange(l,p,nrow = 2)

dev.off()

