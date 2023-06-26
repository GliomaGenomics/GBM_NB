

library(gridExtra)
library(grid)
library(Seurat)
library(cowplot)
library(ggplot2)
set.seed(123)

#get list of genes for protein coding, non ribosomal
pc<-read.table('downloaded_data/gencode/protein_coding_v27_names.txt',header=FALSE)
non_ribo_pc<-pc$V1[grep("^RP[S,L]",pc$V1,invert = TRUE)]

# read in counts and metadata and create seurat object
nb_meta<-read.table('data/nano-biopsy/star_gfp_min0.66_minNmatch3/bam_stats/nb_metrics.txt',sep='\t',header=1)
nb_counts<-ReadMtx(mtx='data/nano-biopsy/star_gfp_min0.66_minNmatch3/output/Gene/raw/nano-biopsy_matrix.mtx',
 features='data/nano-biopsy/star_gfp_min0.66_minNmatch3/output/Gene/raw/nano-biopsy_genes.tsv',
 cells='data/nano-biopsy/star_gfp_min0.66_minNmatch3/output/Gene/raw/nano-biopsy_barcodes.tsv')
row.names(nb_meta)<-nb_meta$sample
# comment out to get prefiltered metrics
nb_counts<-nb_counts[rownames(nb_counts) %in% non_ribo_pc,]
nb_data <- CreateSeuratObject(counts = nb_counts,min.cells=1)
nb_meta<-nb_meta[row.names(nb_data[[]]),]#ensure metadata is in the same order as seurat object
nb_data[[colnames(nb_meta)]]<-nb_meta[colnames(nb_meta)]#add in metadata
Idents(nb_data)<-"category"
nb_data[["percent.mt"]] <- PercentageFeatureSet(nb_data, pattern= "^MT-")

wc_meta<-read.table('data/nano-biopsy_whole-cell/star_gfp_min0.66_minNmatch3/bam_stats/wc_metrics.txt',sep='\t',header=1)
wc_counts<-ReadMtx(mtx='data/nano-biopsy_whole-cell/star_gfp_min0.66_minNmatch3/output/Gene/raw/nano-biopsy_whole-cell_matrix.mtx',
 features='data/nano-biopsy_whole-cell/star_gfp_min0.66_minNmatch3/output/Gene/raw/nano-biopsy_whole-cell_genes.tsv',
 cells='data/nano-biopsy_whole-cell/star_gfp_min0.66_minNmatch3/output/Gene/raw/nano-biopsy_whole-cell_barcodes.tsv')
row.names(wc_meta)<-wc_meta$sample
# comment out to get prefiltered metrics
wc_counts<-wc_counts[rownames(wc_counts) %in% non_ribo_pc,]
wc_data <- CreateSeuratObject(counts = wc_counts)
wc_meta<-wc_meta[row.names(wc_data[[]]),]#ensure metadata is in the same order as seurat object
wc_data[[colnames(wc_meta)]]<-wc_meta[colnames(wc_meta)]#add in metadata
Idents(wc_data)<-"category"
wc_data[["percent.mt"]] <- PercentageFeatureSet(wc_data, pattern= "^MT-")


# filter cells and genes
nb_data[['log_num_mrna']]<-log(nb_data[['mapped']]*nb_data[['PCT_MRNA_BASES']])
nb_data[['proportion_mapped']]<-nb_data[['mapped']]/nb_data[['total']]
wc_data[['log_num_mrna']]<-log(wc_data[['mapped']]*wc_data[['PCT_MRNA_BASES']])
wc_data[['proportion_mapped']]<-wc_data[['mapped']]/wc_data[['total']]
nb_data[['num_mrna']]<-nb_data[['mapped']]*nb_data[['PCT_MRNA_BASES']]
wc_data[['num_mrna']]<-wc_data[['mapped']]*wc_data[['PCT_MRNA_BASES']]

nb_data <- subset(nb_data,idents = c("L_U","L_T","I"))
wc_data <- subset(wc_data,idents = c("M059K_T","M059K_U"))

# # comment out to get prefiltered
nb_data <- subset(nb_data, subset = nFeature_RNA>150 &percent.mt<30 & PCT_RIBOSOMAL_BASES<0.1&num_mrna>0&pass_qc_rep_filter!="N")
nb_data <- subset(nb_data,features=row.names(nb_data[['RNA']]@counts)[rowSums(nb_data[['RNA']]@counts >=4)>=2])
nb_data <- subset(nb_data,features=row.names(nb_data[['RNA']]@counts)[rownames(nb_data[['RNA']]@counts) %in% non_ribo_pc])
wc_data <- subset(wc_data, subset = nFeature_RNA>150 &percent.mt<30 & PCT_RIBOSOMAL_BASES<0.1)#&PCT_MRNA_BASES>0.6&MEDIAN_CV_COVERAGE<1.25
wc_data <- subset(wc_data,features=row.names(wc_data[['RNA']]@counts)[rowSums(wc_data[['RNA']]@counts >=4)>=2])
wc_data <- subset(wc_data,features=row.names(wc_data[['RNA']]@counts)[rownames(wc_data[['RNA']]@counts) %in% non_ribo_pc])

combined<-merge(nb_data,wc_data,add.cell.ids = c("nano-biopsy","whole-cell"))

combined[['category']][combined[['category']]=='I']<-'NB_Day1'
combined[['category']][combined[['category']]=='L_U']<-'NB_U'
combined[['category']][combined[['category']]=='L_T']<-'NB_T'
combined[['category']][combined[['category']]=='M059K_U']<-'WC_U'
combined[['category']][combined[['category']]=='M059K_T']<-'WC_T'


# plot quality metrics
pdf("seurat/seurat1_pre.pdf", width = 6, height = 5)

plots_list <- list()
y_labels <- c('Total reads', '% Reads mapped', 'mRNA bases (log)', '% Mitochondrial counts', '% mRNA bases','Expressed genes')
features <- c('total', 'proportion_mapped', 'log_num_mrna','percent.mt','PCT_MRNA_BASES','nFeature_RNA')

for (i in 1:length(features)) {
  
  # Create a data frame with the values and groups
  df <- data.frame(Value = combined@meta.data[, features[i]], Sample = combined@meta.data[, 'category'])
  
  # Create the ggplot object for the current feature and add it to the list
  p <- ggplot(df, aes(x = Sample, y = Value, fill = Sample)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", color = "black") +
    labs(x = "", y = y_labels[i]) +
    scale_fill_manual(values = c('WC_T' = "#8C804B", 'WC_U' = "#E0C27F", 'NB_T' = "#00876F", 'NB_U' = "#86CABD", 'NB_Initial' = "#d9d9d9")) +  # Customize the fill colors
    theme_classic() +
    theme(legend.position = "none")  # Remove default legend

  plots_list[[length(plots_list) + 1]] <- ggplotGrob(p)
}
# Arrange all the plots in a grid
grid.arrange(grobs = plots_list, ncol = 2)

dev.off()


#process
nb_data<-NormalizeData(nb_data)
nb_data<-ScaleData(nb_data,scale.max = 5,features=rownames(nb_data),vars.to.regress = c("nFeature_RNA","PCT_MRNA_BASES")) #"batch",
wc_data<-NormalizeData(wc_data)
wc_data<-ScaleData(wc_data,scale.max = 5,features=rownames(wc_data),vars.to.regress = c("nFeature_RNA"))

outdata<-as.data.frame(nb_data[['RNA']]@scale.data)
write.table(data.frame(gene = row.names(outdata), outdata), file='seurat/seurat_nb_scaled.txt', quote=FALSE, sep='\t',row.names=FALSE, col.names =TRUE) 

#WC plot
wc_data[['category']][wc_data[['category']]=='M059K_U']<-'WC_U'
wc_data[['category']][wc_data[['category']]=='M059K_T']<-'WC_T'
wc_data<-FindVariableFeatures(wc_data, selection.method = "disp",num.bin=30, nfeatures = 600)
# wc_data<-FindVariableFeatures(wc_data, selection.method = "mvp", dispersion.cutoff = c(0.5, Inf),num.bin=30)
wc_data <- RunPCA(wc_data, verbose = FALSE,n_pcs=30)
wc_data <- RunUMAP(wc_data, reduction = "pca", dims = 1:8, n.neighbors =20)
pdf("seurat/seurat2.pdf", width=4, height=5)
VariableFeaturePlot(wc_data)
ElbowPlot(wc_data, n=30, reduction='pca')
DimPlot(wc_data,reduction = "umap",group.by='category',pt.size=4, cols=c('WC_T' = "#8C804B", 'WC_U' = "#E0C27F"))+ ggtitle("")
dev.off()


# pdf("seurat/seurat3.pdf", width=10, height=10)
# wcp<-DimPlot(wc_data,cols=c("red","royalblue2"),reduction = "umap",group.by='category',pt.size=8)
# p1 <- DimPlot(object = wc_data, reduction="pca",dims=c(1,2), group.by='category') + theme(legend.pos="none")
# p2 <- DimPlot(object = wc_data, reduction="pca",dims=c(3,4), group.by='category') + theme(legend.pos="none")
# p3 <- DimPlot(object = wc_data, reduction="pca",dims=c(6,5), group.by='category') + theme(legend.pos="none")
# p4 <- DimPlot(object = wc_data, reduction="pca",dims=c(7,8), group.by='category') + theme(legend.pos="none")
# p5 <- DimPlot(object = wc_data, reduction="pca",dims=c(9,10), group.by='category') + theme(legend.pos="none")
# p6 <- DimPlot(object = wc_data, reduction="pca",dims=c(11,12), group.by='category') + theme(legend.pos="bottom")
# grid.arrange(p1, p2,p3,p4,p5,p6, ncol=3)
# p1=DimPlot(wc_data,reduction = "umap",group.by='category')
# p2=FeaturePlot(wc_data,reduction = "umap",features='nFeature_RNA')
# p3=FeaturePlot(wc_data,reduction = "umap",features='percent.mt')
# p4=FeaturePlot(wc_data,reduction = "umap",features='nCount_RNA')
# p5=FeaturePlot(wc_data,reduction = "umap",features='total')
# p6=FeaturePlot(wc_data,reduction = "umap",features='unique.total')
# p7=FeaturePlot(wc_data,reduction = "umap",features='PCT_RIBOSOMAL_BASES')
# p8=FeaturePlot(wc_data,reduction = "umap",features='PCT_CODING_BASES')
# p9=FeaturePlot(wc_data,reduction = "umap",features='PCT_MRNA_BASES')
# p10=FeaturePlot(wc_data,reduction = "umap",features='PCT_INTRONIC_BASES')
# p12=FeaturePlot(wc_data,reduction = "umap",features='MEDIAN_CV_COVERAGE')
# p16=FeaturePlot(wc_data,reduction = "umap",split.by='orig.ident',features='percent.mt')
# grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p12,p16, ncol=4)


#NB plot
nb_dataL <- subset(nb_data, idents = c("L_U","L_T"))
nb_dataI<-subset(nb_data, idents = c("I"))
nb_dataR<-nb_data
nb_dataL<-FindVariableFeatures(nb_dataL, selection.method = "disp", nfeatures = 600)
VariableFeaturePlot(nb_dataL)
gen<-VariableFeatures(nb_dataL)
VariableFeatures(nb_dataR)<-gen
nb_dataR <- RunPCA(nb_dataR, verbose = FALSE,npcs = 20)
ElbowPlot(nb_dataR, n=20, reduction='pca')
nb_dataR <- RunUMAP(nb_dataR, reduction = "pca", dims = 1:8, n.neighbors = 50)
DimPlot(nb_dataR,cols=c("orange","royalblue2", "red"),reduction = "umap",group.by='category',cells=colnames(nb_dataR)[nb_dataR@meta.data$category!="I"])

nb_dataR[['category2']]<-nb_dataR[['category']]
nb_dataR[['category2']][nb_dataR[['category2']]=='I']<-'NB_Day1'
nb_dataR[['category2']][nb_dataR[['category2']]=='L_U']<-'NB_U'
nb_dataR[['category2']][nb_dataR[['category2']]=='L_T']<-'NB_T'
dev.off()
pdf("seurat/seurat4.pdf", width=4, height=5)
DimPlot(nb_dataR,cols=c('NB_T' = "#00876F", 'NB_U' = "#86CABD", 'NB_Day1' = "#d9d9d9"),reduction = "umap",group.by='category2',pt.size=4,cells=colnames(nb_dataR)[nb_dataR@meta.data$category!="I"])+ ggtitle("")
dev.off()


pdf("seurat/seurat5.pdf", width=12, height=12)
p1=DimPlot(nb_dataR,reduction = "umap",group.by='category2',cols=c('NB_T' = "#00876F", 'NB_U' = "#86CABD", 'NB_Day1' = "#d9d9d9"))+ ggtitle("Nanobiopsy")
nb_dataR[['batch']][[1]] <- factor(nb_dataR[['batch']][[1]], levels = paste0(0:10))
p2=DimPlot(nb_dataR,reduction = "umap",group.by='batch')+ ggtitle("Batch")
p3=FeaturePlot(nb_dataR,reduction = "umap",features='nFeature_RNA')+ ggtitle("Expressed genes")
p4=FeaturePlot(nb_dataR,reduction = "umap",features='nCount_RNA')+ ggtitle("Total reads")
p5=FeaturePlot(nb_dataR,reduction = "umap",features='proportion_mapped')+ ggtitle("Proportion mapped")
p6=FeaturePlot(nb_dataR,reduction = "umap",features='unique.total')+ ggtitle("Proportion uniquely mapped")
p7=FeaturePlot(nb_dataR,reduction = "umap",features='percent.mt')+ ggtitle("% Mitochondrial counts")
p8=FeaturePlot(nb_dataR,reduction = "umap",features='num_mrna')+ ggtitle("mRNA bases")
p9=FeaturePlot(nb_dataR,reduction = "umap",features='PCT_MRNA_BASES')+ ggtitle("% mRNA bases")
p10=FeaturePlot(nb_dataR,reduction = "umap",features='PCT_CODING_BASES')+ ggtitle("% Coding bases")
p11=FeaturePlot(nb_dataR,reduction = "umap",features='PCT_INTRONIC_BASES')+ ggtitle("% Intronic bases")
p12=FeaturePlot(nb_dataR,reduction = "umap",features='MEDIAN_CV_COVERAGE')+ ggtitle("Exon coverage variance")
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, ncol=4)
dev.off()
pdf("seurat/seurat6.pdf", width=5, height=5)
p1 <- DimPlot(object = nb_dataR, cols=c('NB_T' = "#00876F", 'NB_U' = "#86CABD", 'NB_Day1' = "#d9d9d9"),reduction="pca",dims=c(1,2), group.by='category2',cells=colnames(nb_dataR)[nb_dataR@meta.data$category!="I"]) + ggtitle("")
p2 <- DimPlot(object = nb_dataR, cols=c('NB_T' = "#00876F", 'NB_U' = "#86CABD", 'NB_Day1' = "#d9d9d9"),reduction="pca",dims=c(3,4), group.by='category2',cells=colnames(nb_dataR)[nb_dataR@meta.data$category!="I"]) + ggtitle("")
p3 <- DimPlot(object = nb_dataR, cols=c('NB_T' = "#00876F", 'NB_U' = "#86CABD", 'NB_Day1' = "#d9d9d9"),reduction="pca",dims=c(6,5), group.by='category2',cells=colnames(nb_dataR)[nb_dataR@meta.data$category!="I"]) + ggtitle("")
p4 <- DimPlot(object = nb_dataR, cols=c('NB_T' = "#00876F", 'NB_U' = "#86CABD", 'NB_Day1' = "#d9d9d9"),reduction="pca",dims=c(7,8), group.by='category2',cells=colnames(nb_dataR)[nb_dataR@meta.data$category!="I"]) + ggtitle("")
grid.arrange(p1, p2,p3,p4,ncol=2)
dev.off()




pairs<-read.table('meta/biopsy_pairs.txt',header=TRUE)
vals<- data.frame(matrix(ncol = 5, nrow = 0))
colnames(vals)<-c('i_1','i_2','r_1','r_2','treatment')
lins<- data.frame(matrix(ncol = 4, nrow = 0))
colnames(lins)<-c('xx','yy','long','direction')
initials_paired<-c()
for(i in 1:nrow(pairs)) {
    initial <- pairs[i,2]
    recurrent <- pairs[i,1]
        if (initial %in% colnames(nb_dataR) & recurrent %in% colnames(nb_dataR)){
        initials_paired<-append(initials_paired,initial)  
        i_x<-nb_dataR@ reductions$umap@ cell.embeddings[,1][initial][[1]]
        i_y<-nb_dataR@ reductions$umap@ cell.embeddings[,2][initial][[1]]
        r_x<-nb_dataR@ reductions$umap@ cell.embeddings[,1][recurrent][[1]]
        r_y<-nb_dataR@ reductions$umap@ cell.embeddings[,2][recurrent][[1]]
        tre<-ifelse(pairs[i,3]=="L_U",'untreated','treated')
        lins[nrow(lins) + 1,]<-list(i_x, i_y, recurrent,tre)
        lins[nrow(lins) + 1,]<-list(r_x, r_y, recurrent,tre)
        vals[nrow(vals) + 1,]<-list(i_x, i_y, r_x,r_y,tre)
}
}

pdf("seurat/seurat7.pdf", width=12, height=7)
x_range <- range(nb_dataR@reductions$umap@cell.embeddings[, 1])
y_range <- range(nb_dataR@reductions$umap@cell.embeddings[, 2])
plots_list <- list()
# Loop over untreated samples and create the plot for each sample
data<-vals[vals$treatment=="untreated",]
for (i in 1:nrow(data)) {
  row <- data[i, ]
  p1 <- DimPlot(nb_dataR, dims=c( 1,2 ), cols=c('NB_T' = "#00876F", 'NB_U' = "#86CABD", 'NB_Day1' = "#d9d9d9"), reduction = "umap",
                group.by='category2', pt.size=5, cells=append(colnames(nb_dataR)[nb_dataR@meta.data$category2!="NB_Day1"],initials_paired)) +
    scale_x_continuous(limits = c(x_range[1], x_range[2])) +
    scale_y_continuous(limits = c(y_range[1], y_range[2]))+
     ggtitle("")
  p1 <- p1 + geom_segment(aes(x = row$i_1, y = row$i_2, xend = row$r_1, yend = row$r_2),
                           color = "black", size = 0.25, arrow = arrow(length = unit(0.3, "cm"), type = "closed"))+ 
                           theme(legend.position = "none")
  p1<-p1 + theme(axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
  plots_list[[length(plots_list)+1]] <- ggplot_gtable(ggplot_build(p1))

}

# Arrange plots and the legend
comb<-grid.arrange(grobs = plots_list, ncol = 5,top = textGrob("Untreated",  gp = gpar(fontsize = 25, font = 1)),bottom=textGrob("UMAP1", gp=gpar(fontsize=18)),left = textGrob("UMAP2", rot = 90,gp = gpar(fontsize = 18)))

plots_list <- list()
# Loop over treated samples and create the plot for each sample
data<-vals[vals$treatment=="treated",]
for (i in 1:nrow(data)) {
  row <- data[i, ]
  p1 <- DimPlot(nb_dataR, dims=c( 1,2 ), cols=c('NB_T' = "#00876F", 'NB_U' = "#86CABD", 'NB_Day1' = "#d9d9d9"), reduction = "umap",
                group.by='category2', pt.size=5, cells=append(colnames(nb_dataR)[nb_dataR@meta.data$category2!="NB_Day1"],initials_paired)) +
    scale_x_continuous(limits = c(x_range[1], x_range[2])) +
    scale_y_continuous(limits = c(y_range[1], y_range[2]))+ 
    ggtitle("")
  p1 <- p1 + geom_segment(aes(x = row$i_1, y = row$i_2, xend = row$r_1, yend = row$r_2),
                           color = "black", size = 0.25, arrow = arrow(length = unit(0.3, "cm"), type = "closed"))+
                            theme(legend.position = "none")
p1<-p1 + theme(axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
  plots_list[[length(plots_list)+1]]  <- ggplot_gtable(ggplot_build(p1))
  print(i)
}
p1 <- get_legend(DimPlot(nb_dataR, dims = c(1, 2), cols=c('NB_T' = "#00876F", 'NB_U' = "#86CABD", 'NB_Day1' = "#d9d9d9"), reduction = "umap",
                         group.by = 'category2', pt.size = 8,
                         cells = append(colnames(nb_dataR)[nb_dataR@meta.data$category != "I"], initials_paired)) +
                  theme(legend.text = element_text(size = 15),      # Adjust the legend text size
                        legend.title = element_text(size = 15),
                        legend.key.size = unit(1, "cm")) +
                  guides(color = guide_legend(override.aes = list(size = 6, label.position = "right"))))

plots_list[[length(plots_list)+1]] <- p1
# Arrange all the plots in a grid
comb<-grid.arrange(grobs = plots_list, ncol = 5,top = textGrob("Treated",  gp = gpar(fontsize = 25, font = 1)),bottom=textGrob("UMAP1", gp=gpar(fontsize=18)),left = textGrob("UMAP2", rot = 90,gp = gpar(fontsize = 18)))

dev.off()














