library(GSVA)

nb_data<-as.matrix(read.table('seurat/seurat_nb_scaled.txt',header=TRUE,row.names=1))

#some of the below files have been removed as weren't used
c1_top_100<-as.list(read.table('wang/C1_PC1_top100.txt',header=FALSE))$V1
c1_top_50<-as.list(read.table('wang/C1_PC1_top50.txt',header=FALSE))$V1
sc_top_100<-as.list(read.table('wang/SC_PC1_top100.txt',header=FALSE))$V1
sc_top_50<-as.list(read.table('wang/SC_PC1_top50.txt',header=FALSE))$V1
c1_bottom_100<-as.list(read.table('wang/C1_PC1_bottom100.txt',header=FALSE))$V1
c1_bottom_50<-as.list(read.table('wang/C1_PC1_bottom50.txt',header=FALSE))$V1
sc_bottom_100<-as.list(read.table('wang/SC_PC1_bottom100.txt',header=FALSE))$V1
sc_bottom_50<-as.list(read.table('wang/SC_PC1_bottom50.txt',header=FALSE))$V1
gs<-list('c1_top_100'=c1_top_100,'c1_top_50'=c1_top_50,'sc_top_100'=sc_top_100,'sc_top_50'=sc_top_50,'c1_bottom_100'=c1_bottom_100,'c1_bottom_50'=c1_bottom_50,'sc_bottom_100'=sc_bottom_100,'sc_bottom_50'=sc_bottom_50)

res<-gsva(nb_data,gs)

write.table(t(res), file='wang/gsva_results_scaled.txt')


nb_data2<-as.matrix(read.table('seurat/seurat_nb_scaled.txt',header=TRUE,row.names=1))

Patel_Bernstein_stemness<-as.list(read.table('wang/Patel-Bernstein_stemness.txt',header=FALSE))$V1
Glioma_stemness_PMID_29276782<-as.list(read.table('wang/Glioma_stemness_PMID_29276782.txt',header=FALSE))$V1
Weng_stemness_PMID_30982771<-as.list(read.table('wang/Weng_stemness_PMID_30982771.txt',header=FALSE))$V1
gs2<-list('Patel-Bernstein_stemness'=Patel_Bernstein_stemness,'Glioma_stemness_PMID_29276782'=Glioma_stemness_PMID_29276782,'Weng_stemness_PMID_30982771'=Weng_stemness_PMID_30982771)

res2<-gsva(nb_data2,gs2)

write.table(t(res2), file='wang/gsva_results_scaled_stemness2.txt',sep='\t')



library("ggplot2")
pdf("wang/gsva.pdf", width=5, height=5)
data<-read.table('wang/gsva_paired_results.txt',header=TRUE)

top<-ggplot(data[data$end=='Mes',], aes(x = time, y = GSVA_score)) +
  geom_boxplot(aes(fill=time),outlier.shape = NA,alpha=0.6) +
  geom_line(aes(group=sample,colour=time),show.legend = FALSE) +
  geom_point(aes(group=sample),size=2,shape=21) +
  facet_grid(~category,scales = "free") +
  scale_fill_brewer(palette="GnBu") +
  labs(y = "Mesenchymal score",x='')

bottom<-ggplot(data[data$end=='Pro',], aes(x = time, y = GSVA_score)) +
  geom_boxplot(aes(fill=time),outlier.shape = NA,alpha=0.6) +
  geom_line(aes(group=sample,colour=time),show.legend = FALSE) +
  geom_point(aes(group=sample),size=2,shape=21) +
  facet_grid(~category,scales = "free") +
  scale_fill_brewer(palette="GnBu") +
  labs(y = "Proneural score ",x='')

top
bottom
dev.off()


library("ggplot2")
pdf("wang/gsva_stemness.pdf", width=5, height=5)
data<-read.table('wang/gsva_paired_results_stemness.txt',header=TRUE)

ggplot(data, aes(x = time, y = GSVA_score)) +
  geom_boxplot(aes(fill=time),outlier.shape = NA,alpha=0.6) +
  geom_line(aes(group=sample,colour=time)) +
  geom_point(aes(group=sample),size=2,shape=21) +
  facet_grid(~category,scales = "free") +
  scale_fill_brewer(palette="GnBu") +
  labs(title = "GSVA score for Mes")
