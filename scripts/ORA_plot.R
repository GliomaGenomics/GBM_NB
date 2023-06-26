library(ggplot2)

datain<-read.table('deseq2/nb_UvT_deseq2_over_rep_bp_webgestalt_affinityprop.txt',header=TRUE,sep='\t')
ord<-datain$Description[order(datain$Ratio)]
ord<-as.character(ord)
datain$Description <- factor(datain$Description,levels=ord)
pdf(paste("deseq2/nb_UvT_deseq2_over_rep_bp_webgestalt_affinityprop.pdf", sep=""), width=7, height=3)

ggplot(datain, aes(Ratio, Description, colour=FDR, size=Size))+ 
           geom_point()+ 
           scale_color_gradient2(name='FDR', low='orange', mid='purple',  high='blue',midpoint=0.015)+
           scale_size_area(name='Set size')+
           scale_size(range=c(2,7))+
           theme(legend.key.size = unit(10, 'points'),legend.title = element_text(size=10), legend.text = element_text(size=8),axis.title.y=element_blank())+
           xlab('Enrichment Ratio')+
           scale_x_continuous(breaks=c(1,2,3,4), limits=c(1,4)) 

dev.off()






datain<-read.table('deseq2/nb_T_SvD_deseq2_sig0.05_unique_to_T_over_rep_cus.txt',header=TRUE,sep=' ')
datain$proportion_overlap<-datain$overlap/datain$size
ord<-datain$pathway[order(datain$padj,decreasing = TRUE)]
ord<-as.character(ord)
datain$pathway <- factor(datain$pathway,levels=ord)
datain<-datain[1:30,]

pdf(paste("deseq2/nb_T_SvD_deseq2_sig0.05_unique_to_T_over_rep_cus.pdf", sep=""), width=7, height=5)

ggplot(datain, aes(proportion_overlap, pathway, colour=padj, size=size))+ 
           geom_point()+ 
           scale_color_gradient2(name='FDR', low='orange', mid='purple',  high='blue',midpoint=0.0000002)+
           scale_size_area(name='Set size')+
           scale_size(range=c(2,7))+
           theme(legend.key.size = unit(10, 'points'),legend.title = element_text(size=10), legend.text = element_text(size=8),axis.title.y=element_blank())+
           xlab('Proportion overlap')+
           scale_x_continuous(breaks=c(0,0.25,0.5), limits=c(0,0.5)) 

dev.off()






datain$proportion_overlap<-datain$overlap/datain$size
ord<-datain$pathway[order(datain$padj)]
ord<-as.character(ord)
datain$pathway <- factor(datain$pathway,levels=ord)
datain<-datain[1:10,]
pdf(paste("deseq2/nb_UvT_deseq2_over_rep_bp.pdf", sep=""), width=7, height=5)

ggplot(datain, aes(proportion_overlap, pathway, colour=padj, size=size))+ 
           geom_point()+ 
           scale_color_gradient2(name='FDR', low='red', mid='yellow',  high='blue',midpoint=0.02)+
           scale_size_area(name='Set size')+
           scale_size(range=c(2,7))+
           theme(legend.key.size = unit(10, 'points'),legend.title = element_text(size=10), legend.text = element_text(size=8),axis.title.y=element_blank())+
           xlab('Proportion overlap')+
           scale_x_continuous(breaks=c(1,3,5,7,9,11), limits=c(1,11)) 

dev.off()
