library(fgsea)

conv<-read.table('downloaded_data/gencode.v27_geneidtoname.txt')
con<-as.list(conv$V2)
names(con)<-substr(conv$V1,1,15)

gmt_bp<-gmtPathways('downloaded_data/c5.go.bp.v7.4.symbols.gmt')
# gmt<-gmtPathways('downloaded_data/c5.go.mf.v7.4.symbols.gmt')
gmt_cus<-gmtPathways('downloaded_data/custom_gene_sets.gmt')
for (set in names(gmt_cus)){
    gmt_cus[[set]]<-unlist(unname(con[gmt_cus[set][[1]]]))
    }

#whole-cell
    genes<-read.table('deseq2/wc_UvT_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    write.table(apply(out[,1:5],2,as.character),file='deseq2/wc_UvT_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    write.table(apply(out[,1:5],2,as.character),file='deseq2/wc_UvT_deseq2_over_rep_cus.txt',row.names=FALSE)

#UvT
    genes<-read.table('deseq2/nb_UvT_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_UvT_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_UvT_deseq2_over_rep_custom.txt',row.names=FALSE)

#UvT permutations
    for (i in 1:100){
    genes<-read.table(paste('deseq2/permutations/nb_UvT_deseq2_random_',i,'.txt',sep=''))
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out)) 
    write.table(apply(out[,1:5][out$padj<0.05],2,as.character),file=paste('deseq2/permutations/nb_UvT_deseq2_random_',i,'_over_rep_bp.txt',sep=''),row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out)) 
    write.table(apply(out[,1:5][out$padj<0.05],2,as.character),file=paste('deseq2/permutations/nb_UvT_deseq2_random_',i,'_over_rep_cus.txt',sep=''),row.names=FALSE)
    }

#initials UvT
    genes<-read.table(paste('deseq2/nb_IUvIT_deseq2','.txt',sep=''))
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pval),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IUvIT_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IUvIT_deseq2_over_rep_custom.txt',row.names=FALSE)

#SvD unique to T
    genes<-read.table('deseq2/nb_T_SvD_deseq2_sig0.05_unique_to_T.txt')[[1]]
    background<-read.table('deseq2/nb_TorU_SvD_background.txt')[[1]]
    out<-fora(pathways = gmt_bp, 
        genes = genes, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_T_SvD_deseq2_sig0.05_unique_to_T_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = genes, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_T_SvD_deseq2_sig0.05_unique_to_T_over_rep_cus.txt',row.names=FALSE)

#SvD unique to U
    genes<-read.table('deseq2/nb_U_SvD_deseq2_sig0.05_unique_to_U.txt')[[1]]
    background<-read.table('deseq2/nb_TorU_SvD_background.txt')[[1]]
    out<-fora(pathways = gmt_bp, 
        genes = genes, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_U_SvD_deseq2_sig0.05_unique_to_U_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = genes, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_U_SvD_deseq2_sig0.05_unique_to_U_over_rep_cus.txt',row.names=FALSE)

#SvD shared U and T
    genes<-read.table('deseq2/nb_T_SvD_deseq2_sig0.05_shared_with_U.txt')[[1]]
    background<-read.table('deseq2/nb_TorU_SvD_background.txt')[[1]]
    out<-fora(pathways = gmt_bp, 
        genes = genes, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_T_SvD_deseq2_sig0.05_shared_with_U_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = genes, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_T_SvD_deseq2_sig0.05_shared_with_U_over_rep_cus.txt',row.names=FALSE)


#SvD UvT
    genes<-read.table('deseq2/nb_UvT_SvD_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_UvT_SvD_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_UvT_SvD_deseq2_over_rep_custom.txt',row.names=FALSE)

#SvD T
    genes<-read.table('deseq2/nb_T_SvD_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_T_SvD_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_T_SvD_deseq2_over_rep_custom.txt',row.names=FALSE)

#SvD U
    genes<-read.table('deseq2/nb_U_SvD_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_U_SvD_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_U_SvD_deseq2_over_rep_custom.txt',row.names=FALSE)

#IvL UvT
    genes<-read.table('deseq2/nb_IvL_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvL_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvL_deseq2_over_rep_custom.txt',row.names=FALSE)

#IvL UvT div
    genes<-read.table('deseq2/nb_IvL_div_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvL_div_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvL_div_deseq2_over_rep_custom.txt',row.names=FALSE)

#IvL UvT nondiv
    genes<-read.table('deseq2/nb_IvL_nondiv_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvL_nondiv_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvL_nondiv_deseq2_over_rep_custom.txt',row.names=FALSE)

#UIvL
    genes<-read.table('deseq2/nb_UIvL_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_UIvL_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_UIvL_deseq2_over_rep_custom.txt',row.names=FALSE)

#UIvL div
    genes<-read.table('deseq2/nb_UIvL_div_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_UIvL_div_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_UIvL_div_deseq2_over_rep_custom.txt',row.names=FALSE)

#UIvL nondiv
    genes<-read.table('deseq2/nb_UIvL_nondiv_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_UIvL_nondiv_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_UIvL_nondiv_deseq2_over_rep_custom.txt',row.names=FALSE)

#TIvL
    genes<-read.table('deseq2/nb_TIvL_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_TIvL_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_TIvL_deseq2_over_rep_custom.txt',row.names=FALSE)

#TIvL div
    genes<-read.table('deseq2/nb_TIvL_div_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_TIvL_div_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_TIvL_div_deseq2_over_rep_custom.txt',row.names=FALSE)

#TIvL nondiv
    genes<-read.table('deseq2/nb_TIvL_nondiv_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_TIvL_nondiv_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_TIvL_nondiv_deseq2_over_rep_custom.txt',row.names=FALSE)

#IvU
    genes<-read.table('deseq2/nb_IvU_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvU_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvU_deseq2_over_rep_custom.txt',row.names=FALSE)

#IvU div
    genes<-read.table('deseq2/nb_IvU_div_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvU_div_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvU_div_deseq2_over_rep_custom.txt',row.names=FALSE)

#IvU nondiv
    genes<-read.table('deseq2/nb_IvU_nondiv_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvU_nondiv_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvU_nondiv_deseq2_over_rep_custom.txt',row.names=FALSE)

#IvT
    genes<-read.table('deseq2/nb_IvT_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvT_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvT_deseq2_over_rep_custom.txt',row.names=FALSE)

#IvT div
    genes<-read.table('deseq2/nb_IvT_div_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvT_div_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvT_div_deseq2_over_rep_custom.txt',row.names=FALSE)

#IvT nondiv
    genes<-read.table('deseq2/nb_IvT_nondiv_deseq2.txt')
    degs<-row.names(genes[(genes$padj<0.05 & !is.na(genes$padj)),])
    background<-row.names(genes[!is.na(genes$pvalue),])
    out<-fora(pathways = gmt_bp, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvT_nondiv_deseq2_over_rep_bp.txt',row.names=FALSE)
    out<-fora(pathways = gmt_cus, 
        genes = degs, 
        universe = background, 
        minSize = 10, 
        maxSize = 5000)
    out<-out[,1:5]  
    out$bonf<-out$pval*length(row.names(out))  
    write.table(apply(out[,1:6],2,as.character),file='deseq2/nb_IvT_nondiv_deseq2_over_rep_custom.txt',row.names=FALSE)



