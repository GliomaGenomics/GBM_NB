norm<-read.table('seurat/seurat_nb_norm.txt',header=TRUE)
pairs<-read.table('meta/biopsy_pairs.txt',header=TRUE)
out<-data.frame(row.names = rownames(norm))
for(i in 1:nrow(pairs)) {
    initial <- pairs[i,2]
    recurrent <- pairs[i,1]
    treatment <- pairs[i,3]
        if (initial %in% colnames(norm) & recurrent %in% colnames(norm)){
            out[paste0(initial,'_',recurrent,'_',treatment)]<- 0
            out[paste0(initial,'_',recurrent,'_',treatment)][norm[recurrent]!=0 & norm[initial]!=0,]<-log((norm[recurrent]/norm[initial]), base=10)[norm[recurrent]!=0 & norm[initial]!=0,]
        }
}    
write.table(out,file='seurat/longitudinal_FC.txt')