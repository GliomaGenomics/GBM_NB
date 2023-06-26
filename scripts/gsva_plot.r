library(tidyr)

wang<-read.table('wang/gsva_paired_results.txt',header=T)
wang_tr<-wang[wang$category=='L_T',]
wang_un<-wang[wang$category=='L_U',]
wang_tr<-pivot_wider(wang_tr, names_from = c('time','end'), values_from = 'GSVA_score')
wang_un<-pivot_wider(wang_un, names_from = c('time','end'), values_from = 'GSVA_score')


start_pro_tr <- wang_tr$initial_Pro
end_pro_tr <- wang_tr$longitudinal_Pro
start_mes_tr <- wang_tr$initial_Mes
end_mes_tr <- wang_tr$longitudinal_Mes
samples_tr <- wang_tr$sample
start_pro_un <- wang_un$initial_Pro
end_pro_un <- wang_un$longitudinal_Pro
start_mes_un <- wang_un$initial_Mes
end_mes_un <- wang_un$longitudinal_Mes
samples_un <- wang_un$sample

pdf("seurat/gsva_wang.pdf", width=8, height=8)
# Plotting the lines
par(mfrow = c(2, 1))
# par( mai = c(0.3, 1, 1, 1))
par(mar=c(7,5,2,5))
plot(x = 1:length(samples_un), y = start_pro_un, xlab = "", ylab = "", 
     type = "n", xlim = c(0.5, length(samples_un) + 0.5), ylim = range(c(start_pro_un, end_pro_un,start_mes_un, end_mes_un)), xaxt = "n")
segments(x0 = 1:length(samples_un)-0.3, y0 = start_pro_un, x1 = 1:length(samples_un)+0.3, y1 = end_pro_un, 
         col = "darkorchid4", lwd = 2)
segments(x0 = 1:length(samples_un)-0.3, y0 = start_mes_un, x1 = 1:length(samples_un)+0.3, y1 = end_mes_un, 
         col = "green4", lwd = 2)
         
points(x = 1:length(samples_un)-0.3, y = start_pro_un, pch = 20, col = "orchid2", cex = 2)
points(x = 1:length(samples_un)+0.3, y = end_pro_un, pch = 20, col = "darkorchid", cex = 2)
points(x = 1:length(samples_un)-0.3, y = start_mes_un, pch = 20, col = "palegreen2", cex = 2)
points(x = 1:length(samples_un)+0.3, y = end_mes_un, pch = 20, col = "green4", cex = 2)

title(xlab = "Cell: Day1 -> longitudinal", mgp = c(1, 1, 0))
title(ylab = "GSVA score", mgp = c(3, 1, 0))
axis(side = 1, at = c(1:length(samples_un)), labels=rep("", length(samples_un)), tick=TRUE)

# axis(side = 1, at = sort(c(1:length(samples_un)-0.3,1:length(samples_un)+0.3)), labels = rep(c('initial', 'longitudinal'), times = length(samples_tr)), tick = TRUE, padj = 0.3,las=2)
title('Unreated')
# par( mai = c(0.3, 1, 1, 1))
par(mar=c(7,5,2,5))

plot(x = 1:length(samples_tr), y = start_pro_tr, xlab = "", ylab = "", 
     type = "n", xlim = c(0.5, length(samples_tr) + 0.5), ylim = range(c(start_pro_tr, end_pro_tr,start_mes_tr, end_mes_tr)), xaxt = "n",las=2)
segments(x0 = 1:length(samples_tr)-0.3, y0 = start_pro_tr, x1 = 1:length(samples_tr)+0.3, y1 = end_pro_tr, 
         col = "mediumorchid", lwd = 2)
segments(x0 = 1:length(samples_tr)-0.3, y0 = start_mes_tr, x1 = 1:length(samples_tr)+0.3, y1 = end_mes_tr, 
         col = "green3", lwd = 2)
points(x = 1:length(samples_tr)-0.3, y = start_pro_tr, pch = 20, col = "orchid2", cex = 2)
points(x = 1:length(samples_tr)+0.3, y = end_pro_tr, pch = 20, col = "darkorchid", cex = 2)
points(x = 1:length(samples_tr)-0.3, y = start_mes_tr, pch = 20, col = "palegreen2", cex = 2)
points(x = 1:length(samples_tr)+0.3, y = end_mes_tr, pch = 20, col = "green4", cex = 2)
axis(side = 1, at = c(1:length(samples_tr)), labels=rep("", length(samples_tr)), tick=TRUE)

legend( x=0.1, y=-0.7,xpd=TRUE, bty="n",ncol=2,lty="blank",legend = c("Proneural score - Day1", "Proneural score - Longitudinal", "Mesenchymal score - Day1", "Mesenchymal score - Longitudinal"), 
       pch=c(20,20,20,20),col = c("orchid2", "darkorchid","palegreen2", "green4"), lwd = 2)

# axis(side = 1, at = sort(c(1:length(samples_un)-0.3,1:length(samples_un)+0.3)),labels = rep(c('initial', 'longitudinal'), times = length(samples_un)), tick = TRUE, padj = 0.3,las=2)
title(xlab = "Cell: Day1 -> longitudinal", mgp = c(1, 1, 0))
title(ylab = "GSVA score", mgp = c(3, 1, 0))
title('Treated')

dev.off()
