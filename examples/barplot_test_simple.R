#rm(list=ls())

require(cba)
require(gtable)
require(GSVA)
require(pheatmap)
require(RColorBrewer)
require(SummarizedExperiment)




scores <- 10 * (runif(10) - 0.5)
mat <- matrix(runif(100),ncol=10)

BOR  <- c('#e34a33','#99d8c9','#ffeda0','#636363')
names(BOR) <- c('A','B','C','D')
source('../R/pheatmap.r')

pheatmap(mat,
         main = 'Burstein subtypes',
         scale = 'row',
         cluster_cols = scores,
         color = brewer.pal(11, "RdBu")[11:1],
         cluster_rows = FALSE,
         barplot_legend = BOR, 
         barplot_height = 50,
         fontsize=8)




