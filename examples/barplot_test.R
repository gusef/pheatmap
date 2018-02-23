setwd("C:/Users/GUSENDA1/OneDrive - Novartis Pharma AG/GitHub/pheatmap/examples")
rm(list=ls())

require(cba)
require(GSVA)
require(pheatmap)
require(RColorBrewer)
require(SummarizedExperiment)
require(xlsx)

se <- readRDS('X:/qms/research/analysis/NGDX-P00221/GUSENDA1_rationale_for_combination/derived_data/P218_P220.RDS')

subtypes <- read.xlsx('X:/qms/research/analysis/NGDX-P00221/GUSENDA1_rationale_for_combination/raw_data/Burstein_2014_TNBC_subtypes.xlsx',1)
subtypes$Genes <- sapply(strsplit(as.character(subtypes$Genes),' '),function(x)x[1])
gene_sets <- lapply(levels(subtypes$Subtype),function(x,y)y$Genes[y$Subtype==x],subtypes)
names(gene_sets) <- levels(subtypes$Subtype)

#Using screening samples only
screening <- se[,se$Visit.Name == 'SCREENING']
mat <- assays(screening)$final[unlist(gene_sets),]

hcopt <- function(d, HC=NULL, method = "ward.D", members = NULL)
{
   if ( is.null(HC) ) {
      HC <- hclust(d,method=method,members=members)
   }
   if (ncol(d) > 2){
      ORD <- order.optimal(d,merge=HC$merge)
      HC$merge <- ORD$merge
      HC$order <- ORD$order
   }
   HC
}

#column ordering
hc01.col <- hcopt(dist(t(mat)),method="ward.D")

#reorder each of the sets internally
subsets <- lapply(gene_sets,function(gs,mat)mat[gs,],mat)
subsets <- lapply(subsets,
                  function(x)x[order.dendrogram(as.dendrogram(hcopt(as.dist(1-cor(t(x))),
                                                                    method = "ward.D"))),])
mat <- do.call(rbind,subsets)

#gaps
gaps <- cumsum(sapply(subsets,nrow))
gaps <- gaps[-length(gaps)]

#column annotation
col_annot <- data.frame(colData(screening)[,c('BISRC','treatment','BOR')])
col_annot$BOR <- as.factor(col_annot$BOR)
BOR  <- c('#e34a33','#99d8c9','#ffeda0','#636363')
names(BOR) <- levels(col_annot$BOR)
treatment <- c('#c2a5cf','#40004b')
names(treatment) <- levels(as.factor(col_annot$treatment))
annot_colors <- list(treatment=treatment,
                     BOR=BOR)


scores <- gsva(assays(screening)$final, gene_sets)[1,]


par(mfrow=c(2,1))
barplot(scores)
pheatmap(mat,
         main = 'Burstein subtypes',
         scale = 'row',
         color = brewer.pal(11, "RdBu")[11:1],
         cluster_rows = FALSE,
         annotation_col = col_annot,
         annotation_colors = annot_colors,
         cluster_cols = hc01.col,
         gaps_row = gaps,
         gaps_col = c(8,16,24),
         fontsize=8)



