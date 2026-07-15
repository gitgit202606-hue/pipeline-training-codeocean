# pipeline-training-codeocean

Co-expression Analysis: A co-expression analysis of the prioritized CRC targets (CDH17, GUCY2C, CDCP1, LY6G6D) across Caris datasets, PDX models, and cell lines is certainly feasible. I’ll need to confirm the sources for PDX model and cell line datasets with our colleagues. If you already have specific data in mind, please let me know!
 
CRC Model Selection & Cell Tinder Algorithm: I agree that prioritizing models with robust co-expression and surface antigen correlation is essential. Once the datasets are available, we can explore applying the cell tinder algorithm to assess cell line models for your functional assays.
 
GeoMX Data: Our CO team is launching GeoMX spatial data for multiple cancer indications, and I am currently collecting CART targets for both CRC and prostate cancer. Thank you for sending the prioritized CRC targets! Once the GeoMX data is available, I anticipate it will complement the target surface expression evaluation process I proposed for the rubric ranking scores and be informative for target and model selection.
 
## Loading in RNAseq

#Load
GSE123484.deseq2.stats <- readRDS("GSE123484.deseq2.stats.rds")

#Filtering for Sig Genes
#GSE123484.deseq2.stats <- GSE123484.deseq2.stats[GSE123484.deseq2.stats$pvalue <= 0.05, ]

## Loading in ATACseq

#Load
GSE123486.all.timepoints.unique.peaks <- readRDS("GSE123486.all.timepoints.unique.peaks.rds")

#Filtering for timepoint without IL6 stimulation
GSE123486.all.timepoints.unique.peaks <- GSE123486.all.timepoints.unique.peaks[GSE123486.all.timepoints.unique.peaks$timepoint =="0h", ]

#Exporting peak coordinates for Meme Suite
#write.csv(GSE123486.all.timepoints.unique.peaks,"GSE123486.0h.timepoints.unique.peaks.csv")

## checking overlap..

#creating character objects
expression.genes <- unique(rownames(GSE123484.deseq2.stats))
peaks.genes <- unique(GSE123486.all.timepoints.unique.peaks$SYMBOL)
expression.genes <- unique(na.omit(as.character(expression.genes)))
peaks.genes <- unique(na.omit(as.character(peaks.genes)))

expression.genes
peaks.genes
class(expression.genes)
class(peaks.genes)

#determining overlap
GSE123488.common.genes <- intersect(expression.genes, peaks.genes)
GSE123488.common.genes

library(VennDiagram)
library(grid)

# Create Venn diagram with a slightly different style
venn_object <- draw.pairwise.venn(
  area1 = expression.genes,
  area2 = peaks.genes,
  cross.area = GSE123488.common.genes,
  
  category = c (
    "Ptpn2+/+ vs. Ptpn2+/- \nTreg RNA-seq DEGs",
    "Ptpn2+/+ vs. Ptpn2+/- \nTreg ATAC-seq Genes"
  ),
  
  scaled = F,
  euler.d = F,
  
  fill = c("#c2d6e3", "#f2a7a7"),
  alpha = c(0.7, 0.7),
  lwd = 2,
  col = c("#2f4f61", "#c73c3c"),
  
  label.col = 'black',
  cex = 1.4,
  fontface = 'bold',
  
  cat.cex = 1.1,
  cat.fontface = "bold",
  cat.col = 'black',
  cat.pos = c(190, 10),
  cat.dist = c(0.08,0.08)
)

grid.newpage()
grid.draw(venn.plot)



Error in if (cross.area > area1 | cross.area > area2) { : 
  the condition has length > 1
In addition: Warning messages:
1: In cross.area > area1 :
  longer object length is not a multiple of shorter object length
2: In cross.area > area2 :
  longer object length is not a multiple of shorter object length
3: In cross.area > area1 | cross.area > area2 :
  longer object length is not a multiple of shorter object length


 
