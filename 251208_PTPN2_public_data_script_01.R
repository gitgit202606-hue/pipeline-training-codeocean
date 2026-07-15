#### Libraries ####
setwd("/projects/users/wilsosx11/251208_PTPN2_public_datasets/")
options(bitmapType='cairo')
options(scipen = 999)
options(max.print = 10000)
library(ggplot2)
library(ggpubr)
library(tibble)
library(data.table)
library(tidyr)
library(dplyr)
library(rstatix)
library(ggrepel)
library(tidyverse)
library(Biobase)
library(GSA)
library(reshape)
library(GSVA)
library(pheatmap)
library(gplots)
library(readr)
library(forcats)
library(gridExtra)
library(ComplexHeatmap)
library(vroom)
library(plyr)
library(textshape)
library(Hmisc)
library(ggcorrplot)
library(corrplot)
library(reshape2)
library(genefilter)
library(biomaRt)
library(edgeR)
library(limma)
library(Rtsne)
library(scales)
library(RColorBrewer)
library(topGO)
library(org.Mm.eg.db)
library(ggpmisc)
library(FSA)
library(quantiseqr)
library(dplyr)
library(purrr)
library(readr)
library(patchwork)
library(readr)
library(stringr)
library(GGally)
library(ggnet)
library(ggnetwork)
library(network)
library(DESeq2)
library(rtracklayer)
library(dplyr)
library(readr)

#### Loading in GSE123484: PTPN2 promotes pathogenic loss of FoxP3+ in RORgt+ Tregs [RNA-seq] ####

# Loading in meta data
meta <- read.csv("GSE123484_metadata.csv")
meta$sample <- sub('^[^_]*_(.*)', '\\1', meta$sample)
unique(meta$genotype)

# Loading in counts
GSE123484.counts <- read.csv("GSE123484_all_counts.tsv",sep="\t")

# Formatting for long df format
colnames(GSE123484.counts)
colnames(GSE123484.counts) <- gsub('^X.', '', colnames(GSE123484.counts)) # removing wierd X's from colnames
GSE123484.counts$Feature_name <- as.character(GSE123484.counts$Feature_name)
GSE123484.counts$Feature_name <- gsub("\\$", "", GSE123484.counts$Feature_name) # removing '$' from all values in column
GSE123484.long <- pivot_longer(  GSE123484.counts,  cols = `1_20190_YfpPos_IL17Neg`:`8_23935_YFPNeg_IL17APos`,  names_to = 'sample',  values_to = 'counts')
colnames(GSE123484.long)
GSE123484.long <- GSE123484.long[, c("sample","Feature_name","counts")]
GSE123484.long$sample <- sub('^[^_]*_(.*)', '\\1', GSE123484.long$sample)

# Merging with meta data
GSE123484.long <- merge(GSE123484.long,meta,by=c("sample"))

# Creating a Z-scaled column
GSE123484.long <- GSE123484.long %>%  group_by(strain,genotype, celltype) %>%  mutate(z.scaled = scale(counts))

# Adding a column annotating which mouse strain is in vivo vs. in vitro
unique(GSE123484.long$strain)
GSE123484.long <- GSE123484.long %>%
  mutate(model = case_when(
    strain == 'B6.H2d.SKG' ~ 'In Vivo',
    strain == 'FoxP3eGFP BALB/c' ~ 'In Vitro',
    TRUE ~ NA_character_
  ))

# Changing "Ptpn+/+" and "Ptpnfl/+" to "Ptpn2+/+","Ptpn2+/-"
GSE123484.long <- GSE123484.long %>% mutate(genotype = recode(genotype,'Ptpn+/+' = 'Ptpn2+/+','Ptpnfl/+' = 'Ptpn2+/-'))

#Organizing the data for easier graphing
unique(GSE123484.long$celltype )
GSE123484.long$genotype <- factor(GSE123484.long$genotype, levels = c("Ptpn2+/+","Ptpn2+/-"))
GSE123484.long$celltype <- factor(GSE123484.long$celltype, levels = c("Tregs","ExTregs"))
GSE123484.long$model <- factor(GSE123484.long$model, levels = c("In Vivo","In Vitro"))


# Saving as a rds file
saveRDS(GSE123484.long,"GSE123484_expression_long.rds")

#### DESEQ2 in GSE123484: PTPN2 promotes pathogenic loss of FoxP3+ in RORgt+ Tregs [RNA-seq] ####

# Loading in counts
GSE123484.counts <- read.csv("GSE123484_all_counts.tsv",sep="\t")

# Formatting for long df format
colnames(GSE123484.counts)
colnames(GSE123484.counts) <- gsub('^X.', '', colnames(GSE123484.counts)) # removing wierd X's from colnames
GSE123484.counts$Feature_name <- as.character(GSE123484.counts$Feature_name)
GSE123484.counts$Feature_name <- gsub("\\$", "", GSE123484.counts$Feature_name) # removing '$' from all values in column


# fixing colnames for counts object
colnames(GSE123484.counts) <- sub('^[^_]*_(.*)', '\\1', colnames(GSE123484.counts))

# Loading in meta data
meta <- read.csv("GSE123484_metadata.csv")
meta$sample <- sub('^[^_]*_(.*)', '\\1', meta$sample)
unique(meta$genotype)

#filtering meta for Tregs only
meta <- meta[meta$celltype == 'Tregs', ]
meta <- meta[meta$strain == 'B6.H2d.SKG', ]
meta$strain <- as.factor(meta$strain)
meta$genotype <- as.factor(meta$genotype)

#filtering counts for Treg only samples
meta$sample
GSE123484.counts <- GSE123484.counts[, colnames(GSE123484.counts) %in% c("name","20190_YfpPos_IL17Neg",   "20191_YfpPos_IL17Neg",   "20192_YfpPos_IL17Neg",   "20193_YfpPos_IL17Neg",   "20085_FoxP3Pos_IL17Neg","20086_FoxP3Pos_IL17Neg", "23743_YFPPos_IL17ANeg",  "23745_YFPPos_IL17ANeg",  "23934_YFPPos_IL17ANeg",  "23935_YFPPos_IL17ANeg" )]
rownames(GSE123484.counts) <- GSE123484.counts$name
GSE123484.counts$name <- NULL


# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(
  countData = GSE123484.counts,
  colData = meta,
  design = ~genotype
)

# run analysis
dds <- DESeq(dds)

# extract normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# get results
res <- results(dds)

# convert res to dataframe
GSE123484.deseq2.stats <- as.data.frame(res)
GSE123484.deseq2.stats <- na.omit(GSE123484.deseq2.stats)

# Saving as rds objects 
saveRDS(GSE123484.deseq2.stats, file = "GSE123484.deseq2.stats.rds")

#Loading rds objects
GSE123484.deseq2.stats <- readRDS("GSE123484.deseq2.stats.rds")



#### Plotting Boxplot Genes from GSE123484: PTPN2 promotes pathogenic loss of FoxP3+ in RORgt+ Tregs [RNA-seq] ####

# Load dataframe
GSE123484.long <- readRDS("GSE123484_expression_long.rds")

#filtering out extregs and in vitro
GSE123484.long <- GSE123484.long[GSE123484.long$celltype == 'Tregs', ]
GSE123484.long <- GSE123484.long[GSE123484.long$strain == 'B6.H2d.SKG', ]

#Filtering for Genes of Interest
colnames(GSE123484.long)
GSE123484.filtered <- GSE123484.long[GSE123484.long$Feature_name %in% c("Il6ra","Il2","Cd3d","Cd3e","Cd3g","Lag3","Tox","Tigit","Foxp3","Jak2","Stat5b","Rorc","Ikzf3","Itgae","Ctla4","Il6","Il21","Tbx21","Ikzf4","Ccl20","Ccr6","Gpr83","Tgfb1","Dnmt3a","Prdm1","Stat3","Stat5a","Jak1","Klrc1","Tgfb3","Il1b","Il2ra"), ]


# Plotting Z-scale box plots of Genes 
p <- ggboxplot(GSE123484.filtered, "genotype", "z.scaled", color = "genotype",outliers=F) + theme_classic()
p <- p + facet_wrap(~Feature_name,scales="free_y",ncol=16)
p <- p + labs(x="",y="Z-scale",color="Genotype")
p <- p + theme(axis.text.x = element_text(angle=90,face="bold",size=16,color="black")) 
p <- p + theme(axis.text.y = element_text(face="bold",size=16,color="black")) 
p <- p + theme(strip.text = element_text(face = "bold",size=16,color="black")) 
p <- p + theme(axis.title = element_text(face = "bold",size = 16,color="black"))
p <- p + theme(legend.text = element_text(face = "bold",size=16,color="black"))
p <- p + theme(legend.title = element_text(face = "bold",size=16,color="black"))
p <- p + theme(title = element_text(face = "bold",size=16,color="black"))
p

# Plotting Counts box plots of Genes 
p <- ggboxplot(GSE123484.filtered, "genotype", "counts", color = "genotype",outliers=F) + theme_classic()
p <- p + facet_wrap(~Feature_name,scales="free_y",ncol=7)
p <- p + labs(x="",y="Counts",color="Genotype")
p <- p + theme(axis.text.x = element_text(angle=90,face="bold",size=12,color="black")) 
p <- p + theme(axis.text.y = element_text(face="bold",size=10,color="black")) 
p <- p + theme(strip.text = element_text(face = "bold",size=12,color="black")) 
p <- p + theme(axis.title = element_text(face = "bold",size = 12,color="black"))
p <- p + theme(legend.text = element_text(face = "bold",size=12,color="black"))
p <- p + theme(legend.title = element_text(face = "bold",size=12,color="black"))
p <- p + theme(title = element_text(face = "bold",size=12,color="black"))
p

#### Plotting Heatmap Genes from GSE123484: PTPN2 promotes pathogenic loss of FoxP3+ in RORgt+ Tregs [RNA-seq] ####

# Load dataframe
GSE123484.long <- readRDS("GSE123484_expression_long.rds")

#filtering out extregs and in vitro
GSE123484.long <- GSE123484.long[GSE123484.long$celltype == 'Tregs', ]
GSE123484.long <- GSE123484.long[GSE123484.long$strain == 'B6.H2d.SKG', ]

#Filtering for Genes of Interest
colnames(GSE123484.long)
GSE123484.filtered <- GSE123484.long[GSE123484.long$Feature_name %in% c("Il6ra","Il2","Cd3d","Cd3e","Cd3g","Lag3","Tox","Tigit","Foxp3","Jak2","Stat5b","Rorc","Ikzf3","Itgae","Ctla4","Il6","Il21","Tbx21","Ikzf4","Ccl20","Ccr6","Gpr83","Tgfb1","Dnmt3a","Prdm1","Stat3","Stat5a","Jak1","Klrc1","Tgfb3","Il1b","Il2ra"), ]

#Summarizing the data for mean scores and z-score scaling each gene
GSE123484.filtered <- GSE123484.filtered %>%  group_by(Feature_name,genotype) %>%  dplyr::summarise(average_count = mean(counts, na.rm = TRUE))
GSE123484.filtered <- GSE123484.filtered %>%  group_by(Feature_name) %>%  dplyr::mutate(min_val = min(average_count),max_val = max(average_count),z_scaled_fill = ifelse(min_val == max_val, 0, 2 * (average_count - min_val) / (max_val - min_val) - 1)) %>%dplyr::select(-min_val, -max_val)

#Reorder Genes 
GSE123484.filtered$Feature_name <- factor(GSE123484.filtered$Feature_name, levels = c("Il2","Ctla4","Il21","Ccl20","Tgfb3","Il6ra","Cd3d","Cd3e","Cd3g","Lag3","Tox","Tigit","Foxp3","Jak2","Stat5b","Rorc","Ikzf3","Itgae","Il6","Tbx21","Ikzf4","Ccr6","Gpr83","Tgfb1","Dnmt3a","Prdm1","Stat3","Stat5a","Jak1","Klrc1","Il1b","Il2ra")) 

## GGPLOT Heatmap 
p1 <- ggplot(GSE123484.filtered, aes(x=genotype, y=Feature_name, fill=z_scaled_fill)) + geom_tile() + theme_classic()
p1 <- p1 + labs(x = " ", y ="",fill="Expression")
p1 <- p1 + scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0)
p1 <- p1 + theme(axis.text.x = element_text(face="bold",size=16,color="black")) 
p1 <- p1 + theme(axis.text.y = element_text(face="bold",size=16,color="black")) 
p1 <- p1 + theme(strip.text = element_text(face = "bold",size=16,color="black")) 
p1 <- p1 + theme(axis.title = element_text(face = "bold",size = 16,color="black"))
p1 <- p1 + theme(legend.text = element_text(face = "bold",size=16,color="black"))
p1 <- p1 + theme(legend.title = element_text(face = "bold",size=16,color="black"))
p1

#### AuCell scoring ####

#loading library
library(AUCell)

# Load dataframe
GSE123484.long <- readRDS("GSE123484_expression_long.rds")

#filtering out extregs and in vitro
GSE123484.long <- GSE123484.long[GSE123484.long$celltype == 'Tregs', ]
GSE123484.long <- GSE123484.long[GSE123484.long$strain == 'B6.H2d.SKG', ]

# Convert expression dataframe to wide format
GSE123484.wide <- reshape2::dcast(GSE123484.long, Feature_name ~ sample, value.var="counts", fill=0)
exprMat <- as.matrix(GSE123484.wide[,-1])
rownames(exprMat) <- GSE123484.wide$Feature_name

# Define gene sets (e.g., as a list)
geneSets <- read.csv("/projects/users/wilsosx11/250110_PTPN2_scRNAseq/cherrypicked_treg_activation.csv")

# Calculate AUC scores
cells_rankings <- AUCell_buildRankings(exprMat, plotStats=FALSE, verbose=FALSE)
auc_results <- AUCell_calcAUC(geneSets, cells_rankings)

# View results
auc_matrix <- getAUC(auc_results)
AuCell.scores_long <- melt(auc_matrix)
colnames(AuCell.scores_long) <- c("Barcodes","AuCell.score")


#### Loading in GSE123486: PTPN2 promotes pathogenic loss of FoxP3+ in RORgt+ Tregs [ATAC-seq] ####

#loading in bigwig
bw_files <- list.files(path = "/projects/users/wilsosx11/251208_PTPN2_public_datasets/GSE123486_raw_data/", pattern = "\\.bw$", full.names = TRUE)
bw_list <- lapply(seq_along(bw_files), function(i) {
  gr <- import(bw_files[i])
  mcols(gr)$source_file <- basename(bw_files[i]) # Or use bw_files[i] for the full path
  gr
})

# Combine into one GRanges object
bw_combined <- do.call(c, bw_list)
GSE123486.bw <- bw_combined

#making bigwig dataframe
GSE123486.bw.df <- as.data.frame(GSE123486.bw)

# List all .txt.gz files in your directory
files <- list.files(path = "/projects/users/wilsosx11/251208_PTPN2_public_datasets/GSE123486_raw_data/", pattern = "\\.txt\\.gz$", full.names = TRUE)

# Read and concatenate, adding file name as a column
data_list <- lapply(files, function(f) {
  df <- read_tsv(f, col_names = FALSE)
  df$filename <- basename(f)
  df
})
all_data <- bind_rows(data_list)
GSE123486.peaks <- all_data

#renaming columns within GSE123486.peaks
names(GSE123486.peaks)[1:10] <- c("Chromosome","start","end","name","score","strand","signalValue","pValue(-log10)","qValue(-log10)","peak")

# uploading experiment metadata:
GSE123486.meta <- read_csv("GSE123486_metadata.csv")

# converting filename to sample name 
unique(GSE123486.peaks$filename)
GSE123486.peaks$sample <- GSE123486.peaks$filename
GSE123486.peaks$sample <- sub('_.*', '', GSE123486.peaks$sample)
unique(GSE123486.peaks$sample)

# merging GSE123486.meta with GSE123486.peaks
GSE123486.peaks <- merge(GSE123486.peaks,GSE123486.meta,by=c("sample"))

#creating a sample column 
GSE123486.bw.df$sample <- GSE123486.bw.df$source_file
GSE123486.bw.df$sample <- sub('_.*', '', GSE123486.bw.df$sample)
unique(GSE123486.bw.df$sample)

# merging GSE123486.meta with GSE123486_bw.df
GSE123486.bw.df <- merge(GSE123486.bw.df,GSE123486.meta,by=c("sample"))

# Saving as rds objects 
saveRDS(GSE123486.bw, file = "GSE123486_bigwig_data.rds")
saveRDS(GSE123486.bw.df,file="GSE123486_bigwig_dataframe.rds")
saveRDS(GSE123486.peaks, file = "GSE123486_peaks_data.rds")

#Loading rds objects
GSE123486.bw <- readRDS("GSE123486_bigwig_data.rds")
GSE123486.bw.df <- readRDS("GSE123486_bigwig_dataframe.rds")
GSE123486.peaks <- readRDS("GSE123486_peaks_data.rds")

#### Differential Peaks GSE123486: PTPN2 promotes pathogenic loss of FoxP3+ in RORgt+ Tregs [ATAC-seq] ####

#necessary libraries
library(GenomicRanges)
library(GenomicFeatures)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # Pre-built TxDb for mm10
library(org.Mm.eg.db) # Mouse gene symbols


#Loading rds objects
GSE123486.peaks <- readRDS("GSE123486_peaks_data.rds")

#filtering for only Tregs and 0h
unique(GSE123486.peaks$celltype)
unique(GSE123486.peaks$timepoint)
GSE123486.peaks.filtered <- GSE123486.peaks[GSE123486.peaks$celltype == 'Treg', ]
GSE123486.peaks.filtered <- GSE123486.peaks.filtered[GSE123486.peaks.filtered$timepoint == '0', ]

# creating 0h gr object and identifying unique peaks
gr <- GRanges(seqnames = GSE123486.peaks.filtered$Chromosome,
              ranges = IRanges(start = GSE123486.peaks.filtered$start, end = GSE123486.peaks.filtered$end))
overlaps <- countOverlaps(gr)
non_overlapping <- gr[overlaps == 1]
GSE123486.0h.unique.peaks <- as.data.frame(non_overlapping)

#filtering for only Tregs and 24hh
unique(GSE123486.peaks$celltype)
unique(GSE123486.peaks$timepoint)
GSE123486.peaks.filtered <- GSE123486.peaks[GSE123486.peaks$celltype == 'Treg', ]
GSE123486.peaks.filtered <- GSE123486.peaks.filtered[GSE123486.peaks.filtered$timepoint == '24h', ]

# creating 24h gr object and identifying unique peaks
gr <- GRanges(seqnames = GSE123486.peaks.filtered$Chromosome,
              ranges = IRanges(start = GSE123486.peaks.filtered$start, end = GSE123486.peaks.filtered$end))
overlaps <- countOverlaps(gr)
non_overlapping <- gr[overlaps == 1]
GSE123486.24h.unique.peaks <- as.data.frame(non_overlapping)

#filtering for only Tregs and 48h
unique(GSE123486.peaks$celltype)
unique(GSE123486.peaks$timepoint)
GSE123486.peaks.filtered <- GSE123486.peaks[GSE123486.peaks$celltype == 'Treg', ]
GSE123486.peaks.filtered <- GSE123486.peaks.filtered[GSE123486.peaks.filtered$timepoint == '48h', ]

# creating 48h gr object and identifying unique peaks
gr <- GRanges(seqnames = GSE123486.peaks.filtered$Chromosome,
              ranges = IRanges(start = GSE123486.peaks.filtered$start, end = GSE123486.peaks.filtered$end))
overlaps <- countOverlaps(gr)
non_overlapping <- gr[overlaps == 1]
GSE123486.48h.unique.peaks <- as.data.frame(non_overlapping)

#filtering for only Tregs and 72hh
unique(GSE123486.peaks$celltype)
unique(GSE123486.peaks$timepoint)
GSE123486.peaks.filtered <- GSE123486.peaks[GSE123486.peaks$celltype == 'Treg', ]
GSE123486.peaks.filtered <- GSE123486.peaks.filtered[GSE123486.peaks.filtered$timepoint == '72h', ]

# creating 72h gr object and identifying unique peaks
gr <- GRanges(seqnames = GSE123486.peaks.filtered$Chromosome,
              ranges = IRanges(start = GSE123486.peaks.filtered$start, end = GSE123486.peaks.filtered$end))
overlaps <- countOverlaps(gr)
non_overlapping <- gr[overlaps == 1]
GSE123486.72h.unique.peaks <- as.data.frame(non_overlapping)

##merging unique peaks

#Adding timepoint
GSE123486.0h.unique.peaks$timepoint <- "0h"
GSE123486.24h.unique.peaks$timepoint <- "24h"
GSE123486.48h.unique.peaks$timepoint <- "48h"
GSE123486.72h.unique.peaks$timepoint <- "72h"

#merging data
GSE123486.all.timepoints.unique.peaks <- rbind(GSE123486.0h.unique.peaks,GSE123486.24h.unique.peaks,GSE123486.48h.unique.peaks,GSE123486.72h.unique.peaks)

## adding Gene annotation

# Assuming df_peaks and df_genes are your peak and gene dataframes
gr_peaks <- GRanges(
  seqnames = GSE123486.all.timepoints.unique.peaks$seqnames,
  ranges = IRanges(start = GSE123486.all.timepoints.unique.peaks$start, end = GSE123486.all.timepoints.unique.peaks$end)
)

# Assume gr_peak is your GRanges object with peaks
# Annotate peaks using ChIPseeker
peak_anno <- annotatePeak(
  gr_peaks,
  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
  tssRegion = c(-1000, 1000), # +/- 1kb around TSS
  annoDb = 'org.Mm.eg.db'
)

#creating annotated dataframe
df_peak_anno <- as.data.frame(peak_anno)
df_peak_anno <- df_peak_anno[, c('seqnames', 'start', 'end', 'annotation', 'SYMBOL')]

#merging with GSE123486.all.timepoints.unique.peaks 
GSE123486.all.timepoints.unique.peaks <- merge(GSE123486.all.timepoints.unique.peaks,df_peak_anno,by=c('seqnames', 'start', 'end'))

# Saving as rds objects 
saveRDS(GSE123486.all.timepoints.unique.peaks, file = "GSE123486.all.timepoints.unique.peaks.rds")

#Loading rds objects
GSE123486.all.timepoints.unique.peaks <- readRDS("GSE123486.all.timepoints.unique.peaks.rds")


#### GSE123488 analysis: Overlap of Genes from GSE123484:RNAseq and GSE123486:ATACseq ####

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
  area1 = expression.genes %>% length(),
  area2 = peaks.genes %>% length(),
  cross.area = GSE123488.common.genes %>% length(),
  
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
grid.draw(venn_object)




## Next steps...
# Meme Suite AME analysis of enriched motifs. Check to see which TFs affect Treg function
# Meme Suite FIMO analysis of motifs. Which are the most occuring?
# Enrichr Gene list to see which TF targets are found from common genes
# Enrichr which T cell functions/differentiation are impacted 

#### GSE123488 analysis: Counting the number of binding sites per Treg TF of interest ####

#Loading in FIMO results...update with fixed fimo
GSE123488.fimo <- read.csv("/projects/users/wilsosx11/251208_PTPN2_public_datasets/meme_suite/fimo.tsv",sep="\t")

# Filtering fimo for significant instances...
GSE123488.fimo <- GSE123488.fimo[GSE123488.fimo$p.value <= 0.05, ]

#counting number of instances detected
GSE123488.fimo.summary <- GSE123488.fimo %>% group_by(motif_alt_id) %>% dplyr::summarise(count = n())

# Filtering summary for Treg specific transcription factors
target_motifs <- c(
  'IKZF2', 'IKZF1', 'Runx1', 'BATF3', 'FOXP3', 'Nfatc2', 'TCF7L2', 'IRF4',"STAT3","STAT1","STAT1::STAT2",
  'NR4A1', 'Gata3', 'NR4A2', 'TCF7L1', 'TCF7', 'Stat5a::Stat5b', 'Nfatc1', 'Ikzf3'
)
test <- GSE123488.fimo.summary[GSE123488.fimo.summary$motif_alt_id %in% target_motifs, ]


#### GSE123488 analysis: Creating bed file of TF binding sites for IGV ####

#Loading in FIMO results...update with fixed fimo
GSE123488.fimo <- read.csv("/projects/users/wilsosx11/251208_PTPN2_public_datasets/meme_suite/fimo.tsv",sep="\t")

# Filtering fimo for significant instances...
GSE123488.fimo <- GSE123488.fimo[GSE123488.fimo$p.value <= 0.05, ]

#Subsetting for necessary columns
colnames(GSE123488.fimo)
test <- GSE123488.fimo[c("sequence_name","start","stop","motif_alt_id")]

#sorting the bed file
test <- test[order(test$start), ]
test <- test[order(test$sequence_name), ]

##exporting
#write.table(test, file = "GSE123488_all_tfs.bed", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
# follow up in terminal with bgzip then tabix

## Filtering for only STAT proteins 
test2 <- test[grepl('STAT', test$motif_alt_id, ignore.case = TRUE), ]
write.table(test2, file = "GSE123488_all_STAT_tfs.bed", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
# follow up in terminal with bgzip then tabix

## Filtering for only STAT3 proteins
test3 <- test[test$motif_alt_id == "STAT3", ]
write.table(test3, file = "GSE123488_all_STAT3.bed", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
# follow up in terminal with bgzip then tabix

#### GSE123488 analysis: Gene set enrichment for TF Genes ####

#loading in data
tf.enrichr <- read.csv("260108_GSE123486_TF_geneset_enrichment.csv")

#renaming Dataset column
unique(tf.enrichr$Dataset)
tf.enrichr$Dataset[tf.enrichr$Dataset == "TF_PPI"] <- "PPI"
tf.enrichr$Dataset[tf.enrichr$Dataset == "Jasper_Transfac_PWM"] <- "PWM"
tf.enrichr$Dataset[tf.enrichr$Dataset == "TF_Coexpression"] <- "Coexpression"
tf.enrichr$Dataset[tf.enrichr$Dataset == "TF_Perturbations"] <- "Perturbations"

#removing PWM because it doesn't have significant values
tf.enrichr<- subset(tf.enrichr, Dataset != 'PWM')

#Making a categorical P-value column 
tf.enrichr$P.value.cat <- ifelse(tf.enrichr$P.value < 0.05, 'Significant', 'NS')

# Dotplot 
p <- ggplot(tf.enrichr, aes(x=Odds.Ratio, y=Transcription.Factor,color=P.value.cat)) + geom_point(binaxis='y', stackdir='center',aes(size=Odds.Ratio)) + theme_classic()
p <- p + facet_grid(~Dataset) +xlim(0,5)
p <- p + scale_color_manual(values = c('NS' = 'darkgrey', 'Significant' = 'darkblue'))
p <- p + scale_size(range = c(3,0))
p <- p + theme(axis.text.x = element_text(face="bold",size=12,color="black")) 
p <- p + theme(axis.text.y = element_text(face="bold",size=12,color="black")) 
p <- p + theme(strip.text = element_text(face = "bold",size=12,color="black")) 
p <- p + theme(axis.title = element_text(face = "bold",size = 12,color="black"))
p <- p + theme(legend.text = element_text(face = "bold",size=12,color="black"))
p <- p + theme(legend.title = element_text(face = "bold",size=12,color="black"))
p <- p + theme(title = element_text(face = "bold",size=12,color="black"))
p <- p + labs(x = "Odds Ratio", y = "Transcription Factor", color = "P-value",size="Odds Ratio") 
p 

#### GSE123488 analysis: Gene set enrichment for Treg PTPN2+/- haploinsufficiency ####

#loading in data
path.enrichr <- read.csv("260107_GSE123486_enrichr_pathways.csv")


#Making a categorical P-value column 
path.enrichr$P.value.cat <- ifelse(path.enrichr$P.value < 0.05, 'Significant', 'NS')

# Dotplot 
p <- ggplot(path.enrichr, aes(x=Odds.Ratio, y=reorder(Term,Odds.Ratio),color=P.value)) + geom_point(binaxis='y', stackdir='center',aes(size=Odds.Ratio)) + theme_classic()
p <- p + theme(axis.text.x = element_text(face="bold",size=16,color="black")) 
p <- p + theme(axis.text.y = element_text(face="bold",size=16,color="black")) 
p <- p + theme(strip.text = element_text(face = "bold",size=16,color="black")) 
p <- p + theme(axis.title = element_text(face = "bold",size = 16,color="black"))
p <- p + theme(legend.text = element_text(face = "bold",size=16,color="black"))
p <- p + theme(legend.title = element_text(face = "bold",size=16,color="black"))
p <- p + theme(title = element_text(face = "bold",size=12,color="black"))
p <- p + labs(x = "Odds Ratio", y = "Pathways", color = "P-value",size="Odds Ratio") 
p 

#### The End ####