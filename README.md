# pipeline-training-codeocean

Co-expression Analysis: A co-expression analysis of the prioritized CRC targets (CDH17, GUCY2C, CDCP1, LY6G6D) across Caris datasets, PDX models, and cell lines is certainly feasible. I’ll need to confirm the sources for PDX model and cell line datasets with our colleagues. If you already have specific data in mind, please let me know!
 
CRC Model Selection & Cell Tinder Algorithm: I agree that prioritizing models with robust co-expression and surface antigen correlation is essential. Once the datasets are available, we can explore applying the cell tinder algorithm to assess cell line models for your functional assays.
 
GeoMX Data: Our CO team is launching GeoMX spatial data for multiple cancer indications, and I am currently collecting CART targets for both CRC and prostate cancer. Thank you for sending the prioritized CRC targets! Once the GeoMX data is available, I anticipate it will complement the target surface expression evaluation process I proposed for the rubric ranking scores and be informative for target and model selection.
 
library(Gviz)
library(GenomicRanges)
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#library(org.Mm.eg.db)
#library(dplyr)
#library(rtracklayer)

options(ucscChromosomeNames = F)

target_chr <- 'chr12'
target_start <- 112810000
target_end <- 112820000 ## DNMT3A Exon 1 promotor TSS region

bigwig_dir <- '/projects/users/wilsosx11/251208_PTPN2_public_datasets/GSE123486_raw_data/'
GSE123486.meta <- read.csv('/projects/users/wilsosx11/251208_PTPN2_public_datasets/GSE123486_metadata.csv')

all_bw_files <- list.files(
  path = bigwig_dir,
  pattern = '\\.bw$',
  full.names = F
)

track_mapping <- data.frame(file_name = all_bw_files) %>% 
  mutate(sample = sub('_.*', '',file_name)) %>% 
  inner_join(GSE123486.meta, by = 'sample') %>% 
  filter(celltype == 'Treg') %>% 
  mutate(timepoint_num = as.numeric(gsub('h','',timepoint)),
         genotype_order = ifelse(genotype == 'Ptpn2+/+',1,2)) %>% 
  arrange(timepoint_num,genotype_order)

atac_tracks <- lapply(1:nrow(track_mapping), function(i) {
 
  row <- track_mapping[i, ]
  
  track_color <- ifelse(row$genotype == 'Ptpn2+/+', "darkblue",'red')
  track_label <- paste0(row$genotype, '', row$timepoint, ' Il6 stimulation')
  
  DataTrack(
    range = row$file_path,
    genome = 'mm10',
    name = track_label,
    type = 'mountain',
    col.mountain = track_color,
    fill.mountain = c(track_color,track_color),
    ylim = c(0,5)
  )
  
})

axisTrack <- GenomeAxisTrack()

fimo_ranges = GRanges(
  seqnames = target_chr,
  ranges = IRanges(
    start = c(112813500, 112817200),
    end = c(112813515, 112817215)
  ),
  
  id = c('STAT3 Site 1','STAT3 Site 2')
)

fimoTrack = AnnotationTrack(
  range = fimo_ranges,
  genome = 'mm10',
  name = 'FIMO STAT3',
  fill = 'black',
  col = 'black',
  shape = 'box',
  showFeatureId = T,
  fontcolor.feature = 'black',
  fontsize.feature = 8
)

geneTrack = UcscTrack(
  genome = 'mm10',
  chromosome = target_chr,
  track = 'refGene',
  from = target_start,
  to = target_end,
  trackType =  'GeneRegionTrack',
  rxtRetrieveTranscripts = T,
  showId=T,
  geneSymbol = T,
  name = 'Dnmt3a Locus',
  fill = 'darkgray',
  col = 'black'
)

all_track <- c(
  list(axis_track),
  atac_tracks,
  list(fimoTrack,geneTrack)
)

track_sizes = c(1,rep(2,length(atac_tracks)),1,1.5)

png('~/project_analysis/PTPN2i/20260706_PTPN2i_scRNA_paper_swyz/Panel_D_Dnmt31_Fixed.png',width = 1600,
    height = 750,res = 120)

plotTracks(
  all_track,
  from = target_start,
  to = target_end,
  chromosome = target_chr,
  sizes = track_sizes,
  background.title = 'transparent',
  col.title = 'black',
  col.axis = 'black'

)

dev.off()

> geneTrack = UcscTrack(
+   genome = 'mm10',
+   chromosome = target_chr,
+   track = 'refGene',
+   from = target_start,
+   to = target_end,
+   trackType =  'GeneRegionTrack',
+   rxtRetrieveTranscripts = T,
+   showId=T,
+   geneSymbol = T,
+   name = 'Dnmt3a Locus',
+   fill = 'darkgray',
+   col = 'black'
+ )
Error in if (file == "" || length(file) == 0) stop("empty or no content specified") : 
  missing value where TRUE/FALSE needed






