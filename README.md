# pipeline-training-codeocean

Co-expression Analysis: A co-expression analysis of the prioritized CRC targets (CDH17, GUCY2C, CDCP1, LY6G6D) across Caris datasets, PDX models, and cell lines is certainly feasible. I’ll need to confirm the sources for PDX model and cell line datasets with our colleagues. If you already have specific data in mind, please let me know!
 
CRC Model Selection & Cell Tinder Algorithm: I agree that prioritizing models with robust co-expression and surface antigen correlation is essential. Once the datasets are available, we can explore applying the cell tinder algorithm to assess cell line models for your functional assays.
 
GeoMX Data: Our CO team is launching GeoMX spatial data for multiple cancer indications, and I am currently collecting CART targets for both CRC and prostate cancer. Thank you for sending the prioritized CRC targets! Once the GeoMX data is available, I anticipate it will complement the target surface expression evaluation process I proposed for the rubric ranking scores and be informative for target and model selection.
 
#BiocManager::install(c('Gviz'))

library(Gviz)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

options(ucscChromosomeNames = F)

txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
dnmt3a_gene_id <- select(org.Mm.eg.db,
                         keys = 'Dnmt3a', columns = 'ENTREZID',keytype = 'SYMBOL')$ENTREZID

target_chr <- 'chr12'
target_start <- 112810000
target_end <- 112930000

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

data_tracks <- lapply(seq_len(nrow(track_mapping)), function(i){
  
  row <- track_mapping[i, ]
  
  track_col <- ifelse(row$genotype == 'Ptpn2+/+', "#1f77b4",'#d62728')
  track_name <- paste0(row$genotype, '\n', row$timepoint, ' Il6 stimulation')
  
  DataTrack(
    range = file.path(bigwig_dir, row$file_name),
    genome = 'mm10',
    type = 'mountain',
    name = track_name,
    col.histogram = track_col,
    fill.histogram = track_col,
    ylim = c(0,150),
    background.title = 'white',
    col.axis = 'black',
    col.title = 'black',
    fontcolor.title = 'black',
    fontface.title = 'bold',
    cex.title = 0.6,
    rotate.title = F
  )
  
})


stat_bed_path = 'GSE123488_all_STAT3.bed.gz'

target_region = GRanges(
  seqnames = target_chr,
  ranges = IRanges(start = target_start, end = target_end)
)

Stat_gr = import(stat_bed_path,format = 'BED',which = target_region)

stat_track = AnnotationTrack(
  range = Stat_gr,
  genome = 'mm10',
  chromosome = target_chr,
  name = 'FIMO STAT binding sites',
  stacking = 'dense',
  fill = 'darkblue',
  col = 'darkblue',
  transcriptAnnotation = 'symbol',
  background.title = 'white',
  col.title = 'black',
  fontface.title = 'bold',
  cex.title = 0.6,
  rotate.title = F
)

gene_track = GeneRegionTrack(
  txdb,
  genome = 'mm10',
  chromosome = target_chr,
  start = target_start,
  end = target_end,
  name = 'Dnmt3a loci',
  fill = 'darkblue',
  col = 'darkblue',
  showId=F,
  transcriptAnnotation = 'symbol',
  background.title = 'white',
  col.title = 'black',
  fontface.title = 'bold',
  cex.title = 0.6,
  rotate.title = F
)

axis_track <- GenomeAxisTrack(col='black',fontcolor = 'black')

all_track <- c(
  list(axis_track),
  data_tracks,
  list(stat_track),
  list(gene_track)
)

num_data_tracks = length(data_tracks)

png('~/project_analysis/PTPN2i/20260706_PTPN2i_scRNA_paper_swyz/Panel_D_Dnmt31_Fixed.png',width = 1250,
    height = 800,res = 120)

plotTracks(
  all_track,
  from = target_start,
  to = target_end,
  chromosome = target_chr,
  sizes = c(0.4,rep(1,num_data_tracks),0.5,0.6),
  title.width = 2.8,
  margin = 50
)

dev.off()






