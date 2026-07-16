# pipeline-training-codeocean

Co-expression Analysis: A co-expression analysis of the prioritized CRC targets (CDH17, GUCY2C, CDCP1, LY6G6D) across Caris datasets, PDX models, and cell lines is certainly feasible. I’ll need to confirm the sources for PDX model and cell line datasets with our colleagues. If you already have specific data in mind, please let me know!
 
CRC Model Selection & Cell Tinder Algorithm: I agree that prioritizing models with robust co-expression and surface antigen correlation is essential. Once the datasets are available, we can explore applying the cell tinder algorithm to assess cell line models for your functional assays.
 
GeoMX Data: Our CO team is launching GeoMX spatial data for multiple cancer indications, and I am currently collecting CART targets for both CRC and prostate cancer. Thank you for sending the prioritized CRC targets! Once the GeoMX data is available, I anticipate it will complement the target surface expression evaluation process I proposed for the rubric ranking scores and be informative for target and model selection.
 
/projects/users/wilsosx11/251208_PTPN2_public_datasets/GSE123486_raw_data$ ls
GSM3505003_001_sample_ID_1_S211_L008_R1_001.macs2_peaks.narrowPeak_Q0.01filt.txt.gz
GSM3505003_001_sample_ID_1_S211_L008_R1_001_NormCov.bw
GSM3505004_002_sample_ID_2_S212_L008_R1_001.macs2_peaks.narrowPeak_Q0.01filt.txt.gz
GSM3505004_002_sample_ID_2_S212_L008_R1_001_NormCov.bw
GSM3505005_003_sample_ID_3_S213_L008_R1_001.macs2_peaks.narrowPeak_Q0.01filt.txt.gz
GSM3505005_003_sample_ID_3_S213_L008_R1_001_NormCov.bw
GSM3505006_004_sample_ID_4_S214_L008_R1_001.macs2_peaks.narrowPeak_Q0.01filt.txt.gz
GSM3505006_004_sample_ID_4_S214_L008_R1_001_NormCov.bw
GSM3505007_005_sample_ID_5_S215_L008_R1_001.macs2_peaks.narrowPeak_Q0.01filt.txt.gz
GSM3505007_005_sample_ID_5_S215_L008_R1_001_NormCov.bw
GSM3505008_006_sample_ID_6_S216_L008_R1_001.macs2_peaks.narrowPeak_Q0.01filt.txt.gz
GSM3505008_006_sample_ID_6_S216_L008_R1_001_NormCov.bw
GSM3505009_007_sample_ID_7_S217_L008_R1_001.macs2_peaks.narrowPeak_Q0.01filt.txt.gz
GSM3505009_007_sample_ID_7_S217_L008_R1_001_NormCov.bw
GSM3505010_008_sample_ID_8_S218_L008_R1_001.macs2_peaks.narrowPeak_Q0.01filt.txt.gz
GSM3505010_008_sample_ID_8_S218_L008_R1_001_NormCov.bw
GSM3505011_009_sample_ID_9_S219_L008_R1_001.macs2_peaks.narrowPeak_Q0.01filt.txt.gz
GSM3505011_009_sample_ID_9_S219_L008_R1_001_NormCov.bw
GSM3505012_010_sample_ID_10_S220_L008_R1_001.macs2_peaks.narrowPeak_Q0.01filt.txt.gz
GSM3505012_010_sample_ID_10_S220_L008_R1_001_NormCov.bw
GSM3505013_011_sample_ID_11_S221_L008_R1_001.macs2_peaks.narrowPeak_Q0.01filt.txt.gz
GSM3505013_011_sample_ID_11_S221_L008_R1_001_NormCov.bw
GSM3505014_012_sample_ID_12_S222_L008_R1_001.macs2_peaks.narrowPeak_Q0.01filt.txt.gz
GSM3505014_012_sample_ID_12_S222_L008_R1_001_NormCov.bw

	
sample
title
model
strain
genotype
celltype
timepoint
1
GSM3505003
001_sample_ID_1
In vivo
FoxP3eGFP SKG
Ptpn2+/+
Treg
0
2
GSM3505004
002_sample_ID_2
In vivo
FoxP3eGFP SKG
Ptpn2+/-
Treg
0
3
GSM3505005
003_sample_ID_3
In vivo
FoxP3eGFP SKG
Ptpn2+/+
Th17
0
4
GSM3505006
004_sample_ID_4
In vivo
FoxP3eGFP SKG
Ptpn2+/-
Th17
0
5
GSM3505007
005_sample_ID_5
In vivo
FoxP3eGFP SKG
Ptpn2+/+
Treg
24h
6
GSM3505008
006_sample_ID_6
In vivo
FoxP3eGFP SKG
Ptpn2+/-
Treg
24h
7
GSM3505009
007_sample_ID_7
In vivo
FoxP3eGFP SKG
Ptpn2+/+
Treg
48h
8
GSM3505010
008_sample_ID_8
In vivo
FoxP3eGFP SKG
Ptpn2+/-
Treg
48h
9
GSM3505011
009_sample_ID_9
In vivo
FoxP3eGFP SKG
Ptpn2+/+
Treg
72h
10
GSM3505012
010_sample_ID_10
In vivo
FoxP3eGFP SKG
Ptpn2+/+
exTreg
72h
11
GSM3505013
011_sample_ID_11
In vivo
FoxP3eGFP SKG
Ptpn2+/-
Treg
72h
12
GSM3505014
012_sample_ID_12
In vivo
FoxP3eGFP SKG
Ptpn2+/-
exTreg 
72h




