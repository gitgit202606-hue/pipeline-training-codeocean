# pipeline-training-codeocean

Co-expression Analysis: A co-expression analysis of the prioritized CRC targets (CDH17, GUCY2C, CDCP1, LY6G6D) across Caris datasets, PDX models, and cell lines is certainly feasible. I’ll need to confirm the sources for PDX model and cell line datasets with our colleagues. If you already have specific data in mind, please let me know!
 
CRC Model Selection & Cell Tinder Algorithm: I agree that prioritizing models with robust co-expression and surface antigen correlation is essential. Once the datasets are available, we can explore applying the cell tinder algorithm to assess cell line models for your functional assays.
 
GeoMX Data: Our CO team is launching GeoMX spatial data for multiple cancer indications, and I am currently collecting CART targets for both CRC and prostate cancer. Thank you for sending the prioritized CRC targets! Once the GeoMX data is available, I anticipate it will complement the target surface expression evaluation process I proposed for the rubric ranking scores and be informative for target and model selection.
 
On-Target/Off-Tumor Analysis: Evaluating co-expression in healthy tissues is a key step for safety assessment and de-risking through dual-antigen targeting. I will look into co-expression using GTEx datasets.

 I can run uPAR using the rubric we previously presented for CRC. assess this marker based on LoT or post top in some of our datasets, 

Hi Yao,
 
I wanted to thank you again for a lovely presentation several weeks ago. Your analysis has given us several ideas that we wanted to run by you, as follow up.
 
First, would it be possible to perform a co-expression analysis of some of our current prioritized targets, namely CDH17, GUCY2C, CDCP1, and LY6G6D in CRC separated by Caris, PDX samples, and cell lines?
 
A key aim of ours is to identify CRC PDXs and cell lines that we think may be utilized as workhorse models for functional studies. This would be based on having expression of multiple antigens of interest (e.g. + for CDH17, GUCY2C, and LY6G6D - god willing), not to mention expression of these antigens that correlates well with surface expression levels on tumors from real world patient data. Another bonus would be applying the cell tinder algorithm to identify cell lines that share relevant expression profiles in addition to physiological levels of surface antigen expression.
 
Second, I also heard that we may have Cosmix/GeoMX data for CRC. Do you think that you could analyze those data to inform your rankings?
 
Third, regarding the risk of on target/off tumor tox, is it possible to do co-expression analyses in healthy tissue where we see expression of at least one of these antigens? This might also help us identify de-risking strategies leveraging dual-antigen target approaches.
 
If you have any questions and want to follow up on this email, please feel free to reach out.
 
## %pip install --user upsetplot 
from upsetplot import from_contents, plot

## Define positivity threshold on log1pTPM, an expression > 1.5/2 is standard

thresh = 2.0

possible_targets = ["PLAUR", "GUCY2C", "CDH17","CDCP1","LY6G6D"]

def extract_positive_samples(df,name = 'Dataset'):
    pos_dict = {}

    targets = [gene for gene in possible_targets if gene in df.columns]
    print(f'processing {name}: found {len(targets)} of {len(possible_targets)} targets.')
    
    for gene in targets:
        positive_mask = df[gene] >= thresh

        pos_dict[gene] = df[positive_mask].index.unique().tolist()

    return pos_dict

# format lists for upset plotting

upset_caris = from_contents(extract_positive_samples(caris_long_log,'Caris'))
upset_pdx = from_contents(extract_positive_samples(pdx_long_crc,'PDX'))
upset_ccle = from_contents(extract_positive_samples(ccle_df_crc,'CCLE Lines'))

#plot1 the patient blueprint
fig1 = plt.figure(figsize = (10,6))
plot(upset_caris, fig = fig1, element_size=None,show_counts = True)
plt.suptitle('Fig 1: Caris Patient Multi-Antigen Intersection Niche (CRC)',
            fontsize=14, y = 1.02)
plt.show()
          
#plot2&3, side by side model comparison
fig, axes = plt.subplots(1,2, figsize = (18,6))

plot(upset_pdx, fig= fig, ax = axes[0], show_counts=True)
axes[0].set_title('Fig 2: Internal PDX CRC Overlaps')

plot(upset_ccle, fig= fig, ax = axes[1], show_counts=True)
axes[1].set_title('Fig 3: CCLE Cell Line CRC Overlaps')

plt.tight_layout()
plt.show()




 
