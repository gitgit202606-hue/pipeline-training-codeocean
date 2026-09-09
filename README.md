# pipeline-training-codeocean

Co-expression Analysis: A co-expression analysis of the prioritized CRC targets (CDH17, GUCY2C, CDCP1, LY6G6D) across Caris datasets, PDX models, and cell lines is certainly feasible. I’ll need to confirm the sources for PDX model and cell line datasets with our colleagues. If you already have specific data in mind, please let me know!
 
CRC Model Selection & Cell Tinder Algorithm: I agree that prioritizing models with robust co-expression and surface antigen correlation is essential. Once the datasets are available, we can explore applying the cell tinder algorithm to assess cell line models for your functional assays.
 
GeoMX Data: Our CO team is launching GeoMX spatial data for multiple cancer indications, and I am currently collecting CART targets for both CRC and prostate cancer. Thank you for sending the prioritized CRC targets! Once the GeoMX data is available, I anticipate it will complement the target surface expression evaluation process I proposed for the rubric ranking scores and be informative for target and model selection.



 Try claude code on codeocean vibecoding env to avoid docker container 

2. knowledge graph for car-t target eval. focused on the core abbvie indications, run our car-t target eval rubrix scoroing , 
then in this indication, for these target high expressers, what are the suppressive pathways most likely to happen, 
suggesting the armoring strategy. this is the package we 'd like to bring to odr when we do the target eval work, 
and impact we'd like to bring . maybe extend to toher TA and usecase 

3. co-expression, really data by data. gold standard is protein stain. tell the ODR protein stain is necessary to determine. 
also mention this in cart meeting, ask them to confirm if the target they are interested are included!
<img width="715" height="152" alt="image" src="https://github.com/user-attachments/assets/539cde1f-ec98-4afe-a28b-8c346555ad74" />



Key Findings:

The analysis centered on 3 PC and one CRC dataset to evaluate the TME and tumor-infiltrating lymphocytes (TILs).

T-Cell Phenotypes: Tumor in both PC and CRC exhibited an increase in exhausted and terminal CD8+ T populations. Concurrently, a decrease in CD8 progenitor cells and an increase in effector memory CD+8 cells were observed.

Inflammatory Heterogeneity: Hallmark inflammatory profiling revealed varied TME states; with certain tumors displayed elevated inflammatory signaling, indicative of an inflamed microenvironment.

Homeostatic Disconnect: A signaling starvation state was observed in IL-7 and IL-15. While TILs maintain  high receptor density and the intracellular machinery to respond (as evidenced by high enrichment scores from HCD analysis), the TME is depleted of these essential ligands.

   IL-21 axis evaluation: IL-21R expression is low but preserved in CD8+ effector cells. This suggests that     these cells remain responsive to IL-21, providing a therapeutic window to prevent terminal exhaustion  through targeted signaling.
<img width="1521" height="172" alt="image" src="https://github.com/user-attachments/assets/6622e91d-c0e3-4f22-9912-d21820f9a06a" />

