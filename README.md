Welcome to the code repository for Hochbaum *et al.*, 2024
# MoSeq
The notebook used to generate panels and statistics related to the MoSeq experiment in **$Figure S3$** is located in 'MoSeq_figS3'
## MoSeq data:
The MoSeq related dataframes are located within: https://dataverse.harvard.edu/dataverse/2024_hochbaum_thyroid 

# Reinforcement learning modeling (Q-learning) for 2ABT
To install the Q-learning model for the 2-arm bandit analysis, install the package contained within the Q_learning_2ABT subfolder.
Instructions:
Open a terminal and navigate to the 'Q_learning_2ABT' sub-directory
```
cd Q_learning_2ABT
```
Install the package
```
pip install .
```
## 2ABT modeling notebook:
The notebook used to run Q-learning models is located in '2ABT_run_Qlearning'

## 2ABT data: 
Mouse 2ABT dataframes used for Q-learning modeling and analyses in **$Figures 5/6/S6$** are located within: https://dataverse.harvard.edu/dataverse/2024_hochbaum_thyroid 

# snRNA-seq
Scripts for processing, annotating, and analyzing the single-nucleus RNA sequencing experiments (Figures 2, 3, S4). Raw and filtered data are publicly available as of date of publication on an NCBI GEO data repository (GSE271421) and CELLxGENE repository (https://cellxgene.cziscience.com/collections/c450e15d-321a-42d6-986b-11409d04896d). 

#### 00_scOnline_sourceCode
Source code utilized from Gazestani, V., Kamath, T., Nadaf, N.M., Dougalis, A., Burris, S.J., Rooney, B., Junkkari, A., Vanderburg, C., Pelkonen, A., Gomez-Budia, M., et al. (2023). Early Alzheimer’s disease pathology in human cortex involves transient cell states. Cell 186, 4438–4453.e23. https://doi.org/10.1016/j.cell.2023.08.005.

#### 01_seuratProcessing
Scripts for initial processing of filtered feature barcode matrix data and major cell class annotation.

#### 02_labelTransfer
Example scripts for neuronal and non-neuronal subtype annotation based on published mouse motor cortex reference dataset (Yao, Z., Liu, H., Xie, F., Fischer, S., Adkins, R.S., Aldridge, A.I., Ament, S.A., Bartlett, A., Behrens, M.M., Van den Berge, K., et al. (2021). A transcriptomic and epigenomic cell atlas of the mouse primary motor cortex. Nature 598, 103–110. https://doi.org/10.1038/s41586-021-03500-8). 

#### 03_differentialExpression
Example scripts for the T3 vs. C differential expression (DE) analysis in control animals (Figure 2) and Cre+ vs. Cre- DE analysis in T3+ DN-THR or WT-THR animals (Figure 3). 

#### 04_analyses
Downstream analyses (i.e. ordered gene set enrichment analysis) performed on DE results. 

