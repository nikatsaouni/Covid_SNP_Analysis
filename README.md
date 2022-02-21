# Covid_SNP_Analysis


Functions for the analysis of SNP data.
- summary statistics
- logistic regression
- gene-based test

Steps to run the workflow:
-----------------------------------
1)  src.d01_data_intermediate/load_data.py. ## load data
2)  src.d02_processing/create_plink_inputs.py ## create files for PLINK (.ped, .map, covariate file)
3)  src.d03_modelling/log_regression.py. ## Run PLINK logistic regression
4)  src.d03_modelling/data_gene_based_Test.py. ## Prepare files for gene-based test

5)  Run R script for gene-based test

help functions are saved under src.d00_utils

For  details concerning the layout of the project please refer to layout.txt

Raw data must be saved under /data/01_raw/

