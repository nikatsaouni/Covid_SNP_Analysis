from src.d00_utils.gene_based_test_prepar import generate_patient_snp_matrix
import src.cases.parameters as p
from src.d00_utils.read_process_data import read_csv_file

path_of_SNPs_After_LG = p.Report_Results + "log_Regression_Genes.csv"
list_of_SNPs = read_csv_file(path_of_SNPs_After_LG)
print(list_of_SNPs)
generate_patient_snp_matrix(list_of_SNPs.SNP)
