from src.d00_utils.read_process_data import read_csv_file
from src.d00_utils.plink_functions import PLINK_logistic_regression, PLINK_transform_files
import src.cases.parameters as p

path_of_SNP_data = p.Interm_Results + "SNP_data.csv"

SNP_pat_score = read_csv_file(path_of_SNP_data)
## Run Logistic Regresion with PLINK
PLINK_logistic_regression()

## Transform P-LINK output files to txt
PLINK_transform_files(SNP_pat_score)