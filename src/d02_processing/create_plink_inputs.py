from src.d00_utils.read_process_data import read_csv_file
import src.cases.parameters as p
from src.d00_utils.plink_functions import generate_ped_file, generate_map_file, generate_cov_file


# Read SNP intermediate file
path_of_SNP_data= p.Interm_Results + "SNP_data.csv"
path_of_SNPs_to_Exclude =  p.Interm_Results + "rare_snps.csv"
patient_data_path = p.Interm_Results + "patient_data.csv"


SNP_pat_score = read_csv_file(path_of_SNP_data)
list_of_rare_Snps = read_csv_file(path_of_SNPs_to_Exclude)
patient_data = read_csv_file(patient_data_path)

###Generate ped_file
## SNP_pat_score_GR : df of SNP after excluding the rare varianta(as defined in the parameter file)
SNP_pat_score_GR = generate_ped_file(SNP_pat_score, list_of_rare_Snps, patient_data)
###Generate map_file
generate_map_file(SNP_pat_score, list_of_rare_Snps, SNP_pat_score_GR)

### Generate covar file
generate_cov_file(patient_data)