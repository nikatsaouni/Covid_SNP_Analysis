import src.d00_utils.project_path as path


snp_folder = path.proj_dir + "/data/01_raw/Corrected_csv_files_folder/"

# Directories with data
# snp_folder = "/Users/nikatsaouni/Documents/Disease_Prediction/Covid_Proj/Covid_SNP_Analysis/data/01_raw/Corrected_csv_files_folder/"
patient_data = path.proj_dir + "data/01_raw/Tabelle_Endfassung_040521.xlsx"

# Parameters for calculations
cutoff_SNPs = 4  # set the minimum number of SNPs per clinical score to be accounted for the calculations
clinic_score_split = 2  # How to split the data to high and low clinical score
senario = "Case_1"

# Directories to save results
Interm_Results = path.proj_dir + "data/02_intermediate/" + senario + "/"
Processed_Results = path.proj_dir + "data/03_processed/" + senario + "/"
Model_Results = path.proj_dir + "data/05_model_output/" + senario + "/"
Report_Results = path.proj_dir + "data/06_reporting/" + senario + "/"


