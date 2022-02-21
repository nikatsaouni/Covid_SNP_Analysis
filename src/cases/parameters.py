
# Directories with data
snp_folder = "/Users/nikatsaouni/Documents/Disease_Prediction/Covid_Proj/CovidSNPAnalysis/data/01_raw/Corrected_csv_files_folder/"
patient_data =  "/Users/nikatsaouni/Documents/Disease_Prediction/Covid_Proj/CovidSNPAnalysis/data/01_raw/Tabelle_Endfassung_040521.xlsx"

# Parameters for calculations
cutoff_SNPs = 4  # set the minimum number of SNPs per clinical score to be accounted for the calculations
clinic_score_split = 3  # How to split the data to high and low clinical score
senario = "Case_2"

# Directories to save results
Interm_Results = "/Users/nikatsaouni/Documents/Disease_Prediction/Covid_Proj/CovidSNPAnalysis/data/02_intermediate/" + senario + "/"
Processed_Results = "/Users/nikatsaouni/Documents/Disease_Prediction/Covid_Proj/CovidSNPAnalysis/data/03_processed/" + senario + "/"
Model_Results = "/Users/nikatsaouni/Documents/Disease_Prediction/Covid_Proj/CovidSNPAnalysis/data/05_model_output/" + senario + "/"
Report_Results = "/Users/nikatsaouni/Documents/Disease_Prediction/Covid_Proj/CovidSNPAnalysis/data/06_reporting/" + senario + "/"


