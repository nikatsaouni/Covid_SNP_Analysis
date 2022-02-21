# import parameters as p
import src.cases.parameters as p
import pandas as pd
# disable chained assignments
pd.options.mode.chained_assignment = None
import os
import numpy as np

def read_data():
    SNP_dataframe = read_snp_csv()
    patient_data = read_pat_data()

    ## Keep patient data only for the patients with SNP information
    patient_data = patient_data[patient_data['Lab ID, DK.....'].isin(list(SNP_dataframe.sample_ID))]
    patient_data = patient_data.rename(columns={"Lab ID, DK.....": "sample_ID"}, errors="raise")

    ### Split data to severe and non-severe conditions
    severe_cases = list(patient_data.loc[patient_data['Clinical Score'] > p.clinic_score_split, 'sample_ID'])
    no_severe_cases = list(patient_data.loc[patient_data['Clinical Score'] <= p.clinic_score_split, 'sample_ID'])

    ## Add clinical_score column
    patient_data['Clinical_Score_'] = patient_data.apply(add_clinical_Score, args = (no_severe_cases,), axis=1)
    SNP_dataframe['Clinical_Score'] = SNP_dataframe.apply(add_clinical_Score, args = (no_severe_cases,), axis=1)

    ## SNP data for patients with high clinical Score
    SNP_low_Score = SNP_dataframe[SNP_dataframe['sample_ID'].isin(no_severe_cases)]
    ## SNP data for patients with low clinical Score
    SNP_high_Score = SNP_dataframe[SNP_dataframe['sample_ID'].isin(severe_cases)]
    print('Number of patients with high Clinical Score (>='+ str(p.clinic_score_split+1) + '):', SNP_high_Score['sample_ID'].nunique())
    print('Number of patients with low Clinical Score (<'+ str (p.clinic_score_split+1) + '):', SNP_low_Score['sample_ID'].nunique())

    ### convert snp columns to str
    col_to_str = ["Chromosome", "Genomic"]
    SNP_dataframe[col_to_str] = SNP_dataframe[col_to_str].astype(str)
    SNP_low_Score[col_to_str] = SNP_low_Score[col_to_str].astype(str)
    SNP_high_Score[col_to_str] = SNP_high_Score[col_to_str].astype(str)

    ### add column SNP unique
    SNP_dataframe["SNP"] = SNP_dataframe["Chromosome"] + SNP_dataframe["Genomic"]
    SNP_low_Score["SNP"] = SNP_low_Score["Chromosome"] + SNP_low_Score["Genomic"]
    SNP_high_Score["SNP"] = SNP_high_Score["Chromosome"] + SNP_high_Score["Genomic"]

    ### keep only important information
    SNP_pat_score = SNP_dataframe[['sample_ID', 'Position', 'Chromosome', 'Known', 'Gene', 'Genomic', 'Clinical_Score', "Zygosity", "SNP",
         "Prediction"]]
    total_number_of_patients = SNP_pat_score.sample_ID.unique().size
    print('Total Number of patients:', total_number_of_patients)
    print("--------------------------------------------------")

    ## How many snp occurences per clinical score
    list_of_rare_Snps = summary_snps_per_cs(SNP_pat_score)

    ## add column with allele information
    SNP_pat_score = zygosity(SNP_pat_score)

    ## Save SNP data for all patients to csv file
    path_to_save_SNP_data = p.Interm_Results + "SNP_data.csv"

    SNP_pat_score.to_csv(path_to_save_SNP_data,
                         index=False)
    print("--------------------------------------------------")
    print("SNP data for all patiens saved in: ", path_to_save_SNP_data)

    pd.DataFrame(list_of_rare_Snps).to_csv(p.Interm_Results + "rare_snps.csv",
                         index=False,  header=True)

    patient_data.to_csv(p.Interm_Results + "patient_data.csv",
                         index=False)

    return SNP_pat_score, list_of_rare_Snps



## Calculate rare SNPs and save Summary for SNPs per clinical score
def summary_snps_per_cs(SNP_pat_score, cutoff_SNPs = p.cutoff_SNPs):
    SNP_pat_score_Crostab = pd.crosstab(SNP_pat_score['SNP'], SNP_pat_score['Clinical_Score'])
    SNP_pat_score_Crostab.to_csv(p.Report_Results + 'Summary_SNPs_per_Clinical_Score.csv')
    print(f"Number of SNPs per clinical score saved in: " + p.Interm_Results + 'Summary_SNPs_per_Clinical_Score.csv')
    print("--------------------------------------------------")

    list_of_rare_Snps = list(SNP_pat_score_Crostab[(SNP_pat_score_Crostab[0] < cutoff_SNPs) & (
                SNP_pat_score_Crostab[1] < cutoff_SNPs)].index)
    print(f"Total number of SNPs: " + str(len(SNP_pat_score.SNP.unique())))
    print(f"SNPs with more than " + str(p.cutoff_SNPs) + " occurences for at least one clinical score: " + str(len(SNP_pat_score.SNP.unique())- len(list_of_rare_Snps)))

    return list_of_rare_Snps


## Read csv SNP files and handle formatting issues
def read_snp_csv(snp_folder = p.snp_folder):
    SNP_dataframe = pd.DataFrame() # define dataframe to save all the SNP data

    # for all the csv files
    for filename in os.listdir(snp_folder):
        file_path = snp_folder + filename
        file_df = pd.read_csv(file_path, sep=',', skiprows=[0, 1], usecols=[i for i in range(1, 14)]) ## skip first rwo rows and first col
        file_df.insert(loc=0, column='sample_ID', value=int(filename[6:-4])) #create sample_ID column with the value of the file name
        SNP_dataframe = SNP_dataframe.append(file_df)
    SNP_dataframe = SNP_dataframe.reset_index(drop=True)

    ### Handle frequency issue
    # get all rows by mask
    mask = SNP_dataframe['Frequency'] == 0
    c = ['Gene', 'Prediction', 'Type', 'Known', 'MAF', 'Zygosity', 'Sanger Files']
    # shift columns, but necessary converting to strings
    SNP_dataframe.loc[mask, c] = SNP_dataframe.loc[mask, c].astype(str).shift(-1, axis=1)
    SNP_dataframe = SNP_dataframe.drop(['Sanger Files'], axis=1)
    SNP_dataframe = SNP_dataframe.drop(['Coding', 'Protein', 'Frequency', 'Type'], axis=1)


    return SNP_dataframe

## Add column with clinical score
def add_clinical_Score(row, no_severe_cases):
    if row["sample_ID"] in no_severe_cases:
        val = 0
    else:
        val = 1
    return val

## Read patient data from xlsx file
def read_pat_data(patient_file = p.patient_data):
    # Read xlsx file
    patient_data = pd.read_excel(patient_file)

    return patient_data

## Add zygosity information into the dataframe
def zygosity(SNP_pat_score):
    ### Add allele information
    SNP_pat_score['A1'] = SNP_pat_score['Genomic'].str.split('>', expand=True)[0].str[-1]
    SNP_pat_score['A2'] = SNP_pat_score['A3'] = SNP_pat_score['Genomic'].str.split('>', expand=True)[1]

    heter_dictionary = {
        "A": "G",
        "G": "A",
        "T": "C",
        "C": "T"
    }
    for row in range(SNP_pat_score.shape[0]):
        SNP_pat_score['A3'][row] = heter_dictionary[SNP_pat_score['A1'][row]]

    SNP_pat_score.loc[SNP_pat_score['Zygosity'] == 'Homozygous', 'A1'] = SNP_pat_score["A2"]
    SNP_pat_score['A1_A2'] = SNP_pat_score['A1'] + ' ' + SNP_pat_score['A2']

    SNP_pat_score.A1_A2 = np.where(SNP_pat_score.Zygosity == 'Homozygous',
                                   SNP_pat_score.A1_A2.fillna(SNP_pat_score.A1 + ' ' + SNP_pat_score.A1),
                                   SNP_pat_score.A1_A2)
    SNP_pat_score.A1_A2 = np.where(SNP_pat_score.Zygosity == 'Heterozygous',
                                   SNP_pat_score.A1_A2.fillna(SNP_pat_score.A1 + ' ' + SNP_pat_score.A3),
                                   SNP_pat_score.A1_A2)
    SNP_pat_score.A1_A2 = np.where(SNP_pat_score.Zygosity == 'Homozygous',
                                   SNP_pat_score.A1_A2.fillna(SNP_pat_score.A3 + ' ' + SNP_pat_score.A3),
                                   SNP_pat_score.A1_A2)

    SNP_pat_score = SNP_pat_score.drop(['A1', 'A2', 'A3'], axis=1)

    return SNP_pat_score



def read_csv_file(file_path):
    file_df = pd.read_csv(file_path)

    return file_df



