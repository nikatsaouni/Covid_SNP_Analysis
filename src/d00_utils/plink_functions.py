import os
import pandas as pd
import numpy as np
import src.cases.parameters as p
from statsmodels.stats.multitest import fdrcorrection


#################################
###   FUNCTIONS FOR PLINK #######
#################################
def generate_ped_file(SNP_pat_score, list_of_rare_Snps, patient_data):
    ### SNPs and zygosity per patients
    SNP_pat_score_GR = SNP_pat_score.pivot_table(index='sample_ID', columns='SNP', values='A1_A2', aggfunc=max)
    col_index = 0
    for col in SNP_pat_score_GR.columns:
        if ('>' in col):
            SNP_pat_score_GR[col].fillna(col.split('>')[0][-1] + ' ' + col.split('>')[0][-1], inplace=True)
        else:
            #         SNP_pat_score_GR[col].fillna('0 0', inplace=True)
            SNP_pat_score_GR[col].fillna(SNP_pat_score['Genomic'].str.split('>', expand=True)[0].str[-1] + ' ' +
                                         SNP_pat_score['Genomic'].str.split('>', expand=True)[0].str[-1], inplace=True)
        col_index += 1

    # Uncomment to  exclude rare SNPs


    SNP_pat_score_GR = SNP_pat_score_GR.drop(list(list_of_rare_Snps.iloc[:,0]), axis=1)
    ### SNP_pat_score_GR --> dataframe with SNP information per patient
    ## keep patient ID, Gender and Phenotype inforamtion
    PED_File = patient_data[['sample_ID', 'Sex', 'Clinical_Score_']]

    ## 1:woman 2:man
    PED_File['Sex'] = PED_File['Sex'].replace(['f'], '1')
    PED_File['Sex'] = PED_File['Sex'].replace(['m'], '2')
    ## 1:low clinical score 2:high clinical score
    PED_File['Clinical_Score_'] = PED_File['Clinical_Score_'].replace([1], 2)
    PED_File['Clinical_Score_'] = PED_File['Clinical_Score_'].replace([0], 1)
    ## rename columns
    PED_File.columns = ['Individual_ID', 'Sex', 'Phenotype']
    # PED_File['Family_ID'] = PED_File.index +1
    PED_File['Family_ID'] = PED_File['Individual_ID'].values
    PED_File['Paternal_ID'] = "0"
    PED_File['Maternal_ID'] = "0"
    ## Reorder
    PED_File = PED_File[['Family_ID', 'Individual_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phenotype']]
    # PED_File.index += 1
    PED_File.index = PED_File['Individual_ID'].values

    ##### Add SNP information
    PED_FILE = pd.concat([PED_File, SNP_pat_score_GR], axis=1)
    ## Set all the missing values to 0 according to PLINK requirements
    PED_FILE.fillna(0, inplace=True)

    ## all haploids must be homozygous
    PED_FILE[PED_FILE.columns[pd.Series(PED_FILE.columns).str.startswith('X')]] = 'C C'
    PED_FILE[PED_FILE.columns[pd.Series(PED_FILE.columns).str.startswith('Y')]] = 'C C'

    path_to_save_PED = p.Processed_Results + "covid_data_fin_PED.txt"
    PED_FILE.to_csv(path_to_save_PED, header=None,
                    index=None, sep='\t', mode='a')

    ## save ped file
    os.system('mv {0} {1}'.format(path_to_save_PED, p.Processed_Results + 'covid_data_fin.ped'))

    print("PED file-->" + p.Processed_Results + "covid_data_fin.ped")

    return SNP_pat_score_GR


def generate_map_file(SNP_pat_score, list_of_rare_Snps, SNP_pat_score_GR):
    ### Keep Chromosome name and Position of SNP
    summary_of_SNPs = SNP_pat_score[['Chromosome', 'Known', 'Position', 'SNP']]
    summary_of_SNPs = summary_of_SNPs[~summary_of_SNPs.SNP.duplicated()]
    summary_of_SNPs = summary_of_SNPs.set_index('SNP')

    ### Keep only unique SNPs
    summary_of_SNPs.fillna(0, inplace=True)
    summary_of_SNPs = summary_of_SNPs.drop(list(list_of_rare_Snps.iloc[:,0]), axis='index')

    summary_of_SNPs = summary_of_SNPs.reindex(list(SNP_pat_score_GR.columns))
    summary_of_SNPs['Known'] = summary_of_SNPs.index

    path_to_save_MAP = p.Processed_Results + "covid_data_fin_MAP.txt"
    summary_of_SNPs.to_csv(path_to_save_MAP, header=None,
                    index=None, sep='\t', mode='a')


    os.system('mv {0} {1}'.format(path_to_save_MAP, p.Processed_Results + 'covid_data_fin.map'))

    print("MAP file-->" + p.Processed_Results + "covid_data_fin.map")

def generate_cov_file(patient_data):
    Cov_file = patient_data[['sample_ID', 'Age']]
    Cov_file['Family_ID'] = Cov_file['sample_ID']
    Cov_file = Cov_file[['Family_ID', 'sample_ID', 'Age']]
    Cov_file.fillna('-9', inplace=True)
    Cov_file.columns = ['Family_ID', 'Individual_ID', 'AGE']

    path_to_save_Cov = p.Processed_Results + "cov_age_fin.txt"
    Cov_file.to_csv(path_to_save_Cov, header=None,
                    index=None, sep='\t')
    print("AGE covariate file -->" + path_to_save_Cov)


def PLINK_logistic_regression():
    ## Create bed file
    os.system('./plink --file {0} --allow-no-sex --out {1} --make-bed'.format(p.Processed_Results + "covid_data_fin", p.Processed_Results + "covid_data_fin" ))

    ### Run logistic regression
    os.system('./plink --bfile {0} --sex --allow-no-sex --logistic --covar {1} --out {2}'.format(
        p.Processed_Results + "covid_data_fin", p.Processed_Results + "cov_age_fin.txt", p.Model_Results + "log_regr_plink"
    ))


def PLINK_transform_files(SNP_pat_score):
    import re
    ## TRANSFORM P-LINK output files to txt
    ### Read plink file
    fin = open(p.Model_Results + "log_regr_plink.assoc.logistic", "rt")
    fout = open(p.Report_Results + "log_reg_results.assoc.logistic.txt", "wt")

    for line in fin:
        newline = re.sub('\s+', '\t', line)
        fout.write(newline)
        fout.write('\n')

    fin.close()
    fout.close()
    data = pd.read_csv(p.Report_Results + "log_reg_results.assoc.logistic.txt",
        delimiter="\t", index_col=None)
    data = data[data.columns[1:-1]]
    ## Keep only ADD columns
    data = data[data.TEST == 'ADD']
    data = data.sort_values(by=['P'])
    data = data.dropna()
    data.to_csv(p.Report_Results + 'log_regr_plink.csv')

    ## ADD gene information
    SNP_list_info = SNP_pat_score.drop_duplicates('SNP')
    unique_snps = list(data.SNP)
    SNP_list_info = SNP_list_info[SNP_list_info['SNP'].isin(unique_snps)]
    gene_to_add = []
    for i in range(data.shape[0]):
        gene_to_add.append(SNP_list_info['Gene'][SNP_list_info.SNP == data.SNP.values[i]].values[0])
    data['Gene'] = gene_to_add

    ##  Keep only SNP Gene and p-value column
    data = data[['SNP', 'P', 'Gene']]
    #data_sort = data.sort_values(by=['P'], ascending=True, ignore_index=True)
    data_sort = data.sort_values(by=['P'], ascending=True)


    ## calculate FDR
    rejected, corrected_val = fdrcorrection(data_sort.P)
    print(data_sort.shape)

    ## Uncomment to add multiple testing correction
    # data_sort['FDR'] = corrected_val

    data_sort.to_csv(p.Report_Results + "log_Regression_Genes.csv", index=None)
    print("Log regression results -->" + p.Report_Results + "log_Regression_Genes.csv")

