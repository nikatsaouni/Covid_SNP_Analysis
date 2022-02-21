import pandas as pd
import src.cases.parameters as p

def generate_patient_snp_matrix(SNPs):
    genome_all_pat = pd.read_csv(p.Interm_Results + "SNP_data.csv")


    ### Get all SNPs for each patient
    genome_all_pat_snp = genome_all_pat.groupby(['sample_ID'])['SNP'].apply(lambda x: ','.join(x.astype(str))).reset_index()
    genome_all_pat_snp = genome_all_pat_snp.astype(str)
    genome_all_pat_snp = genome_all_pat_snp.set_index('sample_ID')

    ### create sparse representation
    genome_all_pat_Sparse = genome_all_pat_snp.iloc[:, 0].str.replace(' ', '').str.get_dummies(sep=',')
    genome_all_pat_Sparse = genome_all_pat_Sparse[
        genome_all_pat_Sparse.columns[genome_all_pat_Sparse.columns.isin(SNPs)]]

    #### Add zygosity information 0: aa 1: Aa 2: AA
    for i in range(genome_all_pat_Sparse.shape[0]):
        for j in range(genome_all_pat_Sparse.shape[1]):
            if genome_all_pat_Sparse.iloc[i, j] == 1:

                getzygosity = genome_all_pat[(genome_all_pat.sample_ID == int(genome_all_pat_Sparse.index[i]) & (
                            genome_all_pat.SNP == genome_all_pat_Sparse.columns[j]))].Zygosity.values
                if getzygosity == 'Homozygous':
                    genome_all_pat_Sparse.iloc[i, j] = 2

    ## keep only significant SNPs
    # SNP_attr = genome_all_pat_Sparse[genome_all_pat_Sparse.columns[genome_all_pat_Sparse.columns.isin(SNPs)]]
    # SNP_attr.insert(loc=0, column='sample_ID', value=genome_all_pat['sample_ID'])
    # SNP_attr = SNP_attr.set_index('sample_ID')


    genome_all_pat_Sparse.index = genome_all_pat_Sparse.index.astype(int)
    print('Number of significant SNPs', genome_all_pat_Sparse.shape[1])
    genome_all_pat_Sparse = genome_all_pat_Sparse.reindex(columns = SNPs)
    genome_all_pat_Sparse.to_csv(p.Report_Results + 'rs_correlat.csv', index=False,
        header=True)
    print("Patient-SNP matrix saved in: " + p.Report_Results + "pat_SNP.csv")