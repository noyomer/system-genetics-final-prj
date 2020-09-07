import pandas as pd
import numpy as np
from tqdm import tqdm

def get_associations(reg_results, alpha=0.05):
    num_tests = len(reg_results)
    alpha_star = alpha / num_tests
    reg_results = reg_results[reg_results['p-value'] <= alpha_star]
    return reg_results, num_tests, alpha_star

def cis_trans_annotation(reg_results, mgi_df, delta=2000000):
    """ Annotations of genes to cis-acting or trans-acting genes, with respect to a distnace treshold.
        Default distnace treshold(=delta) is 2Mbp.
        Gebers that are not avaiable in the mgi file would be annotates as 'Unknown'.        
    """
    
    # Drop mgi rows without chromosome data
    mgi_df = mgi_df[~mgi_df['representative genome chromosome'].isna()]

    # Cis and Trans annotation
    reg_results['closeness'] = 'trans'
    gene_list = list(reg_results.gene.unique())
    mgi_genes = list(mgi_df['marker symbol'])

    for gene in tqdm(gene_list):
        if gene not in mgi_genes:
            reg_results['closeness'] = np.where(reg_results['gene'] == gene, 'Unknown', reg_results['closeness'])
            continue

        curr_mgi = mgi_df[mgi_df['marker symbol'] == gene]
        chromosome = curr_mgi['representative genome chromosome'].iloc[0]
        chromosome = 20 if chromosome in ['X', 'Y'] else int(chromosome) 
        start = curr_mgi['representative genome start'].iloc[0]
        end = curr_mgi['representative genome end'].iloc[0]
        cis_conditions = (reg_results['gene'] == gene) & \
                        (reg_results['chromosome'] == chromosome) & \
                        (reg_results["position"] >= start - delta) & \
                        (reg_results["position"] <= end + delta)
        reg_results['closeness'] = np.where(cis_conditions, 'cis', reg_results['closeness'])
    
    return reg_results