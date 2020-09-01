import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import f

def qtl_analysis(phenotype_df, genotype_df, file_prefix=""):
    p_data = phenotype_df.drop(columns=['Phenotype', 'Authors', 'Year', 'Pubmed Id'])
    num_phenotypes = len(phenotype_df)
    reg_results = pd.DataFrame(columns=['SNP', 'chromosome', 'position', 'gene', 'p-value'])

    gene_locus = genotype_df.iloc[:, 0]
    genotypes = genotype_df.drop(columns = ['Locus', 'Chr_Build37', 'Build37_position'])
    
    # Remove heterozygous markers
    genotypes = genotypes.replace({'B': 0, 'b': 0, 'D':2, 'H': np.nan, 'U': np.nan})
     
    for i in tqdm(range(0, num_phenotypes)):
        phenotype = phenotype_df.iloc[i]['Phenotype']
        logpval_list, pval_list = run_regression(genotypes, gene_locus, p_data.iloc[[i]])

        # store results
        values = genotype_df[['Locus','Chr_Build37', 'Build37_position']].copy() 
        values.columns = ['SNP', 'chromosome', 'position']
        values['phenotype'] = phenotype
        values['p-value'] = np.array(pval_list)
        values['minus_log_p-value'] = np.array(logpval_list)

        reg_results = pd.concat([reg_results,values])
        
    reg_results.to_csv(file_prefix + "QTL_results.csv")
    return reg_results

	
def eqtl_analysis(exp_df, genotype_df, file_prefix=""):
    e_data = exp_df.drop(columns='GENE_SYMBOL')
    num_genes = len(exp_df)
    reg_results = pd.DataFrame(columns=['SNP', 'chromosome', 'position', 'gene', 'p-value'])

    gene_locus = genotype_df.iloc[:, 0]
    genotypes = genotype_df.drop(columns = ['Locus', 'Chr_Build37', 'Build37_position'])
    
    # Remove heterozygous markers
    genotypes = genotypes.replace({'B': 0, 'b': 0, 'D':2, 'H': np.nan, 'U': np.nan})
     
    for i in tqdm(range(0, num_genes)):
        gene = exp_df.iloc[i]['GENE_SYMBOL']
        logpval_list, pval_list = run_regression(genotypes, gene_locus, e_data.iloc[[i]])

        # store results
        values = genotype_df[['Locus','Chr_Build37', 'Build37_position']].copy() 
        values.columns = ['SNP', 'chromosome', 'position']
        values['gene'] = gene
        values['p-value'] = np.array(pval_list)
        values['minus_log_p-value'] = np.array(logpval_list)

        reg_results = pd.concat([reg_results,values])
        
    reg_results.to_csv(file_prefix + "reg_results.csv")
    return reg_results

def regression_model(x, y):
    """ Given the samplesm, implement and run a regression model in which
        heterozygous markers are ignored (F-test)
    """
    numerator = 0
    denominator = 0
    sample_size = len(x)
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    
    # 1. calculate beta 1, beta 0
    for i in range(sample_size):
        numerator += (x[i] - x_mean) * (y[i] - y_mean)
        denominator += np.square(x[i] - x_mean) 
        
    beta_1 = numerator/denominator
    beta_0 = np.mean(y) - (beta_1 * np.mean(x))
    
    # 2. calculate Sum of Squares
    SSR = 0
    SSE = 0
    for i in range(sample_size):
        y_hat = beta_0 + (beta_1 * x[i])
        SSR += np.square(y_hat - y_mean)
        SSE += np.square(y_hat - y[i])
    
    SST = SSR + SSE
    R_squared = SSR/SST
    
    # 3.calculate F star and p val
    MSE = SSE/(sample_size-2)
    F = SSR/MSE 
    p_val = (1 - f.cdf(F, 1, sample_size-2))
        
    return F, p_val

	
def run_regression(genotypes, gene_locus, phenotype):   
    """ Run regression, calculate return p-values lists """
    nan_cols = phenotype.columns[phenotype.isna().any()].tolist()
    pval_list = []
    logpval_list = []
    
    for i in range(len(genotypes)):
        SNP = genotypes.iloc[[i]]
        col_to_drop = nan_cols + SNP.columns[SNP.isna().any()].tolist()
        SNP = SNP.drop(columns = col_to_drop)
        adj_phenotype = phenotype.drop(columns = col_to_drop)
        x = SNP.values.tolist()[0]
        y = adj_phenotype.values.tolist()[0]
        F, p_val = regression_model(x, y)
        log_pval = -np.log10(p_val)
        logpval_list.append(log_pval)
        pval_list.append(p_val)

    return logpval_list, pval_list