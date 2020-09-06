import pandas as pd
import numpy as np
import math

def get_L_dependent_dist_params(df, A, l):
    """ Calculate the parameters of the distribution A|L=l.
        L can get 0 or 1
        The event 'A' can be a complex trait (C) or a gene expression (R)
    """
    obs = df[df.L == l][A]
    mu = obs.mean()
    var = obs.var(ddof=0)
    
    df[A + '|L_mean'][df.L == l] = mu
    df[A + '|L_var'][df.L == l] = var
    return df

def get_columns_params(col_a, col_b):
    # column A parameters
    mu_a = col_a.mean()
    var_a = col_a.var(ddof=0)
    # column B parameters
    mu_b = col_b.mean()
    var_b = col_b.var(ddof=0)
    # Correlation coef. between A and B
    corr = col_a.corr(col_b)
    return mu_a, var_a, mu_b, var_b, corr


def get_conditional_dist_params(mu_a, var_a, mu_b, var_b, corr, b_i):
    """ Calculate the parameters of the conditional distribution A|B (in our case it could be C|R or R|C)  
        The calculation is a bit different compared to the conditional probability given L
    """
    std_a = np.sqrt(var_a)
    std_b = np.sqrt(var_b)
    mu = mu_a + corr*std_a/std_b * (b_i - mu_b)
    var = var_a * (1-np.square(corr))
    return mu, var  

def norm_pdf(x, mean, var):
    denom = np.sqrt(2 * math.pi * var)
    num = math.exp(-np.square(x-mean)/(2*var))
    return (num/denom)

def init_models(input_table):
    m1_df = pd.DataFrame(columns=['Individual', 'L', 'R', 'C', 
                                   'R|L_mean', 'R|L_var', 
                                   'C|R_mean', 'C|R_var',
                                   'P(R|L)', 'P(C|R)',
                                   'P(L) * P(C|R) * P(R|L)'])
    m1_df[['Individual', 'L', 'R', 'C']] = input_table

    m2_df = pd.DataFrame(columns=['Individual', 'L', 'R', 'C', 
                                   'C|L_mean', 'C|L_var', 
                                   'R|C_mean', 'R|C_var',
                                   'P(C|L)', 'P(R|C)',
                                   'P(L) * P(R|C) * P(C|L)'])
    m2_df[['Individual', 'L', 'R', 'C']] = input_table

    m3_df = pd.DataFrame(columns=['Individual', 'L', 'R', 'C', 
                                   'R|L_mean', 'R|L_var', 
                                   'C|L_mean', 'C|L_var',
                                   'P(R|L)', 'P(C|L)',
                                   'P(L) * P(R|L) * P(C|L)'])
    m3_df[['Individual', 'L', 'R', 'C']] = input_table
    
    return m1_df, m2_df, m3_df

	
def generate_models_df(input_table):
    # Three models tables init
    m1_df, m2_df, m3_df = init_models(input_table)
    
    # Get and set the parameters of R|L distribution in models M1 and M3
    m1_df = get_L_dependent_dist_params(m1_df, 'R', l=0)
    m1_df = get_L_dependent_dist_params(m1_df, 'R', l=1)
    m3_df = get_L_dependent_dist_params(m3_df, 'R', l=0)
    m3_df = get_L_dependent_dist_params(m3_df, 'R', l=1)

    # Get and set the parameters of C|L distribution models M2 and M3
    m2_df = get_L_dependent_dist_params(m2_df, 'C', l=0)
    m2_df = get_L_dependent_dist_params(m2_df, 'C', l=1)
    m3_df = get_L_dependent_dist_params(m3_df, 'C', l=0)
    m3_df = get_L_dependent_dist_params(m3_df, 'C', l=1)

    # Calculate C|R and R|C parameters
    mu_r, var_r, mu_c, var_c, corr_r_c = get_columns_params(input_table['R'], input_table['C'])
    n = len(m1_df)
    for i in range(0, n): 
        # C|R
        r_i = m1_df['R'].iloc[i]
        mu, var = get_conditional_dist_params(mu_c, var_c, mu_r, var_r, corr_r_c, r_i)
        m1_df.iloc[i, m1_df.columns.get_loc('C|R_mean')] = mu
        m1_df.iloc[i, m1_df.columns.get_loc('C|R_var')] = var
        # R|C
        c_i = m2_df['C'].iloc[i]
        mu, var = get_conditional_dist_params(mu_r, var_r, mu_c, var_c, corr_r_c, c_i) 
        m2_df.iloc[i, m2_df.columns.get_loc('R|C_mean')] = mu
        m2_df.iloc[i, m2_df.columns.get_loc('R|C_var')] = var

    # Calculate normal dist. values for all models
    m1_df['P(R|L)'] = m1_df.apply(lambda x: norm_pdf(x['R'], x['R|L_mean'], x['R|L_var']), axis=1)
    m1_df['P(C|R)'] = m1_df.apply(lambda x: norm_pdf(x['C'], x['C|R_mean'], x['C|R_var']), axis=1)

    m2_df['P(C|L)'] = m2_df.apply(lambda x: norm_pdf(x['C'], x['C|L_mean'], x['C|L_var']), axis=1)
    m2_df['P(R|C)'] = m2_df.apply(lambda x: norm_pdf(x['R'], x['R|C_mean'], x['R|C_var']), axis=1)

    m3_df['P(R|L)'] = m3_df.apply(lambda x: norm_pdf(x['R'], x['R|L_mean'], x['R|L_var']), axis=1)
    m3_df['P(C|L)'] = m3_df.apply(lambda x: norm_pdf(x['C'], x['C|L_mean'], x['C|L_var']), axis=1)

    m1_df['P(L) * P(C|R) * P(R|L)'] = 0.5 * m1_df['P(C|R)'] * m1_df['P(R|L)']
    m2_df['P(L) * P(R|C) * P(C|L)'] = 0.5 * m2_df['P(R|C)'] * m2_df['P(C|L)']
    m3_df['P(L) * P(R|L) * P(C|L)'] = 0.5 * m3_df['P(C|L)'] * m3_df['P(R|L)']
    
    # Save as excel file
    m1_df.to_excel("output/likelihood_formulas_M1.xlsx")
    m2_df.to_excel("output/likelihood_formulas_M2.xlsx")
    m3_df.to_excel("output/likelihood_formulas_M3.xlsx")
    
    return m1_df, m2_df, m3_df

def get_LR(m1_df, m2_df, m3_df):
    L_m1 = m1_df['P(L) * P(C|R) * P(R|L)'].prod()
    L_m2 = m2_df['P(L) * P(R|C) * P(C|L)'].prod()
    L_m3 = m3_df['P(L) * P(R|L) * P(C|L)'].prod()
    likelihood_arr = [L_m1, L_m2, L_m3]

    # Get max likelihoos
    max_model = np.argmax([L_m1, L_m2, L_m3])
    best_model = max_model+1
    others_models = [x for i,x in enumerate(range(0,3)) if i != max_model] 

    # Calculate the likelihood ratio
    LR = likelihood_arr[max_model] / max(likelihood_arr[others_models[0]], likelihood_arr[others_models[1]]) 
    return LR, best_model


def get_shuffled_df(df):
    """ Returns the given dataframe with 2 columns randomly shuffeled (L, R).
        No need to randomize the 3rd column """
    df['C'] = df['C'].sample(frac = 1).values
    df['R'] = df['R'].sample(frac = 1).values
    return df