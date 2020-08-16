import pandas as pd

def data_annotations_merge(df_exp, df_annotations):
    """ Merge gene expression dataframe with annotations file """
    
    # Remove metadata rows to get the raw data only
    df_exp = df_exp[~df_exp['!Sample_title'].str.contains('!Sample_', na=False)]
    df_exp = df_exp[~df_exp['!Sample_title'].str.contains('!series_', na=False)]
    df_exp = df_exp[df_exp['!Sample_title'] != 'ID_REF']

    # Rename ID column to match the annotation matrix
    df_exp = df_exp.rename(columns = {'!Sample_title' : 'ID'})

    # Merge with annotation matrix to get the gene identifier (gene symbol)
    df_annotations = df_annotations[['ID', 'GENE_SYMBOL']]
    input_matrix = df_annotations.merge(df_exp, left_on='ID', right_on='ID')
    
    return input_matrix  