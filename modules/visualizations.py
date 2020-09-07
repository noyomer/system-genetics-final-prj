import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

""" ---------------------------------------------
Visualizations methods for eQTL analysis (part 3)
--------------------------------------------- """

def plot_associated_genes_per_eqtl(assoc_eqtl, genotype_df, tissue_str=""):
	# Count genes per locus
	gene_cnt_df = assoc_eqtl.groupby(["SNP"])["gene"].count().reset_index(name="count")
	gene_cnt_df = gene_cnt_df.rename(columns = {"SNP": "Locus", "count": "num_genes"})

	# Summary of the data to plot
	df = genotype_df[['Locus', 'Chr_Build37']].copy()
	df = df.merge(gene_cnt_df, left_on='Locus', right_on='Locus', how='left')
	df['num_genes'] = df['num_genes'].fillna(value=0) # no associations = 0
	df['index'] = df.index
	df = df.rename(columns={'Chr_Build37': 'chromosome'})
	grp_chr = df.groupby(('chromosome'))

	# plot
	fig = plt.figure()
	ax = fig.add_subplot(111)
	colors = sns.color_palette('coolwarm', n_colors=20)  # a list of RGB tuples
	labels = []
	labels_pos = []

	for num, (chr_num, chr_grp) in enumerate(grp_chr):
		chr_grp.plot(x='index', y='num_genes', color=colors[num % len(colors)],
					 linewidth=3, ax=ax, figsize=(16,5), legend=None)
		labels.append(chr_num)
		curr_pos = (chr_grp['index'].iloc[-1] - (chr_grp['index'].iloc[-1] - chr_grp['index'].iloc[0])/2)
		labels_pos.append(curr_pos)

	ax.set_title(tissue_str + 'Number of genes associated with each eQTL', fontsize = 13)
	ax.set_xlabel('eQTLs Location across the genome', fontsize = 12)
	ax.set_ylabel('Number of genes', fontsize = 12)
	ax.set_xlim([0, len(df)])
	ax.set_xticks(labels_pos)
	ax.set_xticklabels(labels) # replace real index position with chromosome numbers
	sns.set_style("whitegrid", {'axes.grid' : False})
	plt.show()
	return df

	
def plot_association_pval_dist(assoc_eqtl, tissue_str=""):
	# get cis and trans distributrions 
	cis_pval = list(assoc_eqtl[assoc_eqtl["closeness"] == 'cis']["minus_log_p-value"])
	trans_pval = list(assoc_eqtl[assoc_eqtl["closeness"] == 'trans']["minus_log_p-value"])
	# plot distribution
	sns.kdeplot(cis_pval, shade=True, label='cis-acting')
	sns.kdeplot(trans_pval, shade=True, label='trans-acting')
	plt.title(tissue_str + 'Distribution of association P-value scans', fontsize=14)
	plt.xlabel('-log p-value', fontsize=13)
	axes = plt.axes()
	sns.set_style("whitegrid")
	plt.show()

def plot_summary_visualization(reg_results, mgi_df, treshold=0.05, tissue_str=""):
	# Drop mgi rows without chromosome data
	mgi_df = mgi_df[~mgi_df['representative genome chromosome'].isna()]

	# add positions data from mgi file
	df = reg_results
	df = df.rename(columns={'chromosome': 'snp_chromosome', 'position': 'snp_position'})
	df['gene_chromosome'] = np.nan
	df['gene_position'] = np.nan
	gene_list = list(df.gene.unique())
	mgi_genes = list(mgi_df['marker symbol'])
	for gene in gene_list:
		if gene not in mgi_genes:
			df['gene_chromosome'] = np.where(df['gene'] == gene, 'Unknown', df['gene_chromosome'])
			df['gene_position'] = np.where(df['gene'] == gene, 'Unknown', df['gene_position'])
			continue
			
		curr_mgi = mgi_df[mgi_df['marker symbol'] == gene]
		chromosome = curr_mgi['representative genome chromosome'].iloc[0]
		chromosome = 20 if chromosome in ['X', 'Y'] else int(chromosome) 
		start = curr_mgi['representative genome start'].iloc[0]
		
		df['gene_chromosome'] = np.where(df['gene'] == gene, int(chromosome), df['gene_chromosome'])
		df['gene_position'] = np.where(df['gene'] == gene,start, df['gene_position'])
	df = df[df['gene_chromosome'] != 'Unknown']

	# define absolute positions
	delta = df[['gene_position', 'snp_position']].astype(float).min().min()
	max_pos = 200000000
	factor = max_pos - delta
	df['snp_absulute_position'] = df['snp_chromosome'].astype(float) + (df['snp_position'].astype(float)-delta)/factor 
	df['gene_absulute_position'] = df['gene_chromosome'].astype(float)  + (df['gene_position'].astype(float)-delta)/factor 	

	# plot summary visualization
	df2 = df[df['p-value'] <= treshold]
	df2 = df2.astype(str)
	df2[['snp_absulute_position', 'gene_absulute_position']] = df2[['snp_absulute_position', 'gene_absulute_position']].astype(float)
	df2 = df2.sort_values(by=['snp_absulute_position', 'gene_absulute_position'])

	df_grouped = df2.groupby(('snp_chromosome'))

	fig = plt.figure()
	ax = fig.add_subplot(111)
	colors = sns.color_palette('GnBu_d', n_colors=20)  # a list of RGB tuples

	x_labels = []
	x_labels_pos = []
	for num, (name, group) in enumerate(df_grouped):
		group.plot(kind="scatter",figsize=(10,10), x='snp_absulute_position', y='gene_absulute_position',
					s=6, color=colors[num], ax=ax, legend=None)
		x_labels.append(name)
		x_labels_pos.append((group['snp_absulute_position'].iloc[-1] - (group['snp_absulute_position'].iloc[-1] - group['snp_absulute_position'].iloc[0])/2))

	plt.title(tissue_str + 'Summary visualization', fontsize=14)

	ax.set_xticks([i for i in range(1,21)])
	ax.set_yticks([i for i in range(1,21)])
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	ax.set_xlim([0.95,21])
	ax.set_ylim([0.95,21])
	ax.set_ylabel('Gene position (by chromosome)', fontsize = 13)
	ax.set_xlabel('eQTL position (by chromosome)', fontsize = 13)

	ax.grid(False)	
	
	
	
""" ---------------------------------------------
Visualizations methods for QTL analysis (part 4)
--------------------------------------------- """
