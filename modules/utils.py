import pandas as pd

def print_dim(n,m, title=""):
    print(title + "\n* Expression matrix size: %d\n* Genotype matrix size: %d" % (n, m))
    print("* Expected num of tests: {:,}\n".format(n*m))
	
def print_stats(assoc_eqtl, num_tests, title=""):
    bold_s = '\033[1m'
    bold_e = '\033[0m'
    total = len(assoc_eqtl)
    cis = assoc_eqtl[assoc_eqtl['closeness'] == 'cis']
    trans = assoc_eqtl[assoc_eqtl['closeness'] == 'trans']
    print(bold_s + title + bold_e)
    print(bold_s + "Number of tests: " + bold_e, num_tests)
    print("\n" + bold_s + "Number of different significant eQTLs: " + bold_e, len(assoc_eqtl.SNP.unique()))
    print("From which: \n %d - cis-acting \n %d - trans-acting " % (len(cis.SNP.unique()), len(trans.SNP.unique())))
    print("\n" + bold_s + "Number of total significant eQTLs: " + bold_e, len(assoc_eqtl))
    print("From which: \n %d - cis-acting \n %d - trans-acting " % (len(cis), len(trans)))