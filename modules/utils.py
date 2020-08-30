def print_dim(n,m, title=""):
    print(title + "\n* Expression matrix size: %d\n* Genotype matrix size: %d" % (n, m))
    print("* Expected num of tests: {:,}\n".format(n*m))