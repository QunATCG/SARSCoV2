fileIn = open("Covid19_gene_readcount.csv","r")
fileOut = open("Covid19_gene_readcount_new.csv","w+")
for line in fileIn:
	lineTag = line[0:18]
	newline = lineTag + "\t" + line.strip().replace(',','\t')
	fileOut.writelines(newline + "\n")