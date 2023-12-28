fileIn = open("../Data/ExpRNAseq/exp_star_removerRNA_exCovid_20231227/Human_Covid19_removerRNA_removeCovid.counts.txt", "r")
fileOut = open("../Data/ExpRNAseq/exp_star_removerRNA_exCovid_20231227/Human_Covid19_removerRNA_removeCovid.counts", "w+")
for line in fileIn:
    linesets = line.split("\t")
    sql = linesets[0] + "\t" + linesets[6] + "\t" + linesets[7] + "\t" + linesets[8] + "\t" + linesets[9] + "\t" + linesets[10] + "\t" + linesets[11] + "\t" + linesets[12] + "\t" + linesets[13] + "\t" + linesets[14].strip()
    fileOut.writelines(sql + "\n")