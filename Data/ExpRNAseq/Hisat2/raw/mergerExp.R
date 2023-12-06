#Get list of all files in directory
path = "/Users/liqun/Desktop/Projects/Covid19/Data/Code/SARSCoV2/Data/ExpRNAseq/Hisat2/Total/"
files = list.files(path=path)
paths = paste(path, files, sep = "")
#Merge the first to files and store
coNames_file <- c("GeneID", "GeneName", "Reference", "Strand", "Start", "End", "Coverage")
file1 = read.table(paths[1], col.names=c(coNames_file,paste(files[1],"FPKM"), paste(files[1], "TPM")), sep = "\t", header = T)[,c(1,2,8,9)]
file2 = read.table(paths[2], col.names=c(coNames_file,paste(files[2],"FPKM"), paste(files[2], "TPM")), sep = "\t", header = T)[,c(1,2,8,9)]
out.file = merge (file1, file2, by=c("GeneID", "GeneName"))

#For loop to merge contents of remaining files

for(i in 3:length(paths))
{
  file = read.table(paths[i],col.names=c(coNames_file,paste(files[i],"FPKM"), paste(files[i], "TPM")), sep = "\t", header = T)[,c(1,2,8,9)]
  out.file <- merge(out.file, file, by=c("GeneID", "GeneName"))
}

write.table(out.file, file = "exp_all_sample_covid_total.tsv",sep="\t", row.names = FALSE, quote = FALSE)