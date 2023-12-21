###
 # @Descripttion: RNASeq Pipeline
 # @Author: Qun Li
 # @Email: qun.li@ki.se
 # @Date: 2023-11-23 11:25:54
 # @LastEditTime: 2023-12-21 19:35:19
### 

#--------------------------------------------------------------------------------
### Before start
#### you need get reference, annotation files ready.
#### Also, you can download these files from illumina website
#### https://support.illumina.com/sequencing/sequencing_software/igenome.html
#### You can follow these instructions to build your own RNAseq pipeline
#### https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
#### https://github.com/MultiQC/MultiQC
#### https://daehwankimlab.github.io/hisat2/manual/
#### https://subread.sourceforge.net/featureCounts.html
#### https://htseq.readthedocs.io/en/release_0.11.1/count.html
#### https://ccb.jhu.edu/software/stringtie/
#--------------------------------------------------------------------------------


### our analysis pipeline 
#--------------------------------------------------------------------------------
### Default setting
#### samples dir path
SAMPLE_DIR=path_samples
#### genome index path
GENOMEINDEX=path_genome
#### gtf path
GTF=path_GTF
#### FASTQ1 path
FQ1=path_fastq1
#### FASTQ2 path
FQ2=path_fastq2
#### Threads number
THREADS=16
#### TRIMOUT path
TRIM_DIR=path_trimOut
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 0. Remove rRNA
#### Hisat2
#### 0.1 download human rRNA sequence from NCBI. Also, we provide human rRNA fasta in this repository
#### 0.2 build rRNA sequence index (Human_rRNA_Index) using hisat2
#### 0.3 exclude rRNA sequence
hisat2 -p $THREADS -x $Human_rRNA_Index -1 $FQ1 -2 $FQ2 --un-conc-gz ${FQ1}_${FQ2}.rRNA.dep.fastq.gz
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
### 1. QC and Trim
#### software: trim_galore
#### https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
trim_galore -q 30 -stringency 3 -length 20 --phred33 -fastqc \
            -clip_R1 3 -clip_R2 3 -e 0.1 -dont_gzip \
            -paired ${FQ1}_${FQ2}.rRNA.dep.fastq.1.gz ${FQ1}_${FQ2}.rRNA.dep.fastq.2.gz \
            --output_dir $TRIM_DIR
#### https://github.com/MultiQC/MultiQC
multiqc $TRIM_DIR
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 2. Mapping
#### software: Hisat2 
#### https://daehwankimlab.github.io/hisat2/manual/

hisat2 -p $THREADS -x $GENOMEINDEX -1 ${FQ1}_${FQ2}.rRNA.dep.fastq.1.gz_val_1.fq.gz -2 ${FQ1}_${FQ2}.rRNA.dep.fastq.2.gz_val_2.fq.gz \
    --summary-file summary.txt -S ${FQ1}_${FQ2}.sam

samtools view -@ $THREADS -bS ${FQ1}_${FQ2}.sam > ${FQ1}_${FQ2}.bam;
samtools sort -@ $THREADS ${FQ1}_${FQ2}.bam -o ${FQ1}_${FQ2}.sorted.bam;
samtools index -@ $THREADS ${FQ1}_${FQ2}.sorted.bam
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 3. Quantify gene expression
#### Software: Featurecounts, HTseq and stringtie
#### https://subread.sourceforge.net/featureCounts.html
#### https://htseq.readthedocs.io/en/release_0.11.1/count.html
#### https://ccb.jhu.edu/software/stringtie/
featureCounts -T $THREADS -a $GTF -o ${FQ1}_${FQ2}.counts.txt --countReadPairs -p -t exon -g gene_id ${FQ1}_${FQ2}.sorted.bam
htseq-count -f bam -s no -t exon -i gene_id --nonunique=none ${FQ1}_${FQ2}.sorted.bam $GTF > htseq.count
stringtie -e -B -p $THREADS -G $GTF -A ${FQ1}_${FQ2}_gene_abund.tab -o ${FQ1}_${FQ2}.gtf ${FQ1}_${FQ2}.sorted.bam
#--------------------------------------------------------------------------------
