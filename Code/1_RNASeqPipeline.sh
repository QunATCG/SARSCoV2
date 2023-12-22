###
 # @Descripttion: RNASeq Pipeline
 # @Author: Qun Li
 # @Email: qun.li@ki.se
 # @Date: 2023-11-23 11:25:54
 # @LastEditTime: 2023-12-22 14:26:52
### 

#--------------------------------------------------------------------------------
### Before you start
#### you need get reference, annotation files ready.
#### Also, you can download these files from illumina website
#### https://support.illumina.com/sequencing/sequencing_software/igenome.html
#### You can follow these instructions to build your own RNAseq pipeline
#### https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
#### https://github.com/MultiQC/MultiQC
#### https://github.com/alexdobin/STAR
#### https://daehwankimlab.github.io/hisat2/manual/
#### https://github.com/deweylab/RSEM
#### https://subread.sourceforge.net/featureCounts.html
#--------------------------------------------------------------------------------


### Analysis pipeline
### <STAR - RSEM>
### this pipeline was used by some groups, such as
### https://pubmed.ncbi.nlm.nih.gov/38097578/
### https://pubmed.ncbi.nlm.nih.gov/38092806/
### https://pubmed.ncbi.nlm.nih.gov/33870146/
### https://pubmed.ncbi.nlm.nih.gov/34528097/
### https://pubmed.ncbi.nlm.nih.gov/38049398/
#--------------------------------------------------------------------------------
### Default setting
#### samples dir path
SAMPLE_DIR=path_samples
#### genome index path
GENOMEINDEX=path_genome
#### RSEM index
RSEMINDEX=path_RSEM_Index
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
#### https://github.com/MultiQC/MultiQC
trim_galore -q 30 -stringency 3 -length 20 --phred33 -fastqc \
            -clip_R1 3 -clip_R2 3 -e 0.1 -dont_gzip \
            -paired ${FQ1}_${FQ2}.rRNA.dep.fastq.1.gz ${FQ1}_${FQ2}.rRNA.dep.fastq.2.gz \
            --output_dir $TRIM_DIR

multiqc $TRIM_DIR
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 2. Mapping
#### software: STAR
#### https://github.com/alexdobin/STAR
STAR --outSAMtype BAM SortedByCoordinate --runThreadN $THREADS \
        --readFilesCommand zcat --genomeDir $GENOMEINDEX \
        --readFilesIn ${FQ1}_${FQ2}.rRNA.dep.fastq.1.gz_val_1.fq.gz \
        ${FQ1}_${FQ2}.rRNA.dep.fastq.2.gz_val_2.fq.gz \
        --quantMode TranscriptomeSAM GeneCounts \
        --outFileNamePrefix ${FQ1}_${FQ2}
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 3. Quantify gene expression
#### Software: RSEM featureCounts
#### https://subread.sourceforge.net/featureCounts.html
#### https://htseq.readthedocs.io/en/release_0.11.1/count.html
#### https://ccb.jhu.edu/software/stringtie/
rsem-calculate-expression --paired-end -no-bam-output --alignments \
        -p $THREADS ${FQ1}_${FQ2}_Aligned.toTranscriptome.out.bam \
        ${RSEMINDEX} ${FQ1}_${FQ2}


featureCounts -T $THREADS -a $GTF -o ${FQ1}_${FQ2}.counts.txt --countReadPairs -p \ 
    -t exon -g gene_id ${FQ1}_${FQ2}_Aligned.sortedByCoord.out.bam
#--------------------------------------------------------------------------------
