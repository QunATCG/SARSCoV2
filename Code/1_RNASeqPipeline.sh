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
#### https://daehwankimlab.github.io/hisat2/manual/
#### https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
#### https://subread.sourceforge.net/featureCounts.html
#--------------------------------------------------------------------------------


### Analysis pipeline
### <Hisat2 - Stringtie>
#--------------------------------------------------------------------------------
### Default setting
#### samples dir path
SAMPLE_DIR=path_samples
#### Human genome index path
GENOMEINDEX=path_genome
#### Human rRNA genome index path
rRNAINDEX=path_rRNA_index
#### SARS2 genome index path
SARS2INDEX=path_sars_index
#### gtf path
GTF=path_GTF
#### SAMPLE NAME
SAMPLENAME=sampleID
#### FASTQ1 path
FQ1=path_fastq1
#### FASTQ2 path
FQ2=path_fastq2
#### Threads number
THREADS=16
#### TRIMOUT path
TRIM_DIR=path_trimOut
#### mapped out path
MAPPOUTDIR=path_mapped
#### Stringtie
STRINGTIEOUTDIR=path_stringtie
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 0. Trim and QC
#### software: trim_galore
#### https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
#### https://github.com/MultiQC/MultiQC
trim_galore -q 30 -stringency 3 -length 20 --phred33 -fastqc \
	-clip_R1 3 -clip_R2 3 -e 0.1 -dont_gzip \
	-paired ${FQ1} ${FQ2} \
	--output_dir $TRIM_DIR

multiqc $TRIM_DIR
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 1. Remove rRNA
#### Hisat2
#### https://daehwankimlab.github.io/hisat2/manual/
#### 0.1 download human rRNA sequence from NCBI. Also, we provide human rRNA fasta in this repository
#### 0.2 build rRNA sequence index (Human_rRNA_Index) using hisat2
#### 0.3 exclude rRNA sequence
hisat2 -p $THREADS -x $rRNAINDEX -1 $FQ1_trimed -2 $FQ2_trimed --un-conc-gz ${SAMPLE}.rRNA.dep.fastq.gz
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 2. Remove SARS2 reads
#### Hisat2
#### https://daehwankimlab.github.io/hisat2/manual/
#### 1.1 download SARS2 sequence from Ensembl. Also, we provide SARS fasta in this repository
#### 1.2 build SARS2 sequence index (SARS2_Index) using hisat2
#### 1.3 exclude SARS2 sequence
hisat2 -p $THREADS -x $SARS2INDEX -1 ${SAMPLE}.rRNA.dep.fastq.1.gz -2 ${SAMPLE}.rRNA.dep.fastq.2.gz \
	--un-conc-gz ${SAMPLE}.rRNA.sars.dep.fastq.gz
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 3. Mapping
#### 3.1. Mapping reads without rRNA and SARS2
#### software: HISAT2
#### https://daehwankimlab.github.io/hisat2/manual/
hisat2 -p $THREADS -x $GENOMEINDEX -1 $TRIM_DIR/${SAMPLE}.rRNA.sars.dep.fastq.1.gz_val_1.fq.gz \
	-2 $TRIM_DIR/${SAMPLE}.rRNA.sars.dep.fastq.2.gz_val_2.fq.gz \
	--summary-file $MAPPOUTDIR/${SAMPLE}.txt -S $MAPPOUTDIR/${SAMPLE}.sam;

#### 3.2. Mapping reads with SARS2
#### software: HISAT2
#### https://daehwankimlab.github.io/hisat2/manual/
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 4. Quantify gene expression
#### 4.1. Quantify for 3.1 
#### Software: Stringtie featureCounts
#### https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
#### https://subread.sourceforge.net/featureCounts.html
stringtie -e -B -p $THREADS -G $GTF \
	-o ${STRINGTIEOUTDIR}/Sample_${SAMPLE}/${SAMPLE}.gtf \
	$MAPPOUTDIR/${SAMPLE}.sorted.bam

featureCounts -T $THREADS -a $GTF -o ${SAMPLE}.counts.txt --countReadPairs -p \ 
	-t exon -g gene_id $MAPPOUTDIR/${SAMPLE}_Aligned.sortedByCoord.out.bam

#### 4.2. Quantify for 3.2 
#### Software: Stringtie featureCounts
#### https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
#### https://subread.sourceforge.net/featureCounts.html
#--------------------------------------------------------------------------------
