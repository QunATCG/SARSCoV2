###
 # @Descripttion: RNASeq Pipeline
 # @Author: Qun Li
 # @Email: qun.li@ki.se
 # @Date: 2023-11-23 11:25:54
 # @LastEditTime: 2023-11-23 14:16:29
### 

#--------------------------------------------------------------------------------
### Before start
#### you need get reference, annotation files ready.
#### Also, you can download these files from illumina website
#### https://support.illumina.com/sequencing/sequencing_software/igenome.html
#### You can follow these instructions to build index files
#### https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
#### https://github.com/alexdobin/STAR
#### https://deweylab.github.io/RSEM/

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### Default setting
#### samples dir path
SAMPLE_DIR=
#### genome index path
GENOMEINDEX=
#### RSEM index path
RSEMINDEX=
#### gtf path
ANNOTATION=$3
#### FASTQ1 path
FQ1=
#### FASTQ2 path
FQ2=
#### RSEM result path
RSEMOUT=
#### Threads number
THREADS=
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 1. QC and Trim
#### software: trim_galore
#### https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

trim_galore -q 30 -stringency 3 -length 20 --phred33 -fastqc \
            -clip_R1 3 -clip_R2 3 -e 0.1 -dont_gzip \
            -paired $SAMPLE_DIR/$FQ1 $SAMPLE_DIR/$FQ2 \
            --output_dir $TRIM_DIR
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 2. Mapping
#### software: STAR 
#### https://github.com/alexdobin/STAR

STAR --outSAMtype BAM SortedByCoordinate --runThreadN $THREADS \
        --readFilesCommand zcat --genomeDir $GENOMEINDEX \
        --readFilesIn $FQ1 $FQ2 --quantMode TranscriptomeSAM GeneCounts \
        --outFileNamePrefix $STAROUT/${FQ1}_${FQ2}/${FQ1}_${FQ2}
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 3. Quantify gene expression
#### Software: RSEM
#### https://deweylab.github.io/RSEM/

rsem-calculate-expression --paired-end -no-bam-output --alignments \
        -p $THREADS ${FQ1}_${FQ2}Aligned.toTranscriptome.out.bam \
        ${RSEMINDEX} $RSEMOUT/${FQ1}_${FQ2}
#--------------------------------------------------------------------------------

