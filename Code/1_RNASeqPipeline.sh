###
 # @Descripttion: RNASeq Pipeline
 # @Author: Qun Li
 # @Email: qun.li@ki.se
 # @Date: 2023-11-23 11:25:54
 # @LastEditTime: 2023-12-19 13:30:36
### 

#--------------------------------------------------------------------------------
### Before start
#### you need get reference, annotation files ready.
#### Also, you can download these files from illumina website
#### https://support.illumina.com/sequencing/sequencing_software/igenome.html
#### You can follow these instructions to build index files
#### https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
#### https://daehwankimlab.github.io/hisat2/manual/
#### https://subread.sourceforge.net/featureCounts.html
#### https://htseq.readthedocs.io/en/release_0.11.1/count.html
#### https://ccb.jhu.edu/software/stringtie/
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
#### software: Hisat2 
#### https://daehwankimlab.github.io/hisat2/manual/

hisat2 -p $THREADS -x $GENOMEINDEX -1 $FQ1 -2 $FQ2 \
    --summary-file summary.txt -S FQ.sam

samtools view -@ $THREADS -bS FQ.sam > FQ.bam;
samtools sort -@ $THREADS FQ.bam -o FQ.sorted.bam;
samtools index -@ $THREADS FQ.sorted.bam
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
### 3. Quantify gene expression
#### Software: Featurecounts, HTseq and stringtie
#### https://subread.sourceforge.net/featureCounts.html
#### https://htseq.readthedocs.io/en/release_0.11.1/count.html
#### https://ccb.jhu.edu/software/stringtie/
featureCounts -T $THREADS -a $GTF -o counts.txt --countReadPairs -p -t exon -g gene_id *.sorted.bam
htseq-count -f bam -s no -t exon -i gene_id --nonunique=none *.sorted.bam $GTF > htseq.count
stringtie -e -B -p $THREADS -G $GTF -A FQ_gene_abund.tab -o FQ.gtf FQ.sorted.bam
#--------------------------------------------------------------------------------
