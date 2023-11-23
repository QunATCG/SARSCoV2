###
 # @Descripttion: RNASeq Pipeline
 # @Author: Qun Li
 # @Email: qun.li@ki.se
 # @Date: 2023-11-23 11:25:54
 # @LastEditTime: 2023-11-23 11:50:43
### 

### 1. QC and Trim
#### software trim_galore
#### https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

trim_galore -q 30 -stringency 3 -length 20 --phred33 -fastqc -clip_R1 3 -clip_R2 3 -e 0.1 -dont_gzip -paired $SAMPLE_DIR/$R1_FASTQ $SAMPLE_DIR/$R2_FASTQ --output_dir $TRIM_DIR