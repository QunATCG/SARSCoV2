###
 # @Descripttion: Download RNASeq data
 # @Author: LiQun
 # @Email: qun.li@ki.se
 # @Date: 2023-11-23 11:14:05
 # @LastEditTime: 2023-11-23 11:25:16
### 

## download data using wget and 
## we will rename these files when we fished the RNA-seq pipeline analysis
### Option1
wget -c -i ../Data/sampleUrl.txt

### Option2
while read fileLine
do
    wget -c -b ${fileLine}
done < ../Data/sampleUrl.txt
