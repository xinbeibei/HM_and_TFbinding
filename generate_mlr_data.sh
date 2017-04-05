#!/bin/sh

folder=$1
tf=$2
tf_length=$3

bedtools intersect -wa -wb -a $folder/$tf/BSs_regions_for_scanning2.bed -b $folder/$tf/non-BS_10_HM.txt > $folder/$tf/temp.txt

awk '{print $4"\t"0}' $folder/$tf/temp.txt > $folder/$tf/seq_neg.txt
awk '{out=$8; for(i=9; i<=NF; i++){out=out"\t"$i}; print 0"\t"1"\t"out}' $folder/$tf/temp.txt > $folder/$tf/HM_neg.txt

bedtools intersect -wa -wb -a $folder/$tf/BSs_regions_for_scanning1.bed -b $folder/$tf/BS_10_HM.txt > $folder/$tf/temp.txt

awk '{print $4"\t"1}' $folder/$tf/temp.txt > $folder/$tf/seq_pos.txt
awk '{out=$8; for(i=9; i<=NF; i++){out=out"\t"$i}; print 1"\t"1"\t"out}' $folder/$tf/temp.txt > $folder/$tf/HM_pos.txt

awk '{print $_}' $folder/$tf/seq_neg.txt > $folder/$tf/$tf\.txt
awk '{print $_}' $folder/$tf/seq_pos.txt >> $folder/$tf/$tf\.txt

awk '{print $_}' $folder/$tf/HM_neg.txt > $folder/$tf/$tf\.00000000001
awk '{print $_}' $folder/$tf/HM_pos.txt >> $folder/$tf/$tf\.00000000001

rm $folder/$tf/seq_neg.txt $folder/$tf/seq_pos.txt $folder/$tf/temp.txt $folder/$tf/HM_pos.txt $folder/$tf/HM_neg.txt

#use DNAshapeR to generate $tf.fa.MGW $tf.fa.Roll $tf.fa.ProT $tf.fa.HelT
#create data for DNA sequecne model, DNA shape models, make sure DNA shape data are normalized. 
Rscript $folder/encode_custom.R $folder/$tf/$tf\.txt $folder/$tf/encoded $tf_length $folder/feature_list

awk '{out=$3; for(i=4; i<=NF; i++){out=out"\t"$i}; print out}' $folder/$tf/$tf\.00000000001 | paste $folder/$tf/$tf\.10011110000 - > $folder/$tf/$tf\.11000000000

awk '{out=$3; for(i=4; i<=NF; i++){out=out"\t"$i}; print out}' $folder/$tf/$tf\.00000000001 | paste $folder/$tf/$tf\.10000000000 - > $folder/$tf/$tf\.10000000001
