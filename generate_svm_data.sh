#!/bin/sh

tf=$1
tf_length=$2
folder=$3

for a in 1000
do
#mkdir $folder/encoded_$a
bedtools intersect -wa -wb -a $folder/$tf/BSs_regions_for_scanning2.bed -b $folder/$tf/$tf\_h3k9me3_HM_features_neg_$a.txt > $folder/$tf/temp.txt

awk '{print $4"\t"0}' $folder/$tf/temp.txt > $folder/$tf/seq_neg.txt
awk '{out=$8; for(i=9; i<=NF; i++){out=out"\t"$i}; print 0"\t"1"\t"out}' $folder/$tf/temp.txt > $folder/$tf/HM_neg.txt

bedtools intersect -wa -wb -a $folder/$tf/BSs_regions_for_scanning1.bed -b $folder/$tf/$tf\_h3k9me3_HM_features_pos_$a.txt > $folder/$tf/temp.txt

awk '{print $4"\t"1}' $folder/$tf/temp.txt > $folder/$tf/seq_pos.txt
awk '{out=$8; for(i=9; i<=NF; i++){out=out"\t"$i}; print 1"\t"1"\t"out}' $folder/$tf/temp.txt > $folder/$tf/HM_pos.txt

awk '{print $_}' $folder/$tf/seq_neg.txt > $folder/encoded_$a/$tf\.txt
awk '{print $_}' $folder/$tf/seq_pos.txt >> $folder/encoded_$a/$tf\.txt

awk '{print $_}' $folder/$tf/HM_neg.txt > $folder/encoded_$a/$tf\.00000000001
awk '{print $_}' $folder/$tf/HM_pos.txt >> $folder/encoded_$a/$tf\.00000000001

done
#R --no-restore --no-save --args $folder/$tf/$tf\.txt $folder/encoded_$a/ $tf_length $folder/feature_list < $folder/encode_custom.R 2>&1 >/dev/null

#awk '{out=$3; for(i=4; i<=NF; i++){out=out"\t"$i}; print out}' $folder/encoded_$a/$tf\.00000000001 | paste $folder/encoded_$a/$tf\.10011110000 - > $folder/encoded_$a/$tf\.11000000000


#data for training is : (1)HM features only: whole_HM.txt (2) seq / shape/ seq+shape features: whole_seq.00011110000, whole_seq.10000000000, whole_seq.10011110000
#compare the performance of 4 shape features. SS
