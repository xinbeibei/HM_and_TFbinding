#!/bin/sh

folder=$1   #parent directory of tf
tf=$2
pwm_length=$3
flank_length=$4
IsBpWise=$5   # IsBpWise==1 means no need to calculate average HM level among all bases around BS, this option will be used in generate data for Figure 2C. 

if [ -d "$folder/$tf/fimo_out" ]; then
    rm -r $folder/$tf/fimo_out
fi

#convert ChIP-seq peaks into fasta file, since fimo requires .fa for input
awk '{print $1"\t"$2"\t"$3}' $folder/$tf/$tf\.narrowPeak > $folder/$tf/temp.bed
bedtools getfasta -fi $folder/hg19.fa -bed $folder/$tf/temp.bed -fo $folder/$tf/$tf\_ChIP_Seq.fa
rm $folder/$tf/temp.bed

# threshold for fimo is 1e-4
fimo --o $folder/$tf/fimo_out/ $folder/$tf/$tf\_pwm.txt $folder/$tf/$tf\_ChIP_Seq.fa
rm $folder/$tf/$tf\_ChIP_Seq.fa

# remove head of the file which is not used later
sed '1d' $folder/$tf/fimo_out/fimo.txt > $folder/$tf/temp.txt
mv $folder/$tf/temp.txt $folder/$tf/fimo_out/fimo.txt

#sort fimo.txt accoring to column 2 (chromo_start_end) and 7 (p_vale_of_the_motif_occurence)
sort -k2 -k7 $folder/$tf/fimo_out/fimo.txt > $folder/$tf/test.txt

#unique according to the column 2, only save the seq with minimum p_value
sort -uk2,2 $folder/$tf/test.txt > $folder/$tf/fimo_out/unique_fimo_out.txt
rm $folder/$tf/test.txt

#get different BSs_regions_for_scannning1.bed: chr start end seq(add 2 flanking bps), the way I did it is to change fromFimoToBed_consensus1.py
python $folder/fromFimoToBed_consensus1.py $folder/$tf/fimo_out/ $folder/$tf/   #output is BSs_region_for_scanning2.bed which has 4 cols,BSs_regions_for_scanning1.bed has 2 bps flanking on each side 
rm $folder/$tf/BSs_regions_for_getfasta.bed
rm $folder/$tf/BSs_seq_with_2_flanks.fa
mv $folder/$tf/BSs_regions_from_fasta.bed $folder/$tf/BSs_regions_for_scanning1.bed
rm $folder/$tf/BSs_motif_score.bed

awk '{print $1"\t"$2"\t"$3}' $folder/$tf/BSs_regions_for_scanning1.bed  > $folder/$tf/BS.bed #get bed sequences accessible, however, not bound in vivo. 3 cols. at last, find sequences of each region in BSs_regions_for_scanning2.bed

awk -v var="$flank_length" '{print $1"\t"$2-var"\t"$3+var}' $folder/$tf/BS.bed > $folder/$tf/BSs_regions_for_scanning_pos_$flank_length\.bed

#scan for histone modification read coverage (RPM: reads per million)
for histone in h3k4me2 h3k27ac h3k27me3 h3k36me3 h3k4me1 h3k4me3 h3k79me2 h3k9ac h4k20me1
do
    if [ -f $folder/$tf/$histone\_BS_RPM_pos_$flank_length.txt ]; then
        rm $folder/$tf/$histone\_BS_RPM_pos_$flank_length.txt
    fi
    for rep in rep1 rep2
    do
         python $folder/histone_modification_feature.py $folder/10_histone_modification $folder/$tf $folder/10_histone_modification/$histone\_$rep\.bam $folder/$tf/BSs_regions_for_scanning_pos_$flank_length.bed $folder/10_histone_modification/$histone\_$rep\_count.txt $folder/$tf/$histone\_BS_RPM_pos_$flank_length.txt $pwm_length $flank_length $IsBpWise     
    done
    if [ $IsBpWise == 0 ]
    then
    awk '{a[$1"\t"$2"\t"$3"\t"]+=$4}END{for(i in a){print i,a[i]/2}}' $folder/$tf/$histone\_BS_RPM_pos_$flank_length.txt > $folder/$tf/$histone\_avg_BS_RPM_pos_$flank_length.txt
    rm $folder/$tf/$histone\_BS_RPM_pos_$flank_length.txt
    else
    awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"]+=$5}END{for(i in a){print i,a[i]/2}}' $folder/$tf/$histone\_BS_RPM_pos_$flank_length.txt | awk '{a[$1"\t"$2"\t"$3]=a[$1"\t"$2"\t"$3]"\t"$5}END{for(i in a){print i,a[i]}}' - > $folder/$tf/$histone\_avg_BS_RPM_pos_$flank_length_bp_wise.txt
    rm $folder/$tf/$histone\_BS_RPM_pos_$flank_length.txt
    fi
done

if [ -f $folder/$tf/h3k9me3_BS_RPM_pos_$flank_length.txt ]; then
        rm $folder/$tf/h3k9me3_BS_RPM_pos_$flank_length.txt
fi
for rep in rep1 rep2 rep3
    do
         python $folder/histone_modification_feature.py $folder/10_histone_modification $folder/$tf $folder/10_histone_modification/h3k9me3_$rep\.bam $folder/$tf/BSs_regions_for_scanning_pos_$flank_length.bed $folder/10_histone_modification/h3k9me3_$rep\_count.txt $folder/$tf/h3k9me3_BS_RPM_pos_$flank_length.txt $pwm_length $flank_length $IsBpWise        
done
if [ $IsBpWise == 0 ]
    then
    awk '{a[$1"\t"$2"\t"$3"\t"]+=$4}END{for(i in a){print i,a[i]/3}}' $folder/$tf/h3k9me3_BS_RPM_pos_$flank_length.txt > $folder/$tf/h3k9me3_avg_BS_RPM_pos_$flank_length.txt
    rm $folder/$tf/h3k9me3_BS_RPM_pos_$flank_length.txt
else
    awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"]+=$5}END{for(i in a){print i,a[i]/3}}' $folder/$tf/h3k9me3_BS_RPM_pos_$flank_length.txt | awk '{a[$1"\t"$2"\t"$3]=a[$1"\t"$2"\t"$3]"\t"$5}END{for(i in a){print i,a[i]}}' - > $folder/$tf/h3k9me3_avg_BS_RPM_pos_$flank_length_bp_wise.txt
    rm $folder/$tf/h3k9me3_BS_RPM_pos_$flank_length.txt
fi


if [ $IsBpWise == 0 ]
then
sort -k1,3 $folder/$tf/h3k4me2_avg_BS_RPM_pos_$flank_length.txt > $folder/$tf/$tf\_h3k4me2_HM_features_pos_$flank_length.txt

sort -k1,3 $folder/$tf/h3k27ac_avg_BS_RPM_pos_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k4me2_HM_features_pos_$flank_length.txt - > $folder/$tf/$tf\_h3k27ac_HM_features_pos_$flank_length.txt

sort -k1,3 $folder/$tf/h3k27me3_avg_BS_RPM_pos_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k27ac_HM_features_pos_$flank_length.txt - > $folder/$tf/$tf\_h3k27me3_HM_features_pos_$flank_length.txt

sort -k1,3 $folder/$tf/h3k36me3_avg_BS_RPM_pos_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k27me3_HM_features_pos_$flank_length.txt - > $folder/$tf/$tf\_h3k36me3_HM_features_pos_$flank_length.txt

sort -k1,3 $folder/$tf/h3k4me1_avg_BS_RPM_pos_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k36me3_HM_features_pos_$flank_length.txt - > $folder/$tf/$tf\_h3k4me1_HM_features_pos_$flank_length.txt

sort -k1,3 $folder/$tf/h3k4me3_avg_BS_RPM_pos_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k4me1_HM_features_pos_$flank_length.txt - > $folder/$tf/$tf\_h3k4me3_HM_features_pos_$flank_length.txt

sort -k1,3 $folder/$tf/h3k79me2_avg_BS_RPM_pos_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k4me3_HM_features_pos_$flank_length.txt - > $folder/$tf/$tf\_h3k79me2_HM_features_pos_$flank_length.txt

sort -k1,3 $folder/$tf/h3k9ac_avg_BS_RPM_pos_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k79me2_HM_features_pos_$flank_length.txt - > $folder/$tf/$tf\_h3k9ac_HM_features_pos_$flank_length.txt

sort -k1,3 $folder/$tf/h4k20me1_avg_BS_RPM_pos_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k9ac_HM_features_pos_$flank_length.txt - > $folder/$tf/$tf\_h4k20me1_HM_features_pos_$flank_length.txt

sort -k1,3 $folder/$tf/h3k9me3_avg_BS_RPM_pos_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h4k20me1_HM_features_pos_$flank_length.txt - > $folder/$tf/$tf\_h3k9me3_HM_features_pos_$flank_length.txt

cp $folder/$tf/$tf\_h3k9me3_HM_features_pos_$flank_length.txt $folder/$tf/BS_10_HM.txt
fi

for histone in h3k4me2 h3k27ac h3k27me3 h3k36me3 h3k4me1 h3k4me3 h3k79me2 h3k9ac h4k20me1 h3k9me3
do
    rm $folder/$tf/$tf\_$histone\_HM_features_pos_$flank_length.txt
    rm $folder/$tf/$histone\_avg_BS_RPM_pos_$flank_length.txt
done

rm $folder/$tf/BSs_regions_for_scanning_pos_$flank_length\.bed
