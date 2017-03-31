#!/bin/sh

folder=$1   #parent directory of tf
tf=$2
pwm_length=$3
bed_length=$4
flank_length=$5
IsBpWise=$6   # IsBpWise==1 means no need to calculate average HM level among all bases around BS, this option will be used in generate data for Figure 2C. 

if [ -d "$folder/$tf/fimo_out2" ]; then
    rm -r $folder/$tf/fimo_out2
fi

fimo --o $folder/$tf/fimo_out2/ $folder/$tf/$tf\_pwm.txt $folder/DNase_accessible.fa

sed '1d' $folder/$tf/fimo_out2/fimo.txt > $folder/$tf/temp.txt
mv $folder/$tf/temp.txt $folder/$tf/fimo_out2/fimo.txt

sort -k2 -k7 $folder/$tf/fimo_out2/fimo.txt > $folder/$tf/test.txt
sort -uk2,2 $folder/$tf/test.txt > $folder/$tf/fimo_out2/unique_fimo_out.txt
rm $folder/$tf/test.txt

python $folder/fromFimoToBed_consensus1.py $folder/$tf/fimo_out2/ $folder/$tf/ #output is BSs_region_for_scanning2.bed which has 3 cols
rm $folder/$tf/BSs_regions_for_getfasta.bed
rm $folder/$tf/BSs_seq_with_2_flanks.fa
mv $folder/$tf/BSs_regions_from_fasta.bed $folder/$tf/BSs_regions_for_scanning2.bed

bedtools intersect -wa -a $folder/$tf/BSs_regions_for_scanning2.bed -b $folder/$tf/BS.bed -v | awk '{print $1"\t"$2"\t"$3"\t"$4}' - | sort - | uniq > $folder/$tf/test.bed #get bed sequences accessible, however, not bound in vivo. 3 cols. at last, find sequences of each region in BSs_regions_for_scanning2.bed

python $folder/fromBedToFa.py $folder/$tf/test.bed $folder/$tf/test.fa   # make the non-binding sites into fasta format

python $folder/fromBedToFa.py $folder/$tf/BSs_regions_for_scanning1.bed $folder/$tf/BSs_regions_for_scanning1.fa   # make the binding sites into fasta format

python $folder/BiasAway.py g -b $folder/$tf/test.fa -f $folder/$tf/BSs_regions_for_scanning1.fa > $folder/$tf/BS_NO.fa   #randomly select non-binding sites from background to have matched GC content

python $folder/fromFaToBed.py $folder/$tf/BS_NO.fa $folder/$tf/non-BS.bed    # get BS_NO.bed for further analyasis

rm $folder/$tf/test.fa
rm $folder/$tf/test.bed
rm $folder/$tf/BSs_regions_for_scanning1.fa

awk -v var="$flank_length" '{print $1"\t"$2-var"\t"$3+var}' $folder/$tf/non-BS.bed > $folder/$tf/BSs_regions_for_scanning.bed

#scan for histone modification read coverage (RPM: reads per million), data is only for GM12878 cell line, in other cells lines, number of repetitive experiments for each HM may change.
for histone in h3k4me2 h3k27ac h3k27me3 h3k36me3 h3k4me1 h3k4me3 h3k79me2 h3k9ac h4k20me1
do
    if [ -f $folder/$tf/$histone\_BS_RPM_neg_$flank_length.txt ]; then
        rm $folder/$tf/$histone\_BS_RPM_neg_$flank_length.txt
    fi
    for rep in rep1 rep2
    do
         python $folder/histone_modification_feature.py $folder/10_histone_modification $folder/$tf $folder/10_histone_modification/$histone\_$rep\.bam $folder/$tf/BSs_regions_for_scanning_neg_$flank_length $folder/10_histone_modification/$histone\_$rep\_count.txt $folder/$tf/$histone\_BS_RPM_neg_$flank_length.txt $pwm_length

    done
    if [IsBpWise == 0];then
    awk '{a[$1"\t"$2"\t"$3"\t"]+=$4}END{for(i in a){print i,a[i]/2}}' $folder/$tf/$histone\_BS_RPM_neg_$flank_length.txt > $folder/$tf/$histone\_avg_BS_RPM_neg_$flank_length.txt
    rm $folder/$tf/$histone\_BS_RPM_neg_$flank_length.txt
    else
    awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"]+=$5}END{for(i in a){print i,a[i]/2}}' $folder/$tf/$histone\_BS_RPM_neg_$flank_length.txt | awk '{a[$1"\t"$2"\t"$3]=a[$1"\t"$2"\t"$3]"\t"$5}END{for(i in a){print i,a[i]}}' - > $folder/$tf/$histone\_avg_BS_RPM_neg_$flank_length_bp_wise.txt
    rm $folder/$tf/$histone\_BS_RPM_neg_$flank_length.txt
    fi
done

if [ -f $folder/$tf/h3k9me3_BS_RPM_neg_$flank_length.txt ]; then
        rm $folder/$tf/h3k9me3_BS_RPM_neg_$flank_length.txt
fi
for rep in rep1 rep2 rep3
do
         python $folder/histone_modification_feature.py $folder/10_histone_modification $folder/$tf $folder/10_histone_modification/h3k9me3_$rep\.bam $folder/$tf/BSs_regions_for_scanning_neg_$flank_length $folder/10_histone_modification/h3k9me3_$rep\_count.txt $folder/$tf/h3k9me3_BS_RPM_neg_$flank_length.txt $pwm_length

        
done
if [IsBpWise == 0];then
    awk '{a[$1"\t"$2"\t"$3"\t"]+=$4}END{for(i in a){print i,a[i]/3}}' $folder/$tf/h3k9me3_BS_RPM_neg_$flank_length.txt > $folder/$tf/h3k9me3_avg_BS_RPM_neg_$flank_length.txt
    rm $folder/$tf/h3k9me3_BS_RPM_neg_$flank_length.txt
else
    awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"]+=$5}END{for(i in a){print i,a[i]/3}}' $folder/$tf/h3k9me3_BS_RPM_neg_$flank_length.txt | awk '{a[$1"\t"$2"\t"$3]=a[$1"\t"$2"\t"$3]"\t"$5}END{for(i in a){print i,a[i]}}' - > $folder/$tf/h3k9me3_avg_BS_RPM_neg_$flank_length_bp_wise.txt
    rm $folder/$tf/h3k9me3_BS_RPM_neg_$flank_length.txt
fi

if [$IsBpWise == 0];
then
sort -k1,3 $folder/$tf/h3k4me2_avg_BS_RPM_neg_$flank_length.txt > $folder/$tf/$tf\_h3k4me2_HM_features.txt

sort -k1,3 $folder/$tf/h3k27ac_avg_BS_RPM_neg_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k4me2_HM_features.txt - > $folder/$tf/$tf\_h3k27ac_HM_features.txt

sort -k1,3 $folder/$tf/h3k27me3_avg_BS_RPM_neg_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k27ac_HM_features.txt - > $folder/$tf/$tf\_h3k27me3_HM_features.txt

sort -k1,3 $folder/$tf/h3k36me3_avg_BS_RPM_neg_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k27me3_HM_features.txt - > $folder/$tf/$tf\_h3k36me3_HM_features.txt

sort -k1,3 $folder/$tf/h3k4me1_avg_BS_RPM_neg_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k36me3_HM_features.txt - > $folder/$tf/$tf\_h3k4me1_HM_features.txt

sort -k1,3 $folder/$tf/h3k4me3_avg_BS_RPM_neg_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k4me1_HM_features.txt - > $folder/$tf/$tf\_h3k4me3_HM_features.txt

sort -k1,3 $folder/$tf/h3k79me2_avg_BS_RPM_neg_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k4me3_HM_features.txt - > $folder/$tf/$tf\_h3k79me2_HM_features.txt

sort -k1,3 $folder/$tf/h3k9ac_avg_BS_RPM_neg_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k79me2_HM_features.txt - > $folder/$tf/$tf\_h3k9ac_HM_features.txt

sort -k1,3 $folder/$tf/h4k20me1_avg_BS_RPM_neg_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h3k9ac_HM_features.txt - > $folder/$tf/$tf\_h4k20me1_HM_features.txt

sort -k1,3 $folder/$tf/h3k9me3_avg_BS_RPM_neg_$flank_length.txt | awk '{print $4}' - | paste $folder/$tf/$tf\_h4k20me1_HM_features.txt - > $folder/$tf/$tf\_h3k9me3_HM_features.txt
fi

cp $folder/$tf/$tf\_h3k9me3_HM_features.txt $folder/$tf/non-BS_10_HM.txt

for histone in h3k4me2 h3k27ac h3k27me3 h3k36me3 h3k4me1 h3k4me3 h3k79me2 h3k9ac h4k20me1 h3k9me3
do
    rm $folder/$tf/$histone\_avg_BS_RPM_neg_$flank_length.txt
    rm $folder/$tf/$tf\_h3k9ac_HM_features.txt
done


