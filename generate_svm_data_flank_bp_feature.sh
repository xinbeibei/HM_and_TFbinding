#!/bin/sh

tf=$1
tf_length=$2
folder=$3
celltype=$4
shape_flank=$5

#binary code
#00011110000   shape features
#00000000001 	HM features
#10011110000	seq+shape features
#10011110001	seq+shape+HM features

R --no-restore --no-save --args $folder/$tf/$tf.txt $folder/$tf $tf_length $shape_flank < ./encode_custom.R 2>&1 1>/dev/null
#create seq+shape data
awk '{out=$3; for(i=4; i<=NF; i++){out=out"\t"$i}; print out}' $folder/encoded_$a\_10_down/$tf\.00011110000 | paste $folder/encoded_$a\_10_down/$tf\.10000000000 - > $folder/encoded_$a\_10_down/$tf\.10011110000
#create seq+shape+HM data
awk '{out=$3; for(i=4; i<=NF; i++){out=out"\t"$i}; print out}' $folder/encoded_$a\_10_down/$tf\.00000000001 | paste $folder/encoded_$a\_10_down/$tf\.10011110000 - > $folder/encoded_$a\_10_down/$tf\.10011110001

done