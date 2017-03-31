#!/bin/sh

folder=$1
tf=$2
tf_length=$3

R --no-restore --no-save --args $folder/$tf/$tf\.txt $folder/$tf/ $tf_length $folder/feature_list < $folder/encode_custom.R

awk '{out=$3; for(i=4; i<=NF; i++){out=out"\t"$i}; print out}' $folder/$tf/$tf\.00000000001 | paste $folder/$tf/$tf\.10011110000 - > $folder/$tf/$tf\.11000000000

awk '{print $3}' $folder/$tf/$tf\.00000000010 | paste $folder/$tf/$tf\.10011110000 - > $folder/$tf/$tf\.10000000010

awk '{out=$3; for(i=4; i<=NF; i++){out=out"\t"$i}; print out}' $folder/$tf/$tf\.00000000001 | paste $folder/$tf/$tf\.10000000000 - > $folder/$tf/$tf\.10000000001
