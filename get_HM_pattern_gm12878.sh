#!/bin/sh

tf=$1
folder=$2

awk '$1==1{print $_}' $folder/$tf/$tf.00000000001 > $folder/$tf/BS_10_HM.txt
awk '$1==0{print $_}' $folder/$tf/$tf.00000000001 > $folder/$tf/non-BS_10_HM.txt

