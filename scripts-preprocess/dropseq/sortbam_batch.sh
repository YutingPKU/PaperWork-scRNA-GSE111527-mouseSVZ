#!/bin/bash

while read -r line
do
	output="/lustre1/lch3000_pkuhpc/liuyt/SP/public/GSE111527-humanSVZ/velocyto/bamfiles/cellsorted_${line}"
	input="/lustre1/lch3000_pkuhpc/liuyt/SP/public/GSE111527-humanSVZ/velocyto/bamfiles/$line"
	CMD="pkurun-cns 1 20 samtools sort -@ 20 -o $output $input"
	echo $CMD
	eval $CMD
done < bam.ls
