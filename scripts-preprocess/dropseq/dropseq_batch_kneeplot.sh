#!/bin/bash

while read -r line
do
	dirdrop="/lustre1/lch3000_pkuhpc/liuyt/software/Drop-seq_tools-2.3.0"
	output="../dropseq/${line}/cell_readcounts.txt.gz"
	input="../dropseq/${line}/final.bam"
	#mkdir -p "../dropseq/tmp/$line"
	CMD="$dirdrop/BamTagHistogram I=$input O=$output TAG=XC "
	echo $CMD >> knee.plot.sh

done < sra.ls


