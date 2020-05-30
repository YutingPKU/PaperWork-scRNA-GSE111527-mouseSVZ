#!/bin/bash

while read -r line
do
	dirdrop="/lustre1/lch3000_pkuhpc/liuyt/software/Drop-seq_tools-2.3.0"
	output="../dropseq/${line}/dge.txt.gz"
	input="../dropseq/${line}/final.bam"
	#mkdir -p "../dropseq/tmp/$line"
	CMD="$dirdrop/DigitalExpression I=$input O=$output MIN_NUM_GENES_PER_CELL=1 "
	echo $CMD >> DGE.sh

done < sra.ls


