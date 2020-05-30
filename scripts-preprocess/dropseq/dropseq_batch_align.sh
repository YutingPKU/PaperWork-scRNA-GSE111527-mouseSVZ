#!/bin/bash

while read -r line
do
	gn="/lustre1/lch3000_pkuhpc/liuyt/genome/Mm10_STAR_Index"
	ref="/lustre1/lch3000_pkuhpc/liuyt/software/Drop-seq-2.3.0/metadata/mm10/mm10.fasta"
	dirdrop="/lustre1/lch3000_pkuhpc/liuyt/software/Drop-seq_tools-2.3.0"
	outdir="../dropseq/$line"
	mkdir -p $outdir
	input="../cleandata/${line}.unmappedQuerynamesort.bam"
	#mkdir -p "../dropseq/tmp/$line"
	CMD="Drop-seq_alignment.sh -g $gn -r $ref -d $dirdrop -k -o $outdir -t ../dropseq/tmp/$line -e $input "
	echo $CMD
	eval $CMD > ${line}_dropseq.sh
done < sra.ls


