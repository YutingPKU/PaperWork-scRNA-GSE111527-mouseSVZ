#!/bin/bash

while read -r line
do
	f1="../cleandata/${line}_1.fastq"
	f2="../cleandata/${line}_2.fastq"
	picard="/lustre1/lch3000_pkuhpc/liuyt/software/picard-tools-1.140/picard.jar"
	out="../cleandata/${line}.unmappedQuerynamesort.bam"
	CMD="pkurun-cns 1 1  java -jar $picard FastqToSam F1=$f1 F2=$f2 O=$out SM=$line"
	echo $CMD
	eval $CMD
done < sra.ls
