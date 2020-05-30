#!/bin/bash

while read -r line
do
	CMD="pkurun-cns 1 1 fastq-dump --split-files $line"
	echo $CMD
	eval $CMD
done < sra.ls
