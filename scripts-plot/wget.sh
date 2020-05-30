#!/bin/bash

ascp -i /home/liuyf/.aspera/connect/etc/asperaweb_id_dsa.openssh -QT -l2400M anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR681/SRR6814853/SRR6814853.sra . >> SRR6814853.ascp.log 2>&1
	 

#~/.aspera/connect/bin/ascp -i /home/liuyt/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 –T -l200m  anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR165/SRR1658580/SRR1658580.sra . 
