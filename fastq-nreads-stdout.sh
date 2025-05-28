#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Count number of reads and average read length in a FastQ file.

##Input: FastQ file (*.fastq.gz or *.fq.gz).
##Output: Number of reads.

##Usage: 
##fastq-nreads-stdout.sh <FASTQ FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE=$(readlink -f $1)

##Process.

##NREADS=0
##TOTALREADLENGTH=0
##zcat "${INPUTFILE}" | awk 'NR % 4 == 2 { NREADS++; TOTALREADLENGTH += length($0) } END { if (NREADS > 0) { print "Number of reads:", NREADS; print "Average read length:", TOTALREADLENGTH / NREADS } else { print "No reads found." } }'

NREADS=$(echo $(zcat ${INPUTFILE} | wc -l)/4 | bc)
echo "Number of reads: ${NREADS}"
TOTALREADLENGTH=$(zcat ${INPUTFILE} | awk 'NR % 4 == 2' | awk '{TOT += length($0)} END {print TOT}')
AVEREADLENGTH=$((TOTALREADLENGTH / NREADS))
echo "Average read length: ${AVEREADLENGTH}"