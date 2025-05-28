#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Print length of FASTA sequences.

##Input: FASTA file.
##Output: Length of FASTA sequences.

##Usage: 
##fasta-seqlength-stdout.sh <FASTA FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE=$(readlink -f $1)

##Process.

IFS=$'\t'
awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' ${INPUTFILE}
echo "Number of sequences: $(cat ${INPUTFILE} | grep "^>" | wc -l)"
IFS='$ORIGINALIFS'