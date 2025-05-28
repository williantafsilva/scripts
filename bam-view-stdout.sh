#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##View BAM file.
##Requirements: samtools.

##Input: BAM file.
##Output: Print BAM file.

##Usage: 
##bam-view-stdout.sh <BAM FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE=$(readlink -f $1)

##Process.

samtools view -h ${INPUTFILE} | less -S