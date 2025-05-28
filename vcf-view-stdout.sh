#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##View VCF/BCF file.

##Input: VCF/BCF file.
##Output: Print VCF/BCF file.

##Usage: 
##vcf-view-stdout.sh <VCF/BCF FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE=$(readlink -f $1)

##Process.

zcat ${INPUTFILE} | less -S