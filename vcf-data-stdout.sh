#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##View VCF/BCF file data (no header).
##Requirements: bcftools.

##Input: VCF/BCF file.
##Output: VCF/BCF file data (no header).

##Usage: 
##vcf-data-stdout.sh <VCF/BCF FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE=$(readlink -f $1)

##Process.

bcftools view -H ${INPUTFILE} | less -S