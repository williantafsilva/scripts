#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##List samples in a VCF/BCF file.
##Requirements: bcftools.

##Input: VCF/BCF file.
##Output: List of samples.

##Usage: 
##vcf-samplelist-stdout.sh <VCF/BCF FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE=$(readlink -f $1)

##Process.

bcftools query -l ${INPUTFILE}