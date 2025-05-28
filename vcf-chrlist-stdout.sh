#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##View list of chromosomes/scaffolds in a VCF/BCF file.
##Requirements: bcftools.

##Input: Indexed VCF/BCF file.
##Output: List of chromosomes/scaffolds.

##Usage: 
##vcf-chrlist-stdout.sh <VCF/BCF FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE=$(readlink -f $1)

##Process.

bcftools index -s ${INPUTFILE} | cut -f 1