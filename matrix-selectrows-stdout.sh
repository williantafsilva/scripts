#!/bin/bash -l
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Select rows of a matrix in a specific order.

##Input $1: Input tab-separated matrix file.
##Input $2: Input file with list of row names.
##Output: Print matrix subset.

##Usage: 
##matrix-selectrows-stdout.sh <INPUT MATRIX FILE> <INPUT ROWS LIST FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE=$1
LIST_ROWS=$2

##Process.

cat ${LIST_ROWS} | while read ROWNAME ; do
  grep --color=never -P "^${ROWNAME}\t" ${INPUTFILE}
done

echo "------ END ------"