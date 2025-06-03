#!/bin/bash -l
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Transpose a matrix.

##Input $1: Delimiter. Options: comma, tab, space, semicolon.
##Input $2: Matrix file.
##Output: Print transposed matrix.

##Usage: 
##transpose-stdout.sh <DELIMITER> <INPUT FILE>

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