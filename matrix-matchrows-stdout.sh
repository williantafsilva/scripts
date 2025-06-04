#!/bin/bash -l
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Match the rows of two matrices according to the input list of row names or, 
##if such list is not provided, using rows that are common to both input files.

##Input $1: Matrix file 1.
##Input $2: Matrix file 2.
##Input $3: File with list of row names.
##Output: Created output files corresponding to input files.

##Usage: 
##matrix-matchrows-stdout.sh <INPUT FILE 1> <INPUT FILE 2> <INPUT LIST FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE1=$(readlink -f $1)
INPUTFILE2=$(readlink -f $2)
LIST_ROWS=$3

##Output.

OUTPUTFILE1="$(echo "${INPUTFILE1%.*}" | sed 's/\.[^.]*$//').matrixmatchrows.txt"
OUTPUTFILE2="$(echo "${INPUTFILE2%.*}" | sed 's/\.[^.]*$//').matrixmatchrows.txt"

##Process.

if [[ ! -z "${LIST_ROWS}" ]] ; then

  cat ${LIST_ROWS} | while read ROWNAME ; do
    grep --color=never -P "^${ROWNAME}\t" ${INPUTFILE1} >> ${OUTPUTFILE1}
    grep --color=never -P "^${ROWNAME}\t" ${INPUTFILE2} >> ${OUTPUTFILE2}
  done

else

  comm -12 <(cut -f1 ${INPUTFILE1} | sort) <(cut -f1 ${INPUTFILE2} | sort) | while read ROWNAME ; do
    grep --color=never -P "^${ROWNAME}\t" ${INPUTFILE1} >> ${OUTPUTFILE1}
    grep --color=never -P "^${ROWNAME}\t" ${INPUTFILE2} >> ${OUTPUTFILE2}
  done

fi





