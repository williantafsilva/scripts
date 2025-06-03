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

echo "------ END ------"




