#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Select columns of a matrix in a specific order.

##Input $1: Input tab-separated matrix file.
##Input $2: Input file with list of column names.
##Output: Print matrix subset.

##Usage: 
##matrix-selectcolumns-stdout.sh <INPUT MATRIX FILE> <INPUT COLUMNS LIST FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE=$1
LIST_COLS=$2

##Process.

COLVECTOR=""
while read COLNAME ; do
  COLNUMBER=$(head -n1 ${INPUTFILE} | sed -n "s/\(${COLNAME}\(\t\|\$\)\).*/\1/p" | awk '{print NF}' | sort -nu)
  if [[ ! -z "${COLNUMBER}" ]] ; then
    COLVECTOR=$(echo "${COLVECTOR},${COLNUMBER}" | sed 's/^,//')
  fi
done < ${LIST_COLS}

awk -v OFS='\t' -v "COLVECTOR=${COLVECTOR}" '
    BEGIN {split(COLVECTOR, cols, ",")} 
    {for (i in cols) printf("%s\t", $cols[i]); 
     print ""
}' ${INPUTFILE} | sed 's/\t$//g'
