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
##Input $2: Input file with list of column names, or comma-separated list of columns names.
##Output: Print matrix subset.

##Usage: 
##matrix-selectcolumns-stdout.sh <INPUT MATRIX FILE> <INPUT COLUMNS LIST FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE=$1
LIST_COLS=$2

##Process.

if [[ -f "${LIST_COLS}" ]] ; then

  COLVECTOR=""
  while read COLNAME ; do
    COLNUMBER=$(head -n1 ${INPUTFILE} | sed -n "s/\(${COLNAME}\(\t\|\$\)\).*/\1/p" | awk '{print NF}' | sort -nu)
    if [[ ! -z "${COLNUMBER}" ]] ; then
      COLVECTOR=$(echo "${COLVECTOR},${COLNUMBER}" | sed 's/^,//')
    fi
  done < ${LIST_COLS}

else

  COLVECTOR=""
  for COLNAME in $(echo ${LIST_COLS} | tr "," " ") ; do
    COLNUMBER=$(head -n1 ${INPUTFILE} | sed -n "s/\(${COLNAME}\(\t\|\$\)\).*/\1/p" | awk '{print NF}' | sort -nu)
    if [[ ! -z "${COLNUMBER}" ]] ; then
      COLVECTOR=$(echo "${COLVECTOR},${COLNUMBER}" | sed 's/^,//')
    fi
  done

fi

awk -F'\t' -v OFS='\t' -v "COLVECTOR=${COLVECTOR}" '
    BEGIN {
        N = split( COLVECTOR, COLS, ",")
    }{ for (i = 1; i <= N; i++)
            printf "%s%s", $COLS[i], (i < N ? OFS : ORS)
    }' "${INPUTFILE}" | sed 's/\t$//g'

