#!/bin/bash -l
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Match the columns of two matrices according to the input list of column names or, 
##if such list is not provided, using columns that are common to both input files.

##Input $1: Matrix file 1.
##Input $2: Matrix file 2.
##Input $3: File with list of column names.
##Output: Created output files corresponding to input files.

##Usage: 
##matrix-matchcolumns-stdout.sh <INPUT FILE 1> <INPUT FILE 2> <INPUT LIST FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE1=$(readlink -f $1)
INPUTFILE2=$(readlink -f $2)
LIST_COLS=$3

##Output.

OUTPUTFILE1="$(echo "${INPUTFILE1%.*}" | sed 's/\.[^.]*$//').matrixmatchcolumns.txt"
OUTPUTFILE2="$(echo "${INPUTFILE2%.*}" | sed 's/\.[^.]*$//').matrixmatchcolumns.txt"

##Process.

if [[ ! -z "${LIST_COLS}" ]] ; then

  COLVECTOR1=""
  COLVECTOR2=""
  while read COLNAME ; do
    COLNUMBER1=$(head -n1 ${INPUTFILE1} | sed -n "s/\(${COLNAME}\(\t\|\$\)\).*/\1/p" | awk '{print NF}' | sort -nu)
    COLNUMBER2=$(head -n1 ${INPUTFILE2} | sed -n "s/\(${COLNAME}\(\t\|\$\)\).*/\1/p" | awk '{print NF}' | sort -nu)
    COLVECTOR1=$(echo "${COLVECTOR1},${COLNUMBER1}" | sed 's/^,//')
    COLVECTOR2=$(echo "${COLVECTOR2},${COLNUMBER2}" | sed 's/^,//')
  done < ${LIST_COLS}

  awk -v OFS='\t' -v "COLVECTOR=${COLVECTOR1}" '
      BEGIN {split(COLVECTOR, cols, ",")} 
      {for (i in cols) printf("%s\t", $cols[i]); 
       print ""
  }' ${INPUTFILE1} | sed 's/\t$//g' > ${OUTPUTFILE1}

  awk -v OFS='\t' -v "COLVECTOR=${COLVECTOR2}" '
      BEGIN {split(COLVECTOR, cols, ",")} 
      {for (i in cols) printf("%s\t", $cols[i]); 
       print ""
  }' ${INPUTFILE2} | sed 's/\t$//g' > ${OUTPUTFILE2}

else

  COLVECTOR1=""
  COLVECTOR2=""
  while read COLNAME ; do
    COLNUMBER1=$(head -n1 ${INPUTFILE1} | sed -n "s/\(${COLNAME}\(\t\|\$\)\).*/\1/p" | awk '{print NF}' | sort -nu)
    COLNUMBER2=$(head -n1 ${INPUTFILE2} | sed -n "s/\(${COLNAME}\(\t\|\$\)\).*/\1/p" | awk '{print NF}' | sort -nu)
    COLVECTOR1=$(echo "${COLVECTOR1},${COLNUMBER1}" | sed 's/^,//')
    COLVECTOR2=$(echo "${COLVECTOR2},${COLNUMBER2}" | sed 's/^,//')
  done < $(comm -12 <(head -n1 ${INPUTFILE1} | tr '\t' '\n' | sort) <(head -n1 ${INPUTFILE2} | tr '\t' '\n' | sort))

  awk -v OFS='\t' -v "COLVECTOR=${COLVECTOR1}" '
      BEGIN {split(COLVECTOR, cols, ",")} 
      {for (i in cols) printf("%s\t", $cols[i]); 
       print ""
  }' ${INPUTFILE1} | sed 's/\t$//g' > ${OUTPUTFILE1}

  awk -v OFS='\t' -v "COLVECTOR=${COLVECTOR2}" '
      BEGIN {split(COLVECTOR, cols, ",")} 
      {for (i in cols) printf("%s\t", $cols[i]); 
       print ""
  }' ${INPUTFILE2} | sed 's/\t$//g' > ${OUTPUTFILE2}

fi

