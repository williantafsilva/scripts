#!/bin/bash -l
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Given a list of column names, extract those columns from a matrix with a specific delimiter.

##Input $1: Delimiter. Options: comma, tab, space, semicolon.
##Input $2: List of column names.
##Input $3: Matrix file.
##Output: File with matrix with selected columns.

##Usage: 
##extractcolumns-stdout.sh <DELIMITER> <COLUMN LIST> <INPUT MATRIX FILE>

############################################################################
##ACTIONS:

RUNDATE=$(date +"%Y%m%d%H%M%S")

##Input.

DELIMITER=$1
COLUMNLIST=$2
INPUTFILE=$3

##Process.

TMPFILE1="tmp1-extractcolumns-${RUNDATE}.txt"
TMPFILE2="tmp2-extractcolumns-${RUNDATE}.txt"
TMPFILE3="tmp3-extractcolumns-${RUNDATE}.txt"

OUTPUTFILE="output-extractcolumns-${RUNDATE}.txt"

if [[ ${DELIMITER} == comma ]] ; then
  FIRSTCOLUMN=true
  #Create array of column names.
  SAMPLEARRAY=($(head -n1 ${INPUTFILE} | tr ',' '\n'))

  #Select columns.
  cat ${COLUMNLIST} | while read C ; do
    for i in "${!SAMPLEARRAY[@]}"; do
        if [[ "${SAMPLEARRAY[$i]}" = "${C}" ]] ; then
            cat ${INPUTFILE} | cut -d',' -f $((${i}+1)) > "${TMPFILE2}"
            if [[ ${FIRSTCOLUMN} ]] ; then
              "${TMPFILE2}" > "${TMPFILE3}"
              FIRSTCOLUMN=false
            else
              paste -d',' "${TMPFILE1}" "${TMPFILE2}" > "${TMPFILE3}"
            fi
        fi
    done
    cat "${TMPFILE3}" > "${TMPFILE1}"
  done

  cat "${TMPFILE1}" > ${OUTPUTFILE}

  #Delete temporary files.
  rm -f "${TMPFILE1}" "${TMPFILE2}" "${TMPFILE3}"
fi

if [[ ${DELIMITER} == tab ]] ; then
  FIRSTCOLUMN=true
  #Create array of column names.
  SAMPLEARRAY=($(head -n1 ${INPUTFILE} | tr '\t' '\n'))

  #Select columns.
  cat ${COLUMNLIST} | while read C ; do
    for i in "${!SAMPLEARRAY[@]}"; do
        if [[ "${SAMPLEARRAY[$i]}" = "${C}" ]] ; then
            cat ${INPUTFILE} | cut -d'\t' -f $((${i}+1)) > "${TMPFILE2}"
            if [[ ${FIRSTCOLUMN} ]] ; then
              "${TMPFILE2}" > "${TMPFILE3}"
              FIRSTCOLUMN=false
            else
              paste -d'\t' "${TMPFILE1}" "${TMPFILE2}" > "${TMPFILE3}"
            fi
        fi
    done
    cat "${TMPFILE3}" > "${TMPFILE1}"
  done

  cat "${TMPFILE1}" > ${OUTPUTFILE}

  #Delete temporary files.
  rm -f "${TMPFILE1}" "${TMPFILE2}" "${TMPFILE3}"
fi

if [[ ${DELIMITER} == space ]] ; then
  FIRSTCOLUMN=true
  #Create array of column names.
  SAMPLEARRAY=($(head -n1 ${INPUTFILE} | tr ' ' '\n'))

  #Select columns.
  cat ${COLUMNLIST} | while read C ; do
    for i in "${!SAMPLEARRAY[@]}"; do
        if [[ "${SAMPLEARRAY[$i]}" = "${C}" ]] ; then
            cat ${INPUTFILE} | cut -d' ' -f $((${i}+1)) > "${TMPFILE2}"
            if [[ ${FIRSTCOLUMN} ]] ; then
              "${TMPFILE2}" > "${TMPFILE3}"
              FIRSTCOLUMN=false
            else
              paste -d' ' "${TMPFILE1}" "${TMPFILE2}" > "${TMPFILE3}"
            fi
        fi
    done
    cat "${TMPFILE3}" > "${TMPFILE1}"
  done

  cat "${TMPFILE1}" > ${OUTPUTFILE}

  #Delete temporary files.
  rm -f "${TMPFILE1}" "${TMPFILE2}" "${TMPFILE3}"
fi

if [[ ${DELIMITER} == semicolon ]] ; then
  FIRSTCOLUMN=true
  #Create array of column names.
  SAMPLEARRAY=($(head -n1 ${INPUTFILE} | tr ';' '\n'))

  #Select columns.
  cat ${COLUMNLIST} | while read C ; do
    for i in "${!SAMPLEARRAY[@]}"; do
        if [[ "${SAMPLEARRAY[$i]}" = "${C}" ]] ; then
            cat ${INPUTFILE} | cut -d';' -f $((${i}+1)) > "${TMPFILE2}"
            if [[ ${FIRSTCOLUMN} ]] ; then
              "${TMPFILE2}" > "${TMPFILE3}"
              FIRSTCOLUMN=false
            else
              paste -d';' "${TMPFILE1}" "${TMPFILE2}" > "${TMPFILE3}"
            fi
        fi
    done
    cat "${TMPFILE3}" > "${TMPFILE1}"
  done

  cat "${TMPFILE1}" > ${OUTPUTFILE}

  #Delete temporary files.
  rm -f "${TMPFILE1}" "${TMPFILE2}" "${TMPFILE3}"
fi