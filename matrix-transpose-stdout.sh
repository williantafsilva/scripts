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

DELIMITER=$1
INPUTFILE=$2

##Process.

if [[ ${DELIMITER} == comma ]] ; then
    awk -F',' '{
      for (i=1; i<=NF; i++) {
        a[NR,i] = $i;
      }
      max_nf = NF > max_nf ? NF : max_nf;
      max_nr = NR;
    }
    END {
      for (i=1; i<=max_nf; i++) {
        for (j=1; j<=max_nr; j++) {
          printf "%s%s", a[j,i], (j==max_nr ? "\n" : ",");
        }
      }
    }' ${INPUTFILE}
fi

if [[ ${DELIMITER} == tab ]] ; then
    awk -F'\t' '{
      for (i=1; i<=NF; i++) {
        a[NR,i] = $i;
      }
      max_nf = NF > max_nf ? NF : max_nf;
      max_nr = NR;
    }
    END {
      for (i=1; i<=max_nf; i++) {
        for (j=1; j<=max_nr; j++) {
          printf "%s%s", a[j,i], (j==max_nr ? "\n" : "\t");
        }
      }
    }' ${INPUTFILE}
fi

if [[ ${DELIMITER} == space ]] ; then
    awk -F' ' '{
      for (i=1; i<=NF; i++) {
        a[NR,i] = $i;
      }
      max_nf = NF > max_nf ? NF : max_nf;
      max_nr = NR;
    }
    END {
      for (i=1; i<=max_nf; i++) {
        for (j=1; j<=max_nr; j++) {
          printf "%s%s", a[j,i], (j==max_nr ? "\n" : " ");
        }
      }
    }' ${INPUTFILE}
fi

if [[ ${DELIMITER} == semicolon ]] ; then
    awk -F';' '{
      for (i=1; i<=NF; i++) {
        a[NR,i] = $i;
      }
      max_nf = NF > max_nf ? NF : max_nf;
      max_nr = NR;
    }
    END {
      for (i=1; i<=max_nf; i++) {
        for (j=1; j<=max_nr; j++) {
          printf "%s%s", a[j,i], (j==max_nr ? "\n" : ";");
        }
      }
    }' ${INPUTFILE}
fi