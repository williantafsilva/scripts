#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Delete all-X rows from a matrix (skipping column 1).

##Input $1: Matrix file.
##Input $2: X value.
##Input $3: Skip first column? (TRUE/FALSE).
##Output: Print matrix without all-X rows.

##Usage: 
##matrix-removeallXrows-stdout.sh <INPUT FILE> <X> <TRUE/FALSE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE=$1
XVALUE=$2
SKIPCOL1=$3

##Process.
if [[ ${SKIPCOL1} ]]; then

    awk -F'\t' -v XVALUE="${XVALUE}" '{
        all_X = 1
        for (i=2; i<=NF; i++) {
            if ($i != XVALUE) {
                all_X = 0
                break
            }
        }
        if (all_X == 0) {
            print
        }
    }' ${INPUTFILE}

else

    awk -F'\t' -v XVALUE="${XVALUE}" '{
        all_X = 1
        for (i=1; i<=NF; i++) {
            if ($i != XVALUE) {
                all_X = 0
                break
            }
        }
        if (all_X == 0) {
            print
        }
    }' ${INPUTFILE}

fi

