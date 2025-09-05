#!/bin/bash -l
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Delete all-zero rows from a matrix (skipping column 1).

##Input $1: Matrix file.
##Output: Print matrix without all-zero rows.

##Usage: 
##matrix-removeallzerorows-stdout.sh <INPUT FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE=$1

##Process.

awk -F'\t' '{
    all_zero = 1
    for (i=2; i<=NF; i++) {
        if ($i != 0) {
            all_zero = 0
            break
        }
    }
    if (all_zero == 0) {
        print
    }
}' ${INPUTFILE}
