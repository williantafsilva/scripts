#!/bin/bash -l
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Select columns of a matrix in a specific order.

##Input $1: Target list file.
##Input $2: Map file with tab-separated lists of value to replace and corresponding replacement value (VALUE\tREPLACEMENT).
##Output: Print new list.

##Usage: 
##maplist-fromto-stdout.sh <TARGET LIST FILE> <MAP FILE>

############################################################################
##ACTIONS:

##Input.

INPUTLIST=$1
FROMTOLIST=$2

##Process.

awk -v OFS='\t' '
NR == FNR {
   map[$1] = $2
   next
}
{
   for (i=1; i<=NF; ++i)
      $i in map && $i = map[$i]
}
1' ${FROMTOLIST} ${INPUTLIST}

