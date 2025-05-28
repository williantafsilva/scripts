#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Count the number of tab-separated columns in the first line of each file.

##Input: Input files.
##Output: Number of tab-separated columns in the first line in each file.

##Usage: 
##countcolumns-stdout.sh <INPUT FILE>

############################################################################
##ACTIONS:

##Input.

INPUTFILE="$@"

##Process.

IFS=$' '
echo "Number of columns:"
find ${INPUTFILE} -maxdepth 0 | while read F ; do
	FILEX=$(readlink -f ${F})
	echo "-----> ${FILEX}: $(head -n1 ${FILEX} | tr '\t' '\n' | wc -l)"
done
IFS='$ORIGINALIFS'