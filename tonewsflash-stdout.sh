#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Add note to ${PROJHOME}/NEWSFLASH.txt file.

##Input 1: Quoted string to be added to ${PROJHOME}/NEWSFLASH.txt file.
##Output: Note on ${PROJHOME}/NEWSFLASH.txt file.

##Usage: 
##tonewsflash-stdout.sh <QUOTED STRING>

############################################################################
##ACTIONS:

RUNDATE=$(date +"%Y%m%d%H%M%S")

##Input.

ARGS="$@"

#Process.

IFS=$'\t'
echo "########################################################################
User: ${USER}
Date: ${RUNDATE}
${ARGS}
" >> $(echo "${PROJHOME}/NEWSFLASH.txt")
IFS='$ORIGINALIFS'