#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Show similarities and differences between two lists.

##Input 1: File containing a list.
##Input 2: File containing a list.
##Output: Lines that are unique to file 1 (1st column), unique to file 2 (2nd column) and common to both files (3rd column).

##Usage: 
##comparelists-stdout.sh <FILE 1> <FILE 2>

############################################################################
##ACTIONS:

##Input.

#INPUTFILE1=$(readlink -f $1)
#INPUTFILE2=$(readlink -f $2)

##Process.

echo -e "FILE1ONLY\tFILE2ONLY\tBOTH"
echo -e "----------\t----------\t----------"
#comm <(sort ${INPUTFILE1}) <(sort ${INPUTFILE2})
paste -d'\t' <(comm -23 <(sort $1 | uniq) <(sort $2 | uniq)) \
<(comm -13 <(sort $1 | uniq) <(sort $2 | uniq)) \
<(comm -12 <(sort $1 | uniq) <(sort $2 | uniq))
NFILE1ONLY=$(comm -23 <(sort $1 | uniq) <(sort $2 | uniq) | wc -l)
NFILE2ONLY=$(comm -13 <(sort $1 | uniq) <(sort $2 | uniq) | wc -l)
NBOTH=$(comm -12 <(sort $1 | uniq) <(sort $2 | uniq) | wc -l)
echo -e "----------\t----------\t----------"
echo -e "${NFILE1ONLY}\t${NFILE2ONLY}\t${NBOTH}"