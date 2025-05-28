#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Show lines that are common among up to ten lists.

##Input 1: File containing a list.
##Input 2: File containing a list.
##Input 3: File containing a list.
##Input 4: File containing a list.
##Input 5: File containing a list.
##Input 6: File containing a list.
##Input 7: File containing a list.
##Input 8: File containing a list.
##Input 9: File containing a list.
##Input 10: File containing a list.
##Output: Lines that are common to all input files.

##Usage: 
##commonlines-stdout.sh <FILE 1> <FILE 2> <FILE 3> <FILE 4> <FILE 5> <FILE 6> <FILE 7> <FILE 8> <FILE 9> <FILE 10>

############################################################################
##ACTIONS:

##Process.

if [[ "$#" == 2 ]] ; then
    INPUTFILE1=$(readlink -f $1)
	INPUTFILE2=$(readlink -f $2)
	comm -12 <(sort ${INPUTFILE1}) <(sort ${INPUTFILE2})
fi

if [[ "$#" == 3 ]] ; then
    INPUTFILE1=$(readlink -f $1)
	INPUTFILE2=$(readlink -f $2)
	INPUTFILE3=$(readlink -f $3)
	comm -12 <(sort ${INPUTFILE1}) <(sort ${INPUTFILE2}) | \
	comm -12 - <(sort ${INPUTFILE3})
fi

if [[ "$#" == 4 ]] ; then
    INPUTFILE1=$(readlink -f $1)
	INPUTFILE2=$(readlink -f $2)
	INPUTFILE3=$(readlink -f $3)
	INPUTFILE4=$(readlink -f $4)
	comm -12 <(sort ${INPUTFILE1}) <(sort ${INPUTFILE2}) | \
	comm -12 - <(sort ${INPUTFILE3}) | \
	comm -12 - <(sort ${INPUTFILE4})
fi

if [[ "$#" == 5 ]] ; then
    INPUTFILE1=$(readlink -f $1)
	INPUTFILE2=$(readlink -f $2)
	INPUTFILE3=$(readlink -f $3)
	INPUTFILE4=$(readlink -f $4)
	INPUTFILE5=$(readlink -f $5)
	comm -12 <(sort ${INPUTFILE1}) <(sort ${INPUTFILE2}) | \
	comm -12 - <(sort ${INPUTFILE3}) | \
	comm -12 - <(sort ${INPUTFILE4}) | \
	comm -12 - <(sort ${INPUTFILE5})
fi

if [[ "$#" == 6 ]] ; then
    INPUTFILE1=$(readlink -f $1)
	INPUTFILE2=$(readlink -f $2)
	INPUTFILE3=$(readlink -f $3)
	INPUTFILE4=$(readlink -f $4)
	INPUTFILE5=$(readlink -f $5)
	INPUTFILE6=$(readlink -f $6)
	comm -12 <(sort ${INPUTFILE1}) <(sort ${INPUTFILE2}) | \
	comm -12 - <(sort ${INPUTFILE3}) | \
	comm -12 - <(sort ${INPUTFILE4}) | \
	comm -12 - <(sort ${INPUTFILE5}) | \
	comm -12 - <(sort ${INPUTFILE6})
fi

if [[ "$#" == 7 ]] ; then
    INPUTFILE1=$(readlink -f $1)
	INPUTFILE2=$(readlink -f $2)
	INPUTFILE3=$(readlink -f $3)
	INPUTFILE4=$(readlink -f $4)
	INPUTFILE5=$(readlink -f $5)
	INPUTFILE6=$(readlink -f $6)
	INPUTFILE7=$(readlink -f $7)
	comm -12 <(sort ${INPUTFILE1}) <(sort ${INPUTFILE2}) | \
	comm -12 - <(sort ${INPUTFILE3}) | \
	comm -12 - <(sort ${INPUTFILE4}) | \
	comm -12 - <(sort ${INPUTFILE5}) | \
	comm -12 - <(sort ${INPUTFILE6}) | \
	comm -12 - <(sort ${INPUTFILE7})
fi

if [[ "$#" == 8 ]] ; then
    INPUTFILE1=$(readlink -f $1)
	INPUTFILE2=$(readlink -f $2)
	INPUTFILE3=$(readlink -f $3)
	INPUTFILE4=$(readlink -f $4)
	INPUTFILE5=$(readlink -f $5)
	INPUTFILE6=$(readlink -f $6)
	INPUTFILE7=$(readlink -f $7)
	INPUTFILE8=$(readlink -f $8)
	comm -12 <(sort ${INPUTFILE1}) <(sort ${INPUTFILE2}) | \
	comm -12 - <(sort ${INPUTFILE3}) | \
	comm -12 - <(sort ${INPUTFILE4}) | \
	comm -12 - <(sort ${INPUTFILE5}) | \
	comm -12 - <(sort ${INPUTFILE6}) | \
	comm -12 - <(sort ${INPUTFILE7}) | \
	comm -12 - <(sort ${INPUTFILE8})
fi

if [[ "$#" == 9 ]] ; then
    INPUTFILE1=$(readlink -f $1)
	INPUTFILE2=$(readlink -f $2)
	INPUTFILE3=$(readlink -f $3)
	INPUTFILE4=$(readlink -f $4)
	INPUTFILE5=$(readlink -f $5)
	INPUTFILE6=$(readlink -f $6)
	INPUTFILE7=$(readlink -f $7)
	INPUTFILE8=$(readlink -f $8)
	INPUTFILE9=$(readlink -f $9)
	comm -12 <(sort ${INPUTFILE1}) <(sort ${INPUTFILE2}) | \
	comm -12 - <(sort ${INPUTFILE3}) | \
	comm -12 - <(sort ${INPUTFILE4}) | \
	comm -12 - <(sort ${INPUTFILE5}) | \
	comm -12 - <(sort ${INPUTFILE6}) | \
	comm -12 - <(sort ${INPUTFILE7}) | \
	comm -12 - <(sort ${INPUTFILE8}) | \
	comm -12 - <(sort ${INPUTFILE9})
fi

if [[ "$#" == 10 ]] ; then
    INPUTFILE1=$(readlink -f $1)
	INPUTFILE2=$(readlink -f $2)
	INPUTFILE3=$(readlink -f $3)
	INPUTFILE4=$(readlink -f $4)
	INPUTFILE5=$(readlink -f $5)
	INPUTFILE6=$(readlink -f $6)
	INPUTFILE7=$(readlink -f $7)
	INPUTFILE8=$(readlink -f $8)
	INPUTFILE9=$(readlink -f $9)
	INPUTFILE10=$(readlink -f $10)
	comm -12 <(sort ${INPUTFILE1}) <(sort ${INPUTFILE2}) | \
	comm -12 - <(sort ${INPUTFILE3}) | \
	comm -12 - <(sort ${INPUTFILE4}) | \
	comm -12 - <(sort ${INPUTFILE5}) | \
	comm -12 - <(sort ${INPUTFILE6}) | \
	comm -12 - <(sort ${INPUTFILE7}) | \
	comm -12 - <(sort ${INPUTFILE8}) | \
	comm -12 - <(sort ${INPUTFILE9}) | \
	comm -12 - <(sort ${INPUTFILE10})
fi