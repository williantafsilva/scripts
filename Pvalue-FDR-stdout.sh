#!/bin/bash
############################################################################
############################# STDOUT FUNCTION ##############################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Calculate FDR p-values.
##Requirements: R module.

#Input $1: Tab-separated file with header containing a column with p-values.
#Input $2: Index of column that contains p-values.
#Input $3 (optional): Total number of tests (when only significant p-values are provided).
#Output: FDR p-values.

##Usage: 
##Pvalue-FDR-stdout.sh <INPUT FILE> <COLUMN INDEX> <NUMBER OF TESTS>

############################################################################
##ACTIONS:

##Input.
INPUTFILE=$(readlink -f $1)
PVALCOL=$2
NTESTS=$3

if [[ -z "${NTESTS}" ]] ; then

	#Run R to perform FDR correction.
	Rscript --vanilla -e "
#Read input data.
DATA<-read.table('${INPUTFILE}',header=TRUE,sep='\t',stringsAsFactors=FALSE)

#Extract p-values.
PVALUES<-as.numeric(DATA[,${PVALCOL}])

#Perform FDR correction.
FDRPVALUES<-p.adjust(PVALUES,method='fdr')

#Print FDR p-values.
cat(FDRPVALUES,sep='\n')
"

else

    #Run R to perform FDR correction.
	Rscript --vanilla -e "
#Read input data.
DATA<-read.table('${INPUTFILE}',header=TRUE,sep='\t',stringsAsFactors=FALSE)

#Extract p-values.
PVALUES<-as.numeric(DATA[,${PVALCOL}])

#Total number of tests.
TOTALNTESTS<-as.numeric(${NTESTS})

#Perform FDR correction.
FDRPVALUES<-p.adjust(PVALUES,method='fdr',n=TOTALNTESTS)

#Print FDR p-values.
cat(FDRPVALUES,sep='\n')
"

fi