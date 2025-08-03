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
#Output: FDR p-values.

##Usage: 
##Pvalue-FDR-stdout.sh <INPUT FILE> <COLUMN INDEX>

############################################################################
##ACTIONS:

##Input.
INPUTFILE=$(readlink -f $1)
PVALCOL=$2

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