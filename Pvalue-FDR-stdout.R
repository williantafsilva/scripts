#!/usr/bin/env Rscript
rm(list=ls()) #Clear environment.
ARGS=commandArgs(trailingOnly=TRUE)
############################################################################
############################# STDOUT R SCRIPT ##############################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
#SCRIPT DESCRIPTION:

#Description:
#Run the second step of the extended haplotype homozygosity (EHH) analysis on iHH data calculated in the first step of the analysis.

#Input $1: Tab-separated file with header containing a column with p-values.
#Input $2: Index of column that contains p-values.
#Output: FDR p-values.

#Usage: 
#Rscript --vanilla Pvalue-FDR-stdout <INPUT FILE> <COLUMN INDEX>

############################################################################
##ACTIONS:

library(R.utils)
library(tidyverse)

INPUTFILE<-normalizePath(ARGS[1])
PVALUECOLUMN<-ARGS[2]

INPUTDATA<-read.table(INPUTFILE,header=TRUE,sep="\t")

PADJUSTED<-p.adjust(INPUTDATA[,PVALUECOLUMN],method="fdr")

print(PADJUSTED)