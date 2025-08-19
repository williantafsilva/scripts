#!/usr/bin/env Rscript
rm(list=ls()) #Clear environment.
ARGS=commandArgs(trailingOnly=TRUE) 
RUNDATE=format(Sys.time(),"%Y%m%d%H%M%S")
SCRIPTNAME<-"MatrixQTL.R"
JOBID<-ARGS[1]
OUTPUTLOCATION<-normalizePath(ARGS[2])
sink(paste0(OUTPUTLOCATION,"/job",JOBID,".log"),type=c("output","message"))
############################################################################
################################ R SCRIPT ##################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
#SCRIPT DESCRIPTION:

#Description:
#Run MatrixEQTL.

#Input $1: Job ID.
#Input $2: Output location.
#Input $3: Genotype file (tab-separated).
#Input $4: Phenotype file (tab-separated).
#Input $5: Covariates file (tab-separated).
#Input $6: TRUE/FALSE to calculate FDR p-values (requires a lot of memory).
#Output: MatrixEQTL output file (*.txt).

#Usage: 
#Rscript --vanilla MatrixEQTL-otherphenotypes.R <JOB ID> <OUTPUT LOCATION> <GENOTYPE FILE> <PHENOTYPE FILE> <COVARIATES FILE> <FDR?>

############################################################################
#DIAGNOSTICS:

cat(paste0(format(Sys.time(),"%Y-%m-%d @ %H:%M:%S\n")))
cat(paste0("Submission: Rscript --vanilla ",SCRIPTNAME," ",paste(ARGS,collapse=" "),"\n"))
cat(paste0("User: ",Sys.getenv("USER"),"\n"))
cat(paste0("Home: ",Sys.getenv("HOME"),"\n"))
cat(paste0("PATH: ",Sys.getenv("PATH"),"\n"))
cat(paste0("Job ID: ",JOBID,"\n"))
sessionInfo()
cat("\n")

############################################################################
#LOAD TOOLS:

library(MatrixEQTL)

############################################################################
#INPUT:

#Genotype file.
if(ARGS[3]!="NA"){FILE_GENOTYPE<-normalizePath(ARGS[3])}else{FILE_GENOTYPE<-NA}

#Gene expression file.
if(ARGS[4]!="NA"){FILE_PHENOTYPE<-normalizePath(ARGS[4])}else{FILE_PHENOTYPE<-NA}

#Covariates file. In case of no covariates set the variable FILE_COVARIATES to character().
if(ARGS[5]!="NA"){FILE_COVARIATES<-normalizePath(ARGS[5])}else{FILE_COVARIATES<-character()}

#Calculate FDR p-values?
if(ARGS[6]==TRUE){
  SKIPFDRCALC<-FALSE
  OUTPUTNAMESTRING<-"FDR"
}else{
  SKIPFDRCALC<-TRUE
  OUTPUTNAMESTRING<-"noFDR"
}

############################################################################
#OUTPUT:

OUTPUTFILENAME<-paste0("matrixqtl-",OUTPUTNAMESTRING,"-job",JOBID,".txt")
OUTPUTFILE<-paste0(OUTPUTLOCATION,"/",OUTPUTFILENAME)

############################################################################
#ACTIONS:

#Set model.
#Options: 
#modelLINEAR (genotype is assumed to have only additive effect on expression).
#modelANOVA (genotype is assumed to have both additive and dominant effects; genotype data set musts take at most 3 distinct values).
#modelLINEAR_CROSS (test for the significance of the interaction between the genotype and the last covariate; test for equality of effect sizes between two groups of samples).
MODEL<-modelLINEAR #Set model. 

#Set significance threshold for all/distant tests.
#The p-value threshold determines which gene-SNP associations are saved in the output file.
#Note that for larger datasets the threshold should be lower. Setting the threshold to a high 
#value for a large dataset may cause excessively large output files.
if(ARGS[6]==TRUE){
  PTHRESHOLD<-1e-2
}else{
  PTHRESHOLD<-1e-2
}

#Define the covariance matrix for the error term. This parameter is rarely used. 
#If the covariance matrix is a multiple of identity, set it to numeric().
ERRORCOV<-numeric()

#Load files with genotype, gene expression, and covariates. One can set the file delimiter 
#(i.e. tabulation "\t", comma ",", or space " "), the string representation for missing values, 
#the number of rows with column labels, and the number of columns with row labels. 
#Finally, one can change the number of the variables in a slice for the file reading procedure 
#(do not change if not sure).

#Load genotype data.
DATA_GENOTYPE<-SlicedData$new()
DATA_GENOTYPE$fileDelimiter<-"\t" #Delimiter.
DATA_GENOTYPE$fileOmitCharacters<-"NA" #Missing values.
DATA_GENOTYPE$fileSkipRows<-1 #One row of column labels.
DATA_GENOTYPE$fileSkipColumns<-1 #One column of row labels.
DATA_GENOTYPE$fileSliceSize<-1000 #Read file in pieces of 1k rows.
DATA_GENOTYPE$LoadFile(FILE_GENOTYPE)

#Load gene expression data.
DATA_PHENOTYPE<-SlicedData$new()
DATA_PHENOTYPE$fileDelimiter<-"\t"
DATA_PHENOTYPE$fileOmitCharacters<-"NA"
DATA_PHENOTYPE$fileSkipRows<-1
DATA_PHENOTYPE$fileSkipColumns<-1
DATA_PHENOTYPE$fileSliceSize<-1000
DATA_PHENOTYPE$LoadFile(FILE_PHENOTYPE)

#Load covariates data.
DATA_COVARIATES<-SlicedData$new()
DATA_COVARIATES$fileDelimiter<-"\t"
DATA_COVARIATES$fileOmitCharacters<-"NA"
DATA_COVARIATES$fileSkipRows<-1
DATA_COVARIATES$fileSkipColumns<-1
if(length(FILE_COVARIATES)>0){
  DATA_COVARIATES$LoadFile(FILE_COVARIATES)
}

#Run MatrixEQTL.
MATRIXQTLRESULT<-Matrix_eQTL_main(DATA_GENOTYPE,
                                   DATA_PHENOTYPE,
                                   cvrt=DATA_COVARIATES, 
                                   output_file_name=OUTPUTFILE, 
                                   pvOutputThreshold=PTHRESHOLD,
                                   useModel=MODEL, 
                                   errorCovariance=ERRORCOV, 
                                   verbose=TRUE, 
                                   pvalue.hist=FALSE,
                                   noFDRsaveMemory=SKIPFDRCALC)

############################################################################
#SAVE CONTROL FILES:

cat(paste0("############################################################################
Date: ",RUNDATE,"
Job ID: ",JOBID,"
Script: ",SCRIPTNAME,"
Input genotype file: ",FILE_GENOTYPE,"
Input phenotype file: ",FILE_PHENOTYPE,"
Input covariates file: ",FILE_COVARIATES,"
Input FDR request: ",ARGS[6],"
Output file: ",OUTPUTFILE,"

"),file=paste0(OUTPUTLOCATION,"/README.txt"),append=TRUE)

sink()