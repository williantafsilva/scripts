#!/usr/bin/env Rscript
rm(list=ls()) #Clear environment.
ARGS=commandArgs(trailingOnly=TRUE) 
RUNDATE=format(Sys.time(),"%Y%m%d%H%M%S")
SCRIPTNAME<-"MatrixEQTL.R"
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
#Input $4: Expression file (tab-separated).
#Input $5: Covariates file (tab-separated).
#Input $6: SNP location file (tab-separated).
#Input $7: Gene location file (tab-separated).
#Input $8: TRUE/FALSE to calculate FDR p-values (requires a lot of memory).
#Output: MatrixEQTL output files (*.txt).

#Usage: 
#Rscript --vanilla MatrixEQTL.R <JOB ID> <OUTPUT LOCATION> <GENOTYPE FILE> <EXPRESSION FILE> <COVARIATES FILE> <SNP LOCATION FILE> <GENE LOCATION FILE> <FDR?>

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
if(ARGS[4]!="NA"){FILE_EXPRESSION<-normalizePath(ARGS[4])}else{FILE_EXPRESSION<-NA}

#Covariates file. In case of no covariates set the variable FILE_COVARIATES to character().
if(ARGS[5]!="NA"){FILE_COVARIATES<-normalizePath(ARGS[5])}else{FILE_COVARIATES<-character()}

#SNP location file.
if(ARGS[6]!="NA"){
  FILE_SNPLOCATION<-normalizePath(ARGS[6])
  DATA_SNPLOCATION<-read.table(FILE_SNPLOCATION,header=TRUE,stringsAsFactors=FALSE)
}else{
  FILE_SNPLOCATION<-NULL
  DATA_SNPLOCATION<-NULL
}

#Gene location file.
if(ARGS[7]!="NA"){
  FILE_GENELOCATION<-normalizePath(ARGS[7])
  DATA_GENELOCATION<-read.table(FILE_GENELOCATION,header=TRUE,stringsAsFactors=FALSE)
}else{
  FILE_GENELOCATION<-NULL
  DATA_GENELOCATION<-NULL
}

#Calculate FDR p-values?
if(ARGS[8]==TRUE){
  SKIPFDRCALC<-FALSE
  OUTPUTNAMESTRING<-"FDR"
}else{
  SKIPFDRCALC<-TRUE
  OUTPUTNAMESTRING<-"noFDR"
}

############################################################################
#OUTPUT:

OUTPUTFILE1NAME<-paste0("matrixeqtl-transeqtl-",OUTPUTNAMESTRING,"-job",JOBID,".txt")
OUTPUTFILE2NAME<-paste0("matrixeqtl-ciseqtl-",OUTPUTNAMESTRING,"-job",JOBID,".txt")
OUTPUTFILE3NAME<-paste0("matrixeqtl-info-",OUTPUTNAMESTRING,"-job",JOBID,".txt")
OUTPUTFILE1<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE1NAME)
OUTPUTFILE2<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE2NAME)
OUTPUTFILE3<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE3NAME)

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
if(ARGS[8]==TRUE){
  PTHRESHOLD.TRANS<-1e-2
  PTHRESHOLD.CIS<-1e-2
}else{
  PTHRESHOLD.TRANS<-1e-2
  PTHRESHOLD.CIS<-1e-2
}

#Define the covariance matrix for the error term. This parameter is rarely used. 
#If the covariance matrix is a multiple of identity, set it to numeric().
ERRORCOV<-numeric()

#Define cis distance. SNP-gene pairs within this distance are considered local. 
#The distance is measured from the nearest end of the gene. SNPs within a gene are 
#always considered local.
CISDISTANCE<-1e6

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
DATA_EXPRESSION<-SlicedData$new()
DATA_EXPRESSION$fileDelimiter<-"\t"
DATA_EXPRESSION$fileOmitCharacters<-"NA"
DATA_EXPRESSION$fileSkipRows<-1
DATA_EXPRESSION$fileSkipColumns<-1
DATA_EXPRESSION$fileSliceSize<-1000
DATA_EXPRESSION$LoadFile(FILE_EXPRESSION)

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
MATRIXEQTLRESULT<-Matrix_eQTL_main(DATA_GENOTYPE,
                                   DATA_EXPRESSION,
                                   cvrt=DATA_COVARIATES, 
                                   output_file_name=OUTPUTFILE1, 
                                   pvOutputThreshold=PTHRESHOLD.TRANS,
                                   useModel=MODEL, 
                                   errorCovariance=ERRORCOV, 
                                   verbose=TRUE, 
                                   output_file_name.cis=OUTPUTFILE2, 
                                   pvOutputThreshold.cis=PTHRESHOLD.CIS,
                                   snpspos=DATA_SNPLOCATION, 
                                   genepos=DATA_GENELOCATION,
                                   cisDist=CISDISTANCE,
                                   pvalue.hist=FALSE,
                                   noFDRsaveMemory=SKIPFDRCALC)

#Save information in a file.
x<-c(MATRIXEQTLRESULT$time.in.sec[sapply(MATRIXEQTLRESULT$time.in.sec,length)==1],
     MATRIXEQTLRESULT$param[sapply(MATRIXEQTLRESULT$param,length)==1],
     MATRIXEQTLRESULT$all[sapply(MATRIXEQTLRESULT$all,length)==1],
     MATRIXEQTLRESULT$cis[sapply(MATRIXEQTLRESULT$cis,length)==1],
     MATRIXEQTLRESULT$trans[sapply(MATRIXEQTLRESULT$trans,length)==1])

names(x)<-c(paste0("time.in.sec.",names(MATRIXEQTLRESULT$time.in.sec[sapply(MATRIXEQTLRESULT$time.in.sec,length)==1])),
            paste0("param.",names(MATRIXEQTLRESULT$param[sapply(MATRIXEQTLRESULT$param,length)==1])),
            paste0("all.",names(MATRIXEQTLRESULT$all[sapply(MATRIXEQTLRESULT$all,length)==1])),
            paste0("cis.",names(MATRIXEQTLRESULT$cis[sapply(MATRIXEQTLRESULT$cis,length)==1])),
            paste0("trans.",names(MATRIXEQTLRESULT$trans[sapply(MATRIXEQTLRESULT$trans,length)==1])))
INFO<-data.frame(Name=names(x),
                 Value=unlist(x))

write.table(INFO,OUTPUTFILE3,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

############################################################################
#SAVE CONTROL FILES:

cat(paste0("############################################################################
Date: ",RUNDATE,"
Job ID: ",JOBID,"
Script: ",SCRIPTNAME,"
Input genotype file: ",FILE_GENOTYPE,"
Input expression file: ",FILE_EXPRESSION,"
Input covariates file: ",FILE_COVARIATES,"
Input snp location file: ",FILE_SNPLOCATION,"
Input gene location file: ",FILE_GENELOCATION,"
Input FDR request: ",ARGS[8],"
Output file: ",OUTPUTFILE1,"
Output file: ",OUTPUTFILE2,"
Output file: ",OUTPUTFILE3,"

"),file=paste0(OUTPUTLOCATION,"/README.txt"),append=TRUE)

sink()