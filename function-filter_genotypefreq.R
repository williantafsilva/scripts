############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

filter_genotypefreq<-function(MATRIX_GENOTYPES, #Data frame of genotypes with SNPs per row (1st column is SNP ID) and samples per column.
                              MINGENFREQ=0, #Minimum genotype frequency (<1) or count (>1).
                              OUTPUTFILTEREDMATRIX=FALSE, #If true, the function will output a list containing, in addition to the SNP IDs of low frequency genotypes, the filtered genotype matrix.
                              OUTPUTLOWFREQGENOTYPES=FALSE, #If true, the function will output a list containing, in addition to the SNP IDs of low frequency genotypes, the genotype matrix of SNPs with low frequency genotypes.
                              OUTPUTGENCOUNT=FALSE){  #If true, the function will output a list containing, in addition to the SNP IDs of low frequency genotypes, the table of genotype counts.
  
  #Load libraries.
  #library(ggplot2)
  library(tidyverse)
  #library(gtools)
  
  #Transform phased genotypes into unphased genotypes.
  DATA_GENOTYPES<-data.frame(lapply(MATRIX_GENOTYPES,function(x){
    x<-gsub("\\|","/",x)
    x<-gsub("C/A","A/C",x)
    x<-gsub("G/A","A/G",x)
    x<-gsub("T/A","A/T",x)
    x<-gsub("G/C","C/G",x)
    x<-gsub("T/C","C/T",x)
    gsub("T/G","G/T",x)}))
  rownames(DATA_GENOTYPES)<-DATA_GENOTYPES[,1]
  
  #Count samples per genotype and remove SNPs with low genotype counts (<MINGENFREQ per genotype).
  GENOTYPES_COUNTS<-data.frame(SNP=DATA_GENOTYPES[,1],
                               AA=rowSums(DATA_GENOTYPES=="A/A"),
                               AC=rowSums(DATA_GENOTYPES=="A/C"),
                               AG=rowSums(DATA_GENOTYPES=="A/G"),
                               AT=rowSums(DATA_GENOTYPES=="A/T"),
                               CC=rowSums(DATA_GENOTYPES=="C/C"),
                               CG=rowSums(DATA_GENOTYPES=="C/G"),
                               CT=rowSums(DATA_GENOTYPES=="C/T"),
                               GG=rowSums(DATA_GENOTYPES=="G/G"),
                               GT=rowSums(DATA_GENOTYPES=="G/T"),
                               TT=rowSums(DATA_GENOTYPES=="T/T"),
                               MISSING=rowSums(DATA_GENOTYPES=="./."))
  GENOTYPES_COUNTS$N_Samples<-rowSums(GENOTYPES_COUNTS[,2:12])
  GENOTYPES_COUNTS$N_Genotypes<-rowSums(GENOTYPES_COUNTS[,2:12]>0)
  if(MINGENFREQ<1){
    GENOTYPES_COUNTS$N_Genotypes_ltMINGENFREQ<-rowSums(GENOTYPES_COUNTS[,2:12]>0 & (GENOTYPES_COUNTS[,2:12]/GENOTYPES_COUNTS$N_Samples)<MINGENFREQ)
    #SNPs that have at least one genotype with frequency <MINGENFREQ.
    SNPltMINGENFREQ<-GENOTYPES_COUNTS$SNP[which(rowSums(GENOTYPES_COUNTS[,2:12]>0 & (GENOTYPES_COUNTS[,2:12]/GENOTYPES_COUNTS$N_Samples)<MINGENFREQ)>0)]
  }else{
    GENOTYPES_COUNTS$N_Genotypes_ltMINGENFREQ<-rowSums(GENOTYPES_COUNTS[,2:12]>0 & GENOTYPES_COUNTS[,2:12]<MINGENFREQ)
    #SNPs that have at least one genotype with frequency <MINGENFREQ.
    SNPltMINGENFREQ<-GENOTYPES_COUNTS$SNP[which(rowSums(GENOTYPES_COUNTS[,2:12]>0 & GENOTYPES_COUNTS[,2:12]<MINGENFREQ)>0)]
  }
  
  OUTPUT<-list(SNPltMINGENFREQ=SNPltMINGENFREQ)
  
  if(OUTPUTFILTEREDMATRIX){
    #Filter out SNPs that have at least one genotype with frequency <MINGENFREQ.
    DATA_FILTEREDGENOTYPES<-DATA_GENOTYPES %>%
      filter(!(DATA_GENOTYPES[,1] %in% SNPltMINGENFREQ))
    OUTPUT<-c(OUTPUT,list(MATRIXFILTERED=DATA_FILTEREDGENOTYPES))
  }
  
  if(OUTPUTLOWFREQGENOTYPES){
    #Filter out SNPs that have at least one genotype with frequency <MINGENFREQ.
    DATA_GENOTYPESLOWFREQ<-DATA_GENOTYPES %>%
      filter(DATA_GENOTYPES[,1] %in% SNPltMINGENFREQ)
    OUTPUT<-c(OUTPUT,list(MATRIXLOWFREQ=DATA_GENOTYPESLOWFREQ))
  }
  
  if(OUTPUTGENCOUNT){
    OUTPUT<-c(OUTPUT,list(COUNTS=GENOTYPES_COUNTS))
  }
  
  return(OUTPUT)
  
} 

filter_genotypefreq_DT<-function(MATRIX_GENOTYPES, #Data frame of genotypes with SNPs per row (1st column is SNP ID) and samples per column.
                                 MINGENFREQ=0, #Minimum genotype frequency (<1) or count (>1).
                                 OUTPUTFILTEREDMATRIX=FALSE, #If true, the function will output a list containing, in addition to the SNP IDs of low frequency genotypes, the filtered genotype matrix.
                                 OUTPUTLOWFREQGENOTYPES=FALSE, #If true, the function will output a list containing, in addition to the SNP IDs of low frequency genotypes, the genotype matrix of SNPs with low frequency genotypes.
                                 OUTPUTGENCOUNT=FALSE){  #If true, the function will output a list containing, in addition to the SNP IDs of low frequency genotypes, the table of genotype counts.
  
  #Load libraries.
  library(tidyverse)
  library(data.table)

  #Make sure MATRIX_GENOTYPES is a data.table object.
  setDT(MATRIX_GENOTYPES)

  #Function to transform phased genotypes into unphased genotypes.
  UNPHASEGENOTYPE<-function(x){
    x<-gsub("\\|","/",x)
    x<-gsub("C/A","A/C",x)
    x<-gsub("G/A","A/G",x)
    x<-gsub("T/A","A/T",x)
    x<-gsub("G/C","C/G",x)
    x<-gsub("T/C","C/T",x)
    gsub("T/G","G/T",x)}

  #Apply function to MATRIX_GENOTYPES.
  MATRIX_GENOTYPES[,(names(MATRIX_GENOTYPES)):=lapply(.SD,UNPHASEGENOTYPE)]
  setDT(MATRIX_GENOTYPES)
  
  #Define the genotypes to be counted.
  GENOTYPES<-c(AA="A/A",AC="A/C",AG="A/G",AT="A/T", 
               CC="C/C",CG="C/G",CT="C/T", 
               GG="G/G",GT="G/T",
               TT="T/T", 
               MISSING="./.")
  
  #Create the counts table.
  GENOTYPES_COUNTS<-data.table(SNP=MATRIX_GENOTYPES[[1]])
  #Loop through and compute rowSums using .SD (excluding the 1st column).
  for(GENNAME in names(GENOTYPES)){
    GENOTYPES_COUNTS[,(GENNAME):=rowSums(MATRIX_GENOTYPES[,.SD,.SDcols=-1]==GENOTYPES[GENNAME])]
  }
  
  #Add the final summary columns by reference.
  GENOTYPECOLS<-names(GENOTYPES)
  GENOTYPES_COUNTS[,`:=`(N_Samples=rowSums(.SD),
                         N_Genotypes=rowSums(.SD>0)),
                   .SDcols=GENOTYPECOLS]
  setDT(GENOTYPES_COUNTS)

  if(MINGENFREQ<1){
    #Calculate row sums in-place using .SD.
    GENOTYPES_COUNTS[,N_Genotypes_ltMINGENFREQ:=rowSums(.SD>0 & (.SD/N_Samples)<MINGENFREQ),.SDcols=2:12]
  }else{
    #Calculate row sums in-place using .SD.
    GENOTYPES_COUNTS[,N_Genotypes_ltMINGENFREQ:=rowSums(.SD>0 & .SD<MINGENFREQ),.SDcols=2:12]
  }
  #SNPs that have at least one genotype with frequency <MINGENFREQ.
  SNPltMINGENFREQ<-GENOTYPES_COUNTS[N_Genotypes_ltMINGENFREQ>0,SNP]

  #Initialize output list.
  OUTPUT<-list(SNPltMINGENFREQ=SNPltMINGENFREQ)
  
  #Set key on the first column of DATA_GENOTYPES for blistering fast filtering
  SNPCOLNAME<-names(MATRIX_GENOTYPES)[1]
  setkeyv(MATRIX_GENOTYPES,SNPCOLNAME)
  
  if(OUTPUTFILTEREDMATRIX){
    #Filter out SNPs that have at least one genotype with frequency <MINGENFREQ.
    DATA_FILTEREDGENOTYPES<-MATRIX_GENOTYPES[!.(SNPltMINGENFREQ)]
    OUTPUT<-c(OUTPUT,list(MATRIXFILTERED=DATA_FILTEREDGENOTYPES))
  }
  if(OUTPUTLOWFREQGENOTYPES){
    #Filter out SNPs that have at least one genotype with frequency <MINGENFREQ.
    DATA_GENOTYPESLOWFREQ<-MATRIX_GENOTYPES[.(SNPltMINGENFREQ)]
    OUTPUT<-c(OUTPUT,list(MATRIXLOWFREQ=DATA_GENOTYPESLOWFREQ))
  }
  if(OUTPUTGENCOUNT){
    OUTPUT<-c(OUTPUT,list(COUNTS=GENOTYPES_COUNTS))
  }
  
  return(OUTPUT)
  
}
  
  