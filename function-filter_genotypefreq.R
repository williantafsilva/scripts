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
  library(ggplot2)
  library(tidyverse)
  library(gtools)
  
  #Transform phased genotypes into unphased genotypes.
  DATA_GENOTYPES<-data.frame(lapply(DATA_GENOTYPES,function(x){
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
                               TT=rowSums(DATA_GENOTYPES=="T/T"))
  GENOTYPES_COUNTS$N_Samples<-rowSums(GENOTYPES_COUNTS[,2:11])
  GENOTYPES_COUNTS$N_Genotypes<-rowSums(GENOTYPES_COUNTS[,2:11]>0)
  if(MINGENFREQ<1){
    GENOTYPES_COUNTS$N_Genotypes_ltMINGENFREQ<-rowSums(GENOTYPES_COUNTS[,2:11]>0 & (GENOTYPES_COUNTS[,2:11]/GENOTYPES_COUNTS$N_Samples)<MINGENFREQ)
    #SNPs that have at least one genotype with frequency <MINGENFREQ.
    SNPltMINGENFREQ<-GENOTYPES_COUNTS$SNP[which(rowSums(GENOTYPES_COUNTS[,2:11]>0 & (GENOTYPES_COUNTS[,2:11]/GENOTYPES_COUNTS$N_Samples)<MINGENFREQ)>0)]
  }else{
    GENOTYPES_COUNTS$N_Genotypes_ltMINGENFREQ<-rowSums(GENOTYPES_COUNTS[,2:11]>0 & GENOTYPES_COUNTS[,2:11]<MINGENFREQ)
    #SNPs that have at least one genotype with frequency <MINGENFREQ.
    SNPltMINGENFREQ<-GENOTYPES_COUNTS$SNP[which(rowSums(GENOTYPES_COUNTS[,2:11]>0 & GENOTYPES_COUNTS[,2:11]<MINGENFREQ)>0)]
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
  
  