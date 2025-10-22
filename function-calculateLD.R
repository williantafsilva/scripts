############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

#Calculate linkage disequilibrium measurements.
calculateLD<-function(MATRIX_GENOTYPES=data.frame(SNP=c("1:18928702","1:18928916","1:18929231","1:18929569","1:189301351","1:189301455","1:189302827","1:189303295","1:189304348","1:189304652"),
                                                  SAMPLE1=sample(c("A/A","A/T","T/T"),10,replace=TRUE),
                                                  SAMPLE2=sample(c("A/A","A/G","G/G"),10,replace=TRUE),
                                                  SAMPLE3=sample(c("C/C","C/T","T/T"),10,replace=TRUE),
                                                  SAMPLE4=sample(c("G/G","G/T","T/T"),10,replace=TRUE),
                                                  SAMPLE5=sample(c("A/A","A/T","T/T"),10,replace=TRUE),
                                                  SAMPLE6=sample(c("C/C","C/T","T/T"),10,replace=TRUE),
                                                  SAMPLE7=sample(c("C/C","C/G","G/G"),10,replace=TRUE),
                                                  SAMPLE8=sample(c("A/A","A/G","G/G"),10,replace=TRUE),
                                                  SAMPLE9=sample(c("A/A","A/T","T/T"),10,replace=TRUE),
                                                  SAMPLE10=sample(c("A/A","A/C","C/C"),10,replace=TRUE)), #Data frame with biallelic genotypes (first column is SNP ID, other columns are samples).
                      SNP1="1:18929231",
                      SNP2="1:189303295",
                      GENOTYPEFORMAT="TGT"){ #Format of genotypes: TGT (translated genotypes, e.g., "A/T", "G|T"), GT (reference/alternative allele format, e.g., "0/1", "1/1"), NUM (numeric format, e.g., 0, 1, 2).
  
  #Load libraries.
  library(genetics)
  
  #Extract genotypes (excluding SNP_ID column).
  GENOTYPES_SNP1<-as.character(unlist(MATRIX_GENOTYPES[MATRIX_GENOTYPES[,1]==SNP1,-1]))
  GENOTYPES_SNP2<-as.character(unlist(MATRIX_GENOTYPES[MATRIX_GENOTYPES[,1]==SNP2,-1]))
  
  #Reformat genotypes.
  if(GENOTYPEFORMAT=="TGT"){
    #Normalize genotype separators.
    GENOTYPES_SNP1<-unlist(lapply(GENOTYPES_SNP1,function(x){
      x<-gsub("\\|","/",x) #Convert "|" to "/".
      x<-gsub("//","/",x) #Clean double slashes, if any.
      trimws(x)}))
    GENOTYPES_SNP2<-unlist(lapply(GENOTYPES_SNP2,function(x){
      x<-gsub("\\|","/",x)
      x<-gsub("//","/",x)
      trimws(x)}))
  }
  if(GENOTYPEFORMAT=="GT"){
    #Convert GT format (0/0, 0/1, 1/1) to TGT format (A/A, A/B, B/B).
    GENOTYPES_SNP1<-unlist(lapply(GENOTYPES_SNP1,function(x){
      x<-gsub("\\|","/",x) #Convert "|" to "/".
      x<-gsub("//","/",x) #Clean double slashes, if any.
      x<-gsub("0","A",x) #Substitute 0 for A.
      x<-gsub("1","B",x) #Substitute 1 for B.
      trimws(x)}))
    GENOTYPES_SNP2<-unlist(lapply(GENOTYPES_SNP2,function(x){
      x<-gsub("\\|","/",x) #Convert "|" to "/".
      x<-gsub("//","/",x) #Clean double slashes, if any.
      x<-gsub("0","A",x) #Substitute 0 for A.
      x<-gsub("1","B",x) #Substitute 1 for B.
      trimws(x)}))
  }
  if(GENOTYPEFORMAT=="NUM"){
    #Convert NUM format (0,1,2) to TGT format (A/A, A/B, B/B).
    GENOTYPES_SNP1<-unlist(lapply(GENOTYPES_SNP1,function(x){
      ifelse(x==0,"A/A",
             ifelse(x==1,"A/B","B/B"))}))
    GENOTYPES_SNP2<-unlist(lapply(GENOTYPES_SNP2,function(x){
      ifelse(x==0,"A/A",
             ifelse(x==1,"A/B","B/B"))}))
  }
  
  #Convert to "genotype class (genetics package expects format like "A/B").
  GENOTYPES_SNP1<-as.genotype(GENOTYPES_SNP1)
  GENOTYPES_SNP2<-as.genotype(GENOTYPES_SNP2)
  
  #Compute LD (D', r², etc.).
  RESULT_LD<-LD(GENOTYPES_SNP1,GENOTYPES_SNP2)
  
  OUT<-list(D=RESULT_LD$D,
            D_prime=RESULT_LD$`D'`,
            r=RESULT_LD$r,
            R_squared=RESULT_LD$`R^2`,
            N=RESULT_LD$n,
            CHI_squared=RESULT_LD$`X^2`,
            Pvalue=RESULT_LD$`P-value`)
  
  return(OUT)
  
}

#calculateLD()
