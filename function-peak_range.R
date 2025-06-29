############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

peak_range<-function(CHRVECTOR, #Vector of chromosomes.
                     POSVECTOR, #Vector of chromosomal positions.
                     PVALUEVECTOR, #Vector of p-values.
                     SIGTHRESHOLD=0.01, #P-value significance threshold.
                     SEARCHDISTANCE=1e+6, #Maximum distance between adjacent positions within a region.
                     LOCALPEAKTHRESHOLD=0.01){ #Relative threshold of local peak (0.01 represents 2 orders of magnitude).
  
  library(tidyverse)
  library(data.table)
  
  if(missing(CHRVECTOR)){
    CHRVECTOR<-rep("ChromosomeC",length(POSVECTOR))
  }
  
  INPUTDATA<-data.frame(Chromosome=CHRVECTOR,
                        Position=POSVECTOR,
                        Value=PVALUEVECTOR) 
  
  #Base filter: Filter out non-significant points. 
  #This means that non-significant loci can exist within a significant QTL region.
  #The region will be based on the drop in significance among significant loci and their distance,
  #ignoring non-significant loci. This has several advantages: 1. Lower memory requirements and considerably faster 
  #(no need to import complete data set because only significant points are used); 2. This allows for
  #the presence of gaps with non-significant points between two significant points located near each other,
  #without assigning them different regions.
  INPUTDATA<-INPUTDATA %>% filter(Value<=SIGTHRESHOLD) 
  
  #Create output data frame.
  OUTPUT<-data.frame(Chromosome=NA,
                     RegionStart=NA,
                     RegionEnd=NA)
  
  for(C in unique(CHRVECTOR)){
    
    REMAININGDATA<-data.frame(Position=INPUTDATA$Position[INPUTDATA$Chromosome==C],
                              Value=INPUTDATA$Value[INPUTDATA$Chromosome==C]) %>% arrange(Position)
    
    #Continue searching while there are significant positions that haven't been assigned to a range.
    while(nrow(REMAININGDATA)>0){ 
      
      #Find peak value.
      PEAKVALUE<-min(REMAININGDATA$Value) 
      
      #No meaningful peaks left.
      if(PEAKVALUE>SIGTHRESHOLD){break} 
      
      #Local peak dynamic threshold.
      X<-PEAKVALUE/LOCALPEAKTHRESHOLD 
      
      #Find all peak indices.
      PEAKINDICES<-which(REMAININGDATA$Value==PEAKVALUE) 
      
      #Find ranges around each peak.
      RANGES_INDICES<-lapply(PEAKINDICES,function(idx){
        PEAKPOSITION<-REMAININGDATA$Position[idx]
        #Search upstream peak.
        LEFT<-idx
        while(LEFT>1 && REMAININGDATA$Value[LEFT-1]<X && REMAININGDATA$Position[LEFT-1]>=(REMAININGDATA$Position[LEFT]-SEARCHDISTANCE)){
          LEFT<-LEFT-1
        }
        #Search downstream peak.
        RIGHT<-idx
        while(RIGHT<nrow(REMAININGDATA) && REMAININGDATA$Value[RIGHT+1]<X && REMAININGDATA$Position[RIGHT+1]<=(REMAININGDATA$Position[RIGHT]+SEARCHDISTANCE)){
          RIGHT<-RIGHT+1
        }
        return(c(LEFT,RIGHT))
      })
      RANGES_INDICES<-do.call(rbind,RANGES_INDICES)
      RANGES_INDICES<-as.data.frame(RANGES_INDICES)
      colnames(RANGES_INDICES)<-c("Start","End")
      
      #Merge overlapping/adjacent ranges.
      RANGES_INDICES<-RANGES_INDICES[order(RANGES_INDICES$Start,RANGES_INDICES$End),]
      MERGEDRANGES_INDICES<-RANGES_INDICES %>% 
        arrange(Start) %>% 
        group_by(Index=cumsum(cummax(lag(End,default=data.table::first(End)))<Start)) %>% 
        summarise(Start=data.table::first(Start),End=max(End)) %>%
        as.data.frame
      MERGEDRANGES_INDICES<-MERGEDRANGES_INDICES[,c("Start","End")]
      
      #Save ranges.
      for(R in nrow(MERGEDRANGES_INDICES)){
        RANGEPOSITIONS<-REMAININGDATA$Position[MERGEDRANGES_INDICES$Start[R]:MERGEDRANGES_INDICES$End[R]]
        TMP<-data.frame(Chromosome=C,
                        RegionStart=min(RANGEPOSITIONS),
                        RegionEnd=max(RANGEPOSITIONS))
        OUTPUT<-rbind(OUTPUT,TMP)
      }
      
      #Remove obtained range from the data.
      REMOVEINDICES<-c()
      for(i in nrow(MERGEDRANGES_INDICES)){
        REMOVEINDICES<-c(REMOVEINDICES,MERGEDRANGES_INDICES$Start[i]:MERGEDRANGES_INDICES$End[i])
      }
      REMAININGDATA<-REMAININGDATA[-REMOVEINDICES,,drop=FALSE]
      
    }
    
  }
  
  OUTPUT<-OUTPUT[!is.na(OUTPUT$RegionStart),]
  OUTPUT$RegionLength<-OUTPUT$RegionEnd-OUTPUT$RegionStart
  
  return(OUTPUT)
  
}

# DF<-data.frame(Chromosome=rep(1,20),
#                Position=c(c(1,10,1000,2000,10000,100000,200000,500000,700000,900000,
#                           1.8e+6,3.3e+6,3.4e+6,3.6e+6,3.8e+6,4.1e+6,4.3e+6,4.35e+6,4.7e+6,6.1e+6),
#                           c(1,10,1000,2000,10000,100000,200000,500000,700000,900000,
#                             1.8e+6,3.3e+6,3.4e+6,3.6e+6,3.8e+6,4.1e+6,4.3e+6,4.35e+6,4.7e+6,6.1e+6)+5),
#                Pvalue=c(c(5e-7,8e-3,1e-4,5e-7,9.9e-7,3e-8,8e-8,1e-8,9e-7,5e-7,1e-5,5.5e-6,1e-6,1e-5,1e-7,9e-7,5e-7,6e-6,9.8e-6,5e-7),
#                         c(runif(20,min=0.011,max=0.5)))) %>%
#                  arrange(Chromosome,Position)
# 
# p<-ggplot(DF)+
#   geom_point(aes(x=Position,y=-log10(Pvalue)),color="blue")+
#   #geom_line(aes(x=Position,y=-log10(Pvalue)),linetype="dashed")+
#   geom_hline(yintercept=-log10(0.01),linetype="dashed")+
#   labs(x="Position",y=expression("-Log"[10]~"Pvalue"))+
#   theme_bw()
# p
# 
# REGIONS<-peak_range(CHRVECTOR=DF$Chromosome, #Vector of chromosomes.
#                     POSVECTOR=DF$Position, #Vector of chromosomal positions.
#                     PVALUEVECTOR=DF$Pvalue, #Vector of p-values.
#                     SIGTHRESHOLD=0.01, #P-value significance threshold.
#                     SEARCHDISTANCE=1e+6, #Maximum distance from last selected adjacent position.
#                     LOCALPEAKTHRESHOLD=0.01) #Relative threshold of local peak (0.01 represents 2 orders of magnitude).
# REGIONS



