############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

peak_range<-function(CHRVECTOR, #Vector of chromosomes.
                     POSVECTOR, #Vector of chromosomal positions.
                     PVALUEVECTOR, #Vector of p-values.
                     SIGTHRESHOLD=0.01, #P-value significance threshold.
                     SEARCHDISTANCE=1e+6, #Maximum distance from last selected adjacent position.
                     LOCALPEAKTHRESHOLD=0.01){ #Threshold of local peak (orders of magnitude)
  
  library(tidyverse)
  library(data.table)
  
  if(missing(CHRVECTOR)){
    CHRVECTOR<-rep("ChromosomeC",length(POSVECTOR))
  }
  
  INPUTDATA<-data.frame(Chromosome=CHRVECTOR,
                        Position=POSVECTOR,
                        Value=PVALUEVECTOR)
  
  OUTPUT<-data.frame(Chromosome=NA,
                     RegionStart=NA,
                     RegionEnd=NA)
  
  for(C in unique(CHRVECTOR)){
    
    REMAININGDATA<-data.frame(Position=INPUTDATA$Position[INPUTDATA$Chromosome==C],
                              Value=INPUTDATA$Value[INPUTDATA$Chromosome==C]) %>% arrange(Position)
    
    #Continue searching while there are positions that haven't been assigned to a range.
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
  return(OUTPUT)
  
}

DF<-data.frame(Chromosome=sample(1:3,10000,replace=TRUE),
               Position=sample(1:1e+8,10000,replace=FALSE),
               Pvalue=c(sample(seq(from=1e-12,to=1e-8,length.out=10000),5000,replace=FALSE),
                        sample(seq(from=0.01,to=0.1,length.out=10000),5000,replace=FALSE))) %>%
  arrange(Chromosome,Position)

plot_manhattan(CHRVECTOR=DF$Chromosome, #Vector of chromosomes.
               POSVECTOR=DF$Position, #Vector of positions.
               #CHRSET=c(1:39,"Z","W"), #Set of chromosomes to show in the plots.
               VALUES=DF$Pvalue, #Vector of y values (p-values, likelihood values, iHS values, etc).
               #REGIONS, #Data frame (Chromosome,RegionStart,RegionEnd) with chromosomal regions to be shaded in per chromosome plots.
               ABOVETHRESHOLD=Inf, #Threshold above which points are plotted with alpha=1. If both ABOVETHRESHOLD and BELOWTHRESHOLD are set to Inf and -Inf, all points are plotted with alpha=1.
               BELOWTHRESHOLD=-Inf, #Threshold below which points are plotted with alpha=1.
               COLORS=c("deepskyblue","darkblue"), #Alternating colors for adjacent chromosomes.
               STDCHRLENGTH=FALSE, #Use standardized chromosome length?
               PLOTPERCHR=TRUE, #Plot per chromosome? (This options never uses standardized chromosome length).
               XLABEL="Chromosome", #X-axis label.
               YLABEL="P-value", #Y-axis label.
               PLOTTITLE="Manhattan plot")

QTLREGIONS<-peak_range(CHRVECTOR=DF$Chromosome, #Vector of chromosomes.
                       POSVECTOR=DF$Position, #Vector of chromosomal positions.
                       PVALUEVECTOR=DF$Pvalue, #Vector of p-values.
                       SIGTHRESHOLD=0.01, #P-value significance threshold.
                       SEARCHDISTANCE=1e+6, #Maximum distance from last selected adjacent position.
                       LOCALPEAKTHRESHOLD=0.01)
QTLREGIONS$RegionLength<-QTLREGIONS$RegionEnd-QTLREGIONS$RegionStart

max(QTLREGIONS$RegionLength)
