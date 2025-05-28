############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

#Merge intervals that overlap or have a gap of size GAP or smaller between them.
intervals_mergeoverlaps<-function(STARTVECTOR,ENDVECTOR,GAP=0){
  
  library(tidyverse)

  DATA<-data.frame(Start=STARTVECTOR-GAP/2,
                   End=ENDVECTOR+GAP/2)
  DATA<-DATA[order(DATA$Start,DATA$End),]
  MERGED<-DATA %>% 
    arrange(Start) %>% 
    group_by(Index=cumsum(cummax(lag(End,default=first(End)))<Start)) %>% 
    summarise(Start=first(Start),End=max(End)) %>%
    as.data.frame
  MERGED$Start<-MERGED$Start+GAP/2
  MERGED$End<-MERGED$End-GAP/2
  
  return(MERGED[,c("Start","End")])
  
}

#intervals_mergeoverlaps(STARTVECTOR=c(1,10,20,30,40,50),
#                         ENDVECTOR=c(5,15,32,35,51,55),
#                         GAP=0)