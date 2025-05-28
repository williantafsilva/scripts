############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

#This function takes a vector of values and the index of a target value X and finds the range
#of indices around the index of the target value where values fall within the interval 
#[threshold*X,X], with at most ngaps gaps (values that don't fall within the interval) on each
#side of X.
peaktothreshold_indices<-function(values,
                                  targetindex,
                                  threshold=0.5,
                                  ngaps=2){
  L<-targetindex
  R<-targetindex
  SEARCH_LEFT<-TRUE
  SEARCH_RIGHT<-TRUE
  NGAPS_LEFT<-0 
  NGAPS_RIGHT<-0
  REGION<-c(targetindex)
  TRUTHVECTOR<-(values>=(threshold*values[targetindex]))
  while(SEARCH_LEFT==TRUE & SEARCH_RIGHT==TRUE){
    L<-L-1
    R<-R+1
    if(SEARCH_LEFT==TRUE){
      if(L>=1){
        if(TRUTHVECTOR[L]==TRUE){
          REGION<-c(REGION,L)
        }else{
          if(NGAPS_LEFT<=ngaps){
            NGAPS_LEFT<-NGAPS_LEFT+1
          }else{
            SEARCH_LEFT<-FALSE
          }
        }
      }else{
        SEARCH_LEFT<-FALSE
      }
    }
    if(SEARCH_RIGHT==TRUE){
      if(R<=length(values)){
        if(TRUTHVECTOR[R]==TRUE){
          REGION<-c(REGION,R)
        }else{
          if(NGAPS_RIGHT<=ngaps){
            NGAPS_RIGHT<-NGAPS_RIGHT+1
          }else{
            SEARCH_RIGHT<-FALSE
          }
        }
      }else{
        SEARCH_RIGHT<-FALSE
      }
    }
  }
  REGION_INDICES<-sort(REGION)
  return(REGION_INDICES)
}

#DATA<-data.frame(Position=1:100,
#                 Value=rnorm(100,mean=0,sd=1))
#PEAK<-which(DATA$Value==max(DATA$Value))
#THRESHOLD<-0.25
#region_indices<-peaktothreshold_indices(values=DATA$Value,targetindex=PEAK,threshold=THRESHOLD,ngaps=3)
#plot(DATA$Position,DATA$Value,type="l",col="black",lwd=1)
#abline(h=THRESHOLD*max(DATA$Value),lty=2)
#if(length(region_indices)>1){
#  lines(DATA$Position[region_indices],DATA$Value[region_indices],col="blue",lwd=2)
#}else{
#  points(DATA$Position[region_indices],DATA$Value[region_indices],col="blue",pch=16)
#}

peaktothreshold_multiple<-function(values=c(sample(1:100,97,replace=TRUE),c(5,50,80)),
                                   targetvalues=c(5,50,80),
                                   threshold=0.5,
                                   ngaps=2){
  
  OUTPUT<-data.frame(TargetValue=NA,
                     TargetValueIndex=NA,
                     RegionStartIndex=NA,
                     RegionEndIndex=NA)
  
  for(V in targetvalues){
    for(I in which(values==V)){
      L<-I
      R<-I
      SEARCH_LEFT<-TRUE
      SEARCH_RIGHT<-TRUE
      NGAPS_LEFT<-0 
      NGAPS_RIGHT<-0
      REGION<-c(I)
      TRUTHVECTOR<-(values>=(threshold*values[I]))
      while(SEARCH_LEFT==TRUE & SEARCH_RIGHT==TRUE){
        L<-L-1
        R<-R+1
        if(SEARCH_LEFT==TRUE){
          if(L>=1){
            if(TRUTHVECTOR[L]==TRUE){
              REGION<-c(REGION,L)
            }else{
              if(NGAPS_LEFT<=ngaps){
                NGAPS_LEFT<-NGAPS_LEFT+1
              }else{
                SEARCH_LEFT<-FALSE
              }
            }
          }else{
            SEARCH_LEFT<-FALSE
          }
        }
        if(SEARCH_RIGHT==TRUE){
          if(R<=length(values)){
            if(TRUTHVECTOR[R]==TRUE){
              REGION<-c(REGION,R)
            }else{
              if(NGAPS_RIGHT<=ngaps){
                NGAPS_RIGHT<-NGAPS_RIGHT+1
              }else{
                SEARCH_RIGHT<-FALSE
              }
            }
          }else{
            SEARCH_RIGHT<-FALSE
          }
        }
      }
      TMP<-data.frame(TargetValue=V,
                      TargetValueIndex=I,
                      RegionStartIndex=min(REGION),
                      RegionEndIndex=max(REGION))
      OUTPUT<-rbind(OUTPUT,TMP)
    }
  }
  OUTPUT<-OUTPUT[!is.na(OUTPUT$TargetValue),]
  return(OUTPUT)
}
peaktothreshold_multiple()

#DATA<-data.frame(Position=1:100,
#                 Value=rnorm(100,mean=0,sd=1))
#PEAK<-which(DATA$Value==max(DATA$Value))
#REGION_INDICES<-peaktothreshold_multiple(values=DATA$Value,targetvalues=sample(DATA$Value,3),threshold=0.5,ngaps=1)
#REGION_INDICES
#plot(DATA$Position,DATA$Value,type="l",col="black",lwd=1)
#for(i in 1:nrow(REGION_INDICES)){
#  if(REGION_INDICES$RegionStartIndex[i]==REGION_INDICES$RegionEndIndex[i]){
#    points(DATA$Position[REGION_INDICES$RegionStartIndex[i]],
#           DATA$Value[REGION_INDICES$RegionStartIndex[i]],
#           col="blue",pch=16)
#  }else{
#    lines(DATA$Position[REGION_INDICES$RegionStartIndex[i]:REGION_INDICES$RegionEndIndex[i]],
#          DATA$Value[REGION_INDICES$RegionStartIndex[i]:REGION_INDICES$RegionEndIndex[i]],
#          col="blue",lwd=2)
#  }
#}

peaktothreshold_window<-function(values=runif(100),
                                 windowsize=10,
                                 threshold=0.5,
                                 useaverage=FALSE,
                                 ngaps=5){
 
  TMP1<-data.frame(RegionStartIndex=NA,
                   RegionEndIndex=NA)
  
  if(useaverage==TRUE){
    for(S in 1:(length(values)-windowsize)){
      if(mean(values[S:(S+windowsize)])>=threshold){
        TMP2<-data.frame(RegionStartIndex=S,
                         RegionEndIndex=S+windowsize)
        TMP1<-rbind(TMP1,TMP2)
      }
    }
  }else{
    TRUTHVECTOR<-(values>=threshold)
    for(S in 1:(length(values)-windowsize)){
      if(sum(TRUTHVECTOR[S:(S+windowsize)])>(windowsize-ngaps)){
        TMP2<-data.frame(RegionStartIndex=S,
                         RegionEndIndex=S+windowsize)
        TMP1<-rbind(TMP1,TMP2)
      }
    }
  }
  
  TMP1<-TMP1[!is.na(TMP1$RegionStartIndex),]
  if(nrow(TMP1)>1){
    OUTPUT<-TMP1[1,]
    for(i in 2:nrow(TMP1)){
      if(TMP1$RegionStartIndex[i]<=TMP1$RegionEndIndex[i-1]){
        OUTPUT$RegionEndIndex[nrow(OUTPUT)]<-TMP1$RegionEndIndex[i]
      }else{
        OUTPUT<-rbind(OUTPUT,TMP1[i,])
      }
    }
  }
  
  return(OUTPUT)
}

#DATA<-data.frame(Position=1:100,
#                 Value=runif(100))
#REGION_INDICES<-peaktothreshold_window(values=DATA$Value,windowsize=10,threshold=0.5,useaverage=FALSE,ngaps=2)
#REGION_INDICES
#plot(DATA$Position,DATA$Value,type="l",col="black",lwd=1)
#abline(h=0.5,lty=2)
#if(sum(!is.na(REGION_INDICES$RegionStartIndex))>0){
#  for(i in 1:nrow(REGION_INDICES)){
#    lines(c(DATA$Position[REGION_INDICES$RegionStartIndex[i]],DATA$Position[REGION_INDICES$RegionEndIndex[i]]),
#          rep(mean(DATA$Value[REGION_INDICES$RegionStartIndex[i]:REGION_INDICES$RegionEndIndex[i]]),2),
#          col="blue",lwd=2) 
#  }
#}
