############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

intervals_union_intersection<-function(STARTVECTOR=c(100,500,800,150,550,900,200,600,950,250,850,300,500,950), #Vector of start positions.
                                       ENDVECTOR=c(200,700,1000,250,650,1100,300,850,1150,350,1000,400,750,1050), #Vector of end positions.
                                       MINNOVERLAPS=1, #Minimum number of interval overlaps (1=self-overlap).
                                       GAP=0, #Gaps with size <GAP between intervals with be treated as interval overlaps.
                                       GROUPVECTOR){ #Vector of groups for which intervals are separately processed.
  #Merge intervals that are present in at least MINNOVERLAPS groups.
  #Union (at least one groups): MINNOVERLAPS=1
  #Intersection (at least two groups): MINNOVERLAPS=2
  #Intersection (at least three groups): MINNOVERLAPS=3
  
  SETS<-function(START,END,MINNOVERLAPS){
    dd<-rbind(data.frame(pos=START,event=1),data.frame(pos=END,event=-1))
    dd<-aggregate(event~pos,dd,sum)
    dd<-dd[order(dd$pos),]
    dd$open<-cumsum(dd$event)
    r<-rle(dd$open>=MINNOVERLAPS)
    ex<-cumsum(r$lengths-1+rep(1,length(r$lengths))) 
    sx<-ex-r$lengths+1
    cbind(dd$pos[sx[r$values]],dd$pos[ex[r$values]+1])
  } 
  
  if(missing(GROUPVECTOR)){ 
    #If groups are not provided, treat intervals as belonging to the same group.
    DF<-data.frame(START=STARTVECTOR-GAP/2,
                   END=ENDVECTOR+GAP/2)
    
    OUT<-as.data.frame(with(DF,SETS(START=START,END=END,MINNOVERLAPS=MINNOVERLAPS)))
    colnames(OUT)<-c("START","END")
    if(MINNOVERLAPS<=1){
      OUT$START<-OUT$START+GAP/2
      OUT$END<-OUT$END-GAP/2
    }
  }else{
    #If groups are provided, process intervals separately for each group.
    DF<-data.frame(GROUP=GROUPVECTOR,
                   START=STARTVECTOR-GAP/2,
                   END=ENDVECTOR+GAP/2)
    OUT<-data.frame(GROUP=NA,
                    START=NA,
                    END=NA)
    for(G in unique(DF$GROUP)){
      TMP<-as.data.frame(with(DF[DF$GROUP==G,c("START","END")],SETS(START=START,END=END,MINNOVERLAPS=MINNOVERLAPS)))
      if(nrow(TMP)>0){
        colnames(TMP)<-c("START","END")
        TMP$GROUP<-G
        if(MINNOVERLAPS<=1){
          TMP$START<-TMP$START+GAP/2
          TMP$END<-TMP$END-GAP/2
        }
        OUT<-rbind(OUT,TMP)
      }
    }
    OUT<-OUT[!is.na(OUT$GROUP),]
  }

  return(OUT)
  
}

#DF<-data.frame(GROUP=c("A","A","A","B","B","B","C","C","C","D","D","E","E","E"),
#               START=c(100,500,800,150,550,900,200,600,950,250,850,300,500,950),
#               END=c(200,700,1000,250,650,1100,300,850,1150,350,1000,400,750,1050))
#plot(range(DF[,c(2,3)]),c(1,nrow(DF)),type="n",xlab="",ylab="",yaxt="n")
#for(i in 1:nrow(DF)){
#  lines(c(DF[i,2],DF[i,3]),rep(nrow(DF)-i+1,2),col=as.numeric(as.factor(DF$GROUP[i])),lwd=2)
#}
#intervals_union_intersection(STARTVECTOR=DF$START,
#                             ENDVECTOR=DF$END,
#                             MINNOVERLAPS=1,
#                             GAP=0)

