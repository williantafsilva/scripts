############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

#Find the element of a vector with the closest value.
closest_index<-function(myvalue,myvector){
  index<-which(abs(myvalue-myvector)==min(abs(myvalue-myvector)))
  return(index)
}