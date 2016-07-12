minmodelsolve <-
function(ssx,preps){
  modelproportions<-NULL
  #somewhat awkwardly, this assumes that the data consists only of 1 replicate of each component/mix.
  #use a 'summarize' subroutine,elsewhere (collapsereps) to do that.

  for(idx in 2:length(as.character(preps$mixnames))){
    #more generalized...
    form<-paste0(as.character(preps$mixnames[idx]),"~")
    for(cindex in (2:(length(as.character(preps$componentnames))))){
       form<-paste0(form,"I(",as.character(preps$componentnames[cindex]),"*",preps$mfrac[cindex-1],")+")
      } #builds the model formula from the componentnames & mRNA fraction data
    form<-paste0(form,0) #removes the intercept term & closes the formula

    modelproportions<-rbind(modelproportions,c(
      coefficients(lm(data=ssx,form))/
        sum(coefficients(lm(data=ssx,form))),
      as.character(preps$mixnames[idx]))) #solves the model coefficients and normalizes to 1 to make it proportions.
  }
  colnames(modelproportions)<-c(as.character(preps$componentnames)[2:length(preps$componentnames)],"mix")
  modelproportions<-as.data.frame(modelproportions)
  for(I in 1:(length(preps$componentnames)-1)){modelproportions[,I]<-as.numeric(as.character(modelproportions[,I]))}#converts proportions to numeric rather than character.
  return(modelproportions)
}
