minmodelsolve <-
function(ssx=indf,preps){
  modelproportions<-NULL
  #somewhat awkwardly, this assumes that the data consists only of 1 replicate of each component/mix.
  #i suppose i could do some pre-summarizing inside here?  Better: use a 'summarize' subroutine,elsewhere (in normalization?) to do that.

  for(idx in 1:length(preps$mixnames)){
    #more generalized...
    form<-paste0(preps$mixnames[idx],"~")
    for(cindex in (1:(length(preps$componentnames)))){
       form<-paste0(form,"I(",preps$componentnames[cindex],"*",preps$mfrac[cindex],")+")
      }
    form<-paste0(form,0)

    modelproportions<-rbind(modelproportions,c(
      coefficients(lm(data=ssx,form))/
        sum(coefficients(lm(data=ssx,form))),
      preps$mixnames[idx]))
  }
  colnames(modelproportions)<-c(preps$componentnames,"mix")
  modelproportions<-as.data.frame(modelproportions) #converts proportions to numeric rather than character.
  for(I in 1:length(preps$componentnames)){modelproportions[,I]<-as.numeric(as.character(modelproportions[,I]))}
  return(modelproportions)
}
