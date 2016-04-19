magicmodelsolve <-
function(ssx=indf,preps,idvars){
  labround<-NULL
  #awkwardly, this assumes that the data consists only of 1 replicate of each component/mix.

  for(idx in 1:length(preps$mixnames)){
    labround<-rbind(labround,c(
      #Mm=c(signif(coefficients(lm(data=indf,I(mmix*valuemf[2])~I(lmmix*valuemf[3])+I(lmix*valuemf[1])+0))/sum(coefficients(lm(data=indf,I(mmix*valuemf[2])~I(lmmix*valuemf[3])+I(lmix*valuemf[1])+0))),digits=3)[1],0,signif(coefficients(lm(data=indf,I(mmix*valuemf[1])~I(lmmix*valuemf[3])+I(lmix*valuemf[2])+0))/sum(coefficients(lm(data=indf,I(mmix*valuemf[1])~I(lmmix*valuemf[3])+I(lmix*valuemf[2])+0))),digits=3)[2]),

      coefficients(lm(data=ssx,preps$mixnames[idx]~I(preps$componentnames[1]*preps$mfrac[1])+I(preps$componentnames[2]*preps$mfrac[2])+I(preps$componentnames[3]*preps$mfrac[3])+0))/
        sum(coefficients(lm(data=ssx,preps$mixnames[idx]~I(preps$componentnames[1]*preps$mfrac[1])+I(preps$componentnames[2]*preps$mfrac[2])+I(preps$componentnames[3]*preps$mfrac[3])+0))),
      preps$mixnames[idx],idlist))
  } #I could stand to build this out for different numbers of components and mixes, but i'm not sure how.
  return(labround)
}
