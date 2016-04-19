doubleconsensus <-
function(dataf,preps){
  fitmodel<-function(indf,preps){
    fit<-switch(preps$modeltype,
           twomix=minmodelsolve(indf,preps),
           threemix=magicmodelsolve(indf,preps))
    return(fit)
  }
  f2opt<-function(par,indata,preps){prepsn<-c(preps,list(mfrac=c(par,1)));(sum(abs(fitmodel(indata,preps=prepsn)[1,1:length(preps$trueproportions[,1])]-preps$trueproportions[,1]),
                                                     abs(fitmodel(indata,preps=prepsn)[2,1:length(preps$trueproportions[,2])]-preps$trueproportions[,2])))}
  opted<-optim(par=c(1,1),f2opt,indata=dataf,preps=preps)$par
  return(c(opted,1)/sum(c(opted,1)))
}
