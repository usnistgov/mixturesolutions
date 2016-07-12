doubleconsensus <-
function(dataf,preps){
  fitmodel<-function(indf,preps){
    fit<-switch(eval(preps$modeltype),
           twomix=minmodelsolve(indf,preps),
           threemix=magicmodelsolve(indf,preps))
    return(fit)
  }
  f2opt<-function(parx,indata,preps){prepsn<-c(preps,list(mfrac=c(parx,1)));(sum(abs(fitmodel(indata,preps=prepsn)[1,1:length(eval(preps$trueproportions)[,1])]-eval(preps$trueproportions)[,1]),
                                                     abs(fitmodel(indata,preps=prepsn)[2,1:length(eval(preps$trueproportions)[,2])]-eval(preps$trueproportions)[,2])))}
  opted<-optim(par=c(1,1),fn=f2opt,indata=dataf,preps=eval(preps))$par
  return(c(opted,1)/sum(c(opted,1)))
}
