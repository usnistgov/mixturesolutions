unittests <-
function(){
#global testing:
#Test #1:  Does it work with BLMmixdata.
Blmdata=read.csv("Example_2mix_3compBLM.txt")
preptargetplot(mrnatype = "ercc",prenormalized = TRUE,modeltype = "twomix",indf = Blmdata,trueproportions = data.frame(mix1=c(.25,.25,.5),mix2=c(.25,.5,.25)),componentnames=c("bep","lep","mep"),mixnames=c("a1","a2"))
#expected output:  26/26/48 and 29/52/19 BLM for a1 and a2, respectively.
preptargetplot(mrnatype = "internalconsensus",prenormalized = TRUE,modeltype = "twomix",indf = Blmdata,trueproportions = data.frame(mix1=c(.25,.25,.5),mix2=c(.25,.5,.25)),componentnames=c("bep","lep","mep"),mixnames=c("a1","a2"))
#expected output:  25/25/50 and 28/52/20 BLM for a1 and a2, respectively.

#Test #2  Does it work with miRNAmixData:
mirnadata=read.csv("Example_2mix_3comp.txt",sep="\t")
preptargetplot(indf=mirnadata,mrnatype="internalconsensus",prenormalized=FALSE,modeltype="twomix",componentnames=c("Brain","Liver","Placenta"),mixnames=c("Mix1","Mix2"),trueproportions = data.frame(mix1=c(.25,.25,.5),mix2=c(.5,.25,.25)))
##expected output:  25/25/50 and 50/25/25 for Mix1 and Mix2(+/-.006), respectively.
}
