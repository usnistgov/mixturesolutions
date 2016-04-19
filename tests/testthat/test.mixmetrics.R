
test_that("The example BLMmix works properly", {
  data("BLM")
  fval=mixmetrics(mrnatype = "ercc",prenormalized = TRUE,modeltype = "twomix",indf = BLM,
                  trueproportions = data.frame(mix1=c(.25,.25,.5),mix2=c(.25,.5,.25)),componentnames=c("bep","lep","mep"),mixnames=c("a1","a2"))
  expect_equal_to_reference(fval,"testblm.rds")
  #expected output:  26/26/48 and 29/52/19 BLM for a1 and a2, respectively.
  })

test_that("Does it work with miRNAmixData:",{
#mirnadata=read.csv("data/Example_2mix_3comp.txt",sep="\t")
  data("EDRNmir")
  btest=mixmetrics(indf=mirnadata,mrnatype="internalconsensus",prenormalized=FALSE,modeltype="twomix",componentnames=c("Brain","Liver","Placenta"),
                   mixnames=c("Mix1","Mix2"),trueproportions = data.frame(mix1=c(.25,.25,.5),mix2=c(.5,.25,.25)))
expect_equal_to_reference(btest,"testmiRNA.rds")
##expected output:  25/25/50 and 50/25/25 for Mix1 and Mix2(+/-.006), respectively.
})


test_that("Does it work with SEQCData:",{
  data("SEQCagr")
  ctest=mixmetrics(indf=SEQCagr,mrnatype="ercc",prenormalized=FALSE,modeltype="twomix",componentnames=c("A","B"),mixnames=c("C","D"),trueproportions=data.frame(C=c(0.25,0.75),D=c(0.75,0.25)))
  expect_error()
  #eventually will do an expect_equal_to_reference(ctest,"testSEQC.rds")
})