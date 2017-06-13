library(vcfR)
library(parallel)

testdir = '/test-data'
outputdir = '/output'

test.classifier <- function(input.file=file.path(testdir,'sc3_Validation_ClinAnnotations.csv'),
mutationburden=64.51,
mutationFileField='WES_mutationFileMutect',output.file=file.path(outputdir,'predictions.tsv')) {
  clinical_file <- read.csv(file=input.file,header=T,as.is=T,check.names=F)

  getMutationCount <- function(X,testdir) {
    tryCatch({
      curr <- read.vcfR(file.path(testdir,X),verbose=FALSE);
      length(which(grepl(pattern='missense',curr@fix[,8]) & getFIX(curr)[,'FILTER'] == 'PASS'))
    },error=function(e){0})
  }

  mb <- unlist(mclapply(clinical_file[[mutationFileField]],getMutationCount,testdir=testdir,mc.cores=detectCores()))

  output <- data.frame(study=clinical_file$Study,
    patient=clinical_file$Patient,
    predictionscore=mb,
    highriskflag=mb > mutationburden)

  write.table(file=output.file,output,row.names=FALSE,sep="\t",quote=FALSE)
}

test.classifier()
