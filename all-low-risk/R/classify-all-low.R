library(argparse)

testdir = '/test-data'
outputdir = '/output'

message("Reading in clinical file")

clinical_file <- read.csv(file=file.path(testdir,'sc1_Validation_ClinAnnotations.csv'),header=T,as.is=T,check.names=F)


set.seed(1337)

output <- data.frame(study=clinical_file$Study,
  patient=clinical_file$Patient,
  predictionscore=sample(length(clinical_file$Patient)),
  highriskflag=FALSE)

message("Write out prediction file")

write.table(file=file.path(outputdir,"predictions.tsv"),output,row.names=FALSE,sep="\t",quote=FALSE)
