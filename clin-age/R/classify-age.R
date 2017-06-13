library(argparse)

testdir = '/test-data'
outputdir = '/output'

message("Reading in clinical file")

clinical_file <- read.csv(file=file.path(testdir,'sc1_Validation_ClinAnnotations.csv'),header=T,as.is=T,check.names=F)

# Add a rudimentary imputation step
clinical_file$D_Age <- as.numeric(clinical_file$D_Age)
clinical_file$D_Age[is.na(clinical_file$D_Age)] <- median(as.numeric(clinical_file$D_Age),na.rm=T);

output <- data.frame(study=clinical_file$Study,
  patient=clinical_file$Patient,
  predictionscore=clinical_file$D_Age,
  highriskflag=clinical_file$D_Age>65)

message("Write out prediction file")

write.table(file=file.path(outputdir,"predictions.tsv"),output,row.names=FALSE,sep="\t",quote=FALSE)
