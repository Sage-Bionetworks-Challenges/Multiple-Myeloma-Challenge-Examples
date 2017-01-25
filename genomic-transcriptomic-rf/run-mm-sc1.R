#!/usr/bin/env Rscript

## R library dependencies
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))

suppressPackageStartupMessages(library("parallel"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
  options(rf.cores = num.cores - 1, 
          mc.cores = num.cores - 1)  ## Cores for parallel processing  
}

## Read in the model
source("R/classifier-base.R")
source("R/rf.R")

seed <- 1234

option_list <- list(
    make_option(c("--metadata-dir"), action="store",
                default=NULL,
                help="Directory holding any required meta data [which are assumed to be named training-expression-data.tsv and training-anntations.tsv]"),
    make_option(c("--test-dir"), action="store",
                default=NULL,
                help="Directory holding the test name [which is assumed to be named test-expression-data.tsv]"),
    make_option(c("--output-dir"), action="store",
                default=NULL,
                help="Dir in which to store results"))

required.args <- c("metadata-dir", "test-dir", "output-dir")

descr <- "\
   Run the random forest, parameterized by the metadata in the metadatadir, on the test data."

parser <- OptionParser(usage = "%prog [options] file.tsv", option_list=option_list, description=descr)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if ( length(arguments$args) != 0 ) {
  print_help(parser)
  q(status=1)
}

missing.args <- required.args[!(required.args %in% names(opt))]
if( length(missing.args) > 0 ) {
    cat("Need to specify arguments:\n")
    sapply(missing.args, function(x) cat(paste("--", x, "\n", sep="")))
    q(status=-1)
}

## Extract the names of the training, test, and output directories.
metadata.dir <- opt$`metadata-dir`
test.dir <- opt$`test-dir`
output.dir <- opt$`output-dir`

## We assume the test and metadata files has a pre-assigned name:
test.expr.file <- paste0(test.dir, "/", "test-expression-data.Rd")
test.genomic.file <- paste0(test.dir, "/", "test-genomic-data.Rd")
test.anno.file <- paste0(test.dir, "/", "test-annotation-data.Rd")

expr.gene.subset.metadata.file <- paste0(metadata.dir, "/", "expr-gene-subset-metadata.Rd")
genomic.gene.subset.metadata.file <- paste0(metadata.dir, "/", "genomic-gene-subset-metadata.Rd")
rf.metadata.file <- paste0(metadata.dir, "/", "rf-metadata.Rd")

for(file in c(test.expr.file, test.genomic.file, test.anno.file, expr.gene.subset.metadata.file, genomic.gene.subset.metadata.file, rf.metadata.file)) {
  if(!file.exists(file)) {
    stop(paste0("Could not find file ", file, "\n"))
  }
}

load(test.expr.file)
load(test.genomic.file)
load(test.anno.file)

## Read in the metadata file  
load(expr.gene.subset.metadata.file)
load(genomic.gene.subset.metadata.file)
load(rf.metadata.file)

## Create the output directory
dir.create(output.dir)

test.data.binarized <- list(eset = test.eset[expr.gene.subset,], clinical = test.anno[,c("Sex", "ISSstage")], genomic = test.genomic.binarized[genomic.gene.subset,])

set.seed(seed)
test.mat.binarized <- assemble.predictor.matrix(test.data.binarized, already.log2.transformed = TRUE)

## Run the random forest
set.seed(seed)
res <- test.random.forest(rf, X = test.mat.binarized)
rf.res <- merge(as.data.frame(res), test.anno, by="ID")

rf.res$high.risk <- rf.res$high.risk == 1

## Output the results with columns ID, raw.score, and high.risk
output.file <- paste0(output.dir, "/", "predictions.tsv")
write.table(file=output.file, rf.res[,c("ID", "raw.score", "high.risk")], row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

sink("test-sessionInfo.txt")
sessionInfo()
sink()
