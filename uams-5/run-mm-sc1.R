#!/usr/bin/env Rscript

## R library dependencies
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))

## Read in the model
source("R/classifier-base.R")
source("R/uams-5.R")

option_list <- list(
    make_option(c("--test-dir"), action="store",
                default=NULL,
                help="Directory holding the test name [which is assumed to be named test-expression-data.tsv]"),
    make_option(c("--output-dir"), action="store",
                default=NULL,
                help="Dir in which to store results"))

required.args <- c("test-dir", "output-dir")

descr <- "\
   Run the uams-5 model, parameterized by the metadata in the metadatadir, on the test data."

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
test.dir <- opt$`test-dir`
output.dir <- opt$`output-dir`

## We assume the test files has a pre-assigned name:
test.expr.file <- paste0(test.dir, "/", "test-expression-data.tsv")

for(file in c(test.expr.file)) {
  if(!file.exists(file)) {
    stop(paste0("Could not find file ", file, "\n"))
  }
}

## Create the output directory
dir.create(output.dir)

## A little nonsense because fread does not seem to deal well with rows
read.file.with.row.names <- function(file) {
  data <- fread(file, header=TRUE, fill=TRUE)
  data <- as.data.frame(data)
  rownames(data) <- data[,1]
  cols <- colnames(data)[1:(ncol(data)-1)]
  data <- data[,-1]
  colnames(data) <- cols
  data
}

## Read in the metadata file
metadata.file <- paste0("modelstate.tsv")
obj <- read.table(metadata.file, sep="\t", header=TRUE, as.is=TRUE)
threshold <- as.numeric(obj$threshold)

## Read in the test expression data
test.expr.data <- read.file.with.row.names(test.expr.file)

## Run the model on the uams-5 model
test.res <- uams.5.entrez(X = as.matrix(test.expr.data), y = NULL, threshold = threshold)

## Output the results
output.file <- paste0(output.dir, "/", "predictions.tsv")
write.table(file=output.file, test.res$res, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

sink("sessionInfo.txt")
sessionInfo()
sink()