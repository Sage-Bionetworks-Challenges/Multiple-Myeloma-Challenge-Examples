#!/usr/bin/env Rscript

## R library dependencies
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))

## Read in the model
source("R/classifier-base.R")
source("R/uams-5.R")

option_list <- list(
    make_option(c("--training-dir"), action="store",
              default=NULL,
              help="Directory holding the training data [which are assumed to be named training-annotation.tsv and training-expression-data.tsv]"),
    make_option(c("--metadata-dir"), action="store",
                default=NULL,
                help="Directory to which we will output the metadata [in a file metadata.tsv]")
)


required.args <- c("metadata-dir", "training-dir")

descr <- "\
   Train the uams-5 classification model on training data and save the state in the metadata-dir."

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
training.dir <- opt$`training-dir`
metadata.dir <- opt$`metadata-dir`

## We assume the training and test files have pre-assigned names:
training.expr.file <- paste0(training.dir, "/", "training-expression-data.tsv")
training.anno.file <- paste0(training.dir, "/", "training-annotations.tsv")

for(file in c(training.expr.file, training.anno.file)) {
  if(!file.exists(file)) {
    stop(paste0("Could not find file ", file, "\n"))
  }
}

## Create the metadata directory
dir.create(metadata.dir)

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

## Read in the training expression and annotation data
training.expr.data <- read.file.with.row.names(training.expr.file)
training.anno <- read.table(training.anno.file, sep="\t", header=TRUE)

## Subset to samples in both the expression and annotation files
## and but the samples in the same order in both data frames
inters <- intersect(colnames(training.expr.data), rownames(training.anno))
training.expr.data <- training.expr.data[, inters]
training.anno <- training.anno[inters, ]

## Train the uams-5 model (i.e., get the optimal threshold)
training.res <- uams.5.entrez(X = as.matrix(training.expr.data), y = training.anno, threshold = NULL)

metadata.file <- paste0(metadata.dir, "/", "metadata.tsv")
write.table(file=metadata.file, training.res$obj, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

## Get the optimized/trained threshold
threshold <- as.numeric(training.res$obj$threshold)

