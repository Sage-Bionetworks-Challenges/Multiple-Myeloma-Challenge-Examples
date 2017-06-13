#!/usr/bin/env Rscript

## R library dependencies.  These must be installed in the Docker image.
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Matrix.utils"))
suppressPackageStartupMessages(library("plyr"))

use.sequestered.data <- FALSE

if(use.sequestered.data) {
  suppressPackageStartupMessages(library("synapseClient"))
  synapseLogin()
}

## This R script assumes that it runs in the root directory of the Docker image and that
## the test data are mounted in ./test-data,
## the output should be written to ./output,
## and that the entire directory structure of the submitted Docker image
## (e.g., an R object encapsulating trained modeler state) is mounted at ./
docker.image.base <- "./"
test.dir <- "./test-data/"
output.dir <- "./output/"

## Read in the model implementation.
source(paste0(docker.image.base, "R/classifier-base.R"))
source(paste0(docker.image.base, "R/uams-70.R"))

## Read in the trained model state.
model.state.metadata.file <- paste0(docker.image.base, "model-state-metadata.Rd")
load(model.state.metadata.file)
## model.state is assumed defined in model.state.metadata.file.

## It is assumed to be a list with two entries:
## "threshold"--the threshold for our classifier that was optimized on training data.
## "mapping"--the probe to entry mapping
threshold <- unname(model.state[["threshold"]])
mapping <- model.state[["mapping"]]

## This is a kludge to run this without Docker.  Download only the files we are
## interested in from the sequester site (at synId = syn9763945)
chal_data_df <- NULL

if(use.sequestered.data) {
  chal_data_table <- synTableQuery('select id,name from syn9763945')
  chal_data_df <- chal_data_table@values

  ## Download the annotation file
  chal_data_df_anno <- subset(chal_data_df, grepl(name, pattern="sc2_Validation_ClinAnnotations.csv"))
  sapply(chal_data_df_anno$id, function(x) { synGet(x, downloadLocation=test.dir)})
}

## Read in the annotation file
training.validation.file <- paste0(test.dir, "sc2_Validation_ClinAnnotations.csv")
anno.tbl <- read.table(training.validation.file, sep=",", header=TRUE)

## Read in the clinical annotation file that lists the samples and what data
## are available for each in the columns:
data.set.cols <- colnames(anno.tbl)[grepl(pattern="File", colnames(anno.tbl)) & !grepl(pattern="Sampl", colnames(anno.tbl))]

## The mapping between the Patient/row of the clinical annotation file and
## the identifier in each data file is provided by the corresponding ID
## columns of the clinical annotation file.
## These have the same names as the data set columns, but with "SamplId"
## appended.
data.set.patient.id.cols <- unlist(lapply(data.set.cols,
                                          function(str) paste0(str, "SamplId")))

## Restrict to gene-level expression data sets
expression.data.set.flag <- grepl(pattern="geneLevelExp", data.set.cols)
data.set.cols <- data.set.cols[expression.data.set.flag]
data.set.patient.id.cols <- data.set.patient.id.cols[expression.data.set.flag]

if(use.sequestered.data) {
  ## Subset to the files we are interested in.
  all.data.sets <- unique(na.omit(as.vector(as.matrix(anno.tbl[,data.set.cols]))))
  chal_data_df <- subset(chal_data_df, name %in% all.data.sets)
  ## Download the data and put it in the base directory of the docker image,
  ## which is where it would be mounted in the challenge.
  sapply(chal_data_df$id, function(x) { synGet(x, downloadLocation=test.dir)})
}

## Read in each of the data sets and create:
## 1) A list holding each of the data sets.
## 2) A list with an entry for each data set holding a mapping from the
##    sample name in the clinical annotation column and the sample name
##    in the expression data set.
data.sets <- list()
sample.name.mappings <- list()

for(col.indx in 1:length(data.set.cols)) {
    for(data.set in na.omit(unique(anno.tbl[,data.set.cols[col.indx]]))) {
        cat(paste0("Reading in data.set ", data.set, "\n"))
        file <- paste0(test.dir, data.set)
        tbl <- as.data.frame(fread(file, header=TRUE))

        ## The first column is the gene identifier--make it the row
        ## name and drop it from the table
        rownames(tbl) <- tbl[,1]
        tbl <- tbl[,-1]
        data.sets[[data.set]] <- tbl

        ## Extract the sample name mappings from the annotation file
        patient.id.col <- data.set.patient.id.cols[col.indx]
        flag <- !is.na(anno.tbl[,data.set.cols[col.indx]]) & (anno.tbl[,data.set.cols[col.indx]] == data.set)
        map <- anno.tbl[flag,c("Patient", patient.id.col)]
        colnames(map) <- c("patient.id", "sample.id")
        
        sample.name.mappings[[data.set]] <- map
    }
}

## Combine multiple sample columns into a single patient column by
## taking the average of the sample columns
## Columns (i.e., samples) are assumed to have names in sample.to.patient.map$from.
## They will be renamed as sample.to.patient.map$to.
combine_samples_2_patient <- function(expr, sample.to.patient.map) {
    map <- subset(sample.to.patient.map, from %in% colnames(expr))
    rownames(map) <- map$from
    expr <- expr[, colnames(expr) %in% map$from]
    map <- map[colnames(expr),]
    expr <- as.matrix(t(aggregate.Matrix(t(expr), groupings=list(map$to), fun = "mean")))
    expr
}

## Patients listed in the annotation file may have 0, 1, or >1 corresponding
## samples in the data file.  If there are multiple, take the mean across
## samples (NB: since most of these are log expression values, this corresponds
## to the geometric mean in real space).
## Also, subset both the mapping tables and the data to have the
## intersection of patients in both.
for(data.set in names(sample.name.mappings)) {
    tbl <- data.sets[[data.set]]
    map <- sample.name.mappings[[data.set]]

    map <- na.omit(map)

    ## If the mapping from patient to sample id is 1-to-many, split it into multiple
    ## rows, each of which maps the patient to a single sample id.
    map <- ldply(1:nrow(map),
                 .fun = function(i) {
                     sample.id <- unlist(strsplit(as.character(map$sample.id[i]), split=";[ ]*"))
                     df <- data.frame(patient.id = map$patient.id[i], sample.id = sample.id)
                     df
                 })
    
    map <- map[,c("patient.id", "sample.id")]
    colnames(map) <- c("to", "from")
    tbl <- combine_samples_2_patient(tbl, map)

    ## At this point, we have renamed and combined the sample ids columns to patient ids
    ## so that the map is just the identity
    map <- data.frame(patient.id = colnames(tbl), sample.id = colnames(tbl))
    
    data.sets[[data.set]] <- tbl
    sample.name.mappings[[data.set]] <- map
}

map <- do.call("rbind", sample.name.mappings)
rownames(map) <- NULL
map <- unique(map)
if(any(duplicated(map$patient.id))) {
    warning("Was not expecting any patient IDs to be duplicated\n")
}

## The sample IDs are scattered across multiple columns, put them
## in the single "sample.id" column
anno.tbl <- merge(anno.tbl, map, by.x = "Patient", by.y = "patient.id", all = FALSE)

## Subset the expression data sets to only have the genes in the UAMS-70
## data set (and that are in common between all data sets).

## Get the Entrez genes in the model.  This is defined in R/uams-70.R
## model.entrezs <- uams.70.probe.to.entrez.mapping()
## However, we don't have access to the internet (required of biomart in uams.70.probe.to.entrez.mapping).
## We have stored this mapping in the metadata
model.entrezs <- mapping
model.entrezs <- na.omit(model.entrezs$INDEX)

all.genes <- lapply(data.sets, rownames)

common.genes <- Reduce(intersect, all.genes)
common.genes <- intersect(common.genes, model.entrezs)

## Subset all of the data sets
for(data.set in names(data.sets)) {
    data.sets[[data.set]] <- data.sets[[data.set]][as.character(common.genes),]
}

## Z-score each gene, separately for each data set.
## NB: scale normalizes the _columns_ of the matrix
for(data.set in names(data.sets)) {
    data.sets[[data.set]] <- t(scale(t(data.sets[[data.set]]), center = TRUE, scale = TRUE))
}

## Combine all of the expression matrices into a single matrix and
## ensure they are ordered consistently with the survival data.
expr <- do.call("cbind", data.sets)
expr <- expr[, anno.tbl$sample.id]

## Run the (trained) model by passing in the trained/optimized threshold
## The test phase does not provide access to the internet (which would be required for biomart
## as invoked during uams.70.probe.to.entrez.mapping).  So, instead call that during the training
## phase (where we do have access to the internet) and store in the metadata.
## uams.70.entrez <- function(X, y, threshold = NULL) {
##  bm <- uams.70.probe.to.entrez.mapping()
##  eset <- ExpressionSet(as.matrix(X))
##  coefficients <- uams.70.get.coefficients()
##  uams.70.gene(eset, mapping = bm, coefficients = coefficients, anno = y, threshold = threshold)
## }
X <- as.matrix(expr)
y <- NULL
eset <- ExpressionSet(as.matrix(X))
coefficients <- uams.70.get.coefficients()
res <- uams.70.gene(eset, mapping = mapping, coefficients = coefficients, anno = y, threshold = threshold)


## Create the output directory
dir.create(output.dir)

## Output the results with columns ID, raw.score, and high.risk.
output.file <- paste0(output.dir, "/", "predictions.tsv")
predictions.tbl <- res$res[,c("ID", "raw.score", "high.risk")]
colnames(predictions.tbl) <- c("patient", "predictionscore", "highriskflag")

## Merge in the study name, which is required in the prediction output.
predictions.tbl <- merge(predictions.tbl, anno.tbl[,c("Patient", "Study")], by.x = "patient", by.y = "Patient")
predictions.tbl <- predictions.tbl[,c("Study", "patient", "predictionscore", "highriskflag")]
colnames(predictions.tbl) <- c("study", "patient", "predictionscore", "highriskflag")

write.table(file=output.file, predictions.tbl, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

cat("Successfully wrote predictions.\n")

