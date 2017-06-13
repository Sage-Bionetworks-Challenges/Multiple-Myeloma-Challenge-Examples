#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("synapseClient"))
suppressPackageStartupMessages(library("maxstat"))
suppressPackageStartupMessages(library("Matrix.utils"))
suppressPackageStartupMessages(library("plyr"))

synapseLogin()

## docker build -t docker.synapse.org/syn7237116/uams-5:version1 .
## docker images
## docker login docker.synapse.org
## docker push docker.synapse.org/syn7237116/uams-5:version1

base.dir <- "./"

download.dir <- paste0(base.dir, "download/")
dir.create(download.dir)

## Read in the annotation file, which indicates the files
## holding the training data.
chal_data_table <- synTableQuery('select id,name from syn9763946')
chal_data_df <- chal_data_table@values
chal_data_df_anno <- subset(chal_data_df, grepl(name, pattern="sc2_Training_ClinAnnotations.csv"))
sapply(chal_data_df_anno$id, function(x) { synGet(x, downloadLocation=download.dir)})
training.validation.file <- paste0(download.dir, "sc2_Training_ClinAnnotations.csv")
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

## Only download the training data that is pertinent to this subchallenge.
## Subset to the files we are interested in.
all.data.sets <- unique(na.omit(as.vector(as.matrix(anno.tbl[,data.set.cols]))))
chal_data_df <- subset(chal_data_df, name %in% all.data.sets)
## Download the data and put it in the base directory of the docker image,
## which is where it would be mounted in the challenge.
sapply(chal_data_df$id, function(x) { synGet(x, downloadLocation=download.dir)})

## Read in the model implementation.
source(paste0(base.dir, "R/classifier-base.R"))
source(paste0(base.dir, "R/uams-5.R"))

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
        file <- paste0(download.dir, data.set)
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


## Assemble the survival data.
map <- do.call("rbind", sample.name.mappings)
rownames(map) <- NULL
map <- unique(map)
if(any(duplicated(map$patient.id))) {
    warning("Was not expecting any patient IDs to be duplicated\n")
}

## The sample IDs are scattered across multiple columns, put them
## in the single "sample.id" column
anno.tbl <- merge(anno.tbl, map, by.x = "Patient", by.y = "patient.id", all = FALSE)


## Translate Ensembl genes in the MMRF data set to Entrez.
mmrf.data.set <- names(data.sets)[grepl(names(data.sets), pattern="MMRF")]

translate.ensg.to.entrez <- function(ids) {
  # Get a mapping from Ensembl id to entrez id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('ensembl_gene_id', 'entrezgene'), 
              filters = 'ensembl_gene_id', 
              values = ids, 
              mart = ensembl)
  names(bm) <- c("ID", "TRANSLATION")
  bm <- bm[!(bm$TRANSLATION %in% c("")),]
  bm
}

ensg.to.entrez.map <- translate.ensg.to.entrez(rownames(data.sets[[mmrf.data.set]]))
ensg.to.entrez.map <- na.omit(ensg.to.entrez.map)

## Subset the expression data sets to only have the genes in the UAMS-5
## data set (and that are in common between all data sets).

## Get the Entrez genes in the model.  This is defined in R/uams-5.R
model.entrezs <- uams.5.probe.to.entrez.mapping()
model.entrezs <- na.omit(model.entrezs$INDEX)

all.genes <- lapply(data.sets, rownames)
all.genes[[mmrf.data.set]] <- ensg.to.entrez.map$TRANSLATION[ensg.to.entrez.map$ID %in% all.genes[[mmrf.data.set]]]

common.genes <- Reduce(intersect, all.genes)
common.genes <- intersect(common.genes, model.entrezs)

## Finally, translate the MMRF Ensembl genes to Entrez
ensg.to.entrez.map <- subset(ensg.to.entrez.map, TRANSLATION %in% common.genes)
if(any(duplicated(ensg.to.entrez.map$ID))) {
    warning("Ambiguous mapping between Ensembl ID an Entrez TRANSLATION\n")
}

data.sets[[mmrf.data.set]] <- data.sets[[mmrf.data.set]][rownames(data.sets[[mmrf.data.set]]) %in% ensg.to.entrez.map$ID,]
data.sets[[mmrf.data.set]] <- data.sets[[mmrf.data.set]][ensg.to.entrez.map$ID,]
rownames(data.sets[[mmrf.data.set]]) <- ensg.to.entrez.map$TRANSLATION

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

## Train the UAMS-5 model.
y <- anno.tbl[,c("sample.id", "D_PFS", "D_PFS_FLAG")]
colnames(y) <- c("ID", "time", "event")
rownames(y) <- y$ID
y <- na.omit(y)

res <- uams.5.entrez(X = as.matrix(expr), y = y, threshold = NULL)$obj

## Save the trained state of the model, which consistts of
## (1) "threshold" -- the threshold for discretizing the continuous prediction into high vs low risk.
## (2) "mapping" -- a mapping from probe ids to entrez ids.  We will need this mapping during the
## validation phase where we don't have access to the internet.

threshold <- res$threshold
mapping <- uams.5.probe.to.entrez.mapping()
model.state <- list("threshold" = threshold, "mapping" = mapping)

model.state.metadata.file <- paste0(base.dir, "/", "model-state-metadata.Rd")
save(model.state, file=model.state.metadata.file)

cat("Successfully trained model\n")
