#!/usr/bin/env Rscript

## This script:
## (1) sets up the environment for the docker image -- i.e., partitions data into test and training data
##     and puts in format to be used by run-mm-sc1.R
## (2) trains a random forest
## (3) stores the metadata/state of the trained random forest to be used by run-mm-sc1.R

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

source("R/emc-92.R")
source("R/uams-17.R")
source("R/uams-70.R")

seed <- 1234

option_list <- list(
    make_option(c("--metadata-dir"), action="store",
                default=NULL,
                help="Directory holding any required meta data"),
    make_option(c("--training-dir"), action="store",
              default=NULL,
              help="Directory holding the training data"),
    make_option(c("--test-dir"), action="store",
                default=NULL,
                help="Directory holding the test name"))

required.args <- c("metadata-dir", "test-dir", "training-dir")

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

## Extract the names of the training, test, and metadata directories.
metadata.dir <- opt$`metadata-dir`
test.dir <- opt$`test-dir`
training.dir <- opt$`training-dir`

use.synapse <- TRUE

if(use.synapse) {
  library(synapseClient)
  synapseLogin()
}

# This follows http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# counts has genes as rows and samples as columns
get.variable.genes <- function(counts) {
  library(DESeq); library(statmod); library(pcaMethods); library(fastICA)
  require(DESeq)
  lib.size <- estimateSizeFactorsForMatrix(counts)
  ed <- t(t(counts)/lib.size)
  means <- rowMeans(ed)
  vars <- apply(ed,1,var)
  cv2 <- vars/means^2
  pdf("variable.pdf")
  par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
  smoothScatter(log(means),log(cv2))
  
  require(statmod)
  minMeanForFit <- unname( quantile( means[ which( cv2 > .3 ) ], .95 ) )
  minMeanForFit <- 0.001
  useForFit <- means >= minMeanForFit # & spikeins
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  
  print(c(a0, a1))
  # repeat previous plot
  par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2));
  xg <- exp(seq( min(log(means[means>0])), max(log(means[means>0])), length.out=1000 ))
  vfit <- a1/xg + a0
  # add fit line
  lines( log(xg), log(vfit), col="black", lwd=3 )
  df <- ncol(ed) - 1
  # add confidence interval
  lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
  lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")  
  
  afit <- a1/means+a0
  varFitRatio <- vars/(afit*means^2)
  varorder <- order(varFitRatio,decreasing=T)
  oed <- ed[varorder,]

  # repeat previous plot
  par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
  # add top 100 genes

  pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
  adj.pval <- p.adjust(pval,"fdr")
  sigVariedGenes <- adj.pval<1e-3;
  points(log(means[sigVariedGenes]),log(cv2[sigVariedGenes]),col=2)
  print(table(sigVariedGenes))  
  
  d <- dev.off()
  return(names(means)[sigVariedGenes])
}

## Get the non-synonymous variants
ns.variants.data <- NULL
ns.variants.file <- "MMRF_CoMMpass_IA8b_All_Canonical_NS_Variants.txt"
if(!file.exists(ns.variants.file)) {
  if(use.synapse) {
    obj <- synGet(id="syn7450323", downloadFile = TRUE, downloadLocation = ".")
    gz.file <- getFileLocation(obj)
    system(paste0("gunzip ", gz.file))
  }
}
ns.variants.data <- read.table(ns.variants.file, sep="\t", header=TRUE, comment.char="")

data <- NULL
# Get the count data
data.file <- "MMRF_CoMMpass_IA8b_E74GTF_Salmon_Gene_Counts.txt"
if(!file.exists(data.file)) {
  if(use.synapse) {
    obj <- synGet(id="syn7450336", downloadFile = TRUE, downloadLocation = ".")
    gz.file <- getFileLocation(obj)
    system(paste0("gunzip ", gz.file))
  }
}

data <- read.table(data.file, sep="\t", header=TRUE)
# Let's take a single column per patient, favoring "*_1_BM" over "*_1_PB" over "*_2_BM"
# Do this by sorting the columns by name and taking the first non-duplicated.
expr.data <- data[,!(colnames(data) %in% "GENE_ID")]
rownames(expr.data) <- data$GENE_ID

expr.data <- expr.data[,order(colnames(expr.data))]
cols <- colnames(expr.data)
cols <- gsub(cols, pattern="_\\d+_BM", replacement="")
cols <- gsub(cols, pattern="_\\d+_PB", replacement="")

expr.data <- expr.data[,!duplicated(cols, fromLast=FALSE)]
        
colnames(expr.data) <- gsub(colnames(expr.data), pattern="_\\d+_BM", replacement="")
colnames(expr.data) <- gsub(colnames(expr.data), pattern="_\\d+_PB", replacement="")

# Limit the variants to those that are called by at least 2 of 3 callers
flag <- apply(ns.variants.data[,c("SEURAT","MUTECT","STRELKA")], 1, function(row) length(which(row=="true")) >= 2)
ns.variants.data <- ns.variants.data[flag,]

# Let's take a single column per patient, favoring "*_1_BM" over "*_1_PB" over "*_2_BM"
# Do this by sorting the columns by name and taking the first non-duplicated.

ns.variants.data$harmonized.sample <- ns.variants.data$Sample
ns.variants.data$harmonized.sample <- gsub(ns.variants.data$harmonized.sample, pattern="_\\d+_BM", replacement="")
ns.variants.data$harmonized.sample <- gsub(ns.variants.data$harmonized.sample, pattern="_\\d+_PB", replacement="")

## Sort mutations so that those from the patient are ordered with _1_BM first, then _2_BM, _1_PB, etc.
## Then just keep those that are not duplicated--those will be those that are _1_BM or are unique to _2_BM, etc.

## Drop entries that are not mapped to an ENSG gene
## ensg.gene.mutations <- unlist(lapply(ns.variants.data$EFF....GENE, function(x) grepl(x=x, pattern="ENSG")))
## ns.variants.data <- ns.variants.data[ensg.gene.mutations,]

## It appears that are occassionally multiple entries for the same gene/mutation (e.g., based on different)
## transcripts.  Subset to the rows that are unique with respect to the fields we care about
mut.cols <- c("Sample", "X.CHROM", "POS", "REF", "ALT", "GEN.1..AR", "EFF....GENE")
dup <- duplicated(ns.variants.data[,mut.cols])
ns.variants.data <- ns.variants.data[!dup,]

o <- order(do.call(order, ns.variants.data[,mut.cols]))
mut.cols.no.smpl <- c("X.CHROM", "POS", "REF", "ALT", "GEN.1..AR", "EFF....GENE")
## dup <- duplicated(ns.variants.data[,mut.cols.no.smpl], fromLast=TRUE) | duplicated(ns.variants.data[,mut.cols.no.smpl], fromLast=FALSE)
## head(ns.variants.data[dup,])

dup <- duplicated(ns.variants.data[,mut.cols.no.smpl], fromLast=FALSE)
ns.variants.data <- ns.variants.data[!dup,]

## Create a row (ensg gene id) x col (sample) matrix, where the entries are the entries are the highest VAF
## for any mutation in that gene in that sample.
ns.variants.df <- ns.variants.data[,c("harmonized.sample", "EFF....GENE", "GEN.1..AR")]
ns.variants.df <- ns.variants.data[order(ns.variants.data$GEN.1..AR, decreasing=TRUE),]
ns.variants.df <- ns.variants.df[!duplicated(ns.variants.df[,c("harmonized.sample", "EFF....GENE")], fromLast = FALSE),]

suppressPackageStartupMessages(library(reshape2))
gene.by.sample.mutations <- acast(data=ns.variants.df, formula=EFF....GENE ~ harmonized.sample, value.var="GEN.1..AR", fill = 0)

# Pull in the annotation data.  
anno.file <- NULL
if(use.synapse) {
  obj <- synGet(id="syn7488088", downloadFile = TRUE, downloadLocation = ".")
  anno.file <- getFileLocation(obj)
} else {
  anno.file <- "MMRF patient-level clinical and cytogenetic inventory Harmonized.csv"
}

anno.df <- read.table(anno.file, sep=",", header=TRUE)

# This must have columns event, time, and ID.  
# Use PFS annotations to define high risk by setting event = Progression; time = PFStimeMonths; ID = SampleID
anno.df <- anno.df[,c("SampleID", "PFStimeMonths", "Progression", "PFSbasedHighRisk", "Sex", "ISSstage")]
colnames(anno.df) <- c("ID", "time", "event", "PFSbasedHighRisk", "Sex", "ISSstage")

# Exclude any samples that do not have a PFS-based high risk
## cat(paste0("Excluding ", length(which(is.na(anno.df$PFSbasedHighRisk))), " of ", nrow(anno.df), " samples without low/high risk\n"))
anno.df <- anno.df[!is.na(anno.df$PFSbasedHighRisk),]

orig.anno.samples <- nrow(anno.df)
orig.expr.samples <- ncol(expr.data)
orig.genomic.samples <- ncol(gene.by.sample.mutations)
rownames(anno.df) <- anno.df$ID

inters <- intersect(rownames(anno.df), colnames(expr.data))

# Subset the data and annotations to include only samples appearing in both.
anno.df <- anno.df[inters,]
expr.data <- expr.data[,inters]

# Further subset to include only those with genomic data
inters <- intersect(colnames(gene.by.sample.mutations), colnames(expr.data))
anno.df <- anno.df[inters,]
expr.data <- expr.data[, inters]
gene.by.sample.mutations <- gene.by.sample.mutations[, inters]

# Give the genomic (gene) features names that are distinct from the expression (gene) features
rownames(gene.by.sample.mutations) <- unlist(lapply(rownames(gene.by.sample.mutations), function(x) paste0("mut",x)))

# Normalize the data using edgeR
suppressPackageStartupMessages(library(edgeR))

y <- DGEList(counts=expr.data)

## Filter non-expressed genes (unless they are involved in a classifer):
coefficients <- c(emc.92.get.coefficients(), uams.17.get.coefficients(), uams.70.get.coefficients())

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
probesets <- names(coefficients)
bm <- getBM(attributes=c('affy_hg_u133_plus_2', 'ensembl_gene_id'), 
            filters = 'affy_hg_u133_plus_2', 
            values = probesets, 
            mart = ensembl)
names(bm) <- c("PROBE", "GENE")
coefficient.genes <- bm$GENE

A <- aveLogCPM(y)
y <- y[(A>1) | (rownames(y$counts) %in% coefficient.genes),]

# Then normalize and compute log2 counts-per-million with an offset:
# Gordyn uses prior.count = 5 in the above link--I use the default of 0.5
y <- calcNormFactors(y)
expr.tbl <- cpm(y, log=TRUE, prior.count=0.5)
anno.df$ID <- as.character(anno.df$ID)
expr.tbl <- expr.tbl[,anno.df$ID]
eset <- ExpressionSet(as.matrix(expr.tbl))

# Only test those genes that are mutated in some sample
flag <- unlist(apply(gene.by.sample.mutations, 1, function(row) any(row>0)))
genomic.gene.subset <- rownames(gene.by.sample.mutations)[flag]

## Only test genes that are mutated at a rate > 2%
flag <- unlist(apply(gene.by.sample.mutations, 1, function(row) ( length(which(row>0))*1.0/length(row) > 0.02 )))
genomic.gene.subset <- rownames(gene.by.sample.mutations)[flag]
empirically.highly.mutated.genes <- unlist(rownames(gene.by.sample.mutations)[flag], function(x) gsub(x=x, pattern="mut", replacement=""))

translocation.targets <- c("FGFR3", "MMSET", "WHSC1", "CCND3", "CCND1", "MAF", "MAFB")

highly.mutated.genes <- c("CDKN2C", "RB1", "CCND1", "CDKN2A", "NRAS", "KRAS", "BRAF", "MYC", "TRAF3", "CYLD", "DKK1", "FRZB", "DNAH5", "XBP1", "BLIMP1", "PRDM1", "IRF4", "TP53", "MRE11A", "PARP1", "DIS3", "FAM46C", "LRRK2", "KDM6A", "UTX", "MLL", "MMSET", "WHSC1", "HOXA9", "KDM6B", "IGLL5")

sym.to.ensg <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
                     filters = 'hgnc_symbol', 
                     values = c(translocation.targets, highly.mutated.genes, empirically.highly.mutated.genes), 
                     mart = ensembl)
names(sym.to.ensg) <- c("SYMBOL", "GENE")

## When looking at the expression data, check genes that are highly mutated/translocated or are used by one
## of the other classifiers, and are expressed
expr.gene.subset <- intersect(rownames(eset), sym.to.ensg$GENE)
expr.gene.subset <- c(expr.gene.subset, coefficient.genes)
expr.gene.subset <- c(expr.gene.subset, get.variable.genes(expr.data))
expr.gene.subset <- intersect(expr.gene.subset, rownames(eset))
cat(paste0("Number of expression features: ", length(expr.gene.subset), "\n"))

## Let's exclude IGK, IGH, IGLV (but not IGGL5)
exclude.patterns <- c("IGK", "IGH", "IGLV", "MUC")
exclude <- rep(FALSE, length(genomic.gene.subset))
for(pattern in exclude.patterns) {
  exclude <- exclude | grepl(x=genomic.gene.subset, pattern=pattern)
}
genomic.gene.subset <- genomic.gene.subset[!exclude]

genomic.gene.subset <- intersect(genomic.gene.subset, rownames(gene.by.sample.mutations))

# Finally, make sure we don't have any NAs in any of our features.  If we do, drop the corresponding sample.
na.samples <- c()
flag <- unlist(apply(anno.df[,c("Sex", "ISSstage")], 1, function(row) any(is.na(row))))
if(any(flag)) {
  na.samples <- c(na.samples, rownames(anno.df)[flag])
}

flag <- unlist(apply(gene.by.sample.mutations[genomic.gene.subset,], 2, function(col) any(is.na(col))))
if(any(flag)) {
  na.samples <- c(na.samples, colnames(gene.by.sample.mutations)[flag])
}

flag <- unlist(apply(exprs(eset)[expr.gene.subset,], 2, function(col) any(is.na(col))))
if(any(flag)) {
  na.samples <- c(na.samples, colnames(eset)[flag])
}

if(length(na.samples) != 0) {
  ## cat("Dropping samples that were NA in some feature of interest: ", paste0(na.samples, collapse=","), "\n", sep="")
  anno.df <- anno.df[!(rownames(anno.df) %in% na.samples),]
  eset <- eset[,!(colnames(eset) %in% na.samples)]
  gene.by.sample.mutations <- gene.by.sample.mutations[,!(colnames(gene.by.sample.mutations) %in% na.samples)]
}

# Break the MMRF data into training and test data sets, stratifying based on risk.
# Set the seed so that our partition is reproducible
suppressPackageStartupMessages(library(caret))
set.seed(seed)
training.fraction <- 0.7
train_ind <- as.vector(createDataPartition(factor(anno.df$PFSbasedHighRisk), p=training.fraction, list = FALSE))

train.eset <- eset[,train_ind]
train.anno <- anno.df[train_ind,]
test.eset <- eset[,-train_ind]
test.anno <- anno.df[-train_ind,]

if(any(colnames(exprs(train.eset)) != train.anno$ID)) {
  warning("Training data and annotations are not properly aligned\n")
}

if(any(colnames(exprs(test.eset)) != test.anno$ID)) {
  warning("Test data and annotations are not properly aligned\n")
}

rownames(train.anno) <- train.anno$ID
rownames(test.anno) <- test.anno$ID

train.eset <- eset[,train_ind]
train.anno <- anno.df[train_ind,]
train.genomic <- gene.by.sample.mutations[,train_ind]

test.eset <- eset[,-train_ind]
test.anno <- anno.df[-train_ind,]
test.genomic <- gene.by.sample.mutations[,-train_ind]

## Binarize the genomic data
train.genomic.binarized <- train.genomic
test.genomic.binarized <- test.genomic
train.genomic.binarized[train.genomic.binarized > 0] <- 1
test.genomic.binarized[test.genomic.binarized > 0] <- 1

response <- "PFSbasedHighRisk"
train.y <- train.anno[,response]
names(train.y) <- train.anno$ID

## Save the test genomic, expression, and clinical data

## Also, save the metadata expr.gene.subset and genomic.gene.subset

## We assume the test files has a pre-assigned name:
dir.create(test.dir)
test.expr.file <- paste0(test.dir, "/", "test-expression-data.Rd")
test.genomic.file <- paste0(test.dir, "/", "test-genomic-data.Rd")
test.anno.file <- paste0(test.dir, "/", "test-annotation-data.Rd")

save(test.eset, file=test.expr.file)
save(test.genomic.binarized, file=test.genomic.file)
save(test.anno, file=test.anno.file)

dir.create(training.dir)
train.expr.file <- paste0(training.dir, "/", "train-expression-data.Rd")
train.genomic.file <- paste0(training.dir, "/", "train-genomic-data.Rd")
train.anno.file <- paste0(training.dir, "/", "train-annotation-data.Rd")

save(train.eset, file=train.expr.file)
save(train.genomic.binarized, file=train.genomic.file)
save(train.anno, file=train.anno.file)

## Run the random forest with ISS+mut+expr the binarized data
train.data.binarized <- list(eset = train.eset[expr.gene.subset,], clinical = train.anno[,c("Sex", "ISSstage")], genomic = train.genomic.binarized[genomic.gene.subset,])

set.seed(seed)
train.mat.binarized <- assemble.predictor.matrix(train.data.binarized, already.log2.transformed = TRUE)

set.seed(seed)
rf <- train.random.forest(X = train.mat.binarized, y = train.y)

## Save the metadata files
dir.create(metadata.dir)
expr.gene.subset.metadata.file <- paste0(metadata.dir, "/", "expr-gene-subset-metadata.Rd")
genomic.gene.subset.metadata.file <- paste0(metadata.dir, "/", "genomic-gene-subset-metadata.Rd")
rf.metadata.file <- paste0(metadata.dir, "/", "rf-metadata.Rd")

save(expr.gene.subset, file=expr.gene.subset.metadata.file)
save(genomic.gene.subset, file=genomic.gene.subset.metadata.file)
save(rf, file=rf.metadata.file)

sink("training-sessionInfo.txt")
sessionInfo()
sink()
