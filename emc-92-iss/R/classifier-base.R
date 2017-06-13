# These are common functions used by the classifiers.

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

probe.to.ensg.mapping <- function(coefficients) {
  # Get a mapping from affy probe to ensembl id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  probesets <- names(coefficients)
  bm <- getBM(attributes=c('affy_hg_u133_plus_2', 'ensembl_gene_id'), 
              filters = 'affy_hg_u133_plus_2', 
              values = probesets, 
              mart = ensembl)
  names(bm) <- c("PROBE", "INDEX")
  bm <- bm[!(bm$INDEX %in% c("")),]
  bm
}

probe.to.entrez.mapping <- function(coefficients) {
  # Get a mapping from affy probe to entrez id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  probesets <- names(coefficients)
  bm <- getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene'), 
              filters = 'affy_hg_u133_plus_2', 
              values = probesets, 
              mart = ensembl)
  names(bm) <- c("PROBE", "INDEX")
  bm <- bm[!(bm$INDEX %in% c("")),]
  bm
}

## Assume data is gene x samples.  This standardizes the genes so that each row has
## mean 0 and variance 1.
standardize.data <- function(data) {
    means <- rowMeans(data)
    sds <- apply(data,1,sd)
    data <- ((data - means) / sds)
    data
}

## clin.res are clinially-annotated results with columns "time", "event", and "score"
optimize.threshold <- function(clin.res) {
  suppressPackageStartupMessages(library(maxstat))
  suppressPackageStartupMessages(library(survival))
  formula <- as.formula(paste("Surv(time, event)", "score", sep=" ~ "))
  mt <- maxstat.test(formula, data=data.frame(clin.res), smethod="LogRank", pmethod="none")
  return(mt$estimate)
}

# mapping is assumed to have columns "PROBE" and "INDEX"
# this returns a data frame with columns "PROBE", "INDEX", and "coefficient"
translate.probe.coefficients.to.gene.coefficients <- function(coefficients, mapping) {
  coeffs.df <- data.frame(coefficient=coefficients, PROBE=names(coefficients))
  classifier.df <- merge(coeffs.df, mapping)
  return(classifier.df)
}

## Adjust the classifier coefficients to account for a probe that maps
## to multiple genes.  If n genes map to a probe, take the average of the
## n gene expression values.  i.e., if the probe-based score is beta,
## make the gene-based score beta/n for each of the n genes.
## classifier.df is a data frame with columns "PROBE", "INDEX", and "coefficient,"
## as returned by translate.probe.coefficients.to.gene.coefficients.
## This returns a vector of coefficients named according to their respective gene.
## There may be multiple coefficients for a single gene.
adjust.scores.to.reflect.multimappers <- function(classifier.df, eset) {
  
  ## Drop any genes from the coefficients if they are not in our data.
  ## (This should have already happened in the merge in 
  ## translate.probe.coefficients.to.gene.coefficients.)
  classifier.df <- classifier.df[classifier.df$INDEX %in% featureNames(eset),]

  ## Now adjust the score to account for multi-mappers
  
  ## Multiplicative factor to account for multi-mappers
  probe.tbl <- as.data.frame(table(classifier.df$PROBE))
  probe.tbl$Freq <- 1.0 / probe.tbl$Freq
  names(probe.tbl) <- c("PROBE", "factor")
  
  classifier.df <- merge(classifier.df, probe.tbl)
  
  classifier.df$factor <- as.numeric(classifier.df$factor)
  classifier.df$coefficient <- as.numeric(classifier.df$coefficient)    
  classifier.df$adj.coefficient <- apply(classifier.df[,c("factor", "coefficient")], 1, function(row) row[1] * row[2])
  gene.coefficients <- classifier.df$adj.coefficient
  names(gene.coefficients) <- classifier.df$INDEX
  return(gene.coefficients)
}

## If threshold is numeric/non-NULL, use the specified threshold.
## If threshold is NULL and train.eset is non-NULL, train the threshold on the training data.
## In this case, train.anno must be non-NULL and is assumed to have columns
## time, event, and ID (for sample).
## threshold and train.eset should not both be NULL
run.classifier.eset <- function(test.eset, coefficients, test.anno = NULL, already.log2.transformed = FALSE, train.eset = NULL, train.anno = NULL, threshold = NULL, coefficient.prefactors = NULL, standardize = TRUE) {
  if(is.null(threshold) && is.null(train.eset)) { 
    stop("threshold and train.eset should not both be NULL\n")
    return(NULL)
  }

  ## If threshold is NULL and train.eset is non-NULL, train the threshold on the training data.
  if(is.null(threshold)) {
    if(is.null(train.anno)) {
      warning("Training requires that train.anno be specified\n")
      return(NULL)
    }
            
    ## Use the training data to optimize the threshold
    train.data <- exprs(train.eset[intersect(names(coefficients),featureNames(train.eset)),])
    if(!already.log2.transformed) {
      train.data <- log2(train.data + 0.5)
    }

    ## train.data is a gene x sample matrix
    if(standardize) {
      train.data <- standardize.data(train.data)
    }
    train.data <- t(train.data)
    ## data is now a standardized sample x gene sample

    available <- match(names(coefficients),featureNames(train.eset))
    if(any(is.na(available))) {
      warning(paste('The following probesets are missing from training data: ',paste(names(coefficients)[is.na(available)],collapse=', ')))
    }
    available <- na.exclude(available)
    available <- featureNames(train.eset)[available]
  
    raw.score.train <- train.data[,available] %*% coefficients[available]
    score.df <- data.frame(ID=sampleNames(train.eset), score=raw.score.train)
    res <- merge(train.anno, score.df, by="ID")

    threshold <- optimize.threshold(res)
  }

  data <- exprs(test.eset[intersect(names(coefficients),featureNames(test.eset)),])
  
  if(!already.log2.transformed) {
    data <- log2(data + 0.5)
  }

  ## data is a gene x sample matrix
  if(standardize) {
    data <- standardize.data(data)
  }
  data <- t(data)
  ## data is now a standardized sample x gene sample

  available <- match(names(coefficients),featureNames(test.eset))
  if(all(is.na(available))) {
    warning("No probesets could be translated")
    return(NULL)
  }
  if(any(is.na(available))) {
    warning(paste('The following probesets are missing: ',paste(names(coefficients)[is.na(available)],collapse=', ')))
  }
  available <- na.exclude(available)
  available <- featureNames(test.eset)[available]
  
  raw.score <- data[,available] %*% coefficients[available]
  cat(paste0("Dichotomizing raw score with threshold: ", threshold, "\n"))
  obj <- list(threshold=threshold)
  list(obj=obj, res=data.frame(ID=sampleNames(test.eset),raw.score=raw.score,high.risk=(raw.score > threshold)))
}

## Combine any expression, clinical, or genomic data into a single data matrix
## from the "eset", "clinical", or "genomic" slots of data, respectively.
## Additionally, possibly log-transform the expression data, possibly after adding
## an offset to their values.
assemble.predictor.matrix <- function(data, already.log2.transformed = FALSE, offset = 0.05) {

  eset <- data[["eset"]]
  clin <- data[["clinical"]]
  genomic <- data[["genomic"]]

  samples <- NULL
  data.mat <- NULL
  
  if(!is.null(clin)) {
    samples <- rownames(clin)  
  }

  expr <- NULL
  if(!is.null(eset)) {
    expr <- exprs(eset)
  
    if( !already.log2.transformed ) {
      expr <- log2(expr + offset)
    }
    expr <- t(expr)
    if(is.null(samples)) { 
      samples <- rownames(expr)
    } else {
      samples <- intersect(rownames(expr), samples)
    }
  }
  
  if(!is.null(genomic)) {
    genomic <- t(genomic) 
    if(is.null(samples)) { 
      samples <- rownames(genomic)
    } else {
      samples <- intersect(rownames(genomic), samples)
    }
  }

  data.mat <- NULL
  if(!is.null(clin)) {
    clin <- clin[samples,]
    data.mat <- clin
  }
  
  if(!is.null(expr)) {
    expr <- expr[samples,]
    if(is.null(data.mat)) {
      data.mat <- expr
    } else {
      data.mat <- cbind(data.mat, expr)
    }
  }
  
  if(!is.null(genomic)) {
    genomic <- genomic[samples,]
    if(is.null(data.mat)) {
      data.mat <- genomic
    } else {
      data.mat <- cbind(data.mat, genomic)
    }
  }

  data.mat
}
  
