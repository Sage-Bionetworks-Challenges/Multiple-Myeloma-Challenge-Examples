library(Biobase)
library(biomaRt)

## This classifier just takes the average of 5 probes.  So set the score to 1/n_probes,
## where n_probes is the number of probes.  
## This multiplicative factor will be adjusted for multimappers in uams.5.gene
## using adjust.scores.to.reflect.multimappers.
uams.5.get.coefficients <- function() {
  coeff.names <- c('204033_at','200916_at','204023_at','202345_s_at','201231_s_at')
  coefficients <- rep(1.0, length(coeff.names)) / length(coeff.names)
  names(coefficients) <- coeff.names
  coefficients
}

uams.5.probe.to.ensg.mapping <- function() {
  probe.to.ensg.mapping(uams.5.get.coefficients())
}

uams.5.probe.to.entrez.mapping <- function() {
  probe.to.entrez.mapping(uams.5.get.coefficients())
}

# Return a list with "threshold" and "res."
# "res" is a result data frame with columns ID, raw.score, and high.risk.
uams.5.eset <- function(eset, anno = NULL, coefficients = NULL, threshold = NULL) {
  if(is.null(coefficients)) {
    coefficients <- uams.5.get.coefficients()
  }
  if(is.null(anno) && is.null(threshold)) {
    ## Use the published threshold
    threshold <- 10.68
  }
  train.eset <- NULL
  train.anno <- NULL
  if(!is.null(anno) && is.null(threshold)) {
    train.eset <- eset
    train.anno <- anno
  }
  # Don't standardize
  run.classifier.eset(test.eset = eset, test.anno = anno, coefficients = coefficients, already.log2.transformed = TRUE, train.eset = train.eset, train.anno = train.anno, threshold = threshold, standardize = FALSE)
}

# assumes that mapping is a data.frame with 2 columns: INDEX and PROBE
# INDEX = any type of ID (e.g., ENTREZID)
# PROBE = affy probe name (the name of the coefficients)
uams.5.gene <- function(eset, mapping, coefficients, anno = NULL, threshold = NULL) {
  classifier.df <- translate.probe.coefficients.to.gene.coefficients(coefficients, mapping)
  gene.coefficients <- adjust.scores.to.reflect.multimappers(classifier.df, eset)
  uams.5.eset(eset, anno = anno, coefficients = gene.coefficients, threshold = threshold)
}

uams.5.entrez <- function(X, y, threshold = NULL) {
  bm <- uams.5.probe.to.entrez.mapping()
  eset <- ExpressionSet(as.matrix(X))
  coefficients <- uams.5.get.coefficients()
  uams.5.gene(eset, mapping = bm, coefficients = coefficients, anno = y, threshold = threshold)
}

uams.5.ensg <- function(X, y, threshold = NULL) {
  bm <- uams.5.probe.to.ensg.mapping()
  eset <- ExpressionSet(as.matrix(X))
  coefficients <- uams.5.get.coefficients()
  uams.5.gene(eset, mapping = bm, coefficients = coefficients, anno = y, threshold = threshold)
}

