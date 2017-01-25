library(Biobase)
library(biomaRt)

## This classifier subtracts the means of the "down" probes from the        
## means of the "up" probes.  So set the score of the up probes to        
## be 1/n_up and similar for the down probes to be 1/n_down, where           
## n_up = # of up probes and n_down = # of down probes.                   
## Note that this is at the level of probes.  This            
## multiplicative factor will be adjusted for multimappers in uams.70.gene
## using adjust.scores.to.reflect.multimappers.
uams.70.get.coefficients <- function() {
  up<-c("202345_s_at","1555864_s_at","204033_at","206513_at","1555274_a_at","211576_s_at","204016_at","1565951_s_at","219918_s_at","201947_s_at","213535_s_at","204092_s_at","213607_x_at","208117_s_at","210334_x_at","204023_at","201897_s_at","216194_s_at","225834_at","238952_x_at","200634_at","208931_s_at","206332_s_at","220789_s_at","218947_s_at","213310_at","224523_s_at","201231_s_at","217901_at","226936_at","58696_at","200916_at","201614_s_at","200966_x_at","225082_at","242488_at","243011_at","201105_at","224200_s_at","222417_s_at","210460_s_at","200750_s_at","206364_at","201091_s_at","203432_at","221970_s_at","212533_at","213194_at","244686_at","200638_s_at","205235_s_at")
  down<-c("201921_at","227278_at","209740_s_at","227547_at","225582_at","200850_s_at","213628_at","209717_at","222495_at","1557277_a_at","1554736_at","218924_s_at","226954_at","202838_at","230192_at","48106_at","237964_at","202729_s_at","212435_at")
  up.coeffs <- rep(1.0, length(up)) / length(up)
  down.coeffs <- rep(-1.0, length(down)) / length(down)
  coefficients <- c(up.coeffs, down.coeffs)
  names(coefficients) <- c(up, down)
  coefficients
}

uams.70.probe.to.ensg.mapping <- function() {
  probe.to.ensg.mapping(uams.70.get.coefficients())
}

uams.70.probe.to.entrez.mapping <- function() {
  probe.to.entrez.mapping(uams.70.get.coefficients())
}

# Return a list with "threshold" and "res."
# "res" is a result data frame with columns ID, raw.score, and high.risk.
uams.70.eset <- function(eset, anno = NULL, coefficients = NULL, threshold = NULL) {
  if(is.null(coefficients)) {
    coefficients <- uams.70.get.coefficients()
  }
  if(is.null(anno) && is.null(threshold)) {
    ## Use the published threshold
    threshold <- 0.66
  }
  train.eset <- NULL
  train.anno <- NULL
  if(!is.null(anno) && is.null(threshold)) {
    train.eset <- eset
    train.anno <- anno
  }
  run.classifier.eset(test.eset = eset, test.anno = anno, coefficients = coefficients, already.log2.transformed = TRUE, train.eset = train.eset, train.anno = train.anno, threshold = threshold, standardize = FALSE)
}

# assumes that mapping is a data.frame with 2 columns: INDEX and PROBE
# INDEX = any type of ID (e.g., ENTREZID)
# PROBE = affy probe name (the name of the coefficients)
uams.70.gene <- function(eset, mapping, coefficients, anno = NULL, threshold = NULL) {
  classifier.df <- translate.probe.coefficients.to.gene.coefficients(coefficients, mapping)
  gene.coefficients <- adjust.scores.to.reflect.multimappers(classifier.df, eset)
  uams.70.eset(eset, anno = anno, coefficients = gene.coefficients, threshold = threshold)
}

uams.70.entrez <- function(X, y, threshold = NULL) {
  bm <- uams.70.probe.to.entrez.mapping()
  eset <- ExpressionSet(as.matrix(X))
  coefficients <- uams.70.get.coefficients()
  uams.70.gene(eset, mapping = bm, coefficients = coefficients, anno = y, threshold = threshold)
}

uams.70.ensg <- function(X, y, threshold = NULL) {
  bm <- uams.70.probe.to.ensg.mapping()
  eset <- ExpressionSet(as.matrix(X))
  coefficients <- uams.70.get.coefficients()
  uams.70.gene(eset, mapping = bm, coefficients = coefficients, anno = y, threshold = threshold)
}
