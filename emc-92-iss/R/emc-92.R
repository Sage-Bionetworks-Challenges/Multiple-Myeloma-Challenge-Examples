# This is an adaptation of the sourcescript.txt file written by
# Rowan Kuiper.
# The classifier functions themselves are largely unchanged, except
# that they can now be run against (Ensembl) genes instead of Affy probes.

emc.92.get.coefficients <- function() {
  coefficients<-c(0.0594273379801556,-0.110491572253885,-0.108832301005201,-0.062550086484268,-0.0418164824703641,0.0174897923948543,0.00668584689209701,0.0423154071076169,-0.0493203905211073,0.0164925076573208,-0.0344544706411548,0.00866233731636807,0.00455819600852017,-0.0644105573058562,-0.00414107711721288,-0.0183994540236258,-0.0520135932088583,-0.0767532974279582,-0.00063407777663618,0.0489558593887832,-0.0163691988848842,0.0420452994583269,0.0870276088358123,0.0406887458334764,-0.05606504920594,0.0529972408575157,0.0112796922615379,0.0140109440864371,0.00078080860374577,0.075000229909631,0.009265826714752,0.0860546977563988,-0.0996944005531858,0.0558748936987772,0.073011096187736,0.00538232634785429,0.055611094941233,-0.0520445402963303,0.0125985686829356,0.005608194716564,0.0163467078063,-0.0319185908460475,0.0773361612870129,-0.00900009551662273,-0.0575829465263513,0.022119958576764,0.0396248830654306,0.0524561174360331,0.0476703103071683,-0.0422745546301539,-0.0341783399437912,-0.0251675596648738,0.0495596844465701,0.0547508890118559,0.0437224907573062,-0.0371869820707693,0.0205354349665636,0.0445923111233813,-0.0069962048640006,0.025498590737044,0.0208131193381807,0.0685640666167821,-0.0330380997958823,-0.0105781969037762,-0.0585102612125953,0.0128945938676131,-0.0333628118834898,-0.0349253487676344,0.0660971817524709,-0.0617575340318364,-0.0210452256979565,-0.0389953930227488,-0.0873912090180467,-0.0176248165093743,0.0302900524131141,-0.0051520032747514,0.0745893468079515,-0.0322581999457566,0.0200333568025016,0.0115861932788841,-0.00969926642147875,0.00345861530256628,0.0277930303183715,0.0153863421114279,-0.0777964537268526,0.0349480953839344,-0.00022873908620183,-0.0529393125922244,0.0713507977617356,-0.0253569213540041,0.0384170565171902,0.0225490655109955)
  names(coefficients)<-c('204379_s_at','202728_s_at','239054_at','202842_s_at','213002_at','210334_x_at','201795_at','38158_at','208232_x_at','201307_at','226742_at','205046_at','204026_s_at','226218_at','217824_at','233399_x_at','224009_x_at','215177_s_at','202532_s_at','238662_at','212788_x_at','220351_at','202542_s_at','243018_at','209683_at','212282_at','208967_s_at','225366_at','217852_s_at','225601_at','231210_at','214482_at','208942_s_at','219550_at','231989_s_at','202553_s_at','223811_s_at','221041_s_at','221677_s_at','213350_at','200775_s_at','226217_at','217728_at','201930_at','216473_x_at','211714_x_at','221755_at','AFFX-HUMISGF3A/M97935_MA_at','206204_at','217548_at','215181_at','217732_s_at','214612_x_at','202813_at','200875_s_at','201292_at','222680_s_at','233437_at','223381_at','209026_x_at','221606_s_at','231738_at','230034_x_at','213007_at','242180_at','202322_s_at','208904_s_at','214150_x_at','238116_at','208732_at','200701_at','208667_s_at','208747_s_at','218662_s_at','211963_s_at','201555_at','207618_s_at','200933_x_at','221826_at','218355_at','219510_at','218365_s_at','222713_s_at','222154_s_at','228416_at','201102_s_at','203145_at','238780_s_at','202884_s_at','201398_s_at','212055_at','202107_s_at')
  coefficients
}

emc.92.probe.to.ensg.mapping <- function() {
  probe.to.ensg.mapping(emc.92.get.coefficients())
}

emc.92.probe.to.entrez.mapping <- function() {
  probe.to.entrez.mapping(emc.92.get.coefficients())
}

# Return a list with "threshold" and "res."
# "res" is a result data frame with columns ID, raw.score, and high.risk.
emc.92.eset <- function(eset, anno = NULL, coefficients = NULL, threshold = NULL) {
  if(is.null(coefficients)) {
    coefficients <- emc.92.get.coefficients()
  }
  if(is.null(anno) && is.null(threshold)) {
    ## Use the published threshold
    threshold <- 0.827
  }
  train.eset <- NULL
  train.anno <- NULL
  if(!is.null(anno) && is.null(threshold)) {
    train.eset <- eset
    train.anno <- anno
  }
  run.classifier.eset(test.eset = eset, test.anno = anno, coefficients = coefficients, already.log2.transformed = TRUE, train.eset = train.eset, train.anno = train.anno, threshold = threshold)
}

# assumes that mapping is a data.frame with 2 columns: INDEX and PROBE
# INDEX = any type of ID (e.g., ENTREZID)
# PROBE = affy probe name (the name of the coefficients)
emc.92.gene <- function(eset, mapping, coefficients, anno = NULL, threshold = NULL) {
  classifier.df <- translate.probe.coefficients.to.gene.coefficients(coefficients, mapping)
  gene.coefficients <- adjust.scores.to.reflect.multimappers(classifier.df, eset)
  emc.92.eset(eset, anno = anno, coefficients = gene.coefficients, threshold = threshold)
}

# y (if non-NULL) should have a sample identifier column ("ID") and survival related columns
# event and time.
emc.92.ensg <- function(X, y, threshold = NULL) {
  bm <- emc.92.probe.to.ensg.mapping()
  eset <- ExpressionSet(as.matrix(X))
  coefficients <- emc.92.get.coefficients()
  emc.92.gene(eset, mapping = bm, coefficients = coefficients, anno = y, threshold = threshold)
}

train.emc.92.ensg <- function(X, y, threshold = NULL) {
  emc.92.ensg(X, y, threshold)
}

test.emc.92.ensg <- function(obj, X) {
  emc.92.ensg(X, y = NULL, threshold = obj$threshold)
}

# X is assumed to have a column ISSstage
emc.92.ensg.iss <- function(X, y, threshold = NULL) {
  iss.col <- "ISSstage"
  res <- emc.92.ensg(X[,!(colnames(X) %in% iss.col)], y, threshold = NULL)
  res$res$high.risk <- res$res$high.risk | X[,iss.col] %in% c('III', 'II')
  return(list(threshold=res$threshold,res=res$res))
}

train.emc.92.ensg.iss <- function(X, y, threshold = NULL) {
   emc.92.ensg.iss(X, y, threshold)
}

test.emc.92.ensg.iss <- function(obj, X) {
   emc.92.ensg.iss(X, y = NULL, obj$threshold)
}


emc.92.entrez <- function(X, y, threshold = NULL) {
  bm <- emc.92.probe.to.entrez.mapping()
  eset <- ExpressionSet(as.matrix(X))
  coefficients <- emc.92.get.coefficients()
  emc.92.gene(eset, mapping = bm, coefficients = coefficients, anno = y, threshold = threshold)
}

train.emc.92.entrez <- function(X, y, threshold = NULL) {
   emc.92.entrez(X, y, threshold)
}

test.emc.92.entrez <- function(obj, X) {
   emc.92.entrez(X, y = NULL, threshold = obj$threshold)
}

emc.92.entrez.iss <- function(X, y, threshold = NULL) {
  iss.col <- "ISSstage"
  res <- emc.92.entrez(X[,!(colnames(X) %in% iss.col)], y, threshold = NULL)
  res$res$high.risk <- res$res$high.risk | X[,iss.col] %in% c('III', 'II')
  return(list(threshold=res$threshold,res=res$res))
}

train.emc.92.entrez.iss <- function(X, y, threshold = NULL) {
   emc.92.entrez.iss(X, y, threshold)
}

test.emc.92.entrez.iss <- function(obj, X) {
   emc.92.entrez.iss(X, y = NULL, threshold = obj$threshold)
}
