
uams.17.get.coefficients <- function() {
  coefficients <- c(0.283,-0.296,-0.208,+0.314,-0.287,+0.251,+0.193,+0.269,+0.375,+0.158,+0.316,+0.232,-0.251,-0.23,-0.402,+0.191,+0.148)
  names(coefficients) <- c("200638_s_at","1557277_a_at","200850_s_at","201897_s_at","202729_s_at","203432_at","204016_at","205235_s_at","206364_at","206513_at","211576_s_at","213607_x_at","213628_at","218924_s_at","219918_s_at","220789_s_at","242488_at")
  coefficients
}

uams.17.probe.to.ensg.mapping <- function() {
  bm <- probe.to.ensg.mapping(uams.17.get.coefficients())
  # The CKS1BP3 gene is never referred to in the patent or manuscript so
  # we will ignore it from our calculations
  # See: https://www.google.com/patents/US20080187930
  bm <- bm[!(bm$INDEX %in% c("","ENSG00000268942")),]
  bm
}

uams.17.probe.to.entrez.mapping <- function() {
  bm <- probe.to.entrez.mapping(uams.17.get.coefficients())
  # The CKS1BP3 gene is never referred to in the patent or manuscript so
  # we will ignore it from our calculations
  # See: https://www.google.com/patents/US20080187930
  bm <- bm[!(bm$INDEX %in% c("","246715")),]
  bm
}

# Return a list with "threshold" and "res."
# "res" is a result data frame with columns ID, raw.score, and high.risk.
# standardize.using.training.data indicates whether we should normalize the data in eset
# according to the training data used in the original publication.
uams.17.eset <- function(eset, anno = NULL, coefficients = NULL, threshold = NULL, standardize.using.training.data = FALSE) {
  if(is.null(coefficients)) {
    coefficients <- uams.17.get.coefficients()
  }
  
  # Standardize according to the training data
  if(standardize.using.training.data) {
    data <- exprs(eset[intersect(names(coefficients),featureNames(eset)),])
    ## Donwload the GEO training data set ...
    ## training <- log2(exprs(suppressMessages(suppressWarnings(getGEO(GEO='GSE2658', destdir='.')[[1]])))[names(coefficients),])
    ## training <- training[match(rownames(data),rownames(training)),]
    ## trainingMeans <- rowMeans(training)
    ## trainingSDs <- apply(training,1,sd)
    ## data <- t((data - trainingMeans) / trainingSDs)
    
    ## ... but downloading the data is slow and error prone.  Here are the results, hardcoded.
    trainingMeans <- c(12.7122271946531,5.70421405991235,11.2377913501713,10.1350512537729,10.2237252641476,10.1132241529945,6.83684363220325,8.4676803122732,7.89931306956123,11.3257404788118,7.66805660126567,9.07005602886575,10.9986888716969,10.5768811754997,8.62432460604888,7.81197965645517,7.04306601208372)
    names(trainingMeans) <- c('200638_s_at','1557277_a_at','200850_s_at','201897_s_at','202729_s_at','203432_at','204016_at','205235_s_at','206364_at','206513_at','211576_s_at','213607_x_at','213628_at','218924_s_at','219918_s_at','220789_s_at','242488_at')
    trainingSDs <- c(0.566265433826941,1.439493522169,0.497673195366107,0.934992769643516,0.935758059559669,0.684442596437629,1.22889386878789,0.904132681968134,1.01714456624349,0.91353321534587,1.53473985465217,0.942466241663803,0.485836220886129,0.760059867964268,1.75616530808258,1.36275365512788,1.42532285132829)
    names(trainingSDs) <- c('200638_s_at','1557277_a_at','200850_s_at','201897_s_at','202729_s_at','203432_at','204016_at','205235_s_at','206364_at','206513_at','211576_s_at','213607_x_at','213628_at','218924_s_at','219918_s_at','220789_s_at','242488_at')

    # This is a hack.  If the feature names in eset are not probes (i.e., do not match
    # any of the names in the above training data), assume instead that the eset features
    # are entrezs.  In that case, translate the probes in the training data to entrezs.
    use.entrez <- !any(names(trainingMeans) %in% featureNames(eset))
    if(use.entrez) {    
      cat("Assuming data are entrezs")
      bm <- probe.to.entrez.mapping(trainingMeans)
      bm <- bm[!is.na(bm$INDEX),]
      if(any(duplicated(bm$INDEX))) {
        cat("Duplicate genes\n")
      }
      inters <- intersect(names(trainingMeans), bm$PROBE)
      trainingMeans <- trainingMeans[inters]
      trainingSDs <- trainingSDs[inters]
      genes <- bm$INDEX
      names(genes) <- bm$PROBE
      names(trainingSDs) <- as.vector(genes[names(trainingSDs)])
      names(trainingMeans) <- as.vector(genes[names(trainingMeans)])
    }
    
    # Translate these to genes
    trainingMeans <- trainingMeans[match(rownames(data),names(trainingMeans))]
    trainingSDs <- trainingSDs[match(rownames(data),names(trainingSDs))]
    
    data <- t((data - trainingMeans) / trainingSDs)
    exprs(eset) <- t(data)
  }
  
  if(is.null(anno) && is.null(threshold)) {
    ## Use the published threshold
    threshold <- 1.5
  }
  train.eset <- NULL
  train.anno <- NULL
  if(!is.null(anno) && is.null(threshold)) {
    train.eset <- eset
    train.anno <- anno
  }
  if(standardize.using.training.data) {
    cat("Standardizing by training data\n")
    res <- run.classifier.eset(test.eset = eset, test.anno = anno, coefficients = coefficients, already.log2.transformed = TRUE, train.eset = train.eset, train.anno = train.anno, threshold = threshold, standardize = FALSE)
  } else {
    cat("Standardizing by test data\n")
    trainingMeans <- rowMeans(as.matrix(exprs(eset)))
    trainingSDs <- apply(as.matrix(exprs(eset)), 1, sd)
    
    res <- run.classifier.eset(test.eset = eset, test.anno = anno, coefficients = coefficients, already.log2.transformed = TRUE, train.eset = train.eset, train.anno = train.anno, threshold = threshold, standardize = TRUE)
  }
  res
}

# assumes that mapping is a data.frame with 2 columns: INDEX and PROBE
# INDEX = any type of ID (e.g., ENTREZID)
# PROBE = affy probe name (the name of the coefficients)
uams.17.gene <- function(eset, mapping, coefficients, anno = NULL, threshold = NULL, standardize.using.training.data = FALSE) {
  classifier.df <- translate.probe.coefficients.to.gene.coefficients(coefficients, mapping)
  gene.coefficients <- adjust.scores.to.reflect.multimappers(classifier.df, eset)
  uams.17.eset(eset, anno = anno, coefficients = gene.coefficients, threshold = threshold, standardize.using.training.data = standardize.using.training.data)
}

uams.17.entrez <- function(X, y, threshold = NULL, standardize.using.training.data = FALSE) {
  bm <- uams.17.probe.to.entrez.mapping()
  eset <- ExpressionSet(as.matrix(X))
  coefficients <- uams.17.get.coefficients()
  uams.17.gene(eset, mapping = bm, coefficients = coefficients, anno = y, threshold = threshold, standardize.using.training.data = standardize.using.training.data)
}

uams.17.ensg <- function(X, y, threshold = NULL) {
  bm <- uams.17.probe.to.ensg.mapping()
  eset <- ExpressionSet(as.matrix(X))
  coefficients <- uams.17.get.coefficients()
  uams.17.gene(eset, mapping = bm, coefficients = coefficients, anno = y, threshold = threshold)
}

