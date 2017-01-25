suppressPackageStartupMessages(library("randomForest"))
## suppressPackageStartupMessages(library("parallel"))

## num.cores <- detectCores()
## if(!is.na(num.cores) && (num.cores > 1)) {
##   suppressPackageStartupMessages(library("doMC"))
##   cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
##   registerDoMC(cores=(num.cores-1))
## }

train.random.forest <- function(X, y) {

  flag <- apply(X, 1, function(row) any(is.na(row)))
  if(any(flag)) { cat("Excluding NAs during training\n") }
  data.mat <- X[!flag,]

  inters <- rownames(data.mat, names(y))
  data.mat <- data.mat[inters,]
  lbl <- factor(y[inters])
  
  min.factor <- min(table(lbl))
  max.factor <- max(table(lbl))
  num.factors <- length(unique(lbl))
  sampsize <- rep(floor(0.7*min.factor), num.factors)
  upsample <- FALSE
  rf <- NULL
  if(upsample) {
    flag <- lbl == levels(lbl)[1]
    data.1 <- data.mat[flag,]
    data.2 <- data.mat[!flag,]
    lbl.1 <- lbl[flag]
    lbl.2 <- lbl[!flag]
    cat("Upsampling\n")
    sz <- max(floor(1*max.factor), min.factor)
    sz <- 2 * min.factor
    samples1 <- sample.int(n=length(which(flag)), size=sz, replace=TRUE)
    samples2 <- sample.int(n=length(which(!flag)), size=sz, replace=TRUE)
    data.mat <- rbind(data.1[samples1,],data.2[samples2,])
    lbl <- factor(c(lbl.1[samples1],lbl.2[samples2]))
    sampsize <- rep(floor(0.7*sz), 2)
    print(sampsize)
    rf <- randomForest(x=data.mat, y=lbl, keep.forest=TRUE, proximity=FALSE, strata=lbl, sampsize=sampsize, importance=TRUE, replace=FALSE, max.nodes=5)
  } else {
    rf <- randomForest(x=data.mat, y=lbl, keep.forest=TRUE, proximity=FALSE, strata=lbl, sampsize=sampsize, importance=TRUE, replace=FALSE)
  }

  rf
}

test.random.forest <- function(rf, X) {
  flag <- apply(X, 1, function(row) any(is.na(row)))
  if(any(flag)) { cat("Excluding NAs during testing\n") }
  data.mat <- X[!flag,]

  predicted.response <- predict(rf, data.mat, type="response")
  predicted.prob <- predict(rf, data.mat, type="prob")
  list(ID=rownames(data.mat), raw.score=predicted.prob[,"1"], high.risk=predicted.response)
}

