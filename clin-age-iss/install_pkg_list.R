library(utils)


source("https://bioconductor.org/biocLite.R")
installed <- installed.packages();
installed <- as.character(unlist(installed[,1]))
to.install <- as.character(read.table(file='pkg_list.txt')[,1])
to.install <- to.install[!(to.install %in% installed)]
if(length(to.install) > 0) {
    biocLite(to.install)
}
