FROM ubuntu

RUN apt-get update
RUN apt-get -y install r-base
RUN apt-get -y install aptitude libcurl4-openssl-dev  libxml2-dev libxml2-dev

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("Biobase")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("BiocGenerics")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("biomaRt")'
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('optparse')"

COPY run-mm-sc1.R /run-mm-sc1.R
COPY R /R
COPY score_sc1.sh /score_sc1.sh
COPY modelstate.tsv /modelstate.tsv

