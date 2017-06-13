FROM ubuntu
RUN apt-get update
RUN apt-get -y install r-base
RUN apt-get -y install aptitude libcurl4-openssl-dev  libxml2-dev libxml2-dev

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("Biobase")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("biomaRt")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("Matrix.utils")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("plyr")'
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('data.table')"

COPY run-mm-sc2.R /run-mm-sc2.R
COPY R /R
COPY score_sc2.sh /score_sc2.sh
COPY model-state-metadata.Rd /model-state-metadata.Rd
