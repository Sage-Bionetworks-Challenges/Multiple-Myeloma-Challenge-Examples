FROM r-base:3.3.3
MAINTAINER Fadi Towfic <ftowfic@celgene.com>
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2 \
    libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    openssh-client
COPY R /R
COPY score_sc1.sh /score_sc1.sh
RUN chmod +x /score_sc1.sh
COPY . /tmp
RUN cd /tmp && /usr/bin/Rscript install_pkg_list.R
RUN cd /tmp && /usr/bin/Rscript install_synapse_client.R
