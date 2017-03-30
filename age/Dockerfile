FROM ubuntu

RUN apt-get update
RUN apt-get -y install python
RUN apt-get -y install python-pandas

COPY age-model.py /age-model.py
COPY score_sc1.sh /score_sc1.sh