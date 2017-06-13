Description
===========

This is a simple classifier that calls high risk populations based on the ISS
stage (Age > 65 is high risk, others are low risk)

Building
========

To build this container run the following
```
$ docker build --rm -t clin-age:`cat VERSION` .
```

Running
=======

To run this container, invoke the following
```
# Assumes that
# $TESTDATA = <directory housing the clinical file 'sc1_Validation_ClinAnnotations.csv'>
# $OUTPUT = <directory where output file will be written>
$ docker run --rm -v $TESTDATA:/test-data -v $OUTPUT:/output clin-age:`cat VERSION` /score_sc1.sh
```
