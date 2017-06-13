Description
===========

This is a simple classifier that calls high risk populations based on the number
of somatic mutations in a tumor sample (MB = mutation burden)

Building
========

To build this container run the following
```
$ docker build --rm -t dna-mb:`cat VERSION` .
```

Running
=======

To run this container, invoke the following
```
# Assumes that
# $TESTDATA = <directory housing the clinical file 'sc3_Validation_ClinAnnotations.csv'>
# $OUTPUT = <directory where output file will be written>
$ docker run --rm -ti -v $TESTDATA:/test-data -v $OUTPUT:/output dna-mb:`cat VERSION` R
```
