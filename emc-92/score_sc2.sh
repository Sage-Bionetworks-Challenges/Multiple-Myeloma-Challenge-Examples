#!/bin/bash

## Simply invoke the R script that implements our model.
## This R script assumes that it runs in the root directory of the Docker image and that
## the test data are mounted in ./test-data,
## the output should be written to ./output,
## and that the entire directory structure of the submitted Docker image
## (e.g., an R object encapsulating trained modeler state) is mounted at ./
./run-mm-sc2.R
