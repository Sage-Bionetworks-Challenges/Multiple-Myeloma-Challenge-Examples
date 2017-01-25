#!/bin/bash

metadata_dir=$1
test_dir=$2
output_dir=$3

./run-mm-sc1.R --metadata-dir=$metadata_dir --test-dir=$test_dir --output-dir=$output_dir
