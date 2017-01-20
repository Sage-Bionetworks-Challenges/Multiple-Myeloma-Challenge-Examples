## Example model for dockerizing

## This directory contains some supporting files that need _not_ be included in the docker instance.  These are:
README
output/
output/predictions.tsv
sessionInfo.txt
create-model-metadata.R
training-data/
training-data/training-annotations.tsv
training-data/training-expression-data.tsv

## Create metadata state for classifier.
Rscript ./create-model-metadata.R --metadata-dir=metadata --training-dir=training-data/
## The contestant would do the above during the training phase and would just include the metadata
## directory populated by this script in their docker instance.

## Run the classifer.
## These are the only commands that would be included in the docker instance.
## BEGIN RUN CLASSIFIER
metadata_dir=metadata
test_dir=test-data
output_dir=output
./score_sc1.sh $metadata_dir $test_dir $output_dir
## END RUN CLASSIFIER

## The above creates an output file ${output_dir}/predictions.tsv
## The output should have the following format (tsv file with named columns):
## ID: sample name
## raw.score:  a numeric score with higher score corresponding to high-risk samples
## high.risk:  TRUE or FALSE, indicating whether sample is called high risk by the model

## R dependencies that I have explicitly included
## R version 3.3.1 Patched (2016-08-23 r77147)
## optparse (optparse_1.3.2)
## data.table (data.table_1.10.0)
## Biobase (Biobase_2.32.0)
## biomaRt (biomaRt_2.28.0)

## Output of sessionInfo() is in sessionInfo.txt and is:
R version 3.3.1 Patched (2016-08-23 r71147)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 14.04.2 LTS

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  methods   stats     graphics  grDevices utils     datasets 
[8] base     

other attached packages:
[1] Biobase_2.32.0       BiocGenerics_0.18.0  biomaRt_2.28.0      
[4] data.table_1.10.0    optparse_1.3.2       BiocInstaller_1.22.3

loaded via a namespace (and not attached):
 [1] IRanges_2.6.1        DBI_0.5              RCurl_1.95-4.8      
 [4] getopt_1.20.0        AnnotationDbi_1.34.4 RSQLite_1.0.0       
 [7] S4Vectors_0.10.3     stats4_3.3.1         bitops_1.0-6        
[10] XML_3.98-1.4        
