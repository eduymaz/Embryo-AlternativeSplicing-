#!/bin/bash
'''
This project utilizes data from the embryonic stages within the GSE71318 biosample.
Data available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71318.
To proceed, download the metadata file provided on the platform and extract the names listed under the "Run" column. 
Save these names as srrlist.txt and follow the steps outlined below.
'''
# File SRR ID
SRR_LIST="srrlist.txt" # can make a config file like snakefile methodology

# Read and process SRR IDs
while IFS= read -r srr_id; do
    echo "İndiriliyor: $srr_id"
    fastq-dump --gzip "$srr_id" #if you want gunzip :fastq-dump "$srr_id" 
done < "$SRR_LIST"
