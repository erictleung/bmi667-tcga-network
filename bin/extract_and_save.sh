#!/bin/bash

# Description: 
#   Downloads the human protein-protein-interaction pairs in PSI-MITAB format
#   and parses it for protein names in pairwise interactions.
# Output:
#   - homo_sapiens.mitab.interaction.txt (from Reactome website)
#   - prot_interaction.tsv (pairwise prot-prot interactions)
#   - all_unique.txt (list of unique genes in prot-prot interactions)
# Dependencies:
#   - curl (7.43.0 Tested)
#   - gunzip (1.6 Tested)
#   - awk (4.1.3 Tested)
# Versions Tested: 
#   bash 4.3.42
#   bash 3.2.57

# human protein-protein interaction pairs in PSI-MITAB format
DIR=../data/
FILE=homo_sapiens.mitab.interactions.txt
ZIP="$FILE.gz"
ADDRESS=http://www.reactome.org/download/current/$ZIP
PAIR=prot_interaction.tsv
FINAL=all_unique.txt

# download file if file or zipped file not already downloaded
if [ -f $DIR$ZIP ] || [ -f $DIR$FILE ]; then
    echo "File already downloaded."
else
    curl $ADDRESS -o $DIR$ZIP
fi
echo "Successfully downloaded file."

# unzip file
if [ -f $DIR$FILE ]; then
    echo "Found unzipped version of file."
elif [ -f $DIR$ZIP ]; then
    echo "Unzipping file previously downloaded."
    gunzip $DIR$ZIP
else
    echo "Failed to download $FILE. Try again."
    exit 1
fi
echo "Successfully unzipped file."

# parse through file to get pairwise protein-protein interactions
if [ -f $DIR$PAIR ]; then
    echo "Pairwise interactions extracted already."
else
    echo "Extracting pairwise gene interactions."

    # explanation for awk code below
    # 1 extracts only relavent columns
    # 2 removes header
    # 3 extracts only gene names in select columns
    # 4 removes missing gene names represented with '-'
    # 5 alphabetize column values in each row to remove duplicate interactions
    # 6 remove self interactions
    # 7 sort and keep only unique rows
    awk -F '\t' '{ print $3 "\t" $4 }' $DIR$FILE | \
        sed 1d | \
        sed 's/.*:.*_\(.*\)(shortlabel)	.*:.*_\(.*\)(shortlabel)/\1	\2/' | \
        awk '{ if ($1 != "-" && $2 != "-") print $1 "\t" $2 }' | \
        awk '{ if ($1 > $2) print $2 "\t" $1; else print $1 "\t" $2 }' | \
        awk '{ if ($1 != $2) print $1 "\t" $2 }' | \
        sort -u > \
        $DIR$PAIR
fi
echo "Successfully extracted pairwise gene interactions."

# get unique list of genes list of interactions
if [ -f $DIR$PAIR ]; then
    cat $DIR$PAIR | \
        awk '{ print $1 "\n" $2 }' | \
        sort -u > \
        $DIR$FINAL
fi
