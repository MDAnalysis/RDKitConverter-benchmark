#!/bin/bash

IN_FILE="chembl_processed_unique.smi.gz"
DIR=$PWD

# go to data dir
cd $(dirname "$0")/../data

echo Splitting processed file in chunks of 200,000 lines
mkdir -p chunks
zcat ${IN_FILE} | split -d -a 1 --additional-suffix .smi -l 200000 - chunks/part

# counts
cd chunks
wc -l *.smi | grep .smi | awk '{print $1 > ".count_"$2}'

cd $DIR

