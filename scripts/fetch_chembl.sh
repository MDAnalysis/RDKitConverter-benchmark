#!/bin/bash

echo Fetching ChEMBL 33
wget -q --show-progress -P data/ ${CHEMBL_SDF}

echo Counting the number of entries
count=$(zgrep -c "M  END" data/chembl_33.sdf.gz)

echo ${count} > data/.fetched_count
python -c "print(f'Fetched {${count}:,} molecules')"
