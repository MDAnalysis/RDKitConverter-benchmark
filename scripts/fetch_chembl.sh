#!/bin/bash

echo Fetching ChEMBL 30
wget -q --show-progress -P data/ ${CHEMBL_SDF}

echo Counting the number of entries
count=$(zgrep -c "M  END" data/chembl_30.sdf.gz)

echo ${count} > data/.fetched_count
python -c "print(f'Fetched {${count}:,} molecules')"
