#!/bin/bash

IN_FILE="chembl_processed.smi.gz"
OUT_FILE="chembl_processed_unique.smi.gz"

N_THREADS=$(grep -c ^processor /proc/cpuinfo)
if [[ $N_WORKERS -lt 0 ]]; then
    let N_THREADS+=1+$N_WORKERS
else
    N_THREADS=$N_WORKERS
fi

# go to data dir
cd $(dirname "$0")/../data

echo Removing duplicates
# uniquify
zcat "${IN_FILE}" \
| sort -u -t' ' -k3,3 --parallel ${N_THREADS} \
| cut -d' ' -f1,2 \
| gzip > "${OUT_FILE}"

# count
zcat "${OUT_FILE}" | wc -l | tr -d '\n' > .processed_unique_count
COUNT=$(cat .processed_unique_count)
python -c "print(f'Wrote \'data/${OUT_FILE}\' with {${COUNT}:,} unique entries')"

cd - > /dev/null

