#!/bin/bash

DATA_PATH="../data/"
BENCH_BIN_PATH="../build/bin/"
BENCH2D_RSMI="${BENCH_BIN_PATH}bench_rsmi_update"

REAL_DATA_PATH="${DATA_PATH}real/"
SYN_DATA_PATH="${DATA_PATH}synthetic/"
DEFAULT_SYN_DATA_PATH="${DATA_PATH}synthetic/"

RESULT_PATH="../results/update/"

# rmi cannot terminate

# run experiments on FourSquare
data="fs"
for index in "rsmi"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH2D_RSMI} ${index} "${REAL_DATA_PATH}$data" 3680126 "update" > "${RESULT_PATH}${index}_${data}"
done
exit


