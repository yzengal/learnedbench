#!/bin/bash

DATA_PATH="/data/yuxiang"
BENCH_RSMI="../build/bin/bench_rsmi"


REAL_DATA_PATH="${DATA_PATH}real/"
SYN_DATA_PATH="${DATA_PATH}synthetic/"
DEFAULT_SYN_DATA_PATH="${DATA_PATH}synthetic/Default/"

RESULT_PATH="../results/join/"

#mkdir ${RESULT_PATH}

data="fs"
for index in "rsmi" 
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH_RSMI} ${index} "${REAL_DATA_PATH}$data" 3680126 "join" > "${RESULT_PATH}${index}_${data}"
done

data="osm-china"
for index in "rsmi" 
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH_RSMI} ${index} "${REAL_DATA_PATH}$data" 62734869 "join" > "${RESULT_PATH}${index}_${data}"
done
