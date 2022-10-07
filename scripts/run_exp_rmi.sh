#!/bin/bash

DATA_PATH="../data/"
BENCH_BIN_PATH="../build/bin/"
BENCH2D_DEFAULT="${BENCH_BIN_PATH}bench2d_default"
BENCH2D_FS="${BENCH_BIN_PATH}bench2djoin_fs"
BENCH2D_OSM="${BENCH_BIN_PATH}bench2djoin_osm"
BENCH3D_TOR="${BENCH_BIN_PATH}bench3djoin_toronto"

REAL_DATA_PATH="${DATA_PATH}real/"
SYN_DATA_PATH="${DATA_PATH}synthetic/"
DEFAULT_SYN_DATA_PATH="${DATA_PATH}synthetic/"

RESULT_PATH="../results/rmi/"

#mkdir ${RESULT_PATH}

MILLION=1000000

data="fs"
for index in "mli" "flood" "lisa"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH2D_FS} ${index} "${REAL_DATA_PATH}$data" 3680126 "all" > "${RESULT_PATH}${index}_${data}"
done
exit

data="toronto"
for index in "mli" "flood" "lisa"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH3D_TOR} ${index} "${REAL_DATA_PATH}$data" 21567172 "all" > "${RESULT_PATH}${index}_${data}"
done

data="osm-china"
for index in "mli" "flood" "lisa"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH2D_OSM} ${index} "${REAL_DATA_PATH}$data" 62734869 "all" > "${RESULT_PATH}${index}_${data}"
done
