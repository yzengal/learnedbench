#!/bin/bash

DATA_PATH="/data/yuxiang/"
BENCH_BIN_PATH="../build/bin/"
BENCH2D_DEFAULT="${BENCH_BIN_PATH}bench2d_default"
BENCH2D_FS="${BENCH_BIN_PATH}bench2d_fs"
BENCH2D_OSM="${BENCH_BIN_PATH}bench2d_osm"

REAL_DATA_PATH="${DATA_PATH}real/"
SYN_DATA_PATH="${DATA_PATH}synthetic/"
DEFAULT_SYN_DATA_PATH="${DATA_PATH}synthetic/Default/"

RESULT_PATH="../results/qdtree/"

#mkdir ${RESULT_PATH}

MILLION=1000000


data="fs"
for index in "qdtree"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH2D_FS} ${index} "${REAL_DATA_PATH}$data" 3680126 "range" > "${RESULT_PATH}${index}_${data}"
done

data="osm-china"
for index in "qdtree"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH2D_OSM} ${index} "${REAL_DATA_PATH}$data" 62734869 "range" > "${RESULT_PATH}${index}_${data}"
done

for data in "uniform_20m_2_1" "gaussian_20m_2_1" "lognormal_20m_2_1"
do
    for index in "qdtree"
    do
        echo "Benchmark ${index} dataset ${data}"
        ${BENCH2D_DEFAULT} ${index} "${DEFAULT_SYN_DATA_PATH}$data" 20000000 "range" > "${RESULT_PATH}${index}_${data}"
    done
done

for dist in "uniform" "gaussian" "lognormal"
do
    for n in 1 10 50 100
    do
        for index in "qdtree"
        do
            echo "Benchmark ${index} dataset ${SYN_DATA_PATH}${dist}_${n}m_2_1"
            real_n=$[$n * $MILLION]
            ${BENCH2D_DEFAULT} ${index} "${SYN_DATA_PATH}${dist}_${n}m_2_1" $real_n "range" > "${RESULT_PATH}${index}_${dist}_${n}m_2"
        done
    done
done
