#!/bin/bash

DATA_PATH="../data/"
BENCH_BIN_PATH="../build/bin/"
TEST_CASE="octree"
BENCH2D_DEFAULT="${BENCH_BIN_PATH}bench2d${TEST_CASE}_default"
BENCH2D_FS="${BENCH_BIN_PATH}bench2d${TEST_CASE}_fs"
BENCH2D_OSM="${BENCH_BIN_PATH}bench2d${TEST_CASE}_osm"
BENCH3D_TOR="${BENCH_BIN_PATH}bench3d${TEST_CASE}_toronto"

REAL_DATA_PATH="${DATA_PATH}real/"
SYN_DATA_PATH="${DATA_PATH}synthetic/"
DEFAULT_SYN_DATA_PATH="${DATA_PATH}synthetic/Default/"

RESULT_PATH="../results/octree/"

#mkdir ${RESULT_PATH}

MILLION=1000000


data="fs"
for index in "octree"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH2D_FS} ${index} "${REAL_DATA_PATH}$data" 3680126 "all" > "${RESULT_PATH}${index}_${data}"
done

data="toronto"
for index in "octree"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH3D_TOR} ${index} "${REAL_DATA_PATH}$data" 21567172 "all" > "${RESULT_PATH}${index}_${data}"
done

data="osm-china"
for index in "octree"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH2D_OSM} ${index} "${REAL_DATA_PATH}$data" 62734869 "all" > "${RESULT_PATH}${index}_${data}"
done

for data in "uniform_20m_2_1" "gaussian_20m_2_1" "lognormal_20m_2_1"
do
    for index in "octree"
    do
        echo "Benchmark ${index} dataset ${data}"
        ${BENCH2D_DEFAULT} ${index} "${DEFAULT_SYN_DATA_PATH}$data" 20000000 "all" > "${RESULT_PATH}${index}_${data}"
    done
done

for dist in "uniform" "gaussian" "lognormal"
do
    for n in 1 10 50 100
    do
        for index in "octree"
        do
            echo "Benchmark ${index} dataset ${SYN_DATA_PATH}${dist}_${n}m_2_1"
            real_n=$[$n * $MILLION]
            ${BENCH2D_DEFAULT} ${index} "${SYN_DATA_PATH}${dist}_${n}m_2_1" $real_n "all" > "${RESULT_PATH}${index}_${dist}_${n}m_2"
        done
    done
done

for dist in "uniform" "gaussian" "lognormal"
do
    for k in 4 6 8
    do
        for index in "octree"
        do
            echo "Benchmark ${index} dataset ${dist}_${k}"
            "${BENCH_BIN_PATH}bench${k}doctree_default" ${index} "${SYN_DATA_PATH}${dist}_20m_${k}_1" 20000000 "all" > "${RESULT_PATH}${index}_${dist}_20m_${k}"
        done
    done
done
exit
