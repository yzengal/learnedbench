#!/bin/bash

DATA_PATH="../data/"
BENCH_BIN_PATH="../build/bin/"
TEST_CASE="join"
BENCH2D_DEFAULT="${BENCH_BIN_PATH}bench2d${TEST_CASE}_default"
BENCH2D_FS="${BENCH_BIN_PATH}bench2d${TEST_CASE}_fs"
BENCH2D_OSM="${BENCH_BIN_PATH}bench2d${TEST_CASE}_osm"
BENCH3D_TOR="${BENCH_BIN_PATH}bench3d${TEST_CASE}_toronto"

REAL_DATA_PATH="${DATA_PATH}real/"
SYN_DATA_PATH="${DATA_PATH}synthetic/"
DEFAULT_SYN_DATA_PATH="${DATA_PATH}synthetic/"

RESULT_PATH="../results/join/"

#mkdir ${RESULT_PATH}

MILLION=1000000

data="fs"
for index in "rtree" "rstar" "kdtree" "ann" "octree" "fs" "zm" "mli" "lisa"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH2D_FS} ${index} "${REAL_DATA_PATH}$data" 3680126 "join" > "${RESULT_PATH}${index}_${data}"
done

data="toronto"
for index in "rtree" "rstar" "kdtree" "ann" "octree" "fs" "zm" "mli" "lisa"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH3D_TOR} ${index} "${REAL_DATA_PATH}$data" 21567172 "join" > "${RESULT_PATH}${index}_${data}"
done

data="osm-china"
for index in "rtree" "rstar" "kdtree" "ann" "octree" "fs" "zm" "mli" "lisa"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH2D_OSM} ${index} "${REAL_DATA_PATH}$data" 62734869 "join" > "${RESULT_PATH}${index}_${data}"
done
