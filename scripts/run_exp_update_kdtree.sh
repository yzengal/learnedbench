#!/bin/bash

DATA_PATH="../data/"
BENCH2D_FS="../build/bin/bench2dupdate_fs"
BENCH2D_OSM="../build/bin/bench2dupdate_osm"
BENCH3D_TOR="../build/bin/bench3dupdate_toronto"

REAL_DATA_PATH="${DATA_PATH}real/"
SYN_DATA_PATH="${DATA_PATH}synthetic/"
DEFAULT_SYN_DATA_PATH="${DATA_PATH}synthetic/"

RESULT_PATH="../results/update/"

# run experiments on FourSquare
data="fs"
for index in "kdtree"
do
	break
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH2D_FS} ${index} "${REAL_DATA_PATH}$data" 3680126 "update" > "${RESULT_PATH}${index}_${data}"
done

data="osm-china"
for index in "kdtree"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH2D_OSM} ${index} "${REAL_DATA_PATH}$data" 62734869 "update" > "${RESULT_PATH}${index}_${data}"
done

data="toronto"
for index in "kdtree"
do
    echo "Benchmark ${index} dataset ${data}"
    ${BENCH3D_TOR} ${index} "${REAL_DATA_PATH}$data" 21567172 "update" > "${RESULT_PATH}${index}_${data}"
done

