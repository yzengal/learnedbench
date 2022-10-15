#!/bin/bash

DATA_PATH="../data/"
BENCH_BIN_PATH="../build/bin/"
SYN_DATA_PATH="${DATA_PATH}synthetic/"
RESULT_PATH="../results/box/"

# mkdir ${RESULT_PATH}

MILLION=1000000
DEFAULT_N=20
DEFAULT_D=2
DEFAULT_S=1


for s in 0.005 0.01 0.02 0.05 0.1
do
	for index in "rtree" "rstar" "fs" "edg" "ug" "qdtree"
	do
		fname="box_${DEFAULT_N}m_${DEFAULT_D}_${s}"
		echo "Benchmark ${index} dataset ${fname}"
		real_n=$[$DEFAULT_N * $MILLION]
		"${BENCH_BIN_PATH}bench2d_box" ${index} "${SYN_DATA_PATH}${fname}" ${real_n} "range" ${s} > "${RESULT_PATH}${index}_box_${s}"
	done
done
