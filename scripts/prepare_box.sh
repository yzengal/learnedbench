#!/bin/bash

DATA_PATH="../data/"
BENCH_BIN="../build/bin/boxgen"

#mkdir DATA_PATH

REAL_DATA_PATH="${DATA_PATH}real/"
SYN_DATA_PATH="${DATA_PATH}synthetic/"

#mkdir $REAL_DATA_PATH
#mkdir $SYN_DATA_PATH
#mkdir "${SYN_DATA_PATH}Default"

# synthetic data
MILLION=1000000
DEFAULT_N=20
DEFAULT_D=2
DEFAULT_S=1

# default setting
# echo "Generate default data..."
# for dist in "uniform" "gaussian" "lognormal"
# do
	# fname="${dist}_${DEFAULT_N}m_${DEFAULT_D}_${DEFAULT_S}"
	# real_n=$[$DEFAULT_N * $MILLION]
	# $BENCH_BIN -t gen_data -f $fname --dist $dist -n $real_n -d $DEFAULT_D -s $DEFAULT_S
	# mv $fname "${SYN_DATA_PATH}Default"
# done
# exit

# varying dataset size N
echo "Generate data by varying different max_side..."
for s in 0.005 0.01 0.02 0.05 0.1
do
	fname="box_${DEFAULT_N}m_${DEFAULT_D}_${s}"
	real_n=$[$DEFAULT_N * $MILLION] 
	$BENCH_BIN -t gen_box -f $fname -n $real_n -d $DEFAULT_D -s $s
	mv $fname "${SYN_DATA_PATH}"
done
