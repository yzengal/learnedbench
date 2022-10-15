# Evaluating Mtulit-dimensional Learned Indices
This is the source code repo for our emperical study on multi-dimensional learned indices.


## Compared Methods
### Learned Indices
We compare **six** recent multi-dimensional learned indices:
- **ZM-Index** [1]
- **ML-Index** [2] 
- **IF-Index** [3] 
- **RSMI** [4]  
- **LISA** [5] 
- **Flood** [6]

### Non-learned Baselines
- **FullScan**: sequential scan
- **R\*-tree** and **bulk-loading R-tree**: we use the implementation from `boost::geometry`
- **kdtree**: we use a header-only kdtree implementation `nanoflann`  https://github.com/jlblancoc/nanoflann
- **ANN**: another kntree viriant from `ANN` project http://www.cs.umd.edu/~mount/ANN/
- **Quad-tree**: we use the implementation from `GEOS`
- **Grid**: uniform grid (UG) and equal-depth grid (EDG)

## Compilation
### Step 1: Setup Dependencies
- `boost 1.79`: https://www.boost.org/users/history/version_1_79_0.html
- `TPIE`: https://github.com/thomasmoelhave/tpie
- `GEOS`: https://libgeos.org/
- `gperftools`: https://github.com/gperftools/gperftools
- `libtorch`: https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.4.0%2Bcpu.zip
- `numpy` and `matplotlib` for result visualization
- `google-benchmark`: https://github.com/google/benchmark

### Step 2: Build RSMI and ANN
Most of the benchmark and indices (except `RSMI` and `ANN`) are implemented as header-only libraries. 

Compile `RSMI`:
```sh
cd indexes/rsmi
mkdir build && cd build
cmake ..
make
```

Compile `ANN`:
```sh
cd indexes/ann_1.1.2
make
```

### Step 3: Build Benchmark
Modify the following variables in `CMakeLists.txt`:
```
BOOST_ROOT, Boost_INCLUDE_DIR, Boost_LIBRARY_DIR: path to boost
TORCH_PATH: path to libtorch
EXECUTABLE_OUTPUT_PATH: path to compiled benchmark binaries
```

Compile RSMI benchmark:
```sh
mkdir build && cd build
cmake .. -DRSMI=ON
make
```

Compile RSMI benchmark with heap profiling enabled:
```sh
rm -rf * # clear cmake cache
cmake .. -DRSMI=ON -DPROFILE=ON
make
```

Compile benchmark for other indices:
```sh
rm -rf * # clear cmake cache
cmake ..
make
```

Compile benchmark for other indices with heap profiling enabled:
```sh
rm -rf * # clear cmake cache
cmake .. -DPROFILE=ON
make
```

Compile benchmark for complicated shape (rectangle):
```sh
rm -rf * # clear cmake cache
cmake .. -DTESTBOX=ON
make
```

Compile benchmark for data updates:
```sh
rm -rf * # clear cmake cache
cmake .. -DTESTUPDATE=ON
make
```

Compile benchmark for tuning flood with query workload:
```sh
rm -rf * # clear cmake cache
cmake .. -DTESTTUNE=ON
make
```


Compile benchmark for kNN join:
```sh
rm -rf * # clear cmake cache
cmake .. -DTESTJOIN=ON
make
```


Compile benchmark for replacing PGM with RMI:
```sh
rm -rf * # clear cmake cache
cmake .. -DTESTRMI=ON
make
```

Compile benchmark for tuning food with qdtree layout
```sh
rm -rf * # clear cmake cache
cmake .. -DTESTFLOOD=ON
make
```

Compile benchmark for octree
```sh
rm -rf * # clear cmake cache
cmake .. -DTESTOCTREE=ON
make
```

## Run Experiments
We prepare a script to download the real datasets and prepare synthetic datasets:
```sh
cd scripts
bash prepare_data.sh
```

We prepare several scripts to run the experiments.

Run experiments on default settings: `bash run_exp.sh`

Run experiments by varying N: `bash run_exp_n.sh`

Run experiments by varying dim: `bash run_exp_dim.sh`

Run experiments by varying eps: `bash run_exp_eps.sh`

Run experiments of RSMI: `bash rsmi.sh`

Run experiments of complicated shape: `bash run_exp_box.sh`

Run experiments of data updates: `bash run_exp_update.sh` `bash run_exp_update_rsmi.sh`

Run experiments of tuning flood with query workload: `bash run_exp_workload.sh`

Run experiments of kNN join: `bash run_exp_join.sh` `bash run_exp_join_rsmi.sh`

Run experiments of replacing PGM with RMI: `bash run_exp_rmi.sh`

Run experiments of tuning food with qdtree layout: `bash run_exp_rmi.sh`

Run experiments of octre: `bash run_exp_octree.sh`

The results are put in `/project_root/results`, and the figure drawing Jupyter notebooks are put in `/project_root/figures`.

## Reference
[1] Haixin Wang, Xiaoyi Fu, Jianliang Xu, and Hua Lu. 2019. Learned Index for Spatial Queries. In MDM. IEEE, 569–574.

[2] Angjela Davitkova, Evica Milchevski, and Sebastian Michel. 2020. The ML-Index: A Multidimensional, Learned Index for Point, Range, and Nearest-Neighbor Queries. In EDBT. OpenProceedings.org, 407–410.

[3] Ali Hadian, Ankit Kumar, and Thomas Heinis. 2020. Hands-off Model Integration in Spatial Index Structures. In AIDB@VLDB.

[4] Jianzhong Qi, Guanli Liu, Christian S. Jensen, and Lars Kulik. 2020. Effectively Learning Spatial Indices. Proc. VLDB Endow. 13, 11 (2020), 2341–2354.

[5] Pengfei Li, Hua Lu, Qian Zheng, Long Yang, and Gang Pan. 2020. LISA: A Learned Index Structure for Spatial Data. In SIGMOD Conference. ACM, 2119–2133.

[6] Vikram Nathan, Jialin Ding, Mohammad Alizadeh, and Tim Kraska. 2020. Learning Multi-Dimensional Indexes. In SIGMOD Conference. ACM, 985–1000.

