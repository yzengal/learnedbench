cmake_minimum_required(VERSION 3.16)

project(learnedbench)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()


set(CMAKE_CXX_FLAGS "-mbmi2 -O3 -DNDEBUG")
# set(CMAKE_CXX_FLAGS "-mbmi2 -Wall -Wextra")
# set(CMAKE_CXX_FLAGS_DEBUG "-g")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3")

option(PROFILE "whether starting heap profile" OFF)
option(RSMI "whether building RSMI" OFF)


set(BOOST_ROOT "/home/zengyuxiang/boost_1_80_0")
set(Boost_INCLUDE_DIR "/usr/local/include/boost")
set(Boost_LIBRARY_DIR "/usr/local/lib")
set(TORCH_PATH "/home/zengyuxiang/libtorch")
set(RSMI_PATH "/home/zengyuxiang/learnedbench/indexes/rsmi")
set(ANN_PATH "/home/zengyuxiang/learnedbench/indexes/ann_1.1.2")
SET(EXECUTABLE_OUTPUT_PATH "/home/zengyuxiang/learnedbench/build/bin")

# check dependencies
find_package(TPIE REQUIRED)
find_package(Boost 1.79.0 REQUIRED COMPONENTS program_options system filesystem)
find_package(GEOS REQUIRED)

if(PROFILE)
	find_package(benchmark REQUIRED)
	add_definitions(-DHEAP_PROFILE)
endif()

if(RSMI)
	# the RSMI benchmark needs to be compiled using cxx14
	# the dependency libtorch can be downloaded from 
	set(CMAKE_CXX_STANDARD 14)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)

	list(APPEND CMAKE_PREFIX_PATH TORCH_PATH)
	find_package(Torch REQUIRED PATHS ${TORCH_PATH} CONFIG)

	# RSMI library path

	include_directories(${Boost_INCLUDE_DIR} ${TPIE_INCLUDE_DIRS} ${RSMI_PATH})

	if(PROFILE)
		add_executable(bench_rsmi_hp bench/bench_rsmi.cpp)
		target_link_libraries(bench_rsmi_hp ${TPIE_LIBRARIES} ${TORCH_LIBRARIES} Boost::system Boost::filesystem pthread "${RSMI_PATH}/build/librsmi.a" -ltcmalloc -lsnappy)
	else()
		add_executable(bench_rsmi bench/bench_rsmi.cpp)
		target_link_libraries(bench_rsmi ${TPIE_LIBRARIES} ${TORCH_LIBRARIES} Boost::system Boost::filesystem pthread "${RSMI_PATH}/build/librsmi.a" -lsnappy)
	endif()
  
elseif(DATAGEN)
	include_directories(${Boost_INCLUDE_DIR} ${TPIE_INCLUDE_DIRS})

	add_executable(datagen utils/datagen.cpp)
	target_link_libraries(datagen ${TPIE_LIBRARIES} Boost::program_options pthread -lsnappy)

	add_executable(boxgen utils/boxgen.cpp)
	target_link_libraries(boxgen ${TPIE_LIBRARIES} Boost::program_options pthread -lsnappy)

	add_executable(boxshow utils/boxshow.cpp)
	target_link_libraries(boxshow ${TPIE_LIBRARIES} Boost::program_options pthread -lsnappy)
  
elseif(TESTBOX)
	# the other benchmarks are compiled using cxx17 (required by PGM-Index)
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)

	include_directories(${Boost_INCLUDE_DIR} ${TPIE_INCLUDE_DIRS})

	add_executable(bench2d_box bench/bench_box.cpp)
	target_link_libraries(bench2d_box ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread -ltcmalloc -lsnappy)
	target_compile_definitions(bench2d_box PUBLIC PARTITION_NUM=100)
  
elseif(TESTUPDATE)

	# the other benchmarks are compiled using cxx17 (required by PGM-Index)
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
	include_directories(${Boost_INCLUDE_DIR} ${TPIE_INCLUDE_DIRS} ${RSMI_PATH})
	list(APPEND CMAKE_PREFIX_PATH TORCH_PATH)
	find_package(Torch REQUIRED PATHS ${TORCH_PATH} CONFIG)

	add_executable(bench2dupdate_fs bench/bench_update.cpp)
	target_link_libraries(bench2dupdate_fs ${TPIE_LIBRARIES} Boost::program_options pthread -ltcmalloc -lsnappy)
	target_compile_definitions(bench2dupdate_fs PUBLIC PARTITION_NUM=60)

	add_executable(bench2dupdate_osm bench/bench_update.cpp)
	target_link_libraries(bench2dupdate_osm ${TPIE_LIBRARIES} Boost::program_options pthread -ltcmalloc -lsnappy)
	target_compile_definitions(bench2dupdate_osm PUBLIC PARTITION_NUM=250)

	add_executable(bench3dupdate_toronto bench/bench_update.cpp)
	target_link_libraries(bench3dupdate_toronto ${TPIE_LIBRARIES} Boost::program_options pthread -ltcmalloc -lsnappy)
	target_compile_definitions(bench3dupdate_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20)
  
	add_executable(bench_rsmi_update bench/bench_rsmi_update.cpp)
	target_link_libraries(bench_rsmi_update ${TPIE_LIBRARIES} ${TORCH_LIBRARIES} Boost::system Boost::filesystem pthread "${RSMI_PATH}/build/librsmi.a" -lsnappy)
  
  
elseif(TESTTUNE)

	# the other benchmarks are compiled using cxx17 (required by PGM-Index)
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
	include_directories(${Boost_INCLUDE_DIR} ${TPIE_INCLUDE_DIRS})
  
	add_executable(bench2dworkload_fs bench/bench_workload.cpp)
	target_link_libraries(bench2dworkload_fs ${TPIE_LIBRARIES} Boost::program_options pthread -ltcmalloc -lsnappy)
	target_compile_definitions(bench2dworkload_fs PUBLIC PARTITION_NUM=60)

	add_executable(bench2dworkload_osm bench/bench_workload.cpp)
	target_link_libraries(bench2dworkload_osm ${TPIE_LIBRARIES} Boost::program_options pthread -ltcmalloc -lsnappy)
	target_compile_definitions(bench2dworkload_osm PUBLIC PARTITION_NUM=250)

	add_executable(bench3dworkload_toronto bench/bench_workload.cpp)
	target_link_libraries(bench3dworkload_toronto ${TPIE_LIBRARIES} Boost::program_options pthread  -ltcmalloc -lsnappy)
	target_compile_definitions(bench3dworkload_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20)
	
	# add_executable(bench2dworkload_default bench/bench_workload.cpp)
    # target_link_libraries(bench2dworkload_default ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
    # target_compile_definitions(bench2dworkload_default PUBLIC PARTITION_NUM=100)
	
	# add_executable(bench4dworkload_default bench/bench_workload.cpp)
    # target_link_libraries(bench4dworkload_default ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
    # target_compile_definitions(bench4dworkload_default PUBLIC BENCH_DIM=4 PARTITION_NUM=10)

    # add_executable(bench6dworkload_default bench/bench_workload.cpp)
    # target_link_libraries(bench6dworkload_default ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
    # target_compile_definitions(bench6dworkload_default PUBLIC BENCH_DIM=6 PARTITION_NUM=5)

    # add_executable(bench8dworkload_default bench/bench_workload.cpp)
    # target_link_libraries(bench8dworkload_default ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
    # target_compile_definitions(bench8dworkload_default PUBLIC BENCH_DIM=8 PARTITION_NUM=4)
	
  
elseif(TESTJOIN)
	# the other benchmarks are compiled using cxx17 (required by PGM-Index)
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)

	include_directories(${Boost_INCLUDE_DIR} ${TPIE_INCLUDE_DIRS} ${RSMI_PATH} ${ANN_PATH})

	add_executable(bench2djoin_fs bench/bench_join.cpp)
	target_link_libraries(bench2djoin_fs ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
	target_compile_definitions(bench2djoin_fs PUBLIC PARTITION_NUM=60)

	add_executable(bench2djoin_osm bench/bench_join.cpp)
	target_link_libraries(bench2djoin_osm ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
	target_compile_definitions(bench2djoin_osm PUBLIC PARTITION_NUM=250)

	add_executable(bench3djoin_toronto bench/bench_join.cpp)
	target_link_libraries(bench3djoin_toronto ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
	target_compile_definitions(bench3djoin_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20)

	list(APPEND CMAKE_PREFIX_PATH TORCH_PATH)
	find_package(Torch REQUIRED PATHS ${TORCH_PATH} CONFIG)

	include_directories(${Boost_INCLUDE_DIR} ${TPIE_INCLUDE_DIRS} ${RSMI_PATH})
	add_executable(bench_rsmi bench/bench_rsmi.cpp)
	target_link_libraries(bench_rsmi ${TPIE_LIBRARIES} ${TORCH_LIBRARIES} Boost::system Boost::filesystem pthread "${RSMI_PATH}/build/librsmi.a" -lsnappy)
  
elseif(TESTRMI)

	# the other benchmarks are compiled using cxx17 (required by PGM-Index)
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)

	include_directories(${Boost_INCLUDE_DIR} ${TPIE_INCLUDE_DIRS})

	add_executable(bench2drmi_fs bench/bench_rmi.cpp)
	target_link_libraries(bench2drmi_fs ${TPIE_LIBRARIES} Boost::program_options pthread -ltcmalloc -lsnappy)
	target_compile_definitions(bench2drmi_fs PUBLIC PARTITION_NUM=60)

	add_executable(bench2drmi_osm bench/bench_rmi.cpp)
	target_link_libraries(bench2drmi_osm ${TPIE_LIBRARIES} Boost::program_options pthread -ltcmalloc -lsnappy)
	target_compile_definitions(bench2drmi_osm PUBLIC PARTITION_NUM=250)

	add_executable(bench3drmi_toronto bench/bench_rmi.cpp)
	target_link_libraries(bench3drmi_toronto ${TPIE_LIBRARIES} Boost::program_options pthread -ltcmalloc -lsnappy)
	target_compile_definitions(bench3drmi_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20)

	# add_executable(bench2drmi_default bench/bench_rmi.cpp)
	# target_link_libraries(bench2drmi_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
	# target_compile_definitions(bench2drmi_default PUBLIC PARTITION_NUM=100)

	# add_executable(bench4drmi_default bench/bench_rmi.cpp)
	# target_link_libraries(bench4drmi_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
	# target_compile_definitions(bench4drmi_default PUBLIC BENCH_DIM=4 PARTITION_NUM=10)

	# add_executable(bench6drmi_default bench/bench_rmi.cpp)
	# target_link_libraries(bench6drmi_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
	# target_compile_definitions(bench6drmi_default PUBLIC BENCH_DIM=6 PARTITION_NUM=5)

	# add_executable(bench8drmi_default bench/bench_rmi.cpp)
	# target_link_libraries(bench8drmi_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
	# target_compile_definitions(bench8drmi_default PUBLIC BENCH_DIM=8 PARTITION_NUM=4)
	
elseif(TESTFLOOD)
	# the other benchmarks are compiled using cxx17 (required by PGM-Index)
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)

	# ANN library path
	include_directories(${Boost_INCLUDE_DIR} ${TPIE_INCLUDE_DIRS})

	add_executable(bench2dflood_fs bench/bench_flood.cpp)
	target_link_libraries(bench2dflood_fs ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	target_compile_definitions(bench2dflood_fs PUBLIC PARTITION_NUM=60)

	add_executable(bench2dflood_osm bench/bench_flood.cpp)
	target_link_libraries(bench2dflood_osm ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	target_compile_definitions(bench2dflood_osm PUBLIC PARTITION_NUM=250)

	add_executable(bench3dflood_toronto bench/bench_flood.cpp)
	target_link_libraries(bench3dflood_toronto ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	target_compile_definitions(bench3dflood_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20)

	# add_executable(bench2dflood_default bench/bench_flood.cpp)
	# target_link_libraries(bench2dflood_default ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	# target_compile_definitions(bench2dflood_default PUBLIC PARTITION_NUM=100)

	# add_executable(bench4dflood_default bench/bench_flood.cpp)
	# target_link_libraries(bench4dflood_default ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	# target_compile_definitions(bench4dflood_default PUBLIC BENCH_DIM=4 PARTITION_NUM=10)

	# add_executable(bench6dflood_default bench/bench_flood.cpp)
	# target_link_libraries(bench6dflood_default ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	# target_compile_definitions(bench6dflood_default PUBLIC BENCH_DIM=6 PARTITION_NUM=5)

	# add_executable(bench8dflood_default bench/bench_flood.cpp)
	# target_link_libraries(bench8dflood_default ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	# target_compile_definitions(bench8dflood_default PUBLIC BENCH_DIM=8 PARTITION_NUM=4)

elseif(TESTOCTREE)

	# the other benchmarks are compiled using cxx17 (required by PGM-Index)
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)

	# ANN library path
	include_directories(${Boost_INCLUDE_DIR} ${TPIE_INCLUDE_DIRS})

	add_executable(bench2doctree_fs bench/bench_octree.cpp)
	target_link_libraries(bench2doctree_fs ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	target_compile_definitions(bench2doctree_fs PUBLIC PARTITION_NUM=60)

	add_executable(bench2doctree_osm bench/bench_octree.cpp)
	target_link_libraries(bench2doctree_osm ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	target_compile_definitions(bench2doctree_osm PUBLIC PARTITION_NUM=250)

	add_executable(bench3doctree_toronto bench/bench_octree.cpp)
	target_link_libraries(bench3doctree_toronto ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	target_compile_definitions(bench3doctree_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20)

	add_executable(bench2doctree_default bench/bench_octree.cpp)
	target_link_libraries(bench2doctree_default ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	target_compile_definitions(bench2doctree_default PUBLIC PARTITION_NUM=100)

	add_executable(bench4doctree_default bench/bench_octree.cpp)
	target_link_libraries(bench4doctree_default ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	target_compile_definitions(bench4doctree_default PUBLIC BENCH_DIM=4 PARTITION_NUM=10)

	add_executable(bench6doctree_default bench/bench_octree.cpp)
	target_link_libraries(bench6doctree_default ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	target_compile_definitions(bench6doctree_default PUBLIC BENCH_DIM=6 PARTITION_NUM=5)

	add_executable(bench8doctree_default bench/bench_octree.cpp)
	target_link_libraries(bench8doctree_default ${TPIE_LIBRARIES} Boost::program_options  pthread  -ltcmalloc -lsnappy)
	target_compile_definitions(bench8doctree_default PUBLIC BENCH_DIM=8 PARTITION_NUM=4)
		
else()

	# the other benchmarks are compiled using cxx17 (required by PGM-Index)
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)

	include_directories(${Boost_INCLUDE_DIR} ${TPIE_INCLUDE_DIRS} ${RSMI_PATH} ${ANN_PATH})

	if(PROFILE)
		add_executable(bench2d_hp bench/bench.cpp)
		target_link_libraries(bench2d_hp ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
	else()

		add_executable(bench2d_fs bench/bench.cpp)
		target_link_libraries(bench2d_fs ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
		target_compile_definitions(bench2d_fs PUBLIC PARTITION_NUM=60)

		add_executable(bench2d_osm bench/bench.cpp)
		target_link_libraries(bench2d_osm ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
		target_compile_definitions(bench2d_osm PUBLIC PARTITION_NUM=250)

		add_executable(bench2d_default bench/bench.cpp)
		target_link_libraries(bench2d_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
		target_compile_definitions(bench2d_default PUBLIC PARTITION_NUM=100)

		# add_executable(bench3d_toronto bench/bench.cpp)
		# target_link_libraries(bench3d_toronto ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
		# target_compile_definitions(bench3d_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20)

		# add_executable(bench4d_default bench/bench.cpp)
		# target_link_libraries(bench4d_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
		# target_compile_definitions(bench4d_default PUBLIC BENCH_DIM=4 PARTITION_NUM=10)

		# add_executable(bench6d_default bench/bench.cpp)
		# target_link_libraries(bench6d_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
		# target_compile_definitions(bench6d_default PUBLIC BENCH_DIM=6 PARTITION_NUM=5)

		# add_executable(bench8d_default bench/bench.cpp)
		# target_link_libraries(bench8d_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
		# target_compile_definitions(bench8d_default PUBLIC BENCH_DIM=8 PARTITION_NUM=4)

		# add_executable(bench10d_default bench/bench.cpp)
		# target_link_libraries(bench10d_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench10d_default PUBLIC BENCH_DIM=10 PARTITION_NUM=3)

		# add_executable(bench2d_1m_default bench/bench.cpp)
		# target_link_libraries(bench2d_1m_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_1m_default PUBLIC PARTITION_NUM=23)

		# add_executable(bench2d_10m_default bench/bench.cpp)
		# target_link_libraries(bench2d_10m_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_10m_default PUBLIC PARTITION_NUM=70)

		# add_executable(bench2d_50m_default bench/bench.cpp)
		# target_link_libraries(bench2d_50m_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_50m_default PUBLIC PARTITION_NUM=158)

		# add_executable(bench2d_100m_default bench/bench.cpp)
		# target_link_libraries(bench2d_100m_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_100m_default PUBLIC PARTITION_NUM=223)

		# add_executable(bench2d_e4_default bench/bench.cpp)
		# target_link_libraries(bench2d_e4_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_e4_default PUBLIC PARTITION_NUM=100 INDEX_ERROR_THRESHOLD=4)

		# add_executable(bench2d_e16_default bench/bench.cpp)
		# target_link_libraries(bench2d_e16_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_e16_default PUBLIC PARTITION_NUM=100 INDEX_ERROR_THRESHOLD=16)

		# add_executable(bench2d_e64_default bench/bench.cpp)
		# target_link_libraries(bench2d_e64_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_e64_default PUBLIC PARTITION_NUM=100 INDEX_ERROR_THRESHOLD=64)

		# add_executable(bench2d_e256_default bench/bench.cpp)
		# target_link_libraries(bench2d_e256_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_e256_default PUBLIC PARTITION_NUM=100 INDEX_ERROR_THRESHOLD=256)

		# add_executable(bench2d_e1024_default bench/bench.cpp)
		# target_link_libraries(bench2d_e1024_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_e1024_default PUBLIC PARTITION_NUM=100 INDEX_ERROR_THRESHOLD=1024)

		# add_executable(bench2d_e4_fs bench/bench.cpp)
		# target_link_libraries(bench2d_e4_fs ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_e4_fs PUBLIC PARTITION_NUM=60 INDEX_ERROR_THRESHOLD=4)

		# add_executable(bench2d_e16_fs bench/bench.cpp)
		# target_link_libraries(bench2d_e16_fs ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_e16_fs PUBLIC PARTITION_NUM=60 INDEX_ERROR_THRESHOLD=16)

		# add_executable(bench2d_e64_fs bench/bench.cpp)
		# target_link_libraries(bench2d_e64_fs ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_e64_fs PUBLIC PARTITION_NUM=60 INDEX_ERROR_THRESHOLD=64)

		# add_executable(bench2d_e256_fs bench/bench.cpp)
		# target_link_libraries(bench2d_e256_fs ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_e256_fs PUBLIC PARTITION_NUM=60 INDEX_ERROR_THRESHOLD=256)

		# add_executable(bench2d_e1024_fs bench/bench.cpp)
		# target_link_libraries(bench2d_e1024_fs ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench2d_e1024_fs PUBLIC PARTITION_NUM=60 INDEX_ERROR_THRESHOLD=1024)

		# add_executable(bench3d_e4_toronto bench/bench.cpp)
		# target_link_libraries(bench3d_e4_toronto ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench3d_e4_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20 INDEX_ERROR_THRESHOLD=4)

		# add_executable(bench3d_e16_toronto bench/bench.cpp)
		# target_link_libraries(bench3d_e16_toronto ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench3d_e16_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20 INDEX_ERROR_THRESHOLD=16)

		# add_executable(bench3d_e64_toronto bench/bench.cpp)
		# target_link_libraries(bench3d_e64_toronto ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench3d_e64_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20 INDEX_ERROR_THRESHOLD=64)

		# add_executable(bench3d_e256_toronto bench/bench.cpp)
		# target_link_libraries(bench3d_e256_toronto ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench3d_e256_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20 INDEX_ERROR_THRESHOLD=256)

		# add_executable(bench3d_e1024_toronto bench/bench.cpp)
		# target_link_libraries(bench3d_e1024_toronto ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a")
		# target_compile_definitions(bench3d_e1024_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20 INDEX_ERROR_THRESHOLD=1024)

		# add_executable(bench2d_default bench/bench.cpp)
		# target_link_libraries(bench2d_default ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -ltcmalloc -lsnappy)
		# target_compile_definitions(bench2d_default PUBLIC PARTITION_NUM=100)

		# add_executable(bench3d_toronto bench/bench.cpp)
		# target_link_libraries(bench3d_toronto ${TPIE_LIBRARIES} Boost::program_options GEOS::geos pthread "${ANN_PATH}/lib/libANN.a" -lsnappy)
		# target_compile_definitions(bench3d_toronto PUBLIC BENCH_DIM=3 PARTITION_NUM=20)

	endif()
	
endif()

