#include "../utils/datautils.hpp"
#include "../utils/common.hpp"

#include "../indexes/nonlearned/Octree.hpp"

#include "query.hpp"

#include <cstddef>
#include <string>
#include <stdlib.h>

#ifndef BENCH_DIM
#define BENCH_DIM 2
#endif

#ifndef PARTITION_NUM
#define PARTITION_NUM 10
#endif

#ifndef INDEX_ERROR_THRESHOLD 
#define INDEX_ERROR_THRESHOLD 64
#endif

using Point = point_t<BENCH_DIM>;
using Points = std::vector<point_t<BENCH_DIM>>;
using Box = box_t<BENCH_DIM>;

// learned indices
// note RSMI needs to be compiled individually due to the dependency of an old version of libtorch
// see CMakeLists.txt for more details 
using Octree = bench::index::Octree<BENCH_DIM, 128>;     

struct IndexSet {
    Octree*     octree;

    IndexSet() : 
        octree(nullptr) {}

    ~IndexSet() {
        delete octree;
    }
};


static void build_index(IndexSet& idx_set, const std::string& idx_name, Points& points) {
    
    if (idx_name.compare("octree") == 0) {
        idx_set.octree = new Octree(points);
        return;
    }

    std::cout << "index name should be one of [octree]" << std::endl;
    exit(0);
}


int main(int argc, char **argv) {
    assert(argc >= 4);

    std::string index = argv[1]; // index name
    std::string fname = argv[2]; // data file name
    size_t N = std::stoi(argv[3]); // dataset size
    std::string mode = argv[4]; // bench mode {"range", "knn", "all"}

    std::cout << "====================================" << std::endl;
    std::cout << "Load data: " << fname << std::endl;

    Points points;
    bench::utils::read_points(points, fname, N);
    IndexSet idx_set;
    
    auto range_queries = bench::query::sample_range_queries(points);
    auto knn_queries = bench::query::sample_knn_queries(points);

    build_index(idx_set, index, points);

    if (index.compare("octree") == 0) {
        assert(idx_set.octree != nullptr);
        if (mode.compare("range") == 0) {
            bench::query::batch_range_queries(*(idx_set.octree), range_queries);
            return 0;
        }
        if (mode.compare("knn") == 0) {
            bench::query::batch_knn_queries(*(idx_set.octree), knn_queries);
            return 0;
        }
        if (mode.compare("all") == 0) {
            bench::query::batch_range_queries(*(idx_set.octree), range_queries);
            bench::query::batch_knn_queries(*(idx_set.octree), knn_queries);
            return 0;
        }
    }
}
