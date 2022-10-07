#include "../utils/datautils.hpp"
#include "../utils/common.hpp"

#include "../indexes/revision-rmi/learned_index.hpp"

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
using MLI = bench::index::MLIndex<BENCH_DIM, INDEX_ERROR_THRESHOLD>;
using Flood = bench::index::Flood<BENCH_DIM, PARTITION_NUM, INDEX_ERROR_THRESHOLD>;
using Lisa = bench::index::LISA2<BENCH_DIM, PARTITION_NUM, INDEX_ERROR_THRESHOLD>;          

struct IndexSet {
    MLI*       mli;
    Flood*     flood;
    Lisa*      lisa;

    IndexSet() : 
        mli(nullptr),
        flood(nullptr),
        lisa(nullptr) {}

    ~IndexSet() {
        delete mli;
        delete flood;
        delete lisa;
    }
};


static void build_index(IndexSet& idx_set, const std::string& idx_name, Points& points) {
    

    if (idx_name.compare("mli") == 0) {
        idx_set.mli = new MLI(points);
        return;
    }

    if (idx_name.compare("flood") == 0) {
        idx_set.flood = new Flood(points);
        return;
    }

    if (idx_name.compare("lisa") == 0) {
        idx_set.lisa = new Lisa(points);
        return;
    }

    std::cout << "index name should be one of [mli, flood, lisa]" << std::endl;
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
    
    if (index.compare("mli") == 0) {
        assert(idx_set.mli != nullptr);
        if (mode.compare("range") == 0) {
            bench::query::batch_range_queries(*(idx_set.mli), range_queries);
            return 0;
        }
        if (mode.compare("knn") == 0) {
            bench::query::batch_knn_queries(*(idx_set.mli), knn_queries);
            return 0;
        }
        if (mode.compare("all") == 0) {
            bench::query::batch_range_queries(*(idx_set.mli), range_queries);
            bench::query::batch_knn_queries(*(idx_set.mli), knn_queries);
            return 0;
        }
    }

    if (index.compare("flood") == 0) {
        assert(idx_set.flood != nullptr);
        if (mode.compare("range") == 0 || mode.compare("all") == 0) {
            bench::query::batch_range_queries(*(idx_set.flood), range_queries);
            return 0;
        }
    }

    if (index.compare("lisa") == 0) {
        assert(idx_set.lisa != nullptr);
        if (mode.compare("range") == 0) {
            bench::query::batch_range_queries(*(idx_set.lisa), range_queries);
            return 0;
        }
        if (mode.compare("knn") == 0) {
            bench::query::batch_knn_queries(*(idx_set.lisa), knn_queries);
            return 0;
        }
        if (mode.compare("all") == 0) {
            bench::query::batch_range_queries(*(idx_set.lisa), range_queries);
            bench::query::batch_knn_queries(*(idx_set.lisa), knn_queries);
            return 0;
        }
    }
}
