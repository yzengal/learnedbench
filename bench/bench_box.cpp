#include "../utils/datautils.hpp"
#include "../utils/common.hpp"

#include "../indexes/revision-box/nonlearned_index.hpp"
#include "../indexes/revision-box/learned_index.hpp"

#include "query_box.hpp"

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
using Boxes = std::vector<box_t<BENCH_DIM>>;

// non-learned indices
using RTree = bench::index::RTree<BENCH_DIM>;
using RStarTree = bench::index::RStarTree<BENCH_DIM>;
using QDTree = bench::index::OctreeBox<BENCH_DIM>;

// grid indices
using UG = bench::index::UG<BENCH_DIM, PARTITION_NUM>;
using EDG = bench::index::EDG<BENCH_DIM, PARTITION_NUM>;

// linear scan
using FS = bench::index::FullScan<BENCH_DIM>;    

struct IndexSet {
    RTree*     rtree;
    RStarTree* rstartree;
    QDTree*    qdtree;
    UG*        ug;
    EDG*       edg;
    FS*        fs;

    IndexSet() : 
        rtree(nullptr), 
        rstartree(nullptr),
        qdtree(nullptr),
        ug(nullptr),
        edg(nullptr),
        fs(nullptr), {}

    ~IndexSet() {
        delete rtree;
        delete rstartree;
        delete qdtree;
        delete ug;
        delete edg;
        delete fs;
    }
};


static void build_index(IndexSet& idx_set, const std::string& idx_name, Boxes& boxes) {
    if (idx_name.compare("rtree") == 0) {
        idx_set.rtree = new RTree(boxes);
        return;
    }

    if (idx_name.compare("rstar") == 0) {
        idx_set.rstartree = new RStarTree(boxes);
        return;
    }

    if (idx_name.compare("qdtree") == 0) {
        idx_set.qdtree = new QDTree(boxes);
        return;
    }

    if (idx_name.compare("ug") == 0) {
        idx_set.ug = new UG(boxes);
        return;
    }

    if (idx_name.compare("edg") == 0) {
        idx_set.edg = new EDG(boxes);
        return;
    }

    if (idx_name.compare("fs") == 0) {
        idx_set.fs = new FS(boxes);
        return;
    }

    std::cout << "index name should be one of [rtree, rstar, qdtree, ug, edg, fs]" << std::endl;
    exit(0);
}


int main(int argc, char **argv) {
    assert(argc >= 4);

    std::string index = argv[1]; // index name
    std::string fname = argv[2]; // data file name
    size_t N = std::stoi(argv[3]); // dataset size
    std::string mode = argv[4]; // bench mode {"range"}

    std::cout << "====================================" << std::endl;
    std::cout << "Load data: " << fname << std::endl;

    Boxes boxes;
    bench::utils::read_boxes(boxes, fname, N);
    IndexSet idx_set;

#ifdef HEAP_PROFILE
    build_index(idx_set, index, boxes);
    return 0;
#endif

#ifndef HEAP_PROFILE
    
    auto range_queries = bench::query::sample_box_queries();

    build_index(idx_set, index, boxes);

    if (index.compare("rtree") == 0) {
        assert(idx_set.rtree != nullptr);
        if (mode.compare("range") == 0) {
            bench::query::batch_range_queries(*(idx_set.rtree), range_queries);
            return 0;
        }
    }

    if (index.compare("rstar") == 0) {
        assert(idx_set.rstartree != nullptr);
        if (mode.compare("range") == 0) {
            bench::query::batch_range_queries(*(idx_set.rstartree), range_queries);
            return 0;
        }
    }

    if (index.compare("qdtree") == 0) {
        assert(idx_set.qdtree != nullptr);
        if (mode.compare("range") == 0) {
            bench::query::batch_range_queries(*(idx_set.qdtree), range_queries);
            return 0;
        }
    }
    
    if (index.compare("ug") == 0) {
        assert(idx_set.ug != nullptr);
        if (mode.compare("range") == 0) {
            bench::query::batch_range_queries(*(idx_set.ug), range_queries);
            return 0;
        }
    }

    if (index.compare("edg") == 0) {
        assert(idx_set.edg != nullptr);
        if (mode.compare("range") == 0) {
            bench::query::batch_range_queries(*(idx_set.edg), range_queries);
            return 0;
        }
    }

    if (index.compare("fs") == 0) {
        assert(idx_set.fs != nullptr);
        if (mode.compare("range") == 0) {
            bench::query::batch_range_queries(*(idx_set.fs), range_queries);
            return 0;
        }
    }

#endif

	return 0;
}
