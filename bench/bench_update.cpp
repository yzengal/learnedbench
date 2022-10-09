#include "../utils/datautils.hpp"
#include "../utils/common.hpp"

#include "../indexes/revision-update/nonlearned_index.hpp"
#include "../indexes/revision-update/learned_index.hpp"

#include "query_update.hpp"

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

// non-learned indices
using RTree = bench::index::RTree<BENCH_DIM>;
using RStarTree = bench::index::RStarTree<BENCH_DIM>;
using KDTree = bench::index::KDTree<BENCH_DIM>;
using QDTree = bench::index::Octree<BENCH_DIM>;

// grid indices
using UG = bench::index::UG<BENCH_DIM, PARTITION_NUM>;
using EDG = bench::index::EDG<BENCH_DIM, PARTITION_NUM>;

// learned indices
// note RSMI needs to be compiled individually due to the dependency of an old version of libtorch
// see CMakeLists.txt for more details 
using IFI = bench::index::IFIndex<BENCH_DIM>;
using Lisa = bench::index::LISA<BENCH_DIM, PARTITION_NUM, INDEX_ERROR_THRESHOLD>;   
using Lisa2 = bench::index::LISA2<BENCH_DIM, PARTITION_NUM, INDEX_ERROR_THRESHOLD>;         

struct IndexSet {
    RTree*     rtree;
    RStarTree* rstartree;
    KDTree*    kdtree;
    QDTree*    qdtree;
    UG*        ug;
    EDG*       edg;
    IFI*       ifi;
    Lisa*      lisa;
    Lisa2*      lisa2;

    IndexSet() : 
        rtree(nullptr), 
        rstartree(nullptr),
        kdtree(nullptr),
        qdtree(nullptr),
        ug(nullptr),
        edg(nullptr),
        ifi(nullptr),
        lisa(nullptr),
        lisa2(nullptr) {}

    ~IndexSet() {
        delete rtree;
		delete rstartree;
        delete kdtree;
        delete qdtree;
        delete ug;
        delete edg;
        delete ifi;
        delete lisa;
		delete lisa2;
    }
};


static void build_index(IndexSet& idx_set, const std::string& idx_name, Points& points, Points& pointsAll) {
    if (idx_name.compare("rtree") == 0) {
        idx_set.rtree = new RTree(points);
        return;
    }

    if (idx_name.compare("rstar") == 0) {
        idx_set.rstartree = new RStarTree(points);
        return;
    }

    if (idx_name.compare("kdtree") == 0) {
        idx_set.kdtree = new KDTree(points);
        return;
    }

    if (idx_name.compare("qdtree") == 0) {
        idx_set.qdtree = new QDTree(points, pointsAll);
        return;
    }

    if (idx_name.compare("ug") == 0) {
        idx_set.ug = new UG(points);
        return;
    }

    if (idx_name.compare("edg") == 0) {
        idx_set.edg = new EDG(points);
        return;
    }
    if (idx_name.compare("ifi") == 0) {
        idx_set.ifi = new IFI(points);
        return;
    }

    if (idx_name.compare("lisa") == 0) {
        idx_set.lisa = new Lisa(points);
        return;
    }
	
    if (idx_name.compare("lisa2") == 0) {
        idx_set.lisa2 = new Lisa2(points);
        return;
    }

    std::cout << "index name should be one of [rtree, rstar, kdtree, qdtree, ug, edg, ifi, lisa, lisa2]" << std::endl;
    exit(0);
}


int main(int argc, char **argv) {
    assert(argc >= 4);

    std::string index = argv[1]; // index name
    std::string fname = argv[2]; // data file name
    size_t N = std::stoi(argv[3]); // dataset size
    std::string mode = argv[4]; // bench mode {"update", "knn", "all"}

    std::cout << "====================================" << std::endl;
    std::cout << "Load data: " << fname << std::endl;

    Points points;
    bench::utils::read_points(points, fname, N);
    IndexSet idx_set;

	size_t M = N*0.8;
	size_t UPDATE_ROUND = 10;
	
	Points points_backup(points.begin(), points.begin()+M);
    auto update_queries = bench::query::sample_update_queries(M, N, UPDATE_ROUND);

    build_index(idx_set, index, points_backup, points);

    if (index.compare("rtree") == 0) {
        assert(idx_set.rtree != nullptr);
        if (mode.compare("update") == 0) {
            bench::query::batch_update_queries(*(idx_set.rtree), points, points_backup, update_queries);
            return 0;
        }
    }

    if (index.compare("rstar") == 0) {
        assert(idx_set.rstartree != nullptr);
        if (mode.compare("update") == 0) {
            bench::query::batch_update_queries(*(idx_set.rstartree), points, points_backup, update_queries);
            return 0;
        }
    }

    if (index.compare("kdtree") == 0) {
        assert(idx_set.kdtree != nullptr);
        if (mode.compare("update") == 0) {
           bench::query::batch_update_queries_byID(*(idx_set.kdtree), points, points_backup, update_queries);
            return 0;
        }
    }

    if (index.compare("qdtree") == 0) {
        assert(idx_set.qdtree != nullptr);
        if (mode.compare("update") == 0) {
            bench::query::batch_update_queries_byID(*(idx_set.qdtree), points, points_backup, update_queries);
            return 0;
        }
    }
    
    if (index.compare("ug") == 0) {
        assert(idx_set.ug != nullptr);
        if (mode.compare("update") == 0) {
            bench::query::batch_update_queries(*(idx_set.ug), points, points_backup, update_queries);
            return 0;
        }
    }

    if (index.compare("edg") == 0) {
        assert(idx_set.edg != nullptr);
        if (mode.compare("update") == 0) {
            bench::query::batch_update_queries(*(idx_set.edg), points, points_backup, update_queries);
            return 0;
        }
    }

    if (index.compare("ifi") == 0) {
        assert(idx_set.ifi != nullptr);
        if (mode.compare("update") == 0) {
            bench::query::batch_update_queries(*(idx_set.ifi), points, points_backup, update_queries);
            return 0;
        }
    }

    if (index.compare("lisa") == 0) {
        assert(idx_set.lisa != nullptr);
        if (mode.compare("update") == 0) {
            bench::query::batch_update_queries(*(idx_set.lisa), points, points_backup, update_queries);
            return 0;
        }
    }
	
    if (index.compare("lisa2") == 0) {
        assert(idx_set.lisa2 != nullptr);
        if (mode.compare("update") == 0) {
            bench::query::batch_update_queries(*(idx_set.lisa2), points, points_backup, update_queries);
            return 0;
        }
    }

	return 0;
}
