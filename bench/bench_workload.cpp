#include "../utils/datautils.hpp"
#include "../utils/common.hpp"

#include "../indexes/learned/flood.hpp"

#include "query_workload.hpp"

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
using Flood = bench::index::Flood<BENCH_DIM, 128, INDEX_ERROR_THRESHOLD>;
Points points;

static size_t tune_index(Points& points) {
	std::cout << "tune_index" << std::endl;
	
	double minRes = -1.0;
	size_t BEST_DIM = -1;
	Points pointsTmp = points;
	
	for (size_t dimIdx=0; dimIdx<BENCH_DIM; ++dimIdx) {
		const size_t dim = dimIdx;
		std::cout << "dim: " << dim << std::endl;
		for (int i=0; i<pointsTmp.size(); ++i) {
			std::swap(pointsTmp[i][BENCH_DIM-1], pointsTmp[i][dim]);
		}
		std::vector<std::pair<box_t<BENCH_DIM>, size_t>> queriesTmp = bench::query::sample_range_queries<BENCH_DIM>(pointsTmp);
		std::cout << "swap points" << std::endl;
		
		// for (auto& query : queriesTmp) {
			// box_t<BENCH_DIM> box = query.first;
			// Point mn_corner = box.min_corner();
			// Point mx_corner = box.max_corner();
			// std::swap(mn_corner[BENCH_DIM-1], mn_corner[dim]);
			// std::swap(mx_corner[BENCH_DIM-1], mx_corner[dim]);
			// box_t<BENCH_DIM> newBox(mn_corner, mx_corner);
			// query.first = newBox;
		// }
		// std::cout << "swap queries" << std::endl;
			
		Flood flood(pointsTmp);
		std::cout << "build Flood" << std::endl;
		
		double res = bench::query::batch_range_queries(flood, queriesTmp);
		if (minRes<0 || res<minRes) {
			minRes = res;
			BEST_DIM = dim;
		}
		std::cout << "TEST dimension = " << dim << " RESULT = " << std::fixed << std::setprecision(6) << res << " BEST_DIM = " << BEST_DIM << std::endl;
		
		for (int i=0; i<pointsTmp.size(); ++i) {
			std::swap(pointsTmp[i][BENCH_DIM-1], pointsTmp[i][dim]);
		}
		
	}
	
	assert(BEST_DIM>=0 && BEST_DIM<BENCH_DIM);
	return BEST_DIM;
}


static void test_index(Points& points, const size_t dim, std::vector<std::pair<box_t<BENCH_DIM>, size_t>>& range_queries) {
	Points& pointsTmp = points;
	
	std::vector<std::pair<box_t<BENCH_DIM>, size_t>> queriesTmp = range_queries;
	for (int i=0; i<pointsTmp.size(); ++i) {
		std::swap(pointsTmp[i][BENCH_DIM-1], pointsTmp[i][dim]);
	}
	for (auto& query : queriesTmp) {
		box_t<BENCH_DIM> box = query.first;
		Point mn_corner = box.min_corner();
		Point mx_corner = box.max_corner();
		std::swap(mn_corner[BENCH_DIM-1], mn_corner[dim]);
		std::swap(mx_corner[BENCH_DIM-1], mx_corner[dim]);
		box_t<BENCH_DIM> newBox(mn_corner, mx_corner);
		query.first = newBox;
	}
		
	std::cout << "\n\n\n" << std::endl; 
	Flood flood(pointsTmp);
	double res = bench::query::batch_range_queries(flood, queriesTmp);
}

int main(int argc, char **argv) {
    assert(argc >= 4);

    std::string index = argv[1]; // index name
    std::string fname = argv[2]; // data file name
    size_t N = std::stoi(argv[3]); // dataset size
    std::string mode = argv[4]; // bench mode {"range", "all"}

    std::cout << "====================================" << std::endl;
    std::cout << "Load data: " << fname << std::endl;
	std::cout << "Dim: " << BENCH_DIM << std::endl;
	
    bench::utils::read_points<BENCH_DIM>(points, fname, N);
	
	std::cout << "read_points" << std::endl;
    
    std::vector<std::pair<box_t<BENCH_DIM>, size_t>> range_queries = bench::query::sample_range_queries<BENCH_DIM>(points);
	
	std::cout << "sample_range_queries" << std::endl;
	
    // if (index.compare("flood") == 0) {
		size_t bestDim = tune_index(points);
        if (mode.compare("range") == 0 || mode.compare("all") == 0) {
            test_index(points, bestDim, range_queries);
            return 0;
        }
    // }
	
	return -1;
}
