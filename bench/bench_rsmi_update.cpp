#include "../utils/datautils.hpp"
#include "../utils/common.hpp"

#include "../indexes/revision-update/rsmi.hpp"
#include "query_update.hpp"

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#ifndef BENCH_DIM
#define BENCH_DIM 2
#endif

using Point = point_t<BENCH_DIM>;
using Box = box_t<BENCH_DIM>;
using Points = std::vector<point_t<BENCH_DIM>>;

const std::string MODEL_PATH = "/data/yuxiang/model_path/";

// extract filename from a path
std::string get_filename(const std::string& path) {
    auto idx = path.find_last_of('/') + 1;
    return path.substr(idx);
}


int main(int argc, char **argv) {
    assert(argc >= 4);

    std::string index = argv[1]; // index name
    std::string fname = argv[2]; // data file name
    size_t N = std::stoi(argv[3]); // dataset size
    std::string mode = argv[4]; // bench mode

    std::cout << "====================================" << std::endl;
    std::cout << "Load data: " << fname << std::endl;
	
	if (index.compare("rsmi") != 0) {
        return -1;
    }

    Points points;
    bench::utils::read_points(points, fname, N);
	
	size_t M = N*0.8;
	size_t UPDATE_ROUND = 10;
	
	Points points_backup(points.begin(), points.begin()+M);
    auto update_queries = bench::query::sample_update_queries(M, N, UPDATE_ROUND);

    std::string pure_fname = get_filename(fname);
    std::string model_path = MODEL_PATH + pure_fname;
    model_path.push_back('/');

    
    if (! boost::filesystem::is_directory(model_path)) {
        boost::filesystem::create_directory(model_path);
    }
    
    std::cout << "Will save/load model to/from " << model_path << std::endl;

    bench::index::RSMIWrapper<BENCH_DIM> rsmi(points, model_path);

    if (mode.compare("update") == 0) {
		bench::query::batch_update_queries(rsmi, points, points_backup, update_queries);
		return 0;
	}
	
	return -1;
}

