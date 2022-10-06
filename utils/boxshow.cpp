#include <boost/program_options.hpp>
#include <ios>
#include <iostream>
#include <algorithm>
#include <string>
#include "type.hpp"
#include "datautils.hpp"

#ifndef BENCH_DIM
#define BENCH_DIM 2
#endif

namespace po = boost::program_options;
using Point = point_t<BENCH_DIM>;
using Points = std::vector<point_t<BENCH_DIM>>;
using Box = box_t<BENCH_DIM>;
using Boxes = std::vector<Box>;

int main(int argc, char const *argv[]) {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("task,t", po::value<std::string>(), "tasks: show_box")
        ("fname,f", po::value<std::string>(), "file name of generated data")
        ("num,n", po::value<int>(), "number of boxes to be generated")
        ("dim,d", po::value<int>(), "dimension of boxes to be generated")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }
	
	Boxes boxes;

    if (vm.count("task")) {
        std::string task = vm["task"].as<std::string>();

        if (task.compare("show_box") == 0) {
            if (vm.count("fname") && vm.count("num") && vm.count("dim")) {
                std::string fname = vm["fname"].as<std::string>();
                const int n = vm["num"].as<int>();
                const int d = vm["dim"].as<int>();
				
				std::cout << "[BEGIN] read_boxes() " << fname << std::endl;
				bench::utils::read_boxes(boxes, fname, n);
				std::cout << "[END] read_boxes() " << fname << std::endl;	
			
				std::vector<int> L;
				L.reserve(n*d);
				for (auto box : boxes) {
					Point mn_corner = box.min_corner(), mx_corner = box.max_corner();
					for (int i=0; i<d; ++i) {
						L.emplace_back(mx_corner[i]-mn_corner[i]);
					}
				}
				sort(L.begin(), L.end());
				L.erase(unique(L.begin(), L.end()), L.end());
				std::cout << "L.size() = " << L.size() << std::endl;
				for (int i=0,sz=std::min((size_t)10, L.size()); i<sz; ++i)
					std::cout << L[i] << " ";
				std::cout << std::endl;
				
				return 0;
	
            } else {
                std::cout << "Please provide fname, num, and dim." << std::endl;
                return 1;
            }
			
        } else {
            std::cout << "Arg --task is in [show_box]." << std::endl;
            return 1;
        }
		
    } else {
        std::cout << "Specify a task." << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }

    return 0;
}
