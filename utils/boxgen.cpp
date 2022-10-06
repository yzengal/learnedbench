#include <boost/program_options.hpp>
#include <ios>
#include <iostream>
#include <string>
#include <tpie/tpie.h>
#include "type.hpp"
#include "datautils.hpp"

namespace po = boost::program_options;

int main(int argc, char const *argv[]) {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("task,t", po::value<std::string>(), "tasks: gen_box")
        ("fname,f", po::value<std::string>(), "file name of generated data")
        ("num,n", po::value<int>(), "number of boxes to be generated")
        ("dim,d", po::value<int>(), "dimension of boxes to be generated")
        ("size,s", po::value<double>(), "size scale factor of the box")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    tpie::tpie_init();

    if (vm.count("task")) {
        std::string task = vm["task"].as<std::string>();

        if (task.compare("gen_box") == 0) {
            if (vm.count("fname") && vm.count("num") && vm.count("dim") && vm.count("size")) {
                std::string fname = vm["fname"].as<std::string>();
                int n = vm["num"].as<int>();
                int d = vm["dim"].as<int>();
                double s = vm["size"].as<double>();

				// generate box data Size(max_size)
				std::cout << "Generate " << n << " * " << d << "D Boxes " << "SIZE(" << s << ") data to file: " << fname << std::endl;
				bench::utils::gen_box(fname, n, d, s);
				tpie::tpie_finish();
				return 0;
	
            } else {
                std::cout << "Please provide fname, num, dim, and size." << std::endl;
                tpie::tpie_finish();
                return 1;
            }
			
        } else {
            std::cout << "Arg --task is in [gen_box]." << std::endl;
            tpie::tpie_finish();
            return 1;
        }
		
    } else {
        std::cout << "Specify a task." << std::endl;
        std::cout << desc << std::endl;
        tpie::tpie_finish();
        return 1;
    }

    tpie::tpie_finish();
    return 0;
}
