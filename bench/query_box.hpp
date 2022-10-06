#pragma once

#include <cstddef>
#include <random>
#include <map>
#include <algorithm>
#include <cmath>

#include "../utils/type.hpp"
#include "../indexes/nonlearned/fullscan.hpp"


namespace bench { namespace query {

// sample queries from data
// sample point queries
template<size_t dim>
static std::vector<box_t<dim>> sample_box_queries(const int MAXL = 10000, size_t s=100, size_t area=0.01) {
    // seed the generator
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_int_distribution<> dis(0, (int)MAXL);
	const int L = std::sqrt(area) * MAXL;
	

    // generate random indices
    std::vector<box_t<dim>> samples;
    samples.reserve(s);
	point_t<dim> mn_corner, mx_corner;
	
	for (size_t i=0; i<s; ++i) {
		for (int j=0; j<d; ++j) {
			mn_corner[j] = dis(gen);
			mx_corner[j] = mn_corner[j] + L;
		}
		samples.emplace_back(box_t<dim>(mn_corner, mx_corner));
	}

    return samples;
}

template<class Index, size_t Dim>
static void batch_range_queries(Index& index, double& max_side, std::vector<box_t<Dim>>& box_queries) {
    for (auto& box : box_queries) {
		index.range_query(box);
	}
	
	std::cout << "max_side = " << max_side << " Avg. Time: " << index.get_avg_range_time() << " [us]" << std::endl;
	index.reset_timer();
}



}
}


