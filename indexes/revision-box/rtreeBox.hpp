#pragma once

#include "../../utils/type.hpp"
#include "../../utils/common.hpp"
#include "base_index.hpp"

#include <boost/geometry/index/parameters.hpp>
#include <boost/geometry/index/predicates.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <cstddef>
#include <iterator>
#include <chrono>


namespace bgi = boost::geometry::index;

namespace bench { namespace index {

template<size_t dim, size_t MaxElements=128>
class RTree : public BaseIndex {

using Point = point_t<dim>;
using Box = box_t<dim>;
using Boxes = std::vector<Box>;
using rtree_t = bgi::rtree<Box, bgi::linear<MaxElements>>;

public:
RTree(Boxes& boxes) {
    std::cout << "Construct R-tree " << "MaxElements=" << MaxElements << std::endl;

    auto start = std::chrono::steady_clock::now();

    // construct r-tree using packing algorithm
    rtree = new rtree_t(boxes.begin(), boxes.end());

    auto end = std::chrono::steady_clock::now();
    build_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Build Time: " << get_build_time() << " [ms]" << std::endl;
    std::cout << "Index Size: " << index_size() << " Bytes" << std::endl;
}

~RTree() {
    delete this->rtree;
}


Boxes range_query(Box& box) {
    auto start = std::chrono::steady_clock::now();
    
	Boxes results;
    rtree->query(bgi::intersects(box), std::back_inserter(results));
    
	auto end = std::chrono::steady_clock::now();
    range_count++;
    range_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    return results;
}

inline size_t count() {
    return rtree->size();
}

inline size_t index_size() {
    return bench::common::get_boost_rtree_statistics(*rtree);
}

private:
rtree_t* rtree;
};

template<size_t dim, size_t MaxElements=128>
class RStarTree : public BaseIndex {

using Point = point_t<dim>;
using Box = box_t<dim>;
using Boxes = std::vector<Box>;
using rtree_t = bgi::rtree<Box, bgi::rstar<MaxElements>>;

public:
RStarTree(Boxes& boxes) {
    std::cout << "Construct R*-tree " << "MaxElements=" << MaxElements << std::endl;

    auto start = std::chrono::steady_clock::now();

    rtree = new rtree_t(boxes.begin(), boxes.end());

    // for (auto& p : points) {
        // rtree->insert(p);
    // }

    auto end = std::chrono::steady_clock::now();
    build_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Build Time: " << get_build_time() << " [ms]" << std::endl;
    std::cout << "Index Size: " << index_size() << " Bytes" << std::endl;
}

~RStarTree() {
    delete this->rtree;
}

Boxes range_query(Box& box) {
    auto start = std::chrono::steady_clock::now();
	
    Boxes results;
    rtree->query(bgi::intersects(box), std::back_inserter(results));
    
	auto end = std::chrono::steady_clock::now();
    range_count++;
    range_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    return results;
}


inline size_t count() {
    return rtree->size();
}

inline size_t index_size() {
    return bench::common::get_boost_rtree_statistics(*rtree);
}


private:
rtree_t* rtree;

};


}
}
