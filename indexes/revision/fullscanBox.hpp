#pragma once

#include <cstddef>
#include <chrono>
#include "../../utils/type.hpp"
#include "../../utils/common.hpp"
#include "base_index.hpp"

namespace bench { namespace index {


// a naive baseline FullScan
template<size_t dim>
struct FullScanBox : public BaseIndex {
    using Point = point_t<dim>;
    using Box = box_t<dim>;
    using Boxes = std::vector<Box>;

    Boxes& _data;

    FullScanBox(Boxes& boxes) : _data(boxes) {
    }

    inline size_t count() {
        return _data.size();
    }
	
	inline size_t index_size() {
        return 0;
    }

    // linear scan O(N)
    Boxes range_query(const Box& box) {
        auto start = std::chrono::steady_clock::now();
        Boxes results;
        for (auto b : _data) {
            if (bench::common::is_intersect_box(b, box)) {
                results.emplace_back(b);
            }
        }
        auto end = std::chrono::steady_clock::now();
        range_count ++;
        range_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        return results;
    }
}; 

}
}