#pragma once

#include "../../utils/type.hpp"
#include "../../utils/common.hpp"
#include "../base_index.hpp"

#include <algorithm>
#include <array>
#include <boost/geometry/geometries/box.hpp>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <ostream>
#include <streambuf>
#include <tuple>
#include <utility>
#include <vector>
#include <iostream>
#include <chrono>
#include <unordered_set>


namespace bench { namespace index {

// uniform K by K ... grid 
template<size_t dim, size_t K>
class UG : public BaseIndex {

using Point = point_t<dim>;
using Box = box_t<dim>;
using BoxID_t = int;
using Boxes = std::vector<Box>;
using Range = std::pair<size_t, size_t>;

Boxes& _data;

public:
    UG(Boxes& boxes):_data(boxes) {
        std::cout << "Construct Uniform Grid K=" << K << std::endl;
        auto start = std::chrono::steady_clock::now();

        this->num_of_boxes = boxes.size();

        // dimension offsets when computing bucket ID
        for (size_t i=0; i<dim; ++i) {
            this->dim_offset[i] = bench::common::ipow(K, i);
        }

        // boundaries of each dimension
        std::fill(mins.begin(), mins.end(), std::numeric_limits<double>::max());
        std::fill(maxs.begin(), maxs.end(), std::numeric_limits<double>::min());

		for (auto& box : boxes) {
			Box p = box.min_corner(), q = box.max_corner();
			for (size_t i=0; i<dim; ++i) {
				mins[i] = std::min(std::min(p[i], q[i]), mins[i]);
				maxs[i] = std::max(std::max(p[i], q[i]), maxs[i]);
			}
		}

        // widths of each dimension
        for (size_t i=0; i<dim; ++i) {
            widths[i] = (maxs[i] - mins[i]) / K;
        }
        
        // insert points to buckets
        for (BoxID_t idx=boxes.size()-1; idx>=0; --idx) {
			auto gids = compute_ids(boxes[idx]);
			for (auto gid : gids)
				buckets[gid].emplace_back(ix);
        }

        auto end = std::chrono::steady_clock::now();
        build_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Build Time: " << get_build_time() << " [ms]" << std::endl;
        std::cout << "Index Size: " << index_size() << " Bytes" << std::endl;
    }


    Boxes range_query(Box& box) {
        auto start = std::chrono::steady_clock::now();

        // bucket ranges that intersect the query box
        std::vector<Range> ranges;
		unordered_set<BoxID_t> visit;
		
        // search range on the 1-st dimension
        ranges.emplace_back(std::make_pair(get_dim_idx(box.min_corner(), 0), get_dim_idx(box.max_corner(), 0)));
        
        // find all intersect ranges
        for (size_t i=1; i<dim; ++i) {
            auto start_idx = get_dim_idx(box.min_corner(), i);
            auto end_idx = get_dim_idx(box.max_corner(), i);

            std::vector<Range> temp_ranges;
            for (auto idx=start_idx; idx<=end_idx; ++idx) {
                for (size_t j=0; j<ranges.size(); ++j) {
                    temp_ranges.emplace_back(std::make_pair(ranges[j].first + idx*dim_offset[i], ranges[j].second + idx*dim_offset[i]));
                }
            }

            // update the range vector
            ranges = temp_ranges;
        }

        // Points candidates;
        Boxes results;

        // find candidate points
        for (auto range : ranges) {
            auto start_idx = range.first;
            auto end_idx = range.second;

            for (auto idx=start_idx; idx<=end_idx; ++idx) {
                for (BoxID_t candBoxID : this->buckets[idx]) {
					Box& candBox = _data[candBoxID];
                    if (bench::common::is_intersect_box(candBox, box)) {
						if (visit.count(candBoxID) == 0) 
							visit.insert(candBoxID);
							result.emplace_back(candBox);
						}
                    }
                }
            }
        }

        auto end = std::chrono::steady_clock::now();
        range_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        range_count ++;
        
        return results;
    }

    inline size_t count() {
        return this->num_of_boxes;
    }

    inline size_t index_size() {
        size_t ret = dim * (3 * sizeof(double) + sizeof(size_t)); 
		for (int i=buckets.size()-1; i>=0; --i) {
			ret += sizeof(BoxID_t)*buckets[i].size() + sizeof(size_t) + sizeof(void*);
		}
		return ret;
    }


private:
    size_t num_of_boxes;
    std::array<vector<BoxID_t>, common::ipow(K, dim)> buckets;
    std::array<double, dim> mins;
    std::array<double, dim> maxs;
    std::array<double, dim> widths;
    std::array<size_t, dim> dim_offset;

    // compute the index on d-th dimension of a given point
    inline size_t get_dim_idx(const Point& p, const size_t& d) {
        if (p[d] <= mins[d]) {
            return 0;
        } else if (p[d] >= maxs[d]) {
            return K-1;
        } else {
            return (size_t) ((p[d] - mins[d]) / widths[d]);
        }
    }

    // compute the bucket ID of a given point
    inline size_t compute_id(const Point& p) {
        size_t id = 0;

        for (size_t i=0; i<dim; ++i) {
            auto current_idx = get_dim_idx(p, i);
            id += current_idx * dim_offset[i];
        }

        return id;
    }
	
    // compute the bucket ID of a given box
    inline vector<size_t> compute_ids(const Box& box) {
		vector<size_t> ret;
		const Point mxp = box.max_corner(), mnp = box.min_corner();
		Point p;
		
		for (int st=0; st<(1<<dim); ++st) {
			for (int i=0; i<dim; ++i) {
				p[i] = (st & (1<<i)) ? mxp[i] : mnp[i];
			}	
			auto id = compute_id(p);
			bool flag = true;
			for (auto st : ret) {
				if (st == id) {
					flag = false;
					break;
				}
			}
			if (flag) ret.emplace_back(id);
		}
		
        return ret;
    }
};

}
}
