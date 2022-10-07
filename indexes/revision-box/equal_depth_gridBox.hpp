#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <unordered_set>

#include "base_index.hpp"
#include "../../utils/type.hpp"
#include "../../utils/common.hpp"


namespace bench { namespace index {

template<size_t dim, size_t K>
class EDG : public BaseIndex {

	using Point = point_t<dim>;
	using Box = box_t<dim>;
	using BoxID_t = int;
	using Boxes = std::vector<Box>;
	using Range = std::pair<size_t, size_t>;
	using Partition = std::array<double, K>;
	using Partitions = std::array<Partition, dim>;

	public:
	Boxes& _data;
	
	EDG(Boxes& boxes) : _data(boxes) {
		std::cout << "Construct Euqal-Depth Grid K=" << K << std::endl;
		auto start = std::chrono::steady_clock::now();

		// dimension offsets when computing bucket ID
		for (size_t i=0; i<dim; ++i) {
			this->dim_offset[i] = bench::common::ipow(K, i);
		}

		this->num_of_boxes = boxes.size();
		auto bucket_size = num_of_boxes / K;

		// compute equal depth partition boundaries 
		for (size_t i=0; i<dim; ++i) {
			std::vector<double> dim_vector;
			dim_vector.reserve(num_of_boxes);

			for (const auto& box : boxes) {
				double coordinate = (box.min_corner()[i]+box.max_corner()[i]) * 0.5;
				dim_vector.emplace_back(coordinate);
			}
			std::sort(dim_vector.begin(), dim_vector.end());

			for (size_t j=0; j<K; ++j) {
				partitions[i][j] = dim_vector[j * bucket_size];
			}
		}

        // insert boxes to buckets
        for (BoxID_t idx=boxes.size()-1; idx>=0; --idx) {
			auto gids = compute_ids(boxes[idx]);
			for (auto gid : gids)
				buckets[gid].emplace_back(idx);
        }

		auto end = std::chrono::steady_clock::now();
		build_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "Build Time: " << get_build_time() << " [ms]" << std::endl;
		std::cout << "Index Size: " << index_size() << " Bytes" << std::endl;
	}


	Boxes range_query(const Box& box) {
		auto start = std::chrono::steady_clock::now();
		// bucket ranges that intersect the query box
		std::vector<Range> ranges;

		// search range on the 1-st dimension
		ranges.emplace_back(get_dim_idx(box.min_corner(), 0), get_dim_idx(box.max_corner(), 0));
		
		// find all intersect ranges
		for (size_t i=1; i<dim; ++i) {
			auto start_idx = get_dim_idx(box.min_corner(), i);
			auto end_idx = get_dim_idx(box.max_corner(), i);

			std::vector<Range> temp_ranges;
			for (auto idx=start_idx; idx<=end_idx; ++idx) {
				for (size_t j=0; j<ranges.size(); ++j) {
					temp_ranges.emplace_back(ranges[j].first + idx*dim_offset[i], ranges[j].second + idx*dim_offset[i]);
				}
			}

			// update the range vector
			ranges = temp_ranges;
		}

		// Boxes candidates;
		std::unordered_set<BoxID_t> visit;
		Boxes results;

		// find candidate boxes
		for (auto range : ranges) {
			auto start_idx = range.first;
			auto end_idx = range.second;

			for (auto idx=start_idx; idx<=end_idx; ++idx) {
				for (BoxID_t candBoxID : this->buckets[idx]) {
					Box& candBox = _data[candBoxID];
					if (bench::common::is_intersect_box(candBox, box)) {
						if (visit.count(candBoxID) == 0) {
							visit.insert(candBoxID);
							results.emplace_back(candBox);
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
		size_t ret = dim * K * sizeof(double) + dim * sizeof(size_t);
		for (int i=buckets.size()-1; i>=0; --i) {
			ret += sizeof(BoxID_t)*buckets[i].size() + sizeof(size_t) + sizeof(void*);
		}
		return ret;
	}

	void print_partitions() {
		for (auto& partition : partitions) {
			for (auto p : partition) {
				std::cout << p << " ";
			}
			std::cout << std::endl;
		}
	}
		
	private:
	size_t num_of_boxes;
	std::array<std::vector<BoxID_t>, bench::common::ipow(K, dim)> buckets;
	std::array<size_t, dim> dim_offset;
	Partitions partitions; // bucket boundaries on each dimension

	// locate the bucket on d-th dimension using binary search
	inline size_t get_dim_idx(const Point& p, const size_t d) {
		if (p[d] <= *partitions[d].begin()) {
			return 0;
		} else {
			auto upper = std::upper_bound(partitions[d].begin(), partitions[d].end(), p[d]);
			return (size_t) (upper - partitions[d].begin() - 1);
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
    inline std::vector<size_t> compute_ids(const Box& box) {
		std::vector<size_t> ret;
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



