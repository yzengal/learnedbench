#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <tuple>
#include <variant>
#include <boost/multi_array.hpp>
#include <vector>
#include <chrono>
#include <random>
#include <cassert>

#include "../base_index.hpp"
#include "../../utils/type.hpp"
#include "../../utils/common.hpp"
#include "../pgm/pgm_index.hpp"
#include "../pgm/pgm_index_variants.hpp"


namespace bench { namespace index {

// the sort dimension is always the last dimension
template<size_t Dim, size_t K, size_t Eps=64>
class Flood : public BaseIndex {

using Point = point_t<Dim>;
using Points = std::vector<Point>;
using Range = std::pair<size_t, size_t>;
using Box = box_t<Dim>;
using Index = pgm::PGMIndex<double, Eps>;

public:
static size_t SortDim;

class Bucket {
    public:
    Points _local_points;
    // eps for each bucket is fixed to 16 based on a micro benchmark
    pgm::PGMIndex<double, 16>* _local_pgm;

    Bucket() : _local_pgm(nullptr) {}

    ~Bucket() {
        delete this->_local_pgm;
    }

    inline void insert(Point& p) {
        this->_local_points.emplace_back(p);
    }

    inline void build() {
        if (_local_points.size() == 0) {
            return;
        }
        // note points are already sorted by SortDim
        std::vector<double> idx_data;
        idx_data.reserve(_local_points.size());
        for (const auto& p : _local_points) {
            idx_data.emplace_back(p[SortDim]);
        }
		
		for (size_t sz=idx_data.size()-1,i=0; i<sz; ) {
			size_t j = i++;
			while (i<sz && idx_data[i]==idx_data[j]) ++i;
			if (i - j > 1) {
				double delta = (i==sz) ? 1.0 : (idx_data[i]-idx_data[j]);
				delta /= (i - j + 1);
				for (size_t k=j; k<i; ++k) 
					idx_data[k] += delta * (k-j);
			}
		}
		idx_data.erase(std::unique(idx_data.begin(), idx_data.end()), idx_data.end());
		
		try {
			_local_pgm = new pgm::PGMIndex<double, 16>(idx_data);
		} catch (...) {
			for (size_t sz=idx_data.size()-1,i=1; i<sz; ++i) {
				if (idx_data[i] <= idx_data[i-1]) {
					std::cout << i << ": " << std::fixed << std::setprecision(20) << idx_data[i-1] << " " << idx_data[i] << std::endl;
				}
			}
			exit(1);
		}
    }

    inline void search(Points& result, Box& box) {
        if (_local_points.size() == 0) {
            return;
        }
        
		if (_local_pgm == nullptr) {
            for (int i=_local_points.size()-1; i>=0; --i) {
				if (bench::common::is_in_box(_local_points[i], box)) {
					result.emplace_back(this->_local_points[i]);
				}
			}
		} else {
			auto min_key = box.min_corner()[SortDim];
			auto max_key = box.max_corner()[SortDim];
			auto range_lo = this->_local_pgm->search(min_key);
			auto range_hi = this->_local_pgm->search(max_key);

			for (size_t i=range_lo.lo; i<range_hi.hi; ++i) {
				if (bench::common::is_in_box(this->_local_points[i], box)) {
					result.emplace_back(this->_local_points[i]);
				}
			}
		}
    }
};

Flood(Points& points, size_t _SortDim) : _data(points), bucket_size((points.size() + K - 1)/K) {
	this->SortDim = _SortDim;
    std::cout << "Construct Flood " << "K=" << K << " Epsilon=" << Eps << " SortDim=" << SortDim << std::endl;

    auto start = std::chrono::steady_clock::now();

    // dimension offsets when computing bucket ID
    for (size_t i=0,j=0; i<Dim; ++i) {
		if (i == SortDim) continue;
        this->dim_offset[i] = bench::common::ipow(K, j++);
    }

    // sort points by SortDim
    std::sort(_data.begin(), _data.end(), [](auto& p1, auto& p2) {
        return p1[SortDim] < p2[SortDim];
    });

    // boundaries of each dimension
    std::fill(mins.begin(), mins.end(), std::numeric_limits<double>::max());
    std::fill(maxs.begin(), maxs.end(), std::numeric_limits<double>::min());

    // train model on dimension 1 -- Dim-1
    std::vector<double> idx_data;
    idx_data.reserve(points.size());
    for (size_t d=0; d<Dim; ++d) {
		if (d == SortDim) continue;
		
        for (const auto& p : _data) {
            mins[d] = std::min(p[d], mins[d]);
            maxs[d] = std::max(p[d], maxs[d]);
            idx_data.emplace_back(p[d]);
        }
		
		std::sort(idx_data.begin(), idx_data.end());
		for (size_t sz=idx_data.size()-1,i=0; i<sz; ) {
			size_t j = i++;
			while (i<sz && idx_data[i]==idx_data[j]) ++i;
			if (i - j > 1) {
				double delta = (i==sz) ? 1.0 : (idx_data[i]-idx_data[j]);
				delta /= (i - j + 1);
				double initVal = idx_data[j];
				for (size_t k=j; k<i; ++k) 
					idx_data[k] = initVal + delta * (k-j);
			}
		}
		idx_data.erase(std::unique(idx_data.begin(), idx_data.end()), idx_data.end());
		maxs[d] = std::min(maxs[d], *std::max_element(idx_data.begin(), idx_data.end()));
		
		try {
			this->indexes[d] = new Index(idx_data);
		} catch (...) {
			std::cout << "build index failed" << std::endl;
			for (size_t sz=idx_data.size()-1,i=1; i<sz; ++i) {
				if (idx_data[i] <= idx_data[i-1]) {
					std::cout << i << ": " << std::fixed << std::setprecision(20) << idx_data[i-1] << " " << idx_data[i] << std::endl;
				}
			}
			exit(1);
		}
		std::cout << "Index Dim = " << d << std::endl;

        idx_data.clear();
    }

    // note data are sorted by SortDim
    for (auto& p : _data) {
		assert(compute_id(p) < buckets.size());
        buckets[compute_id(p)].insert(p);
    }

    for (auto& b : buckets) {
        b.build();
    }

    auto end = std::chrono::steady_clock::now();
    build_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Build Time: " << get_build_time() << " [ms]" << std::endl;
    std::cout << "Index Size: " << index_size() << " Bytes" << std::endl;
}

Points range_query(Box& box) {
    auto start = std::chrono::steady_clock::now();

    // find all intersected cells
    std::vector<std::pair<size_t, size_t>> ranges;
    find_intersect_ranges(ranges, box);
    
    // search each cell using local models
    Points result;
    for (auto& range : ranges) {
        for (auto idx=range.first; idx<=range.second; ++idx) {
			assert(idx < this->buckets.size());
            this->buckets[idx].search(result, box);
        }
    }

    auto end = std::chrono::steady_clock::now();
    range_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    range_count ++;

    return result;
}

inline size_t count() {
    return _data.size();
}

inline size_t index_size() {
    // size of dimension-level learned index
    size_t cdf_size = 0;
    for (size_t i=0; i<Dim; ++i) {
		if (i == SortDim) continue;
        cdf_size += this->indexes[i]->size_in_bytes();
    }

    // size of models within each bucket
    size_t b_size = 0;
    for (size_t i=0; i<buckets.size(); ++i) {
        if (buckets[i]._local_pgm != nullptr) {
            b_size += buckets[i]._local_pgm->size_in_bytes();
        }
    }

    return cdf_size + b_size + count() * sizeof(size_t);
}


~Flood() {
    for (size_t i=0; i<Dim; ++i) {
		if (i == SortDim) continue;
        delete this->indexes[i];
    }
}


private:
Points& _data;
std::array<Index*, Dim> indexes;
std::array<Bucket, bench::common::ipow(K, Dim-1)> buckets;
std::array<size_t, Dim> dim_offset;

std::array<double, Dim> mins;
std::array<double, Dim> maxs;

const size_t bucket_size;

inline void find_intersect_ranges(std::vector<std::pair<size_t, size_t>>& ranges, Box& qbox) {
	// find all intersect ranges
	for (size_t i=0, j=0; i<Dim; ++i) {
		if (i == SortDim) continue;
		if (j++ == 0) {
			ranges.emplace_back(get_dim_idx(qbox.min_corner(), i), get_dim_idx(qbox.max_corner(), i));
			continue;
		}
		auto start_idx = get_dim_idx(qbox.min_corner(), i);
		auto end_idx = get_dim_idx(qbox.max_corner(), i);

		std::vector<std::pair<size_t, size_t>> temp_ranges;
		for (auto idx=start_idx; idx<=end_idx; ++idx) {
			for (size_t j=0; j<ranges.size(); ++j) {
				temp_ranges.emplace_back(ranges[j].first + idx*dim_offset[i], ranges[j].second + idx*dim_offset[i]);
			}
		}

		// update the range vector
		ranges = temp_ranges;
	}
}

// locate the bucket on d-th dimension using binary search
inline size_t get_dim_idx(Point& p, size_t d) {
    if (p[d] <= this->mins[d]) {
        return 0;
    }
    if (p[d] >= this->maxs[d]) {
        return K-1;
    }
    auto approx_pos = this->indexes[d]->search(p[d]).pos / this->bucket_size;
    return std::min(approx_pos, K-1);
}

inline size_t compute_id(Point& p) {
    size_t id = 0;

    for (size_t i=0; i<Dim; ++i) {
		if (i == SortDim) continue;
        auto current_idx = get_dim_idx(p, i);
        id += current_idx * dim_offset[i];
    }

    return id;
}
};

template<size_t Dim, size_t K, size_t Eps>
size_t Flood<Dim, K, Eps>::SortDim;
}
}
