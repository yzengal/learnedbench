#pragma once

#include <cstddef>
#include <array>
#include <chrono>
#include <algorithm>
#include <unordered_set>
#include "../../utils/type.hpp"
#include "../../utils/common.hpp"
#include "base_index.hpp"
#include "../pgm/pgm_index.hpp"

namespace bench { namespace index {

template<size_t Dim, size_t K=10, size_t Epsilon=64>
class LISA : public BaseIndex {

using Point = point_t<Dim>;
using Box = box_t<Dim>;
using Points = std::vector<Point>;
using Partition = std::array<double, K>;
using Partitions = std::array<Partition, Dim>;
using PGMIdx = pgm::PGMIndex<double, Epsilon>;
	
public:

LISA(Points& points) {
    std::cout << "Construct LISA " << "K=" << K << " Epsilon=" << Epsilon << std::endl;

    auto start = std::chrono::steady_clock::now();
    // dimension offsets when computing bucket ID
    for (size_t i=0; i<Dim; ++i) {
        this->dim_offset[i] = bench::common::ipow(K, i);
    }

    auto bucket_size = points.size() / K;

    this->min_width = std::numeric_limits<double>::max();

    // compute equal depth partition boundaries 
    for (size_t i=0; i<Dim; ++i) {
        std::vector<double> dim_vector;
        dim_vector.reserve(points.size());

        for (const auto& p : points) {
            dim_vector.emplace_back(p[i]);
        }
        std::sort(dim_vector.begin(), dim_vector.end());

        for (size_t j=0; j<K; ++j) {
            if (j == 0) {
                partitions[i][j] = dim_vector[j*bucket_size] - 0.000001; // this is to resolve numerical issues
            } else {
                partitions[i][j] = (dim_vector[j*bucket_size] + dim_vector[j*bucket_size + 1]) / 2.0;
            }
        }

        mins[i] = dim_vector.front();
        maxs[i] = dim_vector.back();

        this->min_width = std::min(this->min_width, (maxs[i] - mins[i]));
    }

    // initialize volumes of grid cells
    std::array<size_t, Dim> comb;
    for (size_t i = 0; i < volumes.size(); ++i) {
        auto I = i;
        for (size_t j = 0; j < Dim; ++j) {
            comb[j] = I % K;
            I /= K;
        }

        auto bucket_id = compute_id(comb);
        double vol = 10.0;
        for (size_t d=0; d<Dim; ++d) {
            if (comb[d] == (K-1)) {
                vol *= (100 * (this->maxs[d] - this->partitions[d][comb[d]]));
            } else {
                vol *= (100 * (this->partitions[d][comb[d]+1] - this->partitions[d][comb[d]]));
            }
        }
        this->volumes[bucket_id] = vol;
    }

    // make sure all volumes are well initialized
    for (auto v : this->volumes) {
        assert(v > 0);
    }

    // compute projections and sort points by projections
    std::vector<std::pair<size_t, double>> pid_and_projection;
    pid_and_projection.reserve(points.size());
    size_t pid = 0;
    for (auto& p : points) {
        pid_and_projection.emplace_back(pid++, project(p));
    }

    std::sort(pid_and_projection.begin(), pid_and_projection.end(), 
        [](auto p1, auto p2) { return std::get<1>(p1) < std::get<1>(p2); });

    std::vector<double> projections;
    projections.reserve(points.size());
    this->_data.reserve(points.size());

    for (auto& pp : pid_and_projection) {
		Points _local_points;
		_local_points.emplace_back(points[std::get<0>(pp)]);
        this->_data.emplace_back(_local_points);
        projections.emplace_back(std::get<1>(pp));
    }

	// make sure there is no duplicate keys
    for (size_t sz=projections.size()-1,i=0; i<sz; ) {
		size_t j = i++;
		while (i<sz && projections[i]==projections[j]) ++i;
		if (i - j > 1) {
			double delta = (i==sz) ? 1.0 : (projections[i]-projections[j]);
			delta /= (i - j);
			for (size_t k=j; k<i; ++k) 
				projections[k] += delta * (k-j);
		}
    }
    
	// train 1-D learned index on projections
	this->_pgm_ptr = new PGMIdx(projections);

    auto end = std::chrono::steady_clock::now();
    build_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Build Time: " << get_build_time() << " [ms]" << std::endl;
    std::cout << "Index Size: " << index_size() << " Bytes" << std::endl;
}

~LISA() {
    delete this->_pgm_ptr;
}

void __insert(const Point& point) {
	auto pid = compute_id(point);
	
	auto range = this->_pgm_ptr->search(pid);
	size_t mid = (range.lo + range.hi) / 2;
	
	if (mid >= this->_data.size()) mid = range.lo;
	this->_data[mid].emplace_back(point);
}

bool __erase(const Point& point) {
	auto pid = compute_id(point);
	auto range = this->_pgm_ptr->search(pid);

    for (size_t i=range.lo; i<=range.hi&&i<this->_data.size(); ++i) {
		if (this->_data[i].empty()) continue;
		for (size_t j=0; j<this->_data[i].size(); ++j) {		
			if (bench::common::is_equal_point(this->_data[i][j], point)) {
				this->_data[i][j] = *this->_data[i].rbegin();
				this->_data[i].pop_back();
				return true;    
			}	
		}
    }
	
	return false;
}

void insert(Point& point) {
	auto start = std::chrono::steady_clock::now();
	
	__insert(point);
	
	auto end = std::chrono::steady_clock::now();
    insert_count ++;
    insert_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

bool erase(Point& point) {
	auto start = std::chrono::steady_clock::now();
	
	bool ret = __erase(point);
	
	auto end = std::chrono::steady_clock::now();
    erase_count ++;
    erase_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	
	return ret;
}

inline size_t count() {
	size_t ret = 0;
	for (size_t i=0; i<this->_data.size(); ++i)
		ret += this->_data[i].size();
    return ret;
}

inline size_t index_size() {
    size_t partition_size = K * Dim * sizeof(double);
    size_t volume_size = this->volumes.size() * sizeof(double);
    size_t rmi_size = this->_pgm_ptr->size_in_bytes();
    return partition_size + volume_size + rmi_size + count() * sizeof(size_t);
}

Points range_query(Box& box) {
    auto start = std::chrono::steady_clock::now();
    std::vector<std::pair<size_t, size_t>> ranges;
    find_intersect_ranges(ranges, box);

    Points result;
	std::unordered_set<unsigned long long> visit;
    for (const auto& range : ranges) {
        search_range(result, static_cast<double>(std::get<0>(range)), static_cast<double>(std::get<1>(range)+1), box, visit);
    }

    auto end = std::chrono::steady_clock::now();
    range_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    range_count ++;

    return result;
}

Points knn_query(Point& point, size_t k) {
    auto start = std::chrono::steady_clock::now();

    Points range_result;

    // initialize search range
    double r = initial_knn_range(point, k);
    if (r == 0.0) {
        r = this->min_width / (2.0 * K);
    }

    while (1) {
        knn_range_helper(range_result, point, r);

        // k results found
        if (range_result.size() >= k) {
            break;
        }

        // increase the search range
        r *= 2;
        range_result.clear();
    }

    auto end = std::chrono::steady_clock::now();
    knn_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    knn_count ++;

    std::sort(range_result.begin(), range_result.end(), 
        [&](auto& p1, auto p2) { 
            return bench::common::eu_dist_square(p1, point) < bench::common::eu_dist_square(p2, point); 
        });

    Points knn_result(range_result.begin(), range_result.begin()+k);

    return knn_result;
}


private:
// dimension offsets in bucket array
std::array<size_t, Dim> dim_offset;

// bucket boundaries on each dimension
Partitions partitions; 

// min corner
std::array<double, Dim> mins;

// max corner
std::array<double, Dim> maxs;

// minimum width
double min_width;

// pre-computed grid volumes
std::array<double, bench::common::ipow(K, Dim)> volumes;

// points ordered by projection function
std::vector<Points>  _data;

// ptr to the underlying 1-d learned index
PGMIdx* _pgm_ptr; 

// find intitial search range
inline double initial_knn_range(Point& q, size_t k) {
    double min_r = std::numeric_limits<double>::max();

    for (size_t i=0; i<Dim; ++i) {
        auto idx = get_dim_idx(q, i);
        if (idx == K-1) {
            min_r = std::min(min_r, (this->maxs[i] - q[i]));
            min_r = std::min(min_r, (q[i] - this->partitions[i][idx]));
        } else {
            min_r = std::min(min_r, (this->partitions[i][idx+1] - q[i]));
            min_r = std::min(min_r, (q[i] - this->partitions[i][idx]));
        }
    }

    return min_r;
}

void knn_range_helper(Points& result_found, Point& q, double r) {
    Point min_p = q;
    Point max_p = q;
    for (size_t i=0; i<Dim; ++i) {
        min_p[i] -= r;
        max_p[i] += r;
    }
    Box qbox(min_p, max_p);

    std::vector<std::pair<size_t, size_t>> ranges;
	std::unordered_set<unsigned long long> visit;
    find_intersect_ranges(ranges, qbox);

    for (const auto& range : ranges) {
        search_range(result_found, static_cast<double>(std::get<0>(range)), static_cast<double>(std::get<1>(range)+1), qbox, visit);
    }
}

// lo and hi are projected values
inline void search_range(Points& result, double lo, double hi, Box& qbox, std::unordered_set<unsigned long long>& visit) {
    auto range_lo = this->_pgm_ptr->search(lo);
    auto range_hi = this->_pgm_ptr->search(hi);

    for (size_t i=range_lo.lo; i<=range_hi.hi&&i<this->_data.size(); ++i) {
		for (size_t j=0; j<this->_data[i].size(); ++j) {		
			if (bench::common::is_in_box(this->_data[i][j], qbox)) {
				unsigned long long keyID = 1ULL * j * this->_data.size() + i;
				if (visit.count(keyID) == 0) {
					result.emplace_back(this->_data[i][j]);
					visit.insert(keyID);
				}    
			}	
		}
    }

    // auto it_lo = this->_data.begin() + range_lo.lo;
    // auto it_hi = this->_data.begin() + range_hi.hi;

    // for (auto it=it_lo; it!=it_hi; ++it) {
    //     if (bench::common::is_in_box(*it, qbox)) {
    //         result.emplace_back(*it);
    //     }
    // }
}

inline void find_intersect_ranges(std::vector<std::pair<size_t, size_t>>& ranges, Box& qbox) {
    // search range on the 1-st dimension
    ranges.emplace_back(get_dim_idx(qbox.min_corner(), 0), get_dim_idx(qbox.max_corner(), 0));
    
    // find all intersect ranges
    for (size_t i=1; i<Dim; ++i) {
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

inline double project(const Point& point) {
    size_t id = compute_id(point);
    double vol = 10.0;
    
    for (size_t i=0; i<Dim; ++i) {
        auto current_idx = get_dim_idx(point, i);
        
        if (point[i] <= this->mins[i]) {
            vol = 0.0;
        } else if (point[i] >= this->maxs[i]) {
            vol *= (100 * (this->maxs[i] - partitions[i][current_idx]));
        } else {
            vol *= (100 * (point[i] - partitions[i][current_idx]));
        }
    }

    return static_cast<double>(id) + vol / this->volumes[id];
}

inline size_t get_dim_idx(const Point& p, size_t d) {
    if (p[d] <= partitions[d][0]) {
        return 0;
    } else {
        auto upper = std::upper_bound(partitions[d].begin(), partitions[d].end(), p[d]);
        return (size_t) (upper - partitions[d].begin() - 1);
    }
}

inline size_t compute_id(const Point& p) {
    size_t id = 0;

    for (size_t i=0; i<Dim; ++i) {
        auto current_idx = get_dim_idx(p, i);
        id += current_idx * dim_offset[i];
    }

    return id;
}

inline size_t compute_id(std::array<size_t, Dim>& idx_comb) {
    size_t id = 0;
    for (size_t i=0; i<Dim; ++i) {
        id += idx_comb[i] * dim_offset[i];
    }
    return id;
}

};

}
}
