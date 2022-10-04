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
#include <unordered_set>

#include "../base_index.hpp"
#include "../../utils/type.hpp"
#include "../../utils/common.hpp"
#include "../rmi/models.hpp"
#include "../rmi/rmi.hpp"


namespace bench { namespace index {

// the sort dimension is always the last dimension
template<size_t Dim, size_t MaxElements=128, size_t Eps=64, size_t SortDim=Dim-1>
class Flood : public BaseIndex {

using Point = point_t<Dim>;
using Points = std::vector<Point>;
using Range = std::pair<size_t, size_t>;
using Box = box_t<Dim>;
using Counts = std::vector<int>;
using layer1_type = rmi::LinearSpline;
using layer2_type = rmi::LinearRegression;
using Index = rmi::RmiLAbs<double, layer1_type, layer2_type>;
	
public:

class OctreeNode {
    public:
	static const int CHILD_SIZE = ((size_t)1)<<dim;
	// Physical position/size. This implicitly defines the bounding 
	// box of this node
	Point origin;         //! The physical center of this node
	Point halfDimension;  //! Half the width/height/depth of this node

	// The tree has up to eight children and can additionally store
	// a point, though in many applications only, the leaves will store data.
	OctreeNode *children[CHILD_SIZE]; //! Pointers to child octants
	
    Points _local_points;
	Counts cnts; // number of same data points
    // eps for each bucket is fixed to 16 based on a micro benchmark
    Index* _local_pgm;

	OctreeNode(const Points& points) {
		int n = points.size();
		int oid = rand() % n;
		origin = points[oid];
		
		halfDimension.fill(0.0);
		for (int i=0; i<n; ++i) {
			for (int j=0; j<dim; ++j) {
				auto delta = fabs(origin[j] - points[i][j]);
				if (delta > halfDimension[j])
					halfDimension[j] = delta;
			}
		}
		
		// Initially, there are no children
		_local_points.clear();
		cnts.clear();
		for(int i=0; i<CHILD_SIZE; ++i) 
			children[i] = NULL;
		
		for (auto point : points) {
			insert(point, 1);
		}
		_local_pgm = NULL;
	}
	
	OctreeNode(const Point& _origin, const Point& _halfDimension) 
		: origin(_origin), halfDimension(_halfDimension) {
			// Initially, there are no children
			_local_points.clear();
			cnts.clear();
			for(int i=0; i<CHILD_SIZE; ++i) 
				children[i] = NULL;
			if (this->_local_pgm != NULL)
				delete this->_local_pgm;
			this->_local_pgm = NULL;
		}


    ~OctreeNode() {
		// Recursively destroy octants
		for(int i=0; i<CHILD_SIZE; ++i) {
			if (children[i] != NULL)
				delete children[i];
		}
		if (this->_local_pgm != NULL)
			delete this->_local_pgm;
    }
	
	inline size_t index_size() {
		size_t ret = sizeof(OctreeNode*)*CHILD_SIZE + sizeof(double)*dim*2;
		
		ret += sizeof(Point*) + sizeof(size_t) + sizeof(double)*dim*_local_points.size()
		ret += sizeof(int*) + sizeof(size_t) + sizeof(int)*cnts.size();
		for(int i=0; i<CHILD_SIZE; ++i) {
			if (children[i] != NULL)
				ret += children[i]->index_size();
		}
		
		return ret;
	}

	// Determine which octant of the tree would contain 'point'
	int getOctantContainingPoint(const Point& point) const {
		int oct = 0;
		for (int i=0; i<dim; ++i) {
			if (point[i] >= origin[i])
				oct |= (1 << i);
		}
		return oct;
	}
	
	bool isLeafNode() const {
		return children[0] == NULL;
	}
	
	void insert(const Point& point, int num=1) {
		
		// If this node doesn't have a data point yet assigned 
		// and it is a leaf, then we're done!
		if(isLeafNode()) {
			
			// check if the new point equals to any current point
			for (int i=0; i<_local_points.size(); ++i) {
				const Point& oldPoint = _local_points[i];
				if (bench::common::is_equal_point<dim>(oldPoint, point)) {
					cnts[i] += num;
					return ;
				}		
			}
			if(ids.size() < MaxElements) {
				cnts.emplace_back(num);
				_local_points.emplace_back(point);
			} else {
				// We're at a leaf, but there's already something here
				// We will split this node so that it has 8 child octants
				// and then insert the old data that was here, along with 
				// this new data point

				// Save this data point that was here for a later re-insert

				// Split the current node and create new empty trees for each
				// child octant.
				Point newHalfDimension = halfDimension;
				for (int i=0; i<dim; ++i)
					newHalfDimension[i] *= .5f;
				for(int j=0; j<CHILD_SIZE; ++j) {
					// Compute new bounding box for this child
					Point newOrigin = origin;
					for (int i=0; i<dim; ++i) {
						newOrigin[i] += halfDimension[i] * ((j&(1<<i)) ? .5f : -.5f);
					}						
					children[j] = new OctreeNode(newOrigin, newHalfDimension);
				}

				// Re-insert the old point, and insert this new point
				// (We wouldn't need to insert from the root, because we already
				// know it's guaranteed to be in this section of the tree)
				for (int i=0; i<_local_points.size(); ++i) {
					const Point& oldPoint = _local_points[i];
					children[getOctantContainingPoint(oldPoint)]->insert(oldPoint, cnts[i]);
				}
				_local_points.clear();
				cnts.clear();
				children[getOctantContainingPoint(point)]->insert(point, num);
			}
		} else {
			// We are at an interior node. Insert recursively into the 
			// appropriate child octant
			int octant = getOctantContainingPoint(point);
			children[octant]->insert(point, num);
		}
	}

    void DFS_build() {
		if(isLeafNode()) {
			if (_local_points.empty())
				return;
			// note points are already sorted by SortDim
			std::vector<double> idx_data;
			idx_data.reserve(_local_points.size());
			for (const auto& p : _local_points) {
				idx_data.emplace_back(std::get<SortDim>(p));
			}

			// _local_pgm = new pgm::PGMIndex<double, 16>(idx_data);
			// train 1-D learned index on projections
			std::size_t index_budget = idx_data.size() * sizeof(double);
			std::size_t layer2_size = (index_budget - 2 * sizeof(double) - 2 * sizeof(std::size_t)) / (2 * sizeof(double));
			if (layer2_size < 8) layer2_size = 8;
			this->_local_pgm = new Index(idx_data, layer2_size);
			
		} else {
			for(int j=0; j<CHILD_SIZE; ++j) {
				if (children[j] != NULL)
					cildren[j]->DFS_build();
			}
		}
    }

    inline void search(Points& result, Box& box, size_t bucketID, std::unordered_set<double>& visit) {
        if (_local_pgm == nullptr) {
            return;
        }
        
        auto min_key = std::get<SortDim>(box.min_corner());
        auto max_key = std::get<SortDim>(box.max_corner());
        auto range_lo = this->_local_pgm->search(min_key);
        auto range_hi = this->_local_pgm->search(max_key);

        for (size_t i=range_lo.lo; i<=range_hi.hi&&i<this->_local_points.size(); ++i) {
            if (bench::common::is_in_box(this->_local_points[i], box)) {
				// double keyID = i*1.0*bench::common::ipow(K, Dim-1) + bucketID;
				// if (visit.count(keyID) == 0) {
					// visit.insert(keyID);
					result.emplace_back(this->_local_points[i]);
				// }
            }
        }
    }
};

Flood(Points& points) : _data(points), bucket_size((points.size() + K - 1)/K) {
    std::cout << "Construct Flood " << "K=" << K << " Epsilon=" << Eps << " SortDim=" << SortDim << std::endl;

    auto start = std::chrono::steady_clock::now();

    // dimension offsets when computing bucket ID
    for (size_t i=0; i<Dim-1; ++i) {
        this->dim_offset[i] = bench::common::ipow(K, i);
    }

    // sort points by SortDim
    std::sort(_data.begin(), _data.end(), [](auto& p1, auto& p2) {
        return std::get<SortDim>(p1) < std::get<SortDim>(p2);
    });

    // boundaries of each dimension
    std::fill(mins.begin(), mins.end(), std::numeric_limits<double>::max());
    std::fill(maxs.begin(), maxs.end(), std::numeric_limits<double>::min());

    // train model on dimension 1 -- Dim-1
    std::vector<double> idx_data;
    idx_data.reserve(points.size());
    for (size_t i=0; i<Dim-1; ++i) {
        for (const auto& p : _data) {
            mins[i] = std::min(p[i], mins[i]);
            maxs[i] = std::max(p[i], maxs[i]);

            idx_data.emplace_back(p[i]);
        }

        std::sort(idx_data.begin(), idx_data.end());
		std::size_t index_budget = idx_data.size() * sizeof(double);
		std::size_t layer2_size = (index_budget - 2 * sizeof(double) - 2 * sizeof(std::size_t)) / (2 * sizeof(double));
		if (layer2_size < 8) layer2_size = 8;
        this->indexes[i] = new Index(idx_data, layer2_size);

        idx_data.clear();
    }


    // note data are sorted by SortDim
    for (auto& p : _data) {
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
	std::unordered_set<double> visit;
    for (auto& range : ranges) {
        for (auto idx=range.first; idx<=range.second; ++idx) {
            this->buckets[idx].search(result, box, idx, visit);
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
    for (size_t i=0; i<Dim-1; ++i) {
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
    for (size_t i=0; i<Dim-1; ++i) {
        delete this->indexes[i];
    }
}


private:
Points& _data;
std::array<Index*, Dim-1> indexes;
std::array<Bucket, bench::common::ipow(K, Dim-1)> buckets;
std::array<size_t, Dim-1> dim_offset;

std::array<double, Dim-1> mins;
std::array<double, Dim-1> maxs;

const size_t bucket_size;

inline void find_intersect_ranges(std::vector<std::pair<size_t, size_t>>& ranges, Box& qbox) {
    if (Dim == 2) {
        ranges.emplace_back(get_dim_idx(qbox.min_corner(), 0), get_dim_idx(qbox.max_corner(), 0));
    } else {
        // search range on the 1-st dimension
        ranges.emplace_back(get_dim_idx(qbox.min_corner(), 0), get_dim_idx(qbox.max_corner(), 0));
        
        // find all intersect ranges
        for (size_t i=1; i<Dim-1; ++i) {
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

    for (size_t i=0; i<Dim-1; ++i) {
        auto current_idx = get_dim_idx(p, i);
        id += current_idx * dim_offset[i];
    }

    return id;
}


};
}
}
