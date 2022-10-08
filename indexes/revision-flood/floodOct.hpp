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
#include <cmath>
#include <unordered_set>

#include "../base_index.hpp"
#include "../../utils/type.hpp"
#include "../../utils/common.hpp"
#include "../pgm/pgm_index.hpp"
#include "../pgm/pgm_index_variants.hpp"


namespace bench { namespace index {

// the sort dimension is always the last dimension
template<size_t Dim, size_t MaxElements=128, size_t Eps=64, size_t SortDim=Dim-1>
class FloodOct : public BaseIndex {

using Point = point_t<Dim>;
using Points = std::vector<Point>;
using Range = std::pair<size_t, size_t>;
using Box = box_t<Dim>;
using Counts = std::vector<int>;
using Index =  pgm::PGMIndex<double, Eps>;
	
public:
static int MAX_DEPTH;

class OctreeNode {
    public:
	static const int LOCAL_DIM = Dim-1;
	using OctreeNodePoint = point_t<LOCAL_DIM>;
	using OctreeNodeBox = box_t<LOCAL_DIM>;
	using PairData_t = std::pair<point_t<Dim>, size_t>;
	static const int CHILD_SIZE = ((size_t)1)<<LOCAL_DIM;
	// Physical position/size. This implicitly defines the bounding 
	// box of this node
	OctreeNodePoint origin;         //! The physical center of this node
	OctreeNodePoint halfDimension;  //! Half the width/height/depth of this node

	// The tree has up to eight children and can additionally store
	// a point, though in many applications only, the leaves will store data.
	OctreeNode *children[CHILD_SIZE]; //! Pointers to child octants
	
    std::vector<PairData_t> _local_data;
    // eps for each bucket is fixed to 16 based on a micro benchmark
    Index* _local_pgm;

	OctreeNode(const Points& points) {
		size_t n = points.size();
		size_t oid = rand() % n;
		
		for (int j=0; j<LOCAL_DIM; ++j)
			origin[j] = points[oid][j];
		
		halfDimension.fill(0.0);
		for (int i=0; i<n; ++i) {
			for (int j=0; j<LOCAL_DIM; ++j) {
				auto delta = fabs(origin[j] - points[i][j]);
				if (delta > halfDimension[j])
					halfDimension[j] = delta;
			}
		}
		
		// Initially, there are no children
		_local_pgm = NULL;
		_local_data.clear();
		for(int i=0; i<CHILD_SIZE; ++i) 
			children[i] = NULL;
		
		for (auto point : points) {
			insert(point, 1);
		}
	}
	
	OctreeNode(const OctreeNodePoint& _origin, const OctreeNodePoint& _halfDimension) 
		: origin(_origin), halfDimension(_halfDimension) {
		// Initially, there are no children
		_local_data.clear();
		for(int i=0; i<CHILD_SIZE; ++i) 
			children[i] = NULL;
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
		size_t ret = sizeof(OctreeNode*)*CHILD_SIZE + sizeof(double)*LOCAL_DIM*2;
		
		// ret += sizeof(Point*) + sizeof(size_t) + sizeof(double)*Dim*_local_data.size();
		ret += sizeof(int*) + sizeof(size_t) + sizeof(int)*_local_data.size();
		for(int i=0; i<CHILD_SIZE; ++i) {
			if (children[i] != NULL)
				ret += children[i]->index_size();
		}
		if (this->_local_pgm != NULL) {
			ret += this->_local_pgm->size_in_bytes();
		}
		
		return ret;
	}

	// Determine which octant of the tree would contain 'point'
	int getOctantContainingPoint(const Point& point) const {
		int oct = 0;
		for (int i=0; i<LOCAL_DIM; ++i) {
			if (point[i] >= origin[i])
				oct |= (1 << i);
		}
		return oct;
	}
	
	bool isLeafNode() const {
		return children[0] == NULL;
	}
	
	void insert(const Point& point, int num=1, int dep=0) {
		
		// If this node doesn't have a data point yet assigned 
		// and it is a leaf, then we're done!
		if(isLeafNode()) {
			
			// check if the new point equals to any current point
			for (int i=0; i<_local_data.size(); ++i) {
				const Point& oldPoint = _local_data[i].first;
				if (bench::common::is_equal_point<Dim>(oldPoint, point)) {
					_local_data[i].second += num;
					return ;
				}		
			}
			if(dep>=MAX_DEPTH || _local_data.size()<MaxElements) {
				_local_data.emplace_back(make_pair(point, num));
			} else {
				// We're at a leaf, but there's already something here
				// We will split this node so that it has 8 child octants
				// and then insert the old data that was here, along with 
				// this new data point

				// Save this data point that was here for a later re-insert

				// Split the current node and create new empty trees for each
				// child octant.
				OctreeNodePoint newHalfDimension = halfDimension;
				for (int i=0; i<LOCAL_DIM; ++i)
					newHalfDimension[i] *= .5f;
				for(int j=0; j<CHILD_SIZE; ++j) {
					// Compute new bounding box for this child
					OctreeNodePoint newOrigin = origin;
					for (int i=0; i<LOCAL_DIM; ++i) {
						newOrigin[i] += halfDimension[i] * ((j&(1<<i)) ? .5f : -.5f);
					}						
					children[j] = new OctreeNode(newOrigin, newHalfDimension);
				}

				// Re-insert the old point, and insert this new point
				// (We wouldn't need to insert from the root, because we already
				// know it's guaranteed to be in this section of the tree)
				for (int i=0; i<_local_data.size(); ++i) {
					const Point& oldPoint = _local_data[i].first;
					const int oldNum = _local_data[i].second;
					children[getOctantContainingPoint(oldPoint)]->insert(oldPoint, oldNum, dep+1);
				}
				_local_data.clear();
				children[getOctantContainingPoint(point)]->insert(point, num, dep+1);
			}
		} else {
			// We are at an interior node. Insert recursively into the 
			// appropriate child octant
			int octant = getOctantContainingPoint(point);
			children[octant]->insert(point, num, dep+1);
		}
	}

    void DFS_build() {
		if(false == isLeafNode()) {
			for(int j=0; j<CHILD_SIZE; ++j) {
				if (children[j] != NULL)
					children[j]->DFS_build();
			}
			return ;
		}
			
		if (_local_data.empty())
			return;
		
		std::sort(_local_data.begin(), _local_data.end(), [](auto& p1, auto& p2) {
			return std::get<SortDim>(p1.first) < std::get<SortDim>(p2.first);
		});
		
		// note points are already sorted by SortDim
		std::vector<double> idx_data;
		idx_data.reserve(_local_data.size());
		for (const auto& p : _local_data) {
			idx_data.emplace_back(std::get<SortDim>(p.first));
		}
		std::sort(idx_data.begin(), idx_data.end());
		
		try {
			_local_pgm = new Index(idx_data);		
		} catch (...) {
			_local_pgm = NULL;
		}
    }

    void range_query(const Box& box, const OctreeNodeBox& boxProject, Points& results) {
		if(false == isLeafNode()) {
			for(int j=0; j<CHILD_SIZE; ++j) {
				if (children[j] == NULL) continue;
				
				// Compute the min/max corners of this child octant
				OctreeNodePoint cmax = children[j]->origin, cmin = children[j]->origin;
				for (int i=0; i<LOCAL_DIM; ++i) {
					cmax[i] += children[j]->halfDimension[i];
					cmin[i] -= children[j]->halfDimension[i];
				}
				
				// If the query rectangle is outside the child's bounding box, 
				// then continue
				OctreeNodeBox cbox(cmin, cmax);
				if (!bench::common::is_intersect_box<LOCAL_DIM>(cbox, boxProject)) continue;
				
				children[j]->range_query(box, boxProject, results);
			} 
			return ;
		}
		
        if (_local_data.empty()) return ;

		if (_local_pgm == nullptr) {
            for (int i=_local_data.size()-1; i>=0; --i) {
				if (bench::common::is_in_box<Dim>(_local_data[i].first, box)) {
					results.insert(results.end(), _local_data[i].second, _local_data[i].first);
				}
			}
        } else {
			auto min_key = std::get<SortDim>(box.min_corner());
			auto max_key = std::get<SortDim>(box.max_corner());
			auto range_lo = this->_local_pgm->search(min_key);
			auto range_hi = this->_local_pgm->search(max_key);

			for (size_t i=range_lo.lo; i<=range_hi.hi&&i<this->_local_data.size(); ++i) {
				if (bench::common::is_in_box<Dim>(_local_data[i].first, box)) {
					results.insert(results.end(), _local_data[i].second, _local_data[i].first);
				}
			}
		}
    }
};

FloodOct(Points& points) {
    std::cout << "Construct Flood-Octree " << "MaxElements=" << MaxElements << " Epsilon=" << Eps << " SortDim=" << SortDim << std::endl;
    auto start = std::chrono::steady_clock::now();
	
	this->num_of_points = points.size();
	MAX_DEPTH = std::ceil(std::log2(this->num_of_points*1.0));
    root = new OctreeNode(points);
	root->DFS_build();

    auto end = std::chrono::steady_clock::now();
    build_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Build Time: " << get_build_time() << " [ms]" << std::endl;
    std::cout << "Index Size: " << index_size() << " Bytes" << std::endl;
}

~FloodOct() {
	if (root != nullptr)
		delete root;
}

Points range_query(const Box& box) {
    auto start = std::chrono::steady_clock::now();

	Points results;
	point_t<Dim-1> mn_corner, mx_corner; 
	for (int i=0; i<Dim-1; ++i) {
		mn_corner[i] = box.min_corner()[i];
		mx_corner[i] = box.max_corner()[i];
	}
	box_t<Dim-1> boxProject(mn_corner, mx_corner);
	if (root != NULL) {
		root->range_query(box, boxProject, results);
	}

    auto end = std::chrono::steady_clock::now();
    range_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    range_count ++;

    return results;
}

inline size_t count() {
    return this->num_of_points;
}

inline size_t index_size() {
	size_t ret = sizeof(size_t) + sizeof(OctreeNode*);
	if (root != NULL)
		ret += root->index_size();
	return ret;
}

private:
OctreeNode *root;
std::size_t num_of_points;

};

template<size_t Dim, size_t MaxElements, size_t Eps, size_t SortDim>
int FloodOct<Dim, MaxElements, Eps, SortDim>::MAX_DEPTH = 20;
}
}
