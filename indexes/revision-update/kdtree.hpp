#pragma once

#include "nanoflann.hpp"
#include "../../utils/type.hpp"
#include "../../utils/common.hpp"
#include "base_index.hpp"

#include <cstddef>
#include <vector>
#include <chrono>
#include <algorithm>

#ifdef HEAP_PROFILE
#include <gperftools/heap-profiler.h>
#endif

// #define LOCAL_DEBUG

namespace bench { namespace index {

// kdtree adapter using nanoflann
template <size_t Dim, size_t MaxSplit=32>
class KDTree : public BaseIndex {

using Point = point_t<Dim>;
using Box = box_t<Dim>;
using Points = std::vector<point_t<Dim>>;


// ===== This example shows how to use nanoflann with these types of containers:
// using my_vector_of_vectors_t = std::vector<std::vector<double> > ;
//
// The next one requires #include <Eigen/Dense>
// using my_vector_of_vectors_t = std::vector<Eigen::VectorXd> ;
// =============================================================================

/** A simple vector-of-vectors adaptor for nanoflann, without duplicating the
 * storage. The i'th vector represents a point in the state space.
 *
 *  \tparam DIM If set to >0, it specifies a compile-time fixed dimensionality
 *      for the points in the data set, allowing more compiler optimizations.
 *  \tparam num_t The type of the point coordinates (typ. double or float).
 *  \tparam Distance The distance metric to use: nanoflann::metric_L1,
 *          nanoflann::metric_L2, nanoflann::metric_L2_Simple, etc.
 *  \tparam IndexType The type for indices in the KD-tree index
 *         (typically, size_t of int)
 */
template <
    class VectorOfVectorsType, typename num_t = double, int DIM = -1,
    class Distance = nanoflann::metric_L2, typename IndexType = size_t>
struct KDTreeVectorOfVectorsAdaptor
{
    using self_t =
        KDTreeVectorOfVectorsAdaptor<VectorOfVectorsType, num_t, DIM, Distance>;
    using metric_t =
        typename Distance::template traits<num_t, self_t>::distance_t;
    using index_t =
        nanoflann::KDTreeSingleIndexDynamicAdaptor<metric_t, self_t, DIM, IndexType>;

    /** The kd-tree index for the user to call its methods as usual with any
     * other FLANN index */
    index_t* index = nullptr;

    /// Constructor: takes a const ref to the vector of vectors object with the
    /// data points
    KDTreeVectorOfVectorsAdaptor(
        const size_t /* dimensionality */, const VectorOfVectorsType& mat,
        const int leaf_max_size = 10)
        : m_data(mat)
    {
        assert(mat.size() != 0 && mat[0].size() != 0);
        const size_t dims = mat[0].size();
        if (DIM > 0 && static_cast<int>(dims) != DIM)
            throw std::runtime_error(
                "Data set dimensionality does not match the 'DIM' template "
                "argument");
        index = new index_t(
            static_cast<int>(dims), *this /* adaptor */,
            nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size));
    }
	
	void buildIndex(IndexType start, IndexType end) {
		index->addPoints(start, end);
	}

    ~KDTreeVectorOfVectorsAdaptor() { delete index; }

    const VectorOfVectorsType& m_data;

    /** Query for the \a num_closest closest points to a given point
     *  (entered as query_point[0:dim-1]).
     *  Note that this is a short-cut method for index->findNeighbors().
     *  The user can also call index->... methods as desired.
     *
     * \note nChecks_IGNORED is ignored but kept for compatibility with
     * the original FLANN interface.
     */
    inline size_t query(
        const num_t* query_point, const size_t num_closest,
        IndexType* out_indices, num_t* out_distances_sq) const
    {
        // nanoflann::KNNResultSet<num_t, IndexType> resultSet(num_closest);
        // resultSet.init(out_indices, out_distances_sq);
        // index->findNeighbors(resultSet, query_point, nanoflann::SearchParams());
		
		nanoflann::KNNResultSet<num_t> resultSet(num_closest);
        resultSet.init(out_indices, out_distances_sq);
        index->findNeighbors(resultSet, query_point, nanoflann::SearchParams());
		
		return resultSet.size();
    }
	
	inline void rangeQuery(
        const num_t* query_point, const num_t radius,
        std::vector<IndexType>& out_indices, std::vector<num_t>& out_distances_sq) const
    {
		std::vector<std::pair<IndexType,num_t> > IndicesDists;
        nanoflann::RadiusResultSet<num_t, IndexType> resultSet(radius, IndicesDists);
		nanoflann::SearchParams params;
        params.sorted = false;
		
		const size_t nFound = index->findNeighbors(resultSet, query_point, params);
		
		out_indices.clear();
		out_distances_sq.clear();
		for (int i=0; i<IndicesDists.size(); ++i) {
			out_indices.push_back(IndicesDists[i].first);
			out_distances_sq.push_back(IndicesDists[i].second);
		}
		
		// std::vector<std::pair<IndexType, num_t> > ret_matches;
        

        // const size_t nFound = index->radiusSearch(
            // &query_point[0], radius, ret_matches, params);
		
		// out_indices.clear();
		// out_distances_sq.clear();
		// for (int i=0; i<nFound; ++i) {
			// out_indices.push_back(ret_matches[i].first);
			// out_distances_sq.push_back(ret_matches[i].second);
		// }
    }

    /** @name Interface expected by KDTreeSingleIndexAdaptor
     * @{ */

    const self_t& derived() const { return *this; }
    self_t&       derived() { return *this; }

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return m_data.size(); }

    // Returns the dim'th component of the idx'th point in the class:
    inline num_t kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        return m_data[idx][dim];
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    // Return true if the BBOX was already computed by the class and returned
    // in "bb" so it can be avoided to redo it again. Look at bb.size() to
    // find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const
    {
        return false;
    }
	
	void addPoint(IndexType id)
	{
		index->addPoints(id, id);
	}
	
	bool removePoint(IndexType id)
	{
		return index->removePoint(id);
	}
 
    /** @} */

};  // end of KDTreeVectorOfVectorsAdaptor


// customized kdtree type
using kdtree_t = KDTreeVectorOfVectorsAdaptor<Points, double, Dim>;

public:
KDTree(Points& points) : _data(points) {
    std::cout << "Construct kd-tree MaxSplit=" << MaxSplit << std::endl;

    auto start = std::chrono::steady_clock::now();

#ifdef HEAP_PROFILE
    HeapProfilerStart("kdtree");
#endif

    kdtree = new kdtree_t(Dim, points, MaxSplit);
    kdtree->buildIndex(0, points.size()-1);

#ifdef HEAP_PROFILE
    HeapProfilerDump("final");
    HeapProfilerStop();
#endif

    auto end = std::chrono::steady_clock::now();
    build_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Build Time: " << get_build_time() << " [ms]" << std::endl;
}

~KDTree() {
    delete kdtree;
}

Points knn_query(Point& q, unsigned int k) {
    const size_t num_of_results = k;
    size_t ret_indexes[k];
    double out_dist_sqr[k];

    auto start = std::chrono::steady_clock::now();
	
    size_t sz = kdtree->query(&q[0], num_of_results, ret_indexes, out_dist_sqr);
	// std::cout << "[FINISH] nano search" << std::endl;
    
	auto end = std::chrono::steady_clock::now();
    knn_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    knn_count ++;

    // final result
    Points result;
    result.reserve(num_of_results);
	for (int i=0; i<sz; ++i) {
		size_t idx = ret_indexes[i];
		assert(idx < _data[idx].size());
        result.emplace_back(_data[idx]);
    }

    return result;
}

Points range_query(Box& box) {
	std::vector<size_t> ret_indexes;
    std::vector<double> out_dist_sqr;
	const double EPS = 1e-5;
	
	auto start = std::chrono::steady_clock::now();
	Point q, mnp = box.min_corner(), mxp = box.max_corner();
	double radius;
	
	for (size_t d=0; d<Dim; ++d) {
		q[d] = (mnp[d]+mxp[d]) * 0.5;
	}
	radius = EPS + std::max(bench::common::eu_dist_square(box.min_corner(), q), 
					bench::common::eu_dist_square(box.max_corner(), q));
    
	#ifdef LOCAL_DEBUG
	cout << "[RangeQuery]: ";
	bench::common::print_point<Dim>(q, false);
	cout << " ";
	bench::common::print_box<Dim>(box);
	#endif
	
    kdtree->rangeQuery(&q[0], radius, ret_indexes, out_dist_sqr);
	
	Points results;
    results.reserve(ret_indexes.size());
    for (auto idx : ret_indexes) {
		 if (bench::common::is_in_box(kdtree->m_data[idx], box)) {
			results.emplace_back(kdtree->m_data[idx]);
		}
    }
    
	auto end = std::chrono::steady_clock::now();
    range_count++;
    range_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    return results;
}


inline size_t count() {
    return kdtree->kdtree_get_point_count();
}

void insert(Point& point) {}
void erase(Point& point) {}

void insert(size_t pid) {
	auto start = std::chrono::steady_clock::now();
	
	kdtree->addPoint(pid);
	
	auto end = std::chrono::steady_clock::now();
    insert_count ++;
    insert_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

bool erase(size_t pid) {
	auto start = std::chrono::steady_clock::now();
	bool ret= false;
	
	// const size_t k = 1;
    // std::vector<size_t> ret_indexes(k);
    // std::vector<double> out_dist_sqr(k);
	// Point q = _data[pid];
	// kdtree->query(&q[0], k, &ret_indexes[0], &out_dist_sqr[0]);
	// if (!ret_indexes.empty()) {
		// int rid = -1;
		// for (int sz=ret_indexes.size(),i=0; i<sz; ++i) {
			// if (out_dist_sqr[i]==0) {
				// rid = i;
				// if (ret_indexes[i] == pid) break;
			// }
		// }
		// if (rid >= 0) {
			// std::cout << "[ERASE] " << ret_indexes[0] << " " << pid << std::endl;
			// kdtree->removePoint(ret_indexes[0]);
			// ret = true;
		// }
	// }
	
	ret = kdtree->removePoint(pid);
	
	auto end = std::chrono::steady_clock::now();
    erase_count ++;
    erase_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	
	return ret;
}

private:
// internal kdtree using nanoflann
kdtree_t* kdtree;
Points& _data;

};

}}

