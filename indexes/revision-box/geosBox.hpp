#pragma once



#include <geos/geom/Coordinate.h>
#include <geos/geom/Envelope.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/index/quadtree/Quadtree.h>
#include <geos/index/strtree/TemplateSTRtree.h>

#ifdef HEAP_PROFILE
#include <gperftools/heap-profiler.h>
#endif

#include <chrono>

#include "base_index.hpp"
#include "../../utils/type.hpp"
#include "../../utils/common.hpp"


namespace bench { namespace index {
    
template<size_t Dim=2>
class QDTree : public BaseIndex {

using QDTree_t = geos::index::quadtree::Quadtree;
using Point = point_t<Dim>;
using Points = std::vector<Point>;
using Box = box_t<Dim>;
using Boxes = std::vector<box_t<Dim>>;

public:
QDTree(Boxes& boxes) : _data(boxes) { 
    std::cout << "Construct QDTree " << "Dim=" << Dim << std::endl;
    auto start = std::chrono::steady_clock::now();
	
	this->N = boxes.size();
    ids = new int[N];
	for (int i=0; i<N; ++i) ids[i] = i;
	
#ifdef HEAP_PROFILE
HeapProfilerStart("qdtree");
#endif

    _qdtree = new QDTree_t();

    for (size_t i=0; i<boxes.size(); ++i) {
		Box box = boxes[i];
		geos::geom::Envelope* enve = new geos::geom::Envelope(box.min_corner()[0], box.max_corner()[0], box.min_corner()[1], box.max_corner()[1]);
        _qdtree->insert(enve, reinterpret_cast<void *>(ids + i));
		delete enve;
    }

#ifdef HEAP_PROFILE
HeapProfilerDump("final");
HeapProfilerStop();
#endif

    auto end = std::chrono::steady_clock::now();
    build_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Build Time: " << get_build_time() << " [ms]" << std::endl;
}


~QDTree() {
    delete _qdtree;
	delete[] ids;
}


Boxes range_query(Box& box) {
    // prepare the query geometry
    geos::geom::GeometryFactory::Ptr gf = geos::geom::GeometryFactory::create();

    geos::geom::Envelope geos_box(box.min_corner()[0], box.max_corner()[0], box.min_corner()[1], box.max_corner()[1]);
    auto query_geometry = gf->toGeometry(&geos_box);

    // query the geometry envolope
    auto start = std::chrono::steady_clock::now();

    std::vector<void*> results;
    _qdtree->query(query_geometry->getEnvelopeInternal(), results);

    auto end = std::chrono::steady_clock::now();
    range_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    range_count ++;

    // collect results
    // the geos libaray returns geometries **may** intersect the query rectangle
    Boxes boxes_results;
    for (size_t i=0; i<results.size(); ++i) {
        auto pg = static_cast<int*>(results[i]);
		int id = *pg;
		assert(id>=0 && id<this->N);
        if (bench::common::is_intersect_box<Dim>(_data[id], box)) {
            boxes_results.emplace_back(_data[id]);
        }
        
    }
    return boxes_results;
}


inline size_t count() {
    return _qdtree->size() + sizeof(int)*N;
}


private:
Boxes& _data;
QDTree_t* _qdtree;
int N;
int* ids;
};


}
}
