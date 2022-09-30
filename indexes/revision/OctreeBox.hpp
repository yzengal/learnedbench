#ifndef Octree_H
#define Octree_H

#include <vector>
#include <cmath>
#include <unordered_set>

#include "../../utils/type.hpp"
#include "../../utils/common.hpp"

// #define LOCAL_DEBUG

namespace bench { namespace index {

	/**!
	 *
	 */
	template<size_t dim=2, size_t MaxElements=1>
	class OctreeBox {
		using Point = point_t<dim>;
		using Box = box_t<dim>;
		using Boxes = std::vector<box_t<dim>>;
		using Counts = std::vector<int>;
		using Indexes = std::vector<int>;
		using BoxID_t = int;
		
		public:
		static const int CHILD_SIZE = ((size_t)1)<<dim;
		static int MAX_DEPTH;
		
		// Physical position/size. This implicitly defines the bounding 
		// box of this node
		Point origin;         //! The physical center of this node
		Point halfDimension;  //! Half the width/height/depth of this node

		// The tree has up to eight children and can additionally store
		// a point, though in many applications only, the leaves will store data.
		OctreeBox *children[CHILD_SIZE]; 	//! Pointers to child octants
		Boxes data;   				 	//! Data point to be stored at a node
		Counts cnts; 					// number of same data points
		Indexes	ids;					// ID of data points

		/*
				Children follow a predictable pattern to make accesses simple.
				Here, - means less than 'origin' in that dimension, + means greater than.
				child:	0 1 2 3 4 5 6 7
				x:      - - - - + + + +
				y:      - - + + - - + +
				z:      - + - + - + - +
		 */

		public:
		OctreeBox(const Boxes& boxes) {
			int n = boxes.size();
			MAX_DEPTH = max(1.0, 1.0+log2(n*1.0));
			
			origin.fill(0.0);
			halfDimension.fill(0.0);
			
			#ifdef LOCAL_DEBUG
			std::cout << "n: " << n << " " << "MAX_DEPTH: " << MAX_DEPTH << " ";
			std::cout << "origin: ";
			bench::common::print_point<dim>(origin, true);
			#endif
			
			for (int i=0; i<n; ++i) {
				Point bound = bench::common::get_boundary_point2box(origin, boxes[i]);
				for (int j=0; j<dim; ++j) {
					halfDimension[j] = max(halfDimension[j], bound[j]);
				}
			}
			
			#ifdef LOCAL_DEBUG
			cout << "halfDimension: ";
			bench::common::print_point<dim>(halfDimension, true);
			#endif
			
			// Initially, there are no children
			data.clear();
			cnts.clear();
			ids.clear();
			for(int i=0; i<CHILD_SIZE; ++i) 
				children[i] = NULL;
			
			int idx = 0;
			for (auto box : boxes) {
				insert(box, idx++);
				#ifdef LOCAL_DEBUG
				std::cout << idx++ << " ";
				bench::common::print_box(box, true);
				std::cout << endl;
				#endif
			}
		}
		
		OctreeBox(const Point& _origin, const Point& _halfDimension) 
			: origin(_origin), halfDimension(_halfDimension) {
				// Initially, there are no children
				data.clear();
				cnts.clear();
				ids.clear();
				for(int i=0; i<CHILD_SIZE; ++i) 
					children[i] = NULL;
			}

		OctreeBox(const OctreeBox& copy)
			: origin(copy.origin), halfDimension(copy.halfDimension), data(copy.data) {
				
			}

		~OctreeBox() {
			// Recursively destroy octants
			for(int i=0; i<CHILD_SIZE; ++i) {
				if (children[i] != NULL)
					delete children[i];
			}
		}

		bool isLeafNode() const {
			// This is correct, but overkill. See below.
			/*
				 for(int i=0; i<CHILD_SIZE; ++i)
				 if(children[i] != NULL) 
				 return false;
				 return true;
			 */

			// We are a leaf iff we have no children. Since we either have none, or 
			// all eight, it is sufficient to just check the first.
			return children[0] == NULL;
		}

		// Determine which octant of the tree would contain 'point'
		std::vector<int> getOctantContainingPoints(const Box& box) {
			std::vector<int> ret;
			Point mxp = box.max_corner(), mnp = box.min_corner();
			Point p;
			
			for (int st=0; st<(1<<dim); ++st) {
				for (int i=0; i<dim; ++i) {
					p[i] = (st & (1<<i)) ? mxp[i] : mnp[i];
				}
				auto oct = getOctantContainingPoint(p);
				bool flag = true;
				for (auto st : ret) {
					if (st == oct) {
						flag = false;
						break;
					}
				}
				if (flag) ret.emplace_back(oct);
			}
			
			return ret;
		}
		
		int getOctantContainingPoint(const Point& point) const {
			int oct = 0;
			for (int i=0; i<dim; ++i) {
				if (point[i] >= origin[i])
					oct |= (1 << i);
			}
			return oct;
		}
		
		void insert(const Box& box, int pid, int num=1, int dep=0) {
			
			#ifdef LOCAL_DEBUG
			for (int i=0; i<dep; ++i) std::cout << " ";
			std::cout << dep << ": ";
			bench::common::print_box<dim>(box, false);
			std::cout << " ";
			bench::common::print_point<dim>(origin, false);
			std::cout << " ";
			bench::common::print_point<dim>(halfDimension, false);
			std::cout << " " << (isLeafNode() ? "Leaf" : "Non-Leaf") << " " << data.size() << endl;	
			#endif			
			// If this node doesn't have a data point yet assigned 
			// and it is a leaf, then we're done!
			if(isLeafNode()) {
				// check if the new point equals to any current point
				for (int i=0; i<data.size(); ++i) {
					if (bench::common::is_equal_box<dim>(data[i], box)) {
						++cnts[i];
						return ;
					}		
				}
				if(dep>=MAX_DEPTH || data.size()<MaxElements) {
					data.emplace_back(box);
					cnts.emplace_back(num);
					ids.emplace_back(pid);
				} else {
					// We're at a leaf, but there's already something here
					// We will split this node so that it has 8 child octants
					// and then insert the old data that was here, along with 
					// this new data point

					// Save this data point that was here for a later re-insert
					// OctreePoint *oldPoint = data;
					// data = NULL;

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
						children[j] = new OctreeBox(newOrigin, newHalfDimension);
					}

					// Re-insert the old point, and insert this new point
					// (We wouldn't need to insert from the root, because we already
					// know it's guaranteed to be in this section of the tree)
					for (int i=0; i<data.size(); ++i) {
						vector<int> octs = getOctantContainingPoints(data[i]);
						for (auto oct : octs) {		
							children[oct]->insert(data[i], ids[i], cnts[i], dep+1);
						}
					}
					data.clear();
					cnts.clear();
					ids.clear();
					vector<int> octs = getOctantContainingPoints(box);
					for (auto oct : octs) {	
						children[oct]->insert(box, pid, num, dep+1);
					}
				}
			} else {
				// We are at an interior node. Insert recursively into the 
				// appropriate child octant
				vector<int> octs = getOctantContainingPoints(box);
				for (auto oct : octs) {
					children[oct]->insert(box, pid, num, dep+1);
				}
			}
		}
		
		Boxes range_query(const Box& box) {
			Boxes results;
			std::unordered_set<int> visit;
			getBoxesInsideBox(box, results, visit);
			return results;
		}
		
		// This is a really simple routine for querying the tree for points
		// within a bounding box defined by min/max points (bmin, bmax)
		// All results are pushed into 'results'
		void getBoxesInsideBox(const Box& box, Boxes& results, std::unordered_set<int>& visit) {
			// If we're at a leaf node, just see if the current data point is inside
			// the query bounding box
			if(isLeafNode()) {
				if(!data.empty()) {
					for (int i=data.size()-1; i>=0; --i) {
						if (bench::common::is_intersect_box<dim>(data[i], box)) {
							if (visit.count(ids[i]) == 0) {
								visit.insert(ids[i]);
								results.insert(results.end(), cnts[i], data[i]);
							}
						}
					}
				}
			} else {
				// We're at an interior node of the tree. We will check to see if
				// the query bounding box lies outside the octants of this node.
				
				for(int j=0; j<CHILD_SIZE; ++j) {
					// Compute the min/max corners of this child octant
					Point cmax = children[j]->origin, cmin = children[j]->origin;
					for (int i=0; i<dim; ++i) {
						cmax[i] += children[i]->halfDimension[i];
						cmin[i] -= children[i]->halfDimension[i];
					}
					
					// If the query rectangle is outside the child's bounding box, 
					// then continue
					Box cbox(cmin, cmax);
					if (!bench::common::is_intersect_box<dim>(cbox, box)) continue;
					
					children[j]->getBoxesInsideBox(box, results, tb);
				} 
			}
		}

	};
	
	template<size_t dim, size_t MaxElements>
	int OctreeBox<dim,MaxElements>::MAX_DEPTH;

}
}
#endif
