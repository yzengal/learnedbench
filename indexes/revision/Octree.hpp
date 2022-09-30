#ifndef Octree_H
#define Octree_H

#include <cmath>
#include <vector>

#include "../../utils/type.hpp"
#include "../../utils/common.hpp"

// #define LOCAL_DEBUG

namespace bench { namespace index {

	/**!
	 *
	 */
	template<size_t dim, size_t MaxElements=1>
	class Octree {
		using Point = point_t<dim>;
		using Box = box_t<dim>;
		using Points = std::vector<point_t<dim> >;
		using Counts = std::vector<int>;
		
		public:
		static const int CHILD_SIZE = ((size_t)1)<<dim;
		
		// Physical position/size. This implicitly defines the bounding 
		// box of this node
		Point origin;         //! The physical center of this node
		Point halfDimension;  //! Half the width/height/depth of this node

		// The tree has up to eight children and can additionally store
		// a point, though in many applications only, the leaves will store data.
		Octree *children[CHILD_SIZE]; //! Pointers to child octants
		Points data;   //! Data point to be stored at a node
		Counts cnts; // number of same data points

		/*
				Children follow a predictable pattern to make accesses simple.
				Here, - means less than 'origin' in that dimension, + means greater than.
				child:	0 1 2 3 4 5 6 7
				x:      - - - - + + + +
				y:      - - + + - - + +
				z:      - + - + - + - +
		 */

		public:
		Octree(const Points& points) {
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
			
			// #ifdef LOCAL_DEBUG
			// cout << "origin: ";
			// bench::common::print_point<dim>(origin, true);
			// cout << "halfDimension: ";
			// bench::common::print_point<dim>(halfDimension, true);
			// #endif
			
			// Initially, there are no children
			data.clear();
			cnts.clear();
			for(int i=0; i<CHILD_SIZE; ++i) 
				children[i] = NULL;
			
			int idx = 0;
			for (auto point : points) {
				insert(point);
				// #ifdef LOCAL_DEBUG
				// cout << idx++ << " ";
				// bench::common::print_point(point, true);
				// cout << endl;
				// #endif
			}
		}
		
		Octree(const Point& _origin, const Point& _halfDimension) 
			: origin(_origin), halfDimension(_halfDimension) {
				// Initially, there are no children
				data.clear();
				cnts.clear();
				for(int i=0; i<CHILD_SIZE; ++i) 
					children[i] = NULL;
			}

		Octree(const Octree& copy)
			: origin(copy.origin), halfDimension(copy.halfDimension), data(copy.data) {
				
			}

		~Octree() {
			// Recursively destroy octants
			for(int i=0; i<CHILD_SIZE; ++i) {
				if (children[i] != NULL)
					delete children[i];
			}
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
		
		#ifdef LOCAL_DEBUG
		void insert(const Point& point, int num=1, int dep=0) {
		#else
		void insert(const Point& point, int num=1) {
		#endif
			
			#ifdef LOCAL_DEBUG
			for (int i=0; i<dep; ++i) std::cout << " ";
			bench::common::print_point<dim>(point, false);
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
					Point& oldPoint = data[i];
					if (bench::common::is_equal_point<dim>(oldPoint, point)) {
						++cnts[i];
						return ;
					}		
				}
				if(data.size() < MaxElements) {
					data.emplace_back(point);
					cnts.emplace_back(num);
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
						children[j] = new Octree(newOrigin, newHalfDimension);
					}

					// Re-insert the old point, and insert this new point
					// (We wouldn't need to insert from the root, because we already
					// know it's guaranteed to be in this section of the tree)
					for (int i=0; i<data.size(); ++i) {
						Point& oldPoint = data[i];
						#ifdef LOCAL_DEBUG
						children[getOctantContainingPoint(oldPoint)]->insert(oldPoint, cnts[i], dep+1);
						#else
						children[getOctantContainingPoint(oldPoint)]->insert(oldPoint, cnts[i]);
						#endif
					}
					data.clear();
					cnts.clear();
					#ifdef LOCAL_DEBUG
					children[getOctantContainingPoint(point)]->insert(point, num, dep+1);
					#else
					children[getOctantContainingPoint(point)]->insert(point, num);
					#endif
				}
			} else {
				// We are at an interior node. Insert recursively into the 
				// appropriate child octant
				int octant = getOctantContainingPoint(point);
				#ifdef LOCAL_DEBUG
				children[octant]->insert(point, num, dep+1);
				#else
				children[octant]->insert(point, num);
				#endif
			}
		}
		
		bool erase(const Point& point) {
			#ifdef LOCAL_DEBUG
			std::cout << "[DELETE]: ";
			bench::common::print_point<dim>(point, true);
			#endif
			bool success = false;
			__erase(point, success);
			#ifdef LOCAL_DEBUG
			std::cout << "\t\t\t" << (success ? "True":"False") << endl;
			#endif
			return success;
		}
		
		int __erase(const Point& point, bool& success) {
			#ifdef LOCAL_DEBUG
			bench::common::print_point<dim>(point, false);
			std::cout << " ";
			bench::common::print_point<dim>(origin, false);
			std::cout << " ";
			bench::common::print_point<dim>(halfDimension, false);
			std::cout << " " << (isLeafNode() ? "Leaf" : "Non-Leaf") << " " << data.size() << endl;	
			#endif
			// If this node doesn't have a data point yet assigned 
			// and it is a leaf, then we're done!
			if(isLeafNode()) {
				assert(data.size() == cnts.size());
				
				#ifdef LOCAL_DEBUG
				std::cout << "\t" << "data.size() = " << data.size() << " ";
				for (auto p : data) {
					bench::common::print_point<dim>(p, false);
					std::cout << " ";
				}
				std::cout << endl;
				#endif
				
				int erased = 0;
				
				for (int i=0; i<data.size(); ++i) {
					#ifdef LOCAL_DEBUG
					std::cout << "\t\t" << i << " " << data.size() << " ";
					#endif
					if (bench::common::is_equal_point<dim>(data[i], point)) {
						if (--cnts[i] == 0) {						
							data[i] = *data.rbegin();
							data.pop_back();
							cnts[i] = *cnts.rbegin();
							cnts.pop_back();
							erased = 1;	
						}
						success = true;
						#ifdef LOCAL_DEBUG
						std::cout << "\t\t" << erased << endl;
						#endif
						break;
					}
					#ifdef LOCAL_DEBUG
					std::cout << "\t\t\t" << (erased ? "erased":"non-erased") << endl;
					#endif
				}
				
				return erased;
				
			} else {
				// We are at an interior node. Insert recursively into the 
				// appropriate child octant
				int octant = getOctantContainingPoint(point);
				int erased = children[octant]->__erase(point, success);
				
				if (0 == erased) return erased;
				
				size_t dataSize = 0;
				for (int j=0; j<CHILD_SIZE; ++j) {
					if (!children[j]->isLeafNode()) 
						return 0;
					else
						dataSize += children[j]->data.size();
				}
				
				// delete all child nodes and merge the child's data into parent node
				if (dataSize > MaxElements) return 0;
				
				for (int j=0; j<CHILD_SIZE; ++j) {
					if (!children[j]->data.empty()) {
						data.insert(data.end(), children[j]->data.begin(), children[j]->data.end());
						cnts.insert(cnts.end(), children[j]->cnts.begin(), children[j]->cnts.end());
					}
					delete children[j];
					children[j] = NULL;
				}
				
				return 1;
			}
			
			return 0;
		}
		
		Points range_query(const Box& box) {
			Points results;
			getPointsInsideBox(box, results);
			return results;
		}
		
		// This is a really simple routine for querying the tree for points
		// within a bounding box defined by min/max points (bmin, bmax)
		// All results are pushed into 'results'
		void getPointsInsideBox(const Box& box, Points& results) {
			// If we're at a leaf node, just see if the current data point is inside
			// the query bounding box
			if(isLeafNode()) {
				if(!data.empty()) {
					for (int i=data.size()-1; i>=0; --i) {
						Point& point = data[i];
						if (bench::common::is_in_box<dim>(point, box)) {
							results.insert(results.end(), cnts[i], point);
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
					
					children[j]->getPointsInsideBox(box, results);
				} 
			}
		}
		
		using PointIterator = typename Points::const_iterator;
		using QueueElement = std::pair<double, PointIterator>;
		using knnQueue_t = std::priority_queue<QueueElement, std::vector<QueueElement>, std::less<QueueElement> > ;
		
		Points knn_query(const Point& q, const size_t k) {
			Point zero;
			zero.fill(0);
			const double EPS = 1e-5;
			double radius = bench::common::eu_dist<dim>(origin, q) + bench::common::eu_dist<dim>(zero, halfDimension);
			
			Point cmin = q, cmax = q;
			for (int i=0; i<dim; ++i) {
				cmin[i] -= radius;
				cmax[i] += radius;
			}
			Box box(cmin, cmax);
			
			Points results;
			knnQueue_t Q;
			getPointsKnn(q, k, box, radius, Q);
			
			while (!Q.empty()) {
				PointIterator iter = Q.top().second;
				results.emplace_back(*iter);
				Q.pop();
			}
			
			return results;
		}

		void getPointsKnn(const Point& q, const size_t k, Box& box, double& radius, knnQueue_t& Q) {
			// If we're at a leaf node, just see if the current data point is inside
			// the query bounding box
			if(isLeafNode()) {
				if(!data.empty()) {
					int idx = 0;
					for (auto iter=data.begin(); iter!=data.end(); ++iter, ++idx) {
						Point& p = *iter;
						double new_dist = bench::common::eu_dist<dim>(p, q);
						for (int j=0; j<cnts[idx]; ++j) {
							if (Q.size() < k) {
								Q.push(std::make_pair(bench::common::eu_dist<dim>(p, q), iter));
							} else {
								if (new_dist < Q.top().first) {
									Q.pop();
									Q.push(std::make_pair(bench::common::eu_dist<dim>(p, q), iter));
								}
							}
						}
					}
					if (radius > Q.top().first) {
						radius = Q.top().first;
						Point cmin = q, cmax = q;
						for (int i=0; i<dim; ++i) {
							cmin[i] -= radius;
							cmax[i] += radius;
						}
						box = Box(cmin, cmax);
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
					
					children[j]->getPointsKnn(q, k, box, radius, Q);
				} 
			}			
		}
	};

}
}
#endif
