#ifndef Octree_H
#define Octree_H

#include <cmath>
#include <vector>
#include <unordered_set>

#include "../../utils/type.hpp"
#include "../../utils/common.hpp"
#include "../base_index.hpp"

// #define LOCAL_DEBUG

namespace bench { namespace index {

	/**!
	 *
	 */
	template<size_t dim, size_t MaxElements=1>
	class OctreeNode {
		using Point = point_t<dim>;
		using Box = box_t<dim>;
		using Points = std::vector<point_t<dim> >;
		using Counts = std::vector<int>;
		using PointID_t = int;
		using Indexes = std::vector<PointID_t>;
		
		public:
		static const int CHILD_SIZE = ((size_t)1)<<dim;
		
		// Physical position/size. This implicitly defines the bounding 
		// box of this node
		Point origin;         //! The physical center of this node
		Point halfDimension;  //! Half the width/height/depth of this node

		// The tree has up to eight children and can additionally store
		// a point, though in many applications only, the leaves will store data.
		OctreeNode *children[CHILD_SIZE]; //! Pointers to child octants
		Counts cnts; // number of same data points
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
			
			#ifdef LOCAL_DEBUG
			cout << "origin: ";
			bench::common::print_point<dim>(origin, true);
			cout << "halfDimension: ";
			bench::common::print_point<dim>(halfDimension, true);
			#endif
			
			// Initially, there are no children
			cnts.clear();
			ids.clear();
			for(int i=0; i<CHILD_SIZE; ++i) 
				children[i] = NULL;
			
			int idx = 0;
			for (auto point : points) {
				insert(points, idx++, 1);
			}
		}
		
		OctreeNode(const Point& _origin, const Point& _halfDimension) 
			: origin(_origin), halfDimension(_halfDimension) {
				// Initially, there are no children
				cnts.clear();
				ids.clear();
				for(int i=0; i<CHILD_SIZE; ++i) 
					children[i] = NULL;
			}

		~OctreeNode() {
			// Recursively destroy octants
			for(int i=0; i<CHILD_SIZE; ++i) {
				if (children[i] != NULL)
					delete children[i];
			}
		}
		
		inline size_t index_size() {
			size_t ret = sizeof(OctreeNode*)*CHILD_SIZE + sizeof(double)*dim*2;
			
			ret += sizeof(int*) + sizeof(int) + sizeof(int)*ids.size();
			ret += sizeof(PointID_t*) + sizeof(int) + sizeof(PointID_t)*ids.size();
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
		void insert(const Points& points, const int pid, int num=1, int dep=0) {
		#else
		void insert(const Points& points, const int pid, int num=1) {
		#endif
			
			const Point& point = points[pid];
			
			#ifdef LOCAL_DEBUG
			for (int i=0; i<dep; ++i) std::cout << " ";
			bench::common::print_point<dim>(point, false);
			std::cout << " ";
			bench::common::print_point<dim>(origin, false);
			std::cout << " ";
			bench::common::print_point<dim>(halfDimension, false);
			std::cout << " " << (isLeafNode() ? "Leaf" : "Non-Leaf") << " " << ids.size() << endl;	
			#endif		
			
			// If this node doesn't have a data point yet assigned 
			// and it is a leaf, then we're done!
			if(isLeafNode()) {
				assert(ids.size() == cnts.size());
				
				// check if the new point equals to any current point
				for (int i=0; i<ids.size(); ++i) {
					const Point& oldPoint = points[ids[i]];
					if (bench::common::is_equal_point<dim>(oldPoint, point)) {
						cnts[i] += num;
						#ifdef LOCAL_DEBUG
						for (int i=0; i<dep; ++i) std::cout << " ";
						bench::common::print_point<dim>(oldPoint, false);
						std::cout << " ";
						std::cout << cnts[i] << endl;
						#endif	
						return ;
					}		
				}
				if(ids.size() < MaxElements) {
					cnts.emplace_back(num);
					ids.emplace_back(pid);
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
					for (int i=0; i<ids.size(); ++i) {
						const Point& oldPoint = points[ids[i]];
						#ifdef LOCAL_DEBUG
						children[getOctantContainingPoint(oldPoint)]->insert(points, ids[i], cnts[i], dep+1);
						#else
						children[getOctantContainingPoint(oldPoint)]->insert(points, ids[i], cnts[i]);
						#endif
					}
					ids.clear();
					cnts.clear();
					#ifdef LOCAL_DEBUG
					children[getOctantContainingPoint(point)]->insert(points, pid, num, dep+1);
					#else
					children[getOctantContainingPoint(point)]->insert(points, pid, num);
					#endif
				}
			} else {
				// We are at an interior node. Insert recursively into the 
				// appropriate child octant
				int octant = getOctantContainingPoint(point);
				#ifdef LOCAL_DEBUG
				children[octant]->insert(points, pid, num, dep+1);
				#else
				children[octant]->insert(points, pid, num);
				#endif
			}
		}
		
		bool erase(const Points& points, const Point& point) {
			#ifdef LOCAL_DEBUG
			std::cout << "[DELETE]: ";
			bench::common::print_point<dim>(point, true);
			#endif
			bool success = false;
			__erase(points, point, success);
			#ifdef LOCAL_DEBUG
			std::cout << "\t\t\t" << (success ? "True":"False") << endl;
			#endif
			return success;
		}
		
		int __erase(const Points& points, const Point& point, bool& success) {
			#ifdef LOCAL_DEBUG
			bench::common::print_point<dim>(point, false);
			std::cout << " ";
			bench::common::print_point<dim>(origin, false);
			std::cout << " ";
			bench::common::print_point<dim>(halfDimension, false);
			std::cout << " " << (isLeafNode() ? "Leaf" : "Non-Leaf") << " " << ids.size() << endl;	
			#endif
			// If this node doesn't have a data point yet assigned 
			// and it is a leaf, then we're done!
			if(isLeafNode()) {
				
				#ifdef LOCAL_DEBUG
				std::cout << "\t" << "ids.size() = " << ids.size() << " ";
				for (auto id : ids) {
					bench::common::print_point<dim>(points[id], false);
					std::cout << " ";
				}
				std::cout << endl;
				#endif
				
				int erased = 0;
				
				for (int i=0; i<ids.size(); ++i) {
					if (bench::common::is_equal_point<dim>(points[ids[i]], point)) {
						--cnts[i];
						if (cnts[i] == 0) {
							cnts[i] = *cnts.rbegin();
							cnts.pop_back();						
							ids[i] = *ids.rbegin();
							ids.pop_back();
							erased = 1;	
						}
						success = true;
						#ifdef LOCAL_DEBUG
						std::cout << "\t\terased = " << erased << " cnts[i] = " << cnts[i] << endl;
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
				int erased = children[octant]->__erase(points, point, success);
				
				if (0 == erased) return erased;
				
				size_t dataSize = 0;
				for (int j=0; j<CHILD_SIZE; ++j) {
					if (!children[j]->isLeafNode()) 
						return 0;
					else
						dataSize += children[j]->ids.size();
				}
				
				// delete all child nodes and merge the child's data into parent node
				if (dataSize > MaxElements) return 0;
				
				for (int j=0; j<CHILD_SIZE; ++j) {
					if (!children[j]->ids.empty()) {
						ids.insert(ids.end(), children[j]->ids.begin(), children[j]->ids.end());
						cnts.insert(cnts.end(), children[j]->cnts.begin(), children[j]->cnts.end());
					}
					delete children[j];
					children[j] = NULL;
				}
				
				return 1;
			}
			
			return 0;
		}
		
		std::vector<PointID_t> range_query(const Points& points, const Box& box) {
			std::vector<PointID_t> results;
			std::unordered_set<PointID_t> visit;
			
			getPointsInsideBox(points, box, results, visit);
			
			return results;
		}
		
		// This is a really simple routine for querying the tree for points
		// within a bounding box defined by min/max points (bmin, bmax)
		// All results are pushed into 'results'
		void getPointsInsideBox(const Points& points, const Box& box, std::vector<PointID_t>& results, std::unordered_set<PointID_t>& visit) {
			// If we're at a leaf node, just see if the current data point is inside
			// the query bounding box
			if(isLeafNode()) {
				for (int i=ids.size()-1; i>=0; --i) {
					if (bench::common::is_in_box<dim>(points[ids[i]], box)) {
						if (visit.count(ids[i]) == 0) {
							visit.insert(ids[i]);
							results.insert(results.end(), cnts[i], ids[i]);
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
					
					children[j]->getPointsInsideBox(points, box, results, visit);
				} 
			}
		}
		
		using QueueElement = std::pair<double, PointID_t>;
		using knnQueue_t = std::priority_queue<QueueElement, std::vector<QueueElement>, std::less<QueueElement> > ;
		
		std::vector<PointID_t> knn_query(const Points& points, const Point& q, const size_t k) {
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
			
			std::vector<PointID_t> results;
			knnQueue_t Q;
			getPointsKnn(points, q, k, box, radius, Q);
			
			while (!Q.empty()) {
				results.emplace_back(Q.top().second);
				Q.pop();
			}
			
			return results;
		}

		void getPointsKnn(const Points& points, const Point& q, const size_t k, Box& box, double& radius, knnQueue_t& Q) {
			// If we're at a leaf node, just see if the current data point is inside
			// the query bounding box
			if(isLeafNode()) {
				if(!ids.empty()) {
					for (int i=0; i<ids.size(); ++i) {
						const Point& p = points[ids[i]];
						double new_dist = bench::common::eu_dist<dim>(p, q);
						for (int j=0; j<cnts[i]; ++j) {
							if (Q.size() < k) {
								Q.push(std::make_pair(bench::common::eu_dist<dim>(p, q), ids[i]));
							} else {
								if (new_dist < Q.top().first) {
									Q.pop();
									Q.push(std::make_pair(bench::common::eu_dist<dim>(p, q), ids[i]));
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
					
					children[j]->getPointsKnn(points, q, k, box, radius, Q);
				} 
			}			
		}
	};

	template<size_t dim=2, size_t MaxElements=1>
	class Octree : public BaseIndex {
		using Point = point_t<dim>;
		using Box = box_t<dim>;
		using Points = std::vector<Point>;
		using PointID_t = int;
		using OctreeNode_t = OctreeNode<dim, MaxElements>;
		
		private:
		size_t num_of_points;
		OctreeNode_t *root;
		Points& _data;
		
		/*
				Children follow a predictable pattern to make accesses simple.
				Here, - means less than 'origin' in that dimension, + means greater than.
				child:	0 1 2 3 4 5 6 7
				x:      - - - - + + + +
				y:      - - + + - - + +
				z:      - + - + - + - +
		 */

		public:
		Octree(Points& points) : _data(points) {
			std::cout << "Construct Octree: MaxElements = " << MaxElements << std::endl;
			auto start = std::chrono::steady_clock::now();

			this->num_of_points = points.size();
			root = new OctreeNode_t(points);
			
			auto end = std::chrono::steady_clock::now();
			build_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			std::cout << "Build Time: " << get_build_time() << " [ms]" << std::endl;
			std::cout << "Index Size: " << index_size() << " Bytes" << std::endl;
		}
		
		~Octree() {
			delete root;
		}

		Points range_query(const Box& box) {
			auto start = std::chrono::steady_clock::now();
			
			Points results;
			if (root != NULL) {
				std::vector<PointID_t> ids = root->range_query(_data, box);
				results.reserve(ids.size());
				for (auto id : ids)
					results.emplace_back(_data[id]);
			}
			
			auto end = std::chrono::steady_clock::now();
			range_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
			range_count ++;
			
			return results;
		}
		
		Points knn_query(const Point& q, const size_t k) {
			auto start = std::chrono::steady_clock::now();
			
			Points results;
			if (root != NULL) {
				std::vector<PointID_t> ids = root->knn_query(_data, q, k);
				results.reserve(ids.size());
				for (auto id : ids)
					results.emplace_back(_data[id]);
			}
			
			auto end = std::chrono::steady_clock::now();
			range_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
			range_count ++;
			
			return results;
		}
		
		void insert(const Point& point) {
			if (root != NULL) {
				int pid = _data.size();
				_data.emplace_back(point);
				root->insert(_data, pid, 1);
				#ifdef LOCAL_DEBUG
				cout << "[INSERT]: pid = " << pid << " ";
				bench::common::print_point(point, true);
				#endif
			}
		}
		
		bool erase(const Point& point) {
			if (root != NULL) {
				return root->erase(_data, point);
			}
			return false;
		}
		
		inline size_t count() {
			return this->num_of_boxes;
		}

		inline size_t index_size() {
			size_t ret = sizeof(size_t) + sizeof(OctreeNode_t*);
			if (root != NULL)
				ret += root->index_size();
			return ret;
		}
	};
}
}
#endif
