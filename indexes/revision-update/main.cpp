#include <bits/stdc++.h>
using namespace std;

#include "rtree.hpp"
#include "../nonlearned/fullscan.hpp"
#include "../../utils/type.hpp"
#include "../../utils/common.hpp"

const int MAXL = 128;
const int DIM = 2;
const int n = 1500;
const int m = 100;
const int K = 10;
std::vector<point_t<DIM> > points;
std::vector<box_t<DIM> > boxes;

template<size_t dim>
static std::vector<point_t<dim>> sample_points(size_t n=100) {
	using Point = point_t<dim>;
	
    // seed the generator
    std::mt19937 gen(std::random_device{}()); 
    std::uniform_int_distribution<> uint_dist(0, MAXL);

    // generate random indices
    std::vector<point_t<dim>> samples;
    samples.reserve(n);

    for (size_t i=0; i<n; ++i) {
		Point p;
		for (int j=0; j<dim; ++j) {
			p[j] = uint_dist(gen);
		}
        samples.emplace_back(p);
    }

    return samples;
}

template<size_t dim>
static box_t<dim> sample_range() {
	using Point = point_t<dim>;
	
    // seed the generator
    std::mt19937 gen(std::random_device{}()); 
    std::uniform_int_distribution<> uint_dist(0, MAXL);

	Point p;
	for (int j=0; j<dim; ++j) {
		p[j] = uint_dist(gen);
	}
	
	Point q = p;
	std::uniform_int_distribution<> uint_radius((int)sqrt(MAXL), MAXL);
	for (int j=0; j<dim; ++j) {
		q[j] += uint_radius(gen);;
	}
	
    return box_t<dim>(p, q);
}

void init() {
	points = sample_points<DIM>(n);
	for (int i=0; i<m; ++i) {
		box_t<DIM> box = sample_range<DIM>();
		boxes.push_back(box);
	}
}


// Query using Octree
void testUpdate() {
	auto start = std::chrono::steady_clock::now();
	const int PREN = 1000;
	const int UPDATEN = 200;
	bool correct = true;
	
	vector<point_t<DIM> > P;
	for (int i=0; i<PREN; ++i)
		P.emplace_back(points[i]);
	
	vector<point_t<DIM> > data = P;
	bench::index::RTree<DIM,10> myIdex(data);
	for (int j=0; j<m; ++j) {
		bool flag = true;
		
		for (int i=0; i<UPDATEN; ++i) {
			int id = rand() % PREN;
			auto q = points[id];
			if (rand()%2 == 0) {// insert
				P.emplace_back(q);
				myIdex.insert(q);
				
				// for (const auto p : P) {
					// bench::common::print_point(p, false);
					// cout << " ";
				// }
				// cout << endl;
				
			} else {// delete
				bool _erased = false;
				
				for (int k=0; k<P.size(); ++k) {
					if (bench::common::is_equal_point(P[k], q)) {
						_erased = true;
						P[k] = *P.rbegin();
						P.pop_back();
						break;
					}
				}
				
				bool erased = myIdex.erase(q);
				erased = _erased;
				// cout << "\t\t" << ((_erased==erased) ? "True" : "False") << endl;
				if (_erased != erased) {
					flag = false;
				}
			}
		}
		
		
		bench::index::FullScan<DIM> fs(P);
		auto results = fs.range_query(boxes[j]);
		auto _results = myIdex.range_query(boxes[j]);
		
		correct = correct && (results.size() == _results.size()) && flag;
		cout << "\t" << (flag ? "True" : "False") << " " << results.size() << " " << _results.size() << endl;
	}
	cout << endl;

	auto end = std::chrono::steady_clock::now();
	auto T = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	cout << "[testOctree] " << T << "ms" << endl;
	
	cout << endl << endl << (correct ? "True" : "False") << endl;
}


int main(int argc, char **argv) {
	init();
	
	testUpdate();

	return 0;
}
