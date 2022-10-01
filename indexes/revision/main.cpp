#include <bits/stdc++.h>
using namespace std;

#include "Octree.hpp"
#include "../nonlearned/fullscan.hpp"
#include "../../utils/type.hpp"
#include "../../utils/common.hpp"

const int MAXL = 128;
const int DIM = 3;
const int n = 150;
const int m = 10;
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

// Query using brute-force
vector<int> testNaive() {
	auto start = std::chrono::steady_clock::now();
	vector<int> ret;
	
	bench::index::FullScan<DIM> fs(points);
	for (int j=0; j<m; ++j) {
		auto results = fs.range_query(boxes[j]);
		cout << results.size() << " ";
		// bench::common::print_box<DIM>(boxes[j]);
		ret.push_back((int) results.size());
	}
	cout << endl;
	
	auto end = std::chrono::steady_clock::now();
	auto T = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	cout << "[testNaive] " << T << "ms" << endl;
	
	return ret;
}

// Query using Octree
vector<int> testOctree() {
	auto start = std::chrono::steady_clock::now();
	vector<int> ret;
	
	bench::index::Octree<DIM,5> octree(points);
	for (int j=0; j<m; ++j) {
		auto results = octree.range_query(boxes[j]);
		cout << results.size() << " ";
		// bench::common::print_box<DIM>(boxes[j]);
		ret.push_back((int) results.size());
	}
	cout << endl;

	auto end = std::chrono::steady_clock::now();
	auto T = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	cout << "[testOctree] " << T << "ms" << endl;
	
	return ret;
}

// Query using brute-force
vector<double> testNaiveKNN() {
	auto start = std::chrono::steady_clock::now();
	vector<double> ret;
	
	bench::index::FullScan<DIM> fs(points);
	for (int j=0; j<m; ++j) {
		auto q = boxes[j].min_corner();
		auto results = fs.knn_query(q, K);
		double knnd = 0.0;
		for (auto point : results) {
			double tmp = bench::common::eu_dist(point, q);
			if (tmp > knnd)
				knnd = tmp;
		}
		cout << knnd << " ";
		ret.push_back(knnd);
	}
	cout << endl;
	
	auto end = std::chrono::steady_clock::now();
	auto T = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	cout << "[testNaive] " << T << "ms" << endl;
	
	return ret;
}

// Query using Octree
vector<double> testOctreeKNN() {
	auto start = std::chrono::steady_clock::now();
	vector<double> ret;
	
	bench::index::Octree<DIM,5> octree(points);
	for (int j=0; j<m; ++j) {
		auto q = boxes[j].min_corner();
		auto results = octree.knn_query(q, K);
		double knnd = 0.0;
		for (auto point : results) {
			double tmp = bench::common::eu_dist(point, q);
			if (tmp > knnd)
				knnd = tmp;
		}
		cout << knnd << " ";
		ret.push_back(knnd);
	}
	cout << endl;

	auto end = std::chrono::steady_clock::now();
	auto T = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	cout << "[testOctree] " << T << "ms" << endl;
	
	return ret;
}

// Query using Octree
void testUpdate() {
	auto start = std::chrono::steady_clock::now();
	const int PREN = 100;
	const int UPDATEN = 100;
	bool correct = true;
	
	vector<point_t<DIM> > P;
	for (int i=0; i<PREN; ++i)
		P.emplace_back(points[i]);
	
	vector<point_t<DIM> > data = P;
	bench::index::Octree<DIM,1> octree(data);
	for (int j=0; j<m; ++j) {
		bool flag = true;
		
		for (int i=0; i<UPDATEN; ++i) {
			int id = rand() % PREN;
			auto q = points[id];
			if (rand()%2 == 1) {// insert
				P.emplace_back(q);
				octree.insert(q);
				
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
				
				bool erased = octree.erase(q);
				cout << "\t\t" << ((_erased==erased) ? "True" : "False") << endl;
				if (_erased != erased) {
					flag = false;
				}
			}
			
			// cout << "\t\t" << data.size() << " " << P.size() << endl;
		}
		
		
		bench::index::FullScan<DIM> fs(P);
		auto results = fs.range_query(boxes[j]);
		auto _results = octree.range_query(boxes[j]);
		
		correct = correct && flag;
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
	// auto va = testNaive();
	// auto vb = testOctree();

	// bool flag = true;
	// for (int i=0; i<m; ++i) {
		// if (va[i] != vb[i]) {
			// flag = false;
			// break;
		// }
	// }
	// cout << (flag ? "True" : "False") << endl;
	
	testUpdate();

	return 0;
}
