#include <bits/stdc++.h>
using namespace std;

#include "floodOct.hpp"
#include "../nonlearned/fullscan.hpp"
#include "../../utils/type.hpp"
#include "../../utils/common.hpp"

const int MAXL = 128;
const int DIM = 2;
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


// Query using KDtree
vector<int> testKDtree() {
	auto start = std::chrono::steady_clock::now();
	vector<int> ret;
	
	bench::index::FloodOct<DIM,24> kdtree(points);
	for (int j=0; j<m; ++j) {
		
		auto results = kdtree.range_query(boxes[j]);
		cout << results.size() << " ";
		ret.push_back((int) results.size());
	}
	cout << endl;

	auto end = std::chrono::steady_clock::now();
	auto T = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	cout << "[testKDtree] " << T << "ms" << endl;
	
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
// vector<double> testOctreeKNN() {
	// auto start = std::chrono::steady_clock::now();
	// vector<double> ret;
	
	// bench::index::Flood<DIM> octree(points);
	// for (int j=0; j<m; ++j) {
		// auto q = boxes[j].min_corner();
		// auto results = octree.knn_query(q, K);
		// double knnd = 0.0;
		// for (auto point : results) {
			// double tmp = bench::common::eu_dist(point, q);
			// if (tmp > knnd)
				// knnd = tmp;
		// }
		// cout << knnd << " ";
		// ret.push_back(knnd);
	// }
	// cout << endl;

	// auto end = std::chrono::steady_clock::now();
	// auto T = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	// cout << "[testOctree] " << T << "ms" << endl;
	
	// return ret;
// }



int main(int argc, char **argv) {
	init();
	auto va = testNaive();
	auto vb = testKDtree();
	
	// auto va = testNaiveKNN();
	// auto vb = testOctreeKNN();

	bool flag = true;
	for (int i=0; i<m; ++i) {
		if (va[i] != vb[i]) {
			flag = false;
			break;
		}
	}
	cout << (flag ? "True" : "False") << endl;

	return 0;
}
