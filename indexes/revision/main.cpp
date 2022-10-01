#include <bits/stdc++.h>
using namespace std;

#include "equal_depth_gridBox.hpp"
#include "fullscanBox.hpp"
#include "../../utils/type.hpp"
#include "../../utils/common.hpp"

const int MAXL = 128;
const int DIM = 2;
const int n = 100;
const int m = 10;
const int K = 10;
std::vector<box_t<DIM> > inputs;
std::vector<box_t<DIM> > boxes;

template<size_t dim>
static box_t<dim> sample_range(const int MAXL) {
	using Point = point_t<dim>;
	
    // seed the generator
    std::mt19937 gen(std::random_device{}()); 
    std::uniform_int_distribution<> uint_dist(0, MAXL);

	Point p;
	for (int j=0; j<dim; ++j) {
		p[j] = uint_dist(gen);
	}
	p.fill(0);
	
	Point q = p;
	std::uniform_int_distribution<> uint_radius((int)sqrt(MAXL), MAXL);
	for (int j=0; j<dim; ++j) {
		q[j] += uint_radius(gen);;
	}
	
    return box_t<dim>(p, q);
}

template<size_t dim>
static vector<box_t<DIM>> sample_boxes(const int n, const int MAXL) {
	using Point = point_t<dim>;
	std::vector<box_t<DIM>> boxes;
	
    // seed the generator
    std::mt19937 gen(std::random_device{}()); 
    std::uniform_int_distribution<> uint_dist(0, MAXL);
	
	for (int i=0; i<n; ++i) {
		Point p;
		for (int j=0; j<dim; ++j) {
			p[j] = uint_dist(gen);
		}
		
		Point q = p;
		std::uniform_int_distribution<> uint_radius((int)sqrt(MAXL), MAXL);
		for (int j=0; j<dim; ++j) {
			q[j] += uint_radius(gen);;
		}
		
		boxes.emplace_back(box_t<dim>(p, q));
	}
	
    return boxes;
}

void init() {
	inputs = sample_boxes<DIM>(n, 1024);
	for (int i=0; i<m; ++i) {
		box_t<DIM> box = sample_range<DIM>(1024);
		boxes.emplace_back(box);
	}
}

// Query using brute-force
vector<int> testNaive() {
	auto start = std::chrono::steady_clock::now();
	vector<int> ret;
	
	bench::index::FullScanBox<DIM> fs(inputs);
	for (int j=0; j<m; ++j) {
		auto results = fs.range_query(boxes[j]);
		cout << results.size() << " ";
		ret.push_back((int) results.size());
	}
	cout << endl;
	
	auto end = std::chrono::steady_clock::now();
	auto T = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	cout << "[testNaive] " << T << "ms" << endl;
	
	return ret;
}

// Query using index
vector<int> testMyIndex() {
	auto start = std::chrono::steady_clock::now();
	vector<int> ret;
	
	bench::index::EDG<DIM,10> myIndex(inputs);
	for (int j=0; j<m; ++j) {
		auto results = myIndex.range_query(boxes[j]);
		cout << results.size() << " ";
		ret.push_back((int) results.size());
	}
	cout << endl;

	auto end = std::chrono::steady_clock::now();
	auto T = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	cout << "[testMyIndex] " << T << "ms" << endl;
	
	return ret;
}



int main(int argc, char **argv) {
	init();
	auto va = testNaive();
	auto vb = testMyIndex();

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
