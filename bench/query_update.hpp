#pragma once

#include <cstddef>
#include <random>
#include <map>
#include <algorithm>

#include "../utils/type.hpp"


namespace bench { namespace query {

// sample queries from data
// sample point queries
template<size_t dim>
static std::vector<point_t<dim>> sample_point_queries(vec_of_point_t<dim>& points, size_t n, size_t s=20) {
	assert(points.size() >= n);
	
    // seed the generator
    std::mt19937 gen(0); 
    std::uniform_int_distribution<> uint_dist(0, n-1);

    // generate random indices
    std::vector<point_t<dim>> samples;
    samples.reserve(s);

    for (size_t i=0; i<s; ++i) {
        auto idx = uint_dist(gen);
        samples.emplace_back(points[idx]);
    }

    return samples;
}


// sample knn queries
// k is in range {1, 10, 100, 500, 1000, 10000}
// each k has s=100 samples
template<size_t dim>
static std::map<size_t, vec_of_point_t<dim>> sample_knn_queries(vec_of_point_t<dim>& points, size_t n, size_t s=20) {
    std::vector<int> ks({1, 10, 100, 500, 1000});
    auto query_points = sample_point_queries(points, n, s);
	std::mt19937 gen(0); 
    std::uniform_int_distribution<> uint_dist(0, ks.size()-1);
	std::map<size_t, vec_of_point_t<dim>> knn_queries;

	for (auto& point : query_points) {
		auto k = ks[uint_dist(gen)];
		knn_queries[k].emplace_back(point);
	}

    return knn_queries;
}

template<size_t dim>
static std::pair<point_t<dim>, point_t<dim>> min_and_max(std::vector<point_t<dim>>& points, size_t n) {
    point_t<dim> mins;
    point_t<dim> maxs;

    // boundaries of each dimension
    std::fill(mins.begin(), mins.end(), std::numeric_limits<double>::max());
    std::fill(maxs.begin(), maxs.end(), std::numeric_limits<double>::min());

    for (size_t i=0; i<dim; ++i) {
        for (size_t j=0; j<n; ++j) {
			point_t<dim>& p = points[j];
            mins[i] = std::min(p[i], mins[i]);
            maxs[i] = std::max(p[i], maxs[i]);
        }
    }

    return std::make_pair(mins, maxs);
}


// sample range queries
// selectivity = range_count(q_box) / N
// for each selectivity we generate s=10 random boxes roughly match the selectivity
template<size_t dim>
static std::vector<box_t<dim>> sample_range_queries(vec_of_point_t<dim>& points, size_t n, size_t s=10) {
    std::vector<double> selectivities({0.001, 0.01, 0.05, 0.1, 0.2});
	std::mt19937 gen(0); 
    std::uniform_int_distribution<> uint_dist(0, selectivities.size()-1);
	
    auto corner_points = sample_point_queries(points, n, s);
    
    std::pair<point_t<dim>, point_t<dim>> min_max = min_and_max(points, n);

    std::vector<box_t<dim>> range_queries;
    range_queries.reserve(s);
    
	for (auto& point : corner_points) {
		point_t<dim> another_corner;
		auto sid = uint_dist(gen);
		for (size_t d=0; d<dim; ++d) {
			// this is a roughly estimation of selectivity by assuming 
			// uniform distribution and dimension independence
			double step = (min_max.second[d] - min_max.first[d]) * std::pow(selectivities[sid], 1.0/dim);
			// make sure the generated box is within the data range
			another_corner[d] = std::min(point[d] + step, min_max.second[d]);
		}
		box_t<dim> box(point, another_corner);
		range_queries.emplace_back(box);
	}
    
    return range_queries;
}

// generate update workload
// ret[i] = <update_type, update_point_id>: -1: delete, 1: insert, 0: test range/knn query
// Input: m: number of points have been inserted, n: total number of points, s: number of test rounds
std::vector<std::pair<int,size_t>> sample_update_queries(const size_t m, const size_t n, const size_t s=10) {
	assert(m+s <= n);
	const size_t delta = (n-m)/s;
	size_t cur = m, next_cur;
	std::vector<std::pair<int,size_t>> update_queries;
	
	for (int i=0; i<s; ++i,cur=next_cur) {
		next_cur = (i==s-1) ? n : (cur+delta);
		// std::cout << cur << " " << next_cur << " " << delta << std::endl;
		// seed the generator
		std::mt19937 gen(i); 
		std::uniform_int_distribution<> uint_type(0, 1);
		std::uniform_int_distribution<> uint_delete(0, cur-1);
		std::vector<std::pair<int,size_t>> vtmp;
		vtmp.reserve(delta+(next_cur-cur)+1);
		size_t delete_n = delta, insert_n = next_cur-cur, insert_id = cur, delete_id = 0;
		
		for (int j=delete_n+insert_n-1; j>=0; --j) {
			auto update_type = uint_type(gen);
			size_t delete_id = uint_delete(gen);
			if (update_type>0 && delete_n>0) {
				--delete_n;
				vtmp.emplace_back(std::make_pair(-1, delete_id));
			} else if (update_type==0 && insert_n>0) {
				--insert_n;
				assert(insert_id < next_cur);
				vtmp.emplace_back(std::make_pair(1, insert_id++));
			} else if (delete_n > 0) {
				--delete_n;
				vtmp.emplace_back(std::make_pair(-1, delete_id));
			} else {
				assert(insert_n > 0);
				--insert_n;
				assert(insert_id < next_cur);
				vtmp.emplace_back(std::make_pair(1, insert_id++));
			}
		}
		vtmp.emplace_back(std::make_pair(0, (size_t)0));
		
		update_queries.insert(update_queries.end(), vtmp.begin(), vtmp.end());
	}
    
    return update_queries;
}


template<class Index, size_t Dim>
static void batch_knn_queries(Index& index, std::map<size_t, vec_of_point_t<Dim>>& knn_queries) {
    
	for (auto iter=knn_queries.begin(); iter!=knn_queries.end(); ++iter) {
		size_t k = iter->first;
		vec_of_point_t<Dim> points = iter->second;
		// std::cout << "k = " << k << " points.size() = " << points.size() << std::endl;
        for (auto q_point : points) {
            index.knn_query(q_point, k);
        }
    }
	
	// std::cout << "[FINISH] knn queries" << std::endl;
}


template<class Index, size_t Dim>
static void batch_range_queries(Index& index, std::vector<box_t<Dim>>& range_queries) {
    
    for (auto& box : range_queries) {
        index.range_query(box);
    }
	
	// std::cout << "[FINISH] range queries" << std::endl;
}



template<class Index, size_t Dim>
static void batch_update_queries(Index& index, vec_of_point_t<Dim>& points, vec_of_point_t<Dim>& points_backup, std::vector<std::pair<int,size_t>>& update_queries) {
    size_t n = 0;
	size_t sz = update_queries.size();
	size_t i = 0;
	
	while (i < sz) {
		size_t j = i;
		while (j<sz && update_queries[j].first!=0) {
			int update_type = update_queries[j].first;
			size_t update_id = update_queries[i].second;
			if (update_type > 0) {
				points_backup.emplace_back(points[update_id]);
			}
			++j;
		}
		
		while (i<sz && i<=j) {
			int update_type = update_queries[i].first;
			size_t update_id = update_queries[i].second;
			if (update_type > 0) {
				index.insert(points[update_id]);
				n = std::max(n, update_id+1);
			} else if (update_type < 0) {
				assert(update_id < points_backup.size());
				index.erase(points[update_id]);
			} else {
				std::vector<box_t<Dim>> range_queries = sample_range_queries<Dim>(points, n, 20);
				std::map<size_t, vec_of_point_t<Dim>> knn_queries = sample_knn_queries<Dim>(points, n, 20);
				batch_range_queries(index, range_queries);
				batch_knn_queries(index, knn_queries);
			}
			++i;
		}
	}
	
	std::cout << std::fixed << std::setprecision(8) << "Insert Avg. Time: " << index.get_avg_insert_time() << " [us] " << index.get_insert_time() << " " << index.get_insert_count() << std::endl;
	std::cout << "Erase Avg. Time: " << index.get_avg_erase_time() << " [us] " << index.get_erase_count() << std::endl;
	std::cout << "Range Queries Avg. Time: " << index.get_avg_range_time() << " [us] " << index.get_range_count() << std::endl;
	std::cout << "Knn Queries Avg. Time: " << index.get_avg_knn_time() << " [us] " << index.get_knn_count() << std::endl;
}

template<class Index, size_t Dim>
static void batch_update_queries_byID(Index& index, vec_of_point_t<Dim>& points, vec_of_point_t<Dim>& points_backup, std::vector<std::pair<int,size_t>>& update_queries) {
    size_t n = 0;
	size_t sz = update_queries.size();
	size_t i = 0;
	
	while (i < sz) {
		size_t j = i;
		while (j<sz && update_queries[j].first!=0) {
			++j;
		}
		
		while (i<sz && i<=j) {
			int update_type = update_queries[i].first;
			size_t update_id = update_queries[i].second;
			if (update_type > 0) {
				// std::cout << "[INSERT]" << std::endl;
				points_backup.emplace_back(points[update_id]);
				index.insert(points_backup.size()-1);
				n = std::max(n, update_id+1);
			} else if (update_type < 0) {
				// std::cout << "[ERASE]" << std::endl;
				assert(update_id < points_backup.size());
				index.erase(update_id);
			} else {
				std::vector<box_t<Dim>> range_queries = sample_range_queries<Dim>(points, n, 20);
				std::map<size_t, vec_of_point_t<Dim>> knn_queries = sample_knn_queries<Dim>(points, n, 20);
				// std::cout << "[RANGE QUERY]" << std::endl;
				batch_range_queries(index, range_queries);
				// std::cout << "[KNN QUERY]" << std::endl;
				batch_knn_queries(index, knn_queries);
			}
			++i;
		}
	}
	
	std::cout << "Insert Avg. Time: " << index.get_avg_insert_time() << " [us] " << index.get_insert_count() << std::endl;
	std::cout << "Erase Avg. Time: " << index.get_avg_erase_time() << " [us] " << index.get_erase_count() << std::endl;
	std::cout << "Range Queries Avg. Time: " << index.get_avg_range_time() << " [us] " << index.get_range_count() << std::endl;
	std::cout << "Knn Queries Avg. Time: " << index.get_avg_knn_time() << " [us] " << index.get_knn_count() << std::endl;
}

}
}


