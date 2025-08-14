//
// Created by xinwei on 9/29/24.
//

#include "filter_duplicates.h"

// Define a function to compare whether two result_t are duplicates
bool is_duplicate(const result_t& res1, const result_t& res2) {
    return res1.contig == res2.contig &&
           res1.gene1 == res2.gene1 &&
           res1.gene2 == res2.gene2 &&
           res1.tstart1 == res2.tstart1 &&
           res1.tstart2 == res2.tstart2 &&
           res1.tstrand1 == res2.tstrand1 &&
           res1.tstrand2 == res2.tstrand2 &&
           res1.chromosome1 == res2.chromosome1 &&
           res1.chromosome2 == res2.chromosome2 &&
           res1.direction1 == res2.direction1 &&
           res1.direction2 == res2.direction2;
}


// Define a function to filter out duplicate results
std::vector<result_t> filter_duplicates(std::vector<result_t>& results) {
    std::vector<result_t> unique_results;

    for (auto& result : results) {
        bool is_dup = false;
        for (const auto& seen : unique_results) {
            if (is_duplicate(result, seen)) {
                is_dup = true;
                break;
            }
        }
        if (is_dup) {
            result.filter_status = "discarded_duplicates";
            unique_results.push_back(result);
        } else {
            result.filter_status = "kept";
            unique_results.push_back(result);
        }
    }

    return unique_results;
}
