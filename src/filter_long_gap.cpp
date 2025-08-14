//
// Created by xinwei on 10/1/24.
//

#include "filter_long_gap.h"

#include <array>


// Define a filter function and return a new filtered result set
std::vector<result_t> filter_long_gap(const std::vector<result_t>& results,
                                      unsigned int long_gap_threshold,
                                      unsigned int short_segment_threshold) {
    std::vector<result_t> filtered_results;

    for (auto result : results) {
        if (result.chromosome1 == result.chromosome2) {
            std::array<unsigned int, 4> positions = {result.tstart1, result.tend1, result.tstart2, result.tend2};
            std::sort(positions.begin(), positions.end());
            unsigned int gap = positions[2] > positions[1] ? positions[2] - positions[1] : 0;

            unsigned int matching_segment_left = result.tstart1 > result.tend1
                                                 ? result.tstart1 - result.tend1
                                                 : result.tend1 - result.tstart1;

            unsigned int matching_segment_right = result.tstart2 > result.tend2
                                                  ? result.tstart2 - result.tend2
                                                  : result.tend2 - result.tstart2;


            if (gap < long_gap_threshold ||
                matching_segment_left <= short_segment_threshold ||
                matching_segment_right <= short_segment_threshold) {
                result.filter_status = "discarded_long_gap";
                filtered_results.push_back(result);
                continue;
                }
        }

        result.filter_status = "kept";
        filtered_results.push_back(result);
    }

    return filtered_results;
}
