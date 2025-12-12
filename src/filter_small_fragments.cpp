//
// Created by xinwei on 1/6/25.
//

#include "filter_small_fragments.h"

#include <vector>
#include <algorithm>
#include <iostream>

std::vector<result_t> filter_small_fragments(const std::vector<result_t>& results, const options_t& options) {
    std::vector<result_t> filtered_results;
    float size_ratio_threshold = options.size_ratio_threshold;

    for (const auto& result : results) {
        int size1 = result.tend1 - result.tstart1;
        int size2 = result.tend2 - result.tstart2;

        int smaller_size = std::min(size1, size2);
        int larger_size = std::max(size1, size2);

        if (static_cast<float>(smaller_size) / larger_size >= size_ratio_threshold) {
            result_t kept = result;
            kept.filter_status = "kept";
            filtered_results.push_back(kept);
        } else {
            result_t discarded = result;
            discarded.filter_status = "discarded_small_fragments";
            filtered_results.push_back(discarded);
        }
    }

    return filtered_results;
}


