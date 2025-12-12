//
// Created by xinwei on 10/1/24.
//

#include "filter_min_support.h"



std::vector<result_t> filter_min_support(std::vector<result_t>& results, unsigned int min_span_reads, unsigned int min_split_reads) {
    std::vector<result_t> kept_results;      // Retained results

    for (auto& result : results) { // Here result is a non-const reference
        if (result.spanReadsCount >= min_span_reads && result.splitReadsCount >= min_split_reads) {
            result.filter_status = "kept";
            kept_results.push_back(result);
        } else {
            result.filter_status = "discarded_min_support";
            kept_results.push_back(result);
        }
    }

    return kept_results; // return retaineed results
}


