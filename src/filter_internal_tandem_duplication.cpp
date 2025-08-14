//
// Created by xinwei on 7/14/25.
//

#include "filter_internal_tandem_duplication.h"

#include <vector>
#include <cmath>
#include <string>
#include <algorithm>




std::vector<result_t> filter_internal_tandem_duplication(
    const std::vector<result_t>& results,
    unsigned int max_itd_length,
    unsigned int min_split_reads,
    float min_itd_fraction
) {
    std::vector<result_t> filtered_results;

    for (auto res : results) {
        // 1. Must be the same gene
        if (res.gene1 != res.gene2) {
            res.filter_status = "kept";
            filtered_results.push_back(res);
            continue;
        }

        // 2. Direction check
        if (!(res.direction1 == "UPSTREAM" && res.direction2 == "DOWNSTREAM")) {
            res.filter_status = "kept";
            filtered_results.push_back(res);
            continue;
        }

        // 3. Breakpoint distance
        unsigned int breakpoint1 = (res.tstart1 + res.tend1) / 2;
        unsigned int breakpoint2 = (res.tstart2 + res.tend2) / 2;
        unsigned int distance = std::abs((int)breakpoint1 - (int)breakpoint2);
        if (distance > max_itd_length) {
            res.filter_status = "kept";
            filtered_results.push_back(res);
            continue;
        }

        // 4. Support read count and fraction
        if (res.splitReadsCount < min_split_reads ||
            (res.coverage > 0 && (float)res.splitReadsCount / res.coverage < min_itd_fraction)) {
            res.filter_status = "discarded_internal_tandem_duplication";
            filtered_results.push_back(res);
            continue;
            }

        // If passed all checks, keep it
        res.filter_status = "kept";
        filtered_results.push_back(res);
    }

    return filtered_results;
}
