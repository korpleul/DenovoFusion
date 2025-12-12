//
// Created by xinwei on 7/23/25.
//

#include "filter_mt.h"

std::vector<result_t> filter_mt(std::vector<result_t>& results) {
    for (auto& result : results) {
        if (result.filter_status == "kept") {
            if (result.chromosome1 == "MT" || result.chromosome2 == "MT" ||
                result.chromosome1 == "chrM" || result.chromosome2 == "chrM") {
                result.filter_status = "discarded_mt";
                }
        }
    }
    return results;
}
