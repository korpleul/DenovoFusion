//
// Created by xinwei on 10/1/24.
//

#ifndef FUSION_DETECTION_2_FILTER_MIN_SUPPORT_H
#define FUSION_DETECTION_2_FILTER_MIN_SUPPORT_H
#include <string>
#include <vector>
#include "output_fusions.h"

std::vector<result_t> filter_min_support(std::vector<result_t>& results, unsigned int min_span_reads, unsigned int min_split_reads);

#endif //FUSION_DETECTION_2_FILTER_MIN_SUPPORT_H
