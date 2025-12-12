//
// Created by xinwei on 10/1/24.
//

#ifndef FUSION_DETECTION_2_FILTER_LONG_GAP_H
#define FUSION_DETECTION_2_FILTER_LONG_GAP_H


#include <vector>
#include "output_fusions.h"

std::vector<result_t> filter_long_gap(const std::vector<result_t>& results, unsigned int long_gap_threshold, unsigned int short_segment_threshold);

#endif //FUSION_DETECTION_2_FILTER_LONG_GAP_H
