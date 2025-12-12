//
// Created by xinwei on 7/14/25.
//
#include "output_fusions.h"

#ifndef FILTER_INTERNAL_TANDEM_DUPLICATION_H
#define FILTER_INTERNAL_TANDEM_DUPLICATION_H

std::vector<result_t> filter_internal_tandem_duplication(
    const std::vector<result_t>& results,
    unsigned int max_itd_length,
    unsigned int min_split_reads,
    float min_itd_fraction
);

#endif //FILTER_INTERNAL_TANDEM_DUPLICATION_H
