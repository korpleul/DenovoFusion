//
// Created by xinwei on 1/16/25.
//

#ifndef FILTER_EDGE_UNALIGNED_H
#define FILTER_EDGE_UNALIGNED_H

#include <vector>
#include <string>
#include "alignment.h"
#include "output_fusions.h"
#include "options.h"


std::vector<result_t> filter_edge_unaligned(const std::vector<result_t>& final_results,
                                             const std::vector<std::pair<alignment_t, alignment_t>>& pairedAlignments,
                                             const options_t& options);

#endif // FILTER_EDGE_UNALIGNED_H
