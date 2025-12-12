//
// Created by xinwei on 1/16/25.
//

#include "filter_edge_unaligned.h"
#include <iostream>
#include <algorithm>


// filter_edge_unaligned
std::vector<result_t> filter_edge_unaligned(const std::vector<result_t>& final_results,
                                             const std::vector<std::pair<alignment_t, alignment_t>>& pairedAlignments,
                                             const options_t& options) {
    std::vector<result_t> filteredResults;


    for (const auto& result : final_results) {
        bool is_unaligned_at_edge = false;


        for (const auto& pair : pairedAlignments) {
            const alignment_t& alignment1 = pair.first;
            const alignment_t& alignment2 = pair.second;


            if (alignment1.query.find(result.contig) != std::string::npos &&
                alignment2.query.find(result.contig) != std::string::npos) {


                int edge_start = std::min(alignment1.qstart, alignment2.qstart);
                int edge_end = std::max(alignment1.qend, alignment2.qend);
                int contig_length = alignment1.query_len;


                if (edge_start > options.edge_unaligned || contig_length - edge_end > options.edge_unaligned) {
                    is_unaligned_at_edge = true;
                    break;
                }
                }
        }


        if (!is_unaligned_at_edge) {
            filteredResults.push_back(result);
        }
    }

    return filteredResults;
}
