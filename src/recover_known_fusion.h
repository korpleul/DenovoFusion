//
// Created by xinwei on 10/1/24.
//

#ifndef FUSION_DETECTION_2_RECOVER_KNOWN_FUSION_H
#define FUSION_DETECTION_2_RECOVER_KNOWN_FUSION_H
#include <string>
#include "output_fusions.h"

// Coordinate structure to hold the parsed coordinates for each gene
struct coordinate_t {
    std::string chromosome;
    int start;
    int end;
    char strand;
};

// Known fusion structure with vectors to hold multiple coordinates
struct known_fusion_t {
    std::string gene1;
    std::string gene2;
    std::vector<coordinate_t> coordinates1;
    std::vector<coordinate_t> coordinates2;
    std::string source;
};

bool parse_coordinates(const std::string& coordinate_str, std::string& chrom, int& start, int& end);

std::vector<known_fusion_t> load_known_fusions(const std::string& filename);

std::vector<result_t> recover_fusions(const std::vector<result_t>& discarded_results, const std::vector<known_fusion_t>& known_fusions);



#endif //FUSION_DETECTION_2_RECOVER_KNOWN_FUSION_H
