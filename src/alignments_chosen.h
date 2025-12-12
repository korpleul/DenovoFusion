//
// Created by xinwei on 6/7/24.
//

#ifndef FUSION_DETECTION_2_ALIGNMENTS_CHOSEN_H
#define FUSION_DETECTION_2_ALIGNMENTS_CHOSEN_H

#include <iostream>
#include <vector>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <mutex>

#include "psl.h"
#include "alignment.h"
#include "options.h"


struct alignment_score_t {
    std::string qname;
    std::vector<std::pair<int, int>> combination;
    double score;
    alignment_t alignment_details;

    alignment_score_t(const std::string& qname, const std::vector<std::pair<int, int>>& combination, double score, const alignment_t& alignment_details)
            : qname(qname), combination(combination), score(score), alignment_details(alignment_details) {}
};




void GetAlignLengthAndFraction(const psl_t& psls, int& align_len, double& align_fract);

void generateCombinations(const std::vector<std::pair<int, int>>& pairs, int start, int k, std::vector<std::pair<int, int>>& current, std::vector<std::vector<std::pair<int, int>>>& allCombinations);

// Calculate overlap score
int calculateOverlapScore(const std::vector<std::pair<int, int>>& combination);

// Calculate inclusion score
int calculateInclusionScore(const std::vector<std::pair<int, int>>& combination);

std::vector<alignment_t> calculate_alignments_score(const std::unordered_map<std::string, std::vector<alignment_t>>& data_by_qname, const options_t& options);


// Parse psl file
void psl_parse(const std::string& filename, std::vector<psl_t>& psls);



#endif //FUSION_DETECTION_2_ALIGNMENTS_CHOSEN_H
