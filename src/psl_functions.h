//
// Created by xinwei on 6/6/24.
//

#ifndef FUSION_DETECTION_2_PSL_FUNCTIONS_H
#define FUSION_DETECTION_2_PSL_FUNCTIONS_H

#include <sstream>
#include <iterator>
#include <string>
#include <algorithm>  // For std::min, std::max, and std::sort
#include <tuple>      // For std::tuple
#include "psl.h"
#include "coord_pair.h"


class AlignBlockCls {
public:
    AlignBlockCls(const std::pair<int, int>& block_coords_tuple);
    int Gap(const AlignBlockCls& other) const;

private:
    int start;
    int end;
    char strand;
    int span;
};

void CalcOverlap(int startA, int endA, int startB, int endB);
void CalcAlignOverlap(const psl_t& align1, const psl_t& align2);

void CalcGap(int startA, int endA, int startB, int endB);
void CalcAlignGap(const psl_t& align1, const psl_t& align2);

std::string FormatChromosomeNameInPsl(const std::string& psl_str, bool use_chr, bool use_mt);
void CompareCoordPairsInGroups(std::unordered_map<std::string, std::vector<CoordPair>>& coord_pairs_q);



#endif //FUSION_DETECTION_2_PSL_FUNCTIONS_H
