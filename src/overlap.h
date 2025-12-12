//
// Created by xinwei on 5/17/24.
//

#ifndef OVERLAP_H
#define OVERLAP_H

#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include "alignment.h"
#include "sam.h"
#include "candidate_group.h"



class OverlapResultCls {
public:
    OverlapResultCls(const std::vector<std::pair<int, int>>& positions, const std::string& baseQueryName);

    std::string toString() const;
    int getOverlapInterval() const;
    int getStart() const;
    int getEnd() const;
    int getContigStart() const;
    int getContigEnd() const;

    std::unordered_map<int, int> calculateOverlapCoverage(const std::vector<sam_t>& reads) const;


    friend std::ostream& operator<<(std::ostream &os, const OverlapResultCls &result);


    std::string query_id_;
    int start_, end_, overlap_interval_;
    int contig_start_, contig_end_;
};

int extractNumber(const std::string& query);
std::string baseQueryName(const std::string& query);

std::vector<std::pair<alignment_t, alignment_t>> pairAlignments(const std::vector<alignment_t>& alignments);
std::vector<OverlapResultCls> processAndMergeAlignments(
    const std::vector<std::pair<alignment_t, alignment_t>>& pairedAlignments,
    const std::unordered_map<std::string, std::string>& mergedsequences);
std::vector<OverlapResultCls> processAlignments(const std::vector<std::pair<alignment_t, alignment_t>>& pairedAlignments);
#endif // OVERLAP_H
