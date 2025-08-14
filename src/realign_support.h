//
// Created by xinwei on 6/19/24.
//

#ifndef FUSION_DETECTION_2_REALIGN_SUPPORT_H
#define FUSION_DETECTION_2_REALIGN_SUPPORT_H

#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include "alignment.h"
#include "overlap.h"
#include "sam.h"
#include "options.h"


extern std::unordered_map<std::string, std::vector<OverlapResultCls>> overlapMap;
extern std::unordered_map<std::string, std::vector<sam_t>> samMap;


void fillOverlapMap(const std::vector<OverlapResultCls>& overlaps);
void fillSamMap(const std::vector<sam_t>& samEntries);
bool isReadSupportingOverlap(const sam_t& read, const OverlapResultCls& overlap, const options_t& options);
bool isSpanningReadPair(const sam_t& read1, const sam_t& read2, const OverlapResultCls& overlap, const options_t& options);

std::unordered_map<std::string, int> countSplitReads(
        const std::unordered_map<std::string, std::vector<OverlapResultCls>>& overlapMap,
        const std::unordered_map<std::string, std::vector<sam_t>>& samMap,
        const options_t& options);
std::unordered_map<std::string, int> countSpanReadPairs(
        const std::unordered_map<std::string, std::vector<OverlapResultCls>>& overlapMap,
        const std::unordered_map<std::string, std::vector<sam_t>>& samMap);

std::unordered_map<std::string, std::vector<sam_t>> collectSplitReads(
    const std::unordered_map<std::string, std::vector<OverlapResultCls>>& overlapMap,
    const std::unordered_map<std::string, std::vector<sam_t>>& samMap,
    const options_t& options);

void writeSplitReadsToFile(
    const std::unordered_map<std::string, std::vector<sam_t>>& splitReadsMap,
    const std::string& outputFilename,
    const options_t& options) ;

std::unordered_set<std::string> filterQueries(
        const std::unordered_map<std::string, int>& splitReadsCount,
        const std::unordered_map<std::string, int>& spanPairsCount);

std::vector<alignment_t> filterAlignmentsByValidBaseQueries(
        const std::vector<alignment_t>& alignments,
        const std::unordered_set<std::string>& validBaseQueries);
std::vector<alignment_t> filterAlignmentsByValidQueries(
        const std::vector<alignment_t>& alignments,
        const std::unordered_set<std::string>& validQueries);

class coordination_t {
public:
    std::string query;
    std::string target;
    int tstart;
    int tend;
    std::string strand;


    coordination_t() : tstart(0), tend(0), strand("") {}



    coordination_t(const std::string& query, const std::string& target, int tstart, int tend, const std::string& strand)
            : query(query), target(target), tstart(tstart), tend(tend), strand(strand) {}

};



std::string extractBaseQueryName(const std::string& query);
std::vector<coordination_t> extractCoordinations(const std::vector<alignment_t>& alignments);
std::vector<coordination_t> mergeAlignments(const std::vector<alignment_t>& alignments);
std::vector<coordination_t> mergeContinuousSegments(const std::vector<coordination_t>& segments);
std::vector<coordination_t> extractSimpleCoordinations(const std::vector<alignment_t>& alignments);

#endif //FUSION_DETECTION_2_REALIGN_SUPPORT_H
