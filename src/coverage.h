//
// Created by xinwei on 7/16/24.
//

#ifndef FUSION_DETECTION_2_COVERAGE_H
#define FUSION_DETECTION_2_COVERAGE_H

#include <vector>
#include <unordered_map>
#include <string>
#include "sam.h"
#include "overlap.h"


std::unordered_map<std::string, std::unordered_map<int, int>> calculateCoverage(const std::vector<sam_t>& reads);


std::unordered_map<int, int> calculateOverlapCoverage(
        const std::unordered_map<int, int>& coverageMap,
        const OverlapResultCls& overlap);



// Function to calculate aligned length from CIGAR string
int calculateAlignedLength(const std::string& cigar);

// Function to calculate per-base coverage within an overlap area for a contig
std::unordered_map<int, int> calculateOverlapCoverage(const std::vector<sam_t>& reads,
                                                      const std::string& contig,
                                                      int overlapStart,
                                                      int overlapEnd);

// Function to calculate average coverage from per-base coverage
float calculateAverageCoverage(const std::unordered_map<int, int>& coverageMap);

// Function to process coverage for all overlaps and reads
void processCoverage(const std::vector<sam_t>& reads,
                     const std::vector<OverlapResultCls>& overlaps,
                     std::unordered_map<std::string, float>& averageCoverageMap,
                     std::unordered_map<std::string, std::unordered_map<int, int>>& perBaseCoverageMap);

// Function to save per-base coverage to TSV file
void savePerBaseCoverageToTSV(
    const std::unordered_map<std::string, std::unordered_map<int, int>>& coverageMap,
    const std::string& filename);

#endif //FUSION_DETECTION_2_COVERAGE_H
