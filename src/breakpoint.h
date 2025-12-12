//
// Created by xinwei on 7/17/24.
//

#ifndef BREAKPOINT_H
#define BREAKPOINT_H

#include <vector>
#include <unordered_map>
#include <utility>
#include <string>
#include "options.h"


// Function declarations for existing factors
std::vector<int> findSignificantIncreases(const std::vector<std::pair<int, float>>& coverage, float coverageDiffer);
std::vector<int> findSignificantDecreases(const std::vector<std::pair<int, float>>& coverage, float coverageDiffer);
std::vector<int> findContinuityChanges(const std::vector<std::pair<int, float>>& coverageRates, float averageRate, float threshold);




std::vector<int> findReadLengthAnomalies(const std::vector<std::pair<int, int>>& readLengths, float threshold, const options_t& options);





// Evaluate breakpoints using all factors
std::unordered_map<std::string, std::vector<int>> evaluateBreakpoints(const std::unordered_map<std::string, std::unordered_map<int, int>>& newPerBaseCoverageMap,
                                                                                  const std::vector<std::pair<int, int>>& readLengths,
                                                                                  float averageRate, float qualityThreshold,
                                                                                  const options_t& options);


#endif // BREAKPOINT_H
