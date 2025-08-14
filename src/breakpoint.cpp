//
// Created by xinwei on 7/17/24.
//

#include "breakpoint.h"
#include "log.h"
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>

#include "utils.h"

// AHC Clustering Function: Groups breakpoints within distance threshold D
std::vector<std::vector<int>> ahcClustering(const std::vector<int>& positions, int D) {
    std::vector<std::vector<int>> clusters;
    if (positions.empty()) return clusters;

    // Sort positions for hierarchical merging
    std::vector<int> sortedPositions = positions;
    std::sort(sortedPositions.begin(), sortedPositions.end());

    // Initialize clusters (each position starts as its own cluster)
    for (int pos : sortedPositions) {
        clusters.push_back({pos});
    }

    // Merge clusters iteratively until distance condition is met
    while (clusters.size() > 1) {
        int minDist = INT_MAX;
        int mergeIdx1 = -1, mergeIdx2 = -1;

        // Find the closest clusters
        for (size_t i = 0; i < clusters.size() - 1; ++i) {
            int dist = std::abs(clusters[i].back() - clusters[i + 1].front());
            if (dist < minDist) {
                minDist = dist;
                mergeIdx1 = i;
                mergeIdx2 = i + 1;
            }
        }

        // Stop merging if the minimum distance exceeds the threshold D
        if (minDist > D) break;

        // Merge clusters
        clusters[mergeIdx1].insert(clusters[mergeIdx1].end(), clusters[mergeIdx2].begin(), clusters[mergeIdx2].end());
        clusters.erase(clusters.begin() + mergeIdx2);
    }

    return clusters;
}

// Select a representative breakpoint from each cluster
std::vector<int> selectBreakpointsFromClusters(
    const std::vector<std::vector<int>>& clusters,
    const std::unordered_map<int, int>& contigScores) {

    std::vector<int> finalBreakpoints;

    for (const auto& cluster : clusters) {
        int bestPos = cluster[0];
        int bestScore = contigScores.at(bestPos);

        // Choose the highest-scoring position in the cluster
        for (int pos : cluster) {
            if (contigScores.at(pos) > bestScore) {
                bestPos = pos;
                bestScore = contigScores.at(pos);
            }
        }

        finalBreakpoints.push_back(bestPos);
    }

    return finalBreakpoints;
}


// Implement existing functions
std::vector<int> findSignificantIncreases(const std::vector<std::pair<int, float>>& coverage, float coverageDiffer) {
    std::vector<int> significantPositions;
    for (size_t i = 1; i < coverage.size(); ++i) {
        float diff = (coverage[i].second - coverage[i - 1].second) / coverage[i - 1].second;
        if (diff >= coverageDiffer) {
            significantPositions.push_back(coverage[i].first);
        }
    }
    return significantPositions;
}

std::vector<int> findSignificantDecreases(const std::vector<std::pair<int, float>>& coverage, float coverageDiffer) {
    std::vector<int> significantPositions;
    for (size_t i = 1; i < coverage.size(); ++i) {
        float diff = (coverage[i - 1].second - coverage[i].second) / coverage[i - 1].second;
        if (diff >= coverageDiffer) {
            significantPositions.push_back(coverage[i].first);
        }
    }
    return significantPositions;
}



std::vector<int> findContinuityChanges(const std::vector<std::pair<int, float>>& coverage, float averageRate, float threshold) {
    std::vector<int> continuityChanges;
    int count = 0;
    for (int i = 0; i < coverage.size(); ++i) {
        if (coverage[i].second < averageRate * (1 - threshold)) {
            count++;
        } else {
            if (count >= 3) {
                continuityChanges.push_back(coverage[i-1].first);
            }
            count = 0;
        }
    }
    return continuityChanges;
}




std::vector<int> findReadLengthAnomalies(const std::vector<std::pair<int, int>>& readLengths, float threshold, const options_t& options) {
    std::vector<int> anomalyPoints;
    int readLength = options.read_length;

    // Simulate detection logic (if the read length is always the same, you can directly return an empty result)
    for (const auto& length : readLengths) {
        if (abs(length.second - readLength) > threshold) {
            anomalyPoints.push_back(length.first);
        }
    }
    return anomalyPoints;
}



// Main breakpoint evaluation function
std::unordered_map<std::string, std::vector<int>> evaluateBreakpoints(
    const std::unordered_map<std::string, std::unordered_map<int, int>>& newPerBaseCoverageMap,
    const std::vector<std::pair<int, int>>& readLengths,
    float averageRate, float qualityThreshold,
    const options_t& options) {

    std::unordered_map<std::string, std::vector<int>> finalBreakpointsMap;

    for (const auto& [contig, coverageMap] : newPerBaseCoverageMap) {
        Logger::Info( get_time_string() + "Evaluating breakpoints for contig: " + contig);

        if (coverageMap.size() == 1) {
            for (const auto& [pos, count] : coverageMap) {
                finalBreakpointsMap[contig].push_back(pos);
                Logger::Debug(get_time_string() + "Single base contig detected at position " +
                              std::to_string(pos) + ", marked as breakpoint\n");
            }
            continue;
        }

        std::vector<std::pair<int, float>> coverage;
        for (const auto& [pos, cov] : coverageMap) {
            coverage.emplace_back(pos, static_cast<float>(cov));
        }
        std::sort(coverage.begin(), coverage.end());

        std::unordered_map<int, int> contigScores;
        for (const auto& [pos, _] : coverageMap) {
            contigScores[pos] = 0;
        }

        // Identify significant changes based on coverage
        auto significantIncreases = findSignificantIncreases(coverage, options.coverage_differ);
        auto significantDecreases = findSignificantDecreases(coverage, options.coverage_differ);

        for (int pos : significantIncreases) contigScores[pos] += 8;
        for (int pos : significantDecreases) contigScores[pos] += 8;

        auto continuityChanges = findContinuityChanges(coverage, averageRate, options.coverage_differ);
        auto readLengthAnomalies = findReadLengthAnomalies(readLengths, 50, options);

        for (int pos : continuityChanges) contigScores[pos] += 2;
        for (int pos : readLengthAnomalies) contigScores[pos] += 2;

        // Extract all potential breakpoint positions
        std::vector<int> candidatePositions;
        for (const auto& [pos, score] : contigScores) {
            if (score > 0) { // Only include positions with evidence
                candidatePositions.push_back(pos);
            }
        }

        // Apply AHC clustering
        int clusteringThreshold = 10;  // Set based on expected breakpoint precision
        auto clusters = ahcClustering(candidatePositions, clusteringThreshold);

        // Select representative breakpoints
        auto finalBreakpoints = selectBreakpointsFromClusters(clusters, contigScores);
        finalBreakpointsMap[contig] = finalBreakpoints;
    }

    return finalBreakpointsMap;
}
