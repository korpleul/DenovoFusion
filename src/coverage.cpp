//
// Created by xinwei on 7/16/24.
//

#include "coverage.h"
#include <regex>
#include <algorithm>

// Calculate aligned length from CIGAR string (for operations that affect reference positions)
int calculateAlignedLength(const std::string& cigar) {
    std::regex cigarRegex(R"((\d+)([MIDNSHP=X]))");
    std::sregex_iterator begin(cigar.begin(), cigar.end(), cigarRegex);
    std::sregex_iterator end;

    int alignedLength = 0;

    for (auto it = begin; it != end; ++it) {
        int length = std::stoi((*it)[1].str());
        char type = (*it)[2].str()[0];

        // Include operations that consume reference positions
        if (type == 'M' || type == 'D' || type == 'N' || type == '=' || type == 'X') {
            alignedLength += length;
        }
    }

    return alignedLength;
}

// Calculate per-base coverage for an overlap region
std::unordered_map<int, int> calculateOverlapCoverage(const std::vector<sam_t>& reads,
                                                      const std::string& contig,
                                                      int overlapStart,
                                                      int overlapEnd) {
    std::unordered_map<int, int> coverageMap;

    for (const auto& read : reads) {
        if (read.rname != contig) {
            continue; // Skip reads that do not belong to this contig
        }

        int start = std::max(read.pos, overlapStart);
        int end = std::min(read.pos + calculateAlignedLength(read.cigar), overlapEnd);

        // Handle case where start == end (single base)
        if (start == end) {
            coverageMap[start]++;
        }
        // Handle general case where start < end (multiple bases)
        else {
            for (int i = start; i < end; ++i) {
                coverageMap[i]++;
            }
        }
    }

    return coverageMap;
}


// Calculate average coverage from the coverage map
float calculateAverageCoverage(const std::unordered_map<int, int>& coverageMap) {
    int totalCoverage = 0;
    int basesCount = coverageMap.size(); // Number of bases with coverage information

    for (const auto& [position, count] : coverageMap) {
        totalCoverage += count;
    }

    return basesCount > 0 ? static_cast<float>(totalCoverage) / basesCount : 0.0f;
}

// Process coverage for all overlaps and store results in maps
void processCoverage(const std::vector<sam_t>& reads,
                     const std::vector<OverlapResultCls>& overlaps,
                     std::unordered_map<std::string, float>& averageCoverageMap,
                     std::unordered_map<std::string, std::unordered_map<int, int>>& perBaseCoverageMap) {
    for (const auto& overlap : overlaps) {
        std::string contig = overlap.query_id_;
        int start = overlap.getStart();
        int end = overlap.getEnd();

        // Calculate per-base coverage for the contig's overlap area
        auto overlapCoverage = calculateOverlapCoverage(reads, contig, start, end);

        // Store the per-base coverage for the contig
        perBaseCoverageMap[contig] = overlapCoverage;

        // Calculate average coverage for the contig's overlap area
        float average = calculateAverageCoverage(overlapCoverage);
        averageCoverageMap[contig] = average;
    }
}





void savePerBaseCoverageToTSV(
    const std::unordered_map<std::string, std::unordered_map<int, int>>& coverageMap,
    const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    // Write the header
    outFile << "Contig\tPosition\tCoverage\n";

    // Write the data
    for (const auto& [contig, positions] : coverageMap) {
        for (const auto& [position, count] : positions) {
            outFile << contig << '\t' << position << '\t' << count << '\n';
        }
    }

    outFile.close();
}
