//
// Created by xinwei on 5/17/24.
//

#include "overlap.h"
#include "coverage.h"

#include <iostream>
#include <sstream>
#include <algorithm>



OverlapResultCls::OverlapResultCls(const std::vector<std::pair<int, int>>& positions, const std::string& baseQueryName)
        : query_id_(baseQueryName) {
    std::vector<int> adjusted_positions;
    for (const auto& pos : positions) {
        adjusted_positions.push_back(pos.first + (pos.second - 1) * 1000);
    }
    std::sort(adjusted_positions.begin(), adjusted_positions.end());


    // Compute overlap range if enough data points are available
    if (adjusted_positions.size() >= 4) {
        start_ = adjusted_positions[1];
        end_ = adjusted_positions[2];
        overlap_interval_ = end_ - start_ + 1;
    } else {
        // Fallback values for insufficient data
        start_ = end_ = overlap_interval_ = 0;
    }
}


std::string OverlapResultCls::toString() const {
    std::ostringstream oss;
    oss << query_id_
        << ": " << start_ << "-" << end_
        << " (" << overlap_interval_ << "), Contig: ["
        << contig_start_ << "-" << contig_end_ << "]";
    return oss.str();
}

int OverlapResultCls::getOverlapInterval() const {
    return overlap_interval_;
}

int OverlapResultCls::getStart() const {
    return start_;
}

int OverlapResultCls::getEnd() const {
    return end_;
}

int OverlapResultCls::getContigStart() const {
    return contig_start_;
}

int OverlapResultCls::getContigEnd() const {
    return contig_end_;
}


std::vector<std::pair<alignment_t, alignment_t>> pairAlignments(const std::vector<alignment_t>& alignments) {
    std::map<std::string, std::vector<alignment_t>> groupedAlignments;
    for (const auto& alignment : alignments) {
        groupedAlignments[alignment.query].push_back(alignment);
    }

    std::vector<std::pair<alignment_t, alignment_t>> pairedAlignments;
    for (const auto& entry : groupedAlignments) {
        if (entry.second.size() == 2) {
            pairedAlignments.push_back(std::make_pair(entry.second[0], entry.second[1]));
        }
    }
    return pairedAlignments;
}



std::vector<OverlapResultCls> processAndMergeAlignments(
    const std::vector<std::pair<alignment_t, alignment_t>>& pairedAlignments,
    const std::unordered_map<std::string, std::string>& mergedsequences)
{
    std::vector<OverlapResultCls> overlapResults;

    for (const auto& pair : pairedAlignments) {

        bool is_single = false;

        if (pair.second.query.empty()) {
            is_single = true;
        }

        if (is_single) {

            std::vector<std::pair<int, int>> positions = {
                {pair.first.qstart, extractNumber(pair.first.query)},
                {pair.first.qend, extractNumber(pair.first.query)}
            };

            OverlapResultCls result(positions, baseQueryName(pair.first.query));

            std::string contig_name = baseQueryName(pair.first.query);
            if (mergedsequences.find(contig_name) != mergedsequences.end()) {
                int contig_length = static_cast<int>(mergedsequences.at(contig_name).length());
                result.contig_start_ = 0;
                result.contig_end_ = contig_length - 1;
            }

            overlapResults.push_back(std::move(result));
        } else {

            std::vector<std::pair<int, int>> positions = {
                {pair.first.qstart, extractNumber(pair.first.query)},
                {pair.first.qend, extractNumber(pair.first.query)},
                {pair.second.qstart, extractNumber(pair.second.query)},
                {pair.second.qend, extractNumber(pair.second.query)}
            };

            OverlapResultCls result(positions, baseQueryName(pair.first.query));

            std::string contig_name = baseQueryName(pair.first.query);
            if (mergedsequences.find(contig_name) != mergedsequences.end()) {
                int contig_length = static_cast<int>(mergedsequences.at(contig_name).length());
                result.contig_start_ = 0;
                result.contig_end_ = contig_length - 1;
            }

            overlapResults.push_back(std::move(result));
        }
    }

    return overlapResults;
}




// Process alignments into OverlapResultCls objects, handling both overlaps and gaps
std::vector<OverlapResultCls> processAlignments(const std::vector<std::pair<alignment_t, alignment_t>>& pairedAlignments) {
    std::vector<OverlapResultCls> overlapResults;

    for (const auto& pair : pairedAlignments) {
        // Calculate the start and end points of the overlap or gap region
        int overlapStart = std::max(pair.first.qstart, pair.second.qstart);
        int overlapEnd = std::min(pair.first.qend, pair.second.qend);

        // Initialize variables for gap detection
        bool isGap = false;
        int gapStart = 0, gapEnd = 0;

        if (overlapStart <= overlapEnd) {
            // There is an overlap
            std::vector<std::pair<int, int>> positions = {
                {overlapStart, overlapEnd}
            };

            // Create an OverlapResultCls object for overlap
            OverlapResultCls result(positions, pair.first.query);
            result.start_ = overlapStart;
            result.end_ = overlapEnd;
            result.overlap_interval_ = overlapEnd - overlapStart;
            result.contig_start_ = std::min(pair.first.qstart, pair.second.qstart);
            result.contig_end_ = std::max(pair.first.qend, pair.second.qend);

            overlapResults.push_back(result);
        } else {
            // There is a gap
            isGap = true;
            gapStart = pair.first.qend;
            gapEnd = pair.second.qstart;

            std::vector<std::pair<int, int>> gapPositions = {
                {gapStart, gapEnd}
            };

            // Create an OverlapResultCls object for gap
            OverlapResultCls result(gapPositions, pair.first.query);
            result.start_ = gapStart;
            result.end_ = gapEnd;
            result.overlap_interval_ = gapEnd - gapStart; // Negative or 0 means no overlap
            result.contig_start_ = std::min(pair.first.qstart, pair.second.qstart);
            result.contig_end_ = std::max(pair.first.qend, pair.second.qend);

            overlapResults.push_back(result);
        }
    }

    return overlapResults;
}


// Helper function to extract number from query
int extractNumber(const std::string& query) {
    size_t pos = query.rfind('_');
    if (pos != std::string::npos) {
        std::string last_part = query.substr(pos + 1);
        if (last_part.size() <= 3 &&
            std::all_of(last_part.begin(), last_part.end(), ::isdigit)) {
            return std::stoi(last_part);
            }
    }
    return 1;
}


// Extract base query name
std::string baseQueryName(const std::string& query) {
    size_t pos = query.rfind('_');
    if (pos != std::string::npos) {
        std::string last_part = query.substr(pos + 1);
        if (last_part.size() <= 3 &&
            std::all_of(last_part.begin(), last_part.end(), ::isdigit)) {
            return query.substr(0, pos);
            }
    }
    return query;
}
