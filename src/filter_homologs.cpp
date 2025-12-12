//
// Created by xinwei on 1/15/25.
//

#include "filter_homologs.h"
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <iostream>

// Constructor
FilterHomologs::FilterHomologs(const std::vector<result_t>& final_results,
                               const std::vector<std::pair<alignment_t, alignment_t>>& pairedAlignments,
                               const std::vector<alignment_t>& alignments)
    : final_results_(final_results), pairedAlignments_(pairedAlignments), alignments_(alignments) {}

// Function to apply the homolog filter
std::vector<result_t> FilterHomologs::filter_homologs() {
    std::vector<result_t> filteredResults;

    // Iterate through each result in final_results
    for (auto result : final_results_) {
        bool is_repeated = false;

        // Extract contig name from result
        std::string contig_name = result.contig;

        // Get matching alignments from pairedAlignments
        std::vector<std::pair<alignment_t, alignment_t>> matchingPairedAlignments;
        for (const auto& pair : pairedAlignments_) {
            const alignment_t& alignment1 = pair.first;
            const alignment_t& alignment2 = pair.second;

            // Match contig name
            if (alignment1.query.find(contig_name + "_") == 0 &&
                alignment2.query.find(contig_name + "_") == 0) {
                matchingPairedAlignments.push_back(pair);
            }
        }

        // Get matching alignments from alignments
        std::vector<alignment_t> matchingAlignments;
        for (const auto& alignment : alignments_) {
            if (alignment.query.find(contig_name + "_") == 0) {
                matchingAlignments.push_back(alignment);
            }
        }

        // Iterate through matching paired alignments and check for repeated coordinates
        for (const auto& pair : matchingPairedAlignments) {
            const alignment_t& alignment1 = pair.first;
            const alignment_t& alignment2 = pair.second;

            // Check if (qstart, qend) of alignment1 or alignment2 is repeated in matchingAlignments
            int count1 = countCoordinatePairs(alignment1.qstart, alignment1.qend, matchingAlignments);
            int count2 = countCoordinatePairs(alignment2.qstart, alignment2.qend, matchingAlignments);

            // If any coordinate pair is repeated, mark as repeated and break
            if (count1 > 1 || count2 > 1) {
                is_repeated = true;
                break;
            }
        }

        // Set filter status
        if (!is_repeated) {
            result.filter_status = "kept";
            filteredResults.push_back(result);
        } else {
            result.filter_status = "discarded_homologs";
            filteredResults.push_back(result);
        }
    }

    return filteredResults;
}

// Helper function to count occurrences of a coordinate pair (qstart, qend) in alignments
int FilterHomologs::countCoordinatePairs(int qstart, int qend, const std::vector<alignment_t>& alignments) {
    int count = 0;

    // Count the number of times (qstart, qend) appears in alignments
    for (const auto& alignment : alignments) {
        if (alignment.qstart == qstart && alignment.qend == qend) {
            ++count;
        }
    }

    return count;
}
