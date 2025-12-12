//
// Created by xinwei on 6/7/24.
//

#include "alignments_chosen.h"

void GetAlignLengthAndFraction(const psl_t& psls, int& align_len, double& align_fract) {
    align_len = (psls.qEnd - psls.qStart) + 1;  // Calculate alignment length
    align_fract = static_cast<double>(align_len) / static_cast<double>(psls.qSize);  // Calculate alignment fraction

}

using namespace std;

// Generate all possible combinations
    void generateCombinations(const vector<pair<int, int>>& pairs, int start, int k, vector<pair<int, int>>& current, vector<vector<pair<int, int>>>& allCombinations) {
        if (k == 0) {
            allCombinations.push_back(current);
            return;
        }

        for (int i = start; i <= pairs.size() - k; ++i) {
            current.push_back(pairs[i]);
            generateCombinations(pairs, i + 1, k - 1, current, allCombinations);
            current.pop_back();
        }
    }

// Calculating overlap score
    int calculateOverlapScore(const vector<pair<int, int>>& combination) {
        unordered_map<int, int> frequencyMap;
        for (const auto& p : combination) {
            for (int i = p.first; i <= p.second; ++i) {
                frequencyMap[i]++;
            }
        }

        int score = 0;
        for (const auto& entry : frequencyMap) {
            if (entry.second > 1) {
                score++;
            }
        }

        return score;
    }

// Calculates the value contained
int calculateInclusionScore(const std::vector<std::pair<int, int>>& combination) {
    std::unordered_set<int> coveredValues;
    for (const auto& p : combination) {
        for (int i = p.first; i <= p.second; ++i) {
            coveredValues.insert(i);
        }
    }

    return coveredValues.size();
}




std::vector<alignment_t> calculate_alignments_score(const std::unordered_map<std::string, std::vector<alignment_t>>& data_by_qname, const options_t& options) {
    std::vector<alignment_t> best_alignment;

    // Traverse the grouped data
    for (const auto& [qname, alignments_group] : data_by_qname) {
        std::vector<std::pair<int, int>> pairs;
        int qSize = alignments_group[0].query_len;  // Assume that all qSizes under the same qName are the same

        // Extract pairs from alignments_group
        for (const auto& alignment : alignments_group) {
            pairs.emplace_back(alignment.qstart, alignment.qend);
        }

        // Variables are used to track the best score for each qName and its information
        double best_score = -std::numeric_limits<double>::infinity();
        std::vector<std::pair<int, int>> best_combination;
        std::vector<alignment_t> current_best_alignment; // Vector storing the best alignment data

        // Used to mark whether a perfect alignment has been found
        bool perfect_alignment_found = false;

        // Traverse all possible combinations
        for (int k = 1; k <= std::min(static_cast<int>(pairs.size()), options.max_pair_combination); ++k) {
            std::vector<std::vector<std::pair<int, int>>> allCombinations;
            std::vector<std::pair<int, int>> current;
            generateCombinations(pairs, 0, k, current, allCombinations);

            // Process each combination
            for (const auto& combination : allCombinations) {
                // Calculate the identity sum
                double identity_sum = 0.0;
                std::vector<alignment_t> current_alignment; // Temporarily store the alignment data of the current combination
                for (const auto& p : combination) {
                    // Find the corresponding alignment_t object
                    for (const auto& alignment : alignments_group) {
                        if (alignment.qstart == p.first && alignment.qend == p.second) {
                            identity_sum += alignment.identity;
                            current_alignment.push_back(alignment); // Stores the alignment data for the current combination
                            break;
                        }
                    }
                }

                // Calculating overlap and inclusion scores
                int overlap = calculateOverlapScore(combination);
                float overlap_fraction = static_cast<float>(overlap) / qSize;

                int inclusion = calculateInclusionScore(combination);
                float inclusion_fraction = static_cast<float>(inclusion) / qSize;

                // calculate final score for all single or combinations alignments
                double score = identity_sum + options. inclusion_fraction_weight * inclusion_fraction - options.inclusion_fraction_weight * overlap_fraction - options.size_weight * combination.size();

                // Compare and update the best scores and their information
                if (score > best_score) {
                    best_score = score;
                    best_combination = combination;
                    current_best_alignment = current_alignment; // Update best alignment data
                }

                // Check if perfect alignment is found
                if (current_alignment.size() == 1 &&
                    (current_alignment[0].matches + current_alignment[0].mismatch) == qSize) {
                    perfect_alignment_found = true;
                    break;  // Stop Iteration
                }
            }

            // If a perfect alignment is found, stop iterating.
            if (perfect_alignment_found) {
                break;
            }
        }

        // Store the best alignment found into best_alignment
        best_alignment.insert(best_alignment.end(), current_best_alignment.begin(), current_best_alignment.end());
    }

    return best_alignment;
}
