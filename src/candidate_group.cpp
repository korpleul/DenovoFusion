//
// Created by xinwei on 6/7/24.
//

#include "candidate_group.h"



std::vector<alignment_t> single_alignments(const std::vector<alignment_t>& alignments) {
    std::vector<alignment_t> single_alignments;
    std::unordered_map<std::string, int> queryCount;

    // Count the number of times each query name occurs
    for (const auto& alignment : alignments) {
        queryCount[alignment.query]++;
    }

    // Filter alignments where the block number is 1 and the query name appears uniquely in all alignments
    for (const auto& alignment : alignments) {
        if (alignment.blockcount == 1 && queryCount[alignment.query] == 1) {
            single_alignments.push_back(alignment);
        }
    }
    return single_alignments;
}

std::vector<alignment_t> gap_alignments(const std::vector<alignment_t>& alignments) {
    std::vector<alignment_t> gap_alignments;
    std::unordered_map<std::string, int> queryCount;

    // Count the number of times each query name occurs
    for (const auto& alignment : alignments) {
        queryCount[alignment.query]++;
    }

    // Filter alignments where the block number is 1 and the query name is not the only one in all alignments
    for (const auto& alignment : alignments) {
        if (alignment.blockcount > 1 && queryCount[alignment.query] == 1) {
            gap_alignments.push_back(alignment);
        }
    }
    return gap_alignments;
}

std::vector<std::pair<alignment_t, alignment_t>> pair_alignments(const std::vector<alignment_t>& alignments) {
    std::unordered_map<std::string, std::vector<alignment_t>> queryMap;

    // Group each alignment by query name
    for (const auto& alignment : alignments) {
        queryMap[alignment.query].push_back(alignment);
    }

    std::vector<std::pair<alignment_t, alignment_t>> pair_alignments;

    // Filter out queries from the map that have exactly two alignments, and make sure they form a pair
    for (const auto& pair : queryMap) {
        if (pair.second.size() == 2) {
            pair_alignments.emplace_back(pair.second[0], pair.second[1]);
        }
    }

    return pair_alignments;
}

std::vector<alignment_t> multiple_alignments(const std::vector<alignment_t>& alignments) {
    std::vector<alignment_t> multiple_alignments;
    std::unordered_map<std::string, std::vector<alignment_t>> queryAlignments;

    // Store the alignment corresponding to each query name in a vector
    for (const auto& alignment : alignments) {
        queryAlignments[alignment.query].push_back(alignment);
    }

    // Only add all alignments with a query name whose number of alignments is greater than 2
    for (const auto& query_pair : queryAlignments) {
        if (query_pair.second.size() > 2) {
            multiple_alignments.insert(multiple_alignments.end(), query_pair.second.begin(), query_pair.second.end());
        }
    }

    return multiple_alignments;
}

void group_alignments(const std::vector<alignment_t>& alignments) {
    auto singles = single_alignments(alignments);

    auto gaps = gap_alignments(alignments);

    auto pairs = pair_alignments(alignments);

    auto multiples = multiple_alignments(alignments);
}


std::string AlignmentCategoryToString(AlignmentCategory category) {
    static const std::unordered_map<AlignmentCategory, std::string> categoryToString = {
            {CONTAINS_SAME_STRAND, "CONTAINS_SAME_STRAND"},
            {CONTAINS_DIFFERENT_STRAND, "CONTAINS_DIFFERENT_STRAND"},
            {OVERLAPS_SAME_STRAND, "OVERLAPS_SAME_STRAND"},
            {OVERLAPS_DIFFERENT_STRAND, "OVERLAPS_DIFFERENT_STRAND"},
            {GAP_SAME_STRAND, "GAP_SAME_STRAND"},
            {GAP_DIFFERENT_STRAND, "GAP_DIFFERENT_STRAND"}
    };

    auto it = categoryToString.find(category);
    if (it != categoryToString.end()) {
        return it->second;
    } else {
        return "UNKNOWN_CATEGORY";
    }
}

std::unordered_map<AlignmentCategory, std::vector<std::pair<alignment_t, alignment_t>>> classify_alignments(
        const std::vector<std::pair<alignment_t, alignment_t>>& paired_alignments, const options_t& options) {

    std::unordered_map<AlignmentCategory, std::vector<std::pair<alignment_t, alignment_t>>> classified_alignments;

    for (const auto& alignment_pair : paired_alignments) {
        const alignment_t& align1 = alignment_pair.first;
        const alignment_t& align2 = alignment_pair.second;

        CoordPair coord1(align1.qstart, align1.qend, align1.query_strand == '+', align1.query);
        CoordPair coord2(align2.qstart, align2.qend, align2.query_strand == '+', align2.query);

        bool same_strand = coord1.pos_strand == coord2.pos_strand;
        int overlap_size = std::max(0, std::min(coord1.end, coord2.end) - std::max(coord1.start, coord2.start));
        int gap_size = std::max(0, std::max(coord1.start, coord2.start) - std::min(coord1.end, coord2.end));

        // Check for containment
        if (coord1.Contains(coord2) || coord2.Contains(coord1)) {
            if (same_strand) {
                classified_alignments[CONTAINS_SAME_STRAND].push_back(alignment_pair);
            } else {
                classified_alignments[CONTAINS_DIFFERENT_STRAND].push_back(alignment_pair);
            }
        }
            // Check for overlaps
        else if (coord1.Overlaps(coord2)) {
            if (same_strand && overlap_size <= options.max_overlap_size) {
                classified_alignments[OVERLAPS_SAME_STRAND].push_back(alignment_pair);
            } else if (!same_strand && overlap_size <= options.max_overlap_size) {
                classified_alignments[OVERLAPS_DIFFERENT_STRAND].push_back(alignment_pair);
            }
        }
            // Check for gaps
        else if (coord1.Gap(coord2)) {
            if (same_strand && gap_size <= options.max_gap_size) {
                classified_alignments[GAP_SAME_STRAND].push_back(alignment_pair);
            } else if (!same_strand && gap_size <= options.max_gap_size) {
                classified_alignments[GAP_DIFFERENT_STRAND].push_back(alignment_pair);
            }
        }
    }

    return classified_alignments;
}

// Add filter_and_combine_same_strand function implementation
std::vector<alignment_t> filter_and_combine_same_strand(
        const std::unordered_map<AlignmentCategory, std::vector<std::pair<alignment_t, alignment_t>>>& classified_alignments,
        AlignmentCategory category, int max_value, bool is_gap) {

    std::vector<alignment_t> result;

    if (classified_alignments.find(category) == classified_alignments.end()) {
        return result;  // If there is no such category in the classification, return an empty result directly
    }

    for (const auto& alignment_pair : classified_alignments.at(category)) {
        const alignment_t& align1 = alignment_pair.first;
        const alignment_t& align2 = alignment_pair.second;

        CoordPair coord1(align1.qstart, align1.qend, align1.query_strand == '+', align1.query);
        CoordPair coord2(align2.qstart, align2.qend, align2.query_strand == '+', align2.query);

        int value;
        if (is_gap) {
            value = std::max(0, coord2.start - coord1.end);  // 计算 gap 的大小
        } else {
            value = std::max(0, std::min(coord1.end, coord2.end) - std::max(coord1.start, coord2.start));  // 计算 overlap 的大小
        }

        if (value <= max_value) {
            result.push_back(align1);
            result.push_back(align2);
        }
    }

    return result;
}

std::vector<alignment_t> filter_and_combine_diff_strand(
        const std::unordered_map<AlignmentCategory, std::vector<std::pair<alignment_t, alignment_t>>>& classified_alignments,
        AlignmentCategory category, int max_value, bool is_gap) {

    std::vector<alignment_t> result;

    // If the category is not in the classification, return an empty result directly
    if (classified_alignments.find(category) == classified_alignments.end()) {
        return result;
    }

    for (const auto& alignment_pair : classified_alignments.at(category)) {
        const alignment_t& align1 = alignment_pair.first;
        const alignment_t& align2 = alignment_pair.second;

        // create the coordinate pairs for each alignment
        CoordPair coord1(align1.qstart, align1.qend, align1.query_strand == '+', align1.query);
        CoordPair coord2(align2.qstart, align2.qend, align2.query_strand == '+', align2.query);

        // calculate the overlap or gap size
        int value;
        if (is_gap) {
            value = std::max(0, coord2.start - coord1.end);
        } else {
            value = std::max(0, std::min(coord1.end, coord2.end) - std::max(coord1.start, coord2.start));  // 计算重叠大小
        }

        // if the overlap or gap size is less than or equal to the maximum value, add the alignments to the result list
        if (value <= max_value) {
            result.push_back(align1);
            result.push_back(align2);
        }
    }

    return result;
}


// Extract the base part of the query name from the alignment data and avoid duplication
std::unordered_set<std::string> extractBaseQueryNames(const std::vector<alignment_t>& alignments) {
    std::unordered_set<std::string> base_query_names;
    for (const auto& alignment : alignments) {
        size_t second_underscore = alignment.query.find('_', alignment.query.find('_') + 1);
        if (second_underscore != std::string::npos) {
            std::string base_query = alignment.query.substr(0, second_underscore);
            base_query_names.insert(base_query);
        } else {
            // If there is no second underscore, just add the entire query name
            base_query_names.insert(alignment.query);
        }
    }
    return base_query_names;
}

// Filter alignments based on base query name
std::vector<alignment_t> filterAlignmentsByBaseQuery(
    const std::vector<alignment_t>& alignments,
    const std::unordered_set<std::string>& base_query_names
) {
    std::vector<alignment_t> filtered_alignments;
    for (const auto& alignment : alignments) {
        // ✅ Added: also match the full query name directly
        if (base_query_names.find(alignment.query) != base_query_names.end()) {
            filtered_alignments.push_back(alignment);
            continue; // no need to parse underscores
        }

        size_t second_underscore = alignment.query.find('_', alignment.query.find('_') + 1);
        if (second_underscore != std::string::npos) {
            std::string base_query = alignment.query.substr(0, second_underscore);
            if (base_query_names.find(base_query) != base_query_names.end()) {
                filtered_alignments.push_back(alignment);
            }
        }
    }
    return filtered_alignments;
}

// Check and remove discontinuous alignments of the target chromosome
std::vector<alignment_t> removeDiscontinuousChromosomes(const std::vector<alignment_t>& alignments) {
    std::unordered_map<std::string, std::set<std::string>> chromosomeMap;
    std::vector<alignment_t> filtered_fragment_alignments;

    // First, collect the target chromosome for each base query name
    for (const auto& alignment : alignments) {
        std::string baseName = alignment.query.substr(0, alignment.query.find('_', alignment.query.find('_') + 1));
        chromosomeMap[baseName].insert(alignment.target);
    }

    // Then, only add the alignment to the results list if the number of target chromosomes is less than 3
    for (const auto& alignment : alignments) {
        std::string baseName = alignment.query.substr(0, alignment.query.find('_', alignment.query.find('_') + 1));
        if (chromosomeMap[baseName].size() < 3) {
            filtered_fragment_alignments.push_back(alignment);
        }
    }

    return filtered_fragment_alignments;
}


// Parsing query name and number with safe check
std::pair<std::string, int> parseQuery(const std::string& query) {
    size_t lastUnder = query.rfind('_');
    if (lastUnder != std::string::npos) {
        std::string last_part = query.substr(lastUnder + 1);
        if (last_part.size() <= 3 &&
            std::all_of(last_part.begin(), last_part.end(), ::isdigit)) {
            std::string baseName = query.substr(0, lastUnder);
            int number = std::stoi(last_part);
            return {baseName, number};
            }
    }
    return {query, -1};
}


// compare
bool compareAlignments(const alignment_t& a, const alignment_t& b) {
    auto [baseA, numA] = parseQuery(a.query);
    auto [baseB, numB] = parseQuery(b.query);
    if (baseA != baseB) return baseA < baseB;
    if (numA != numB) return numA < numB;
    return a.qstart < b.qstart;
}

// Extract the query name from the alignment data and avoid duplication
std::unordered_set<std::string> extractQueryNames(const std::vector<alignment_t>& alignments) {
    std::unordered_set<std::string> base_query_names;
    for (const auto& alignment : alignments) {
        base_query_names.insert(alignment.query);
    }
    return base_query_names;
}

// Filter alignments based on base query name
std::vector<alignment_t> filterAlignmentsByQuery(const std::vector<alignment_t>& alignments, const std::unordered_set<std::string>& query_names) {
    std::vector<alignment_t> filtered_alignments;

    for (const auto& alignment : alignments) {
        // Directly check if the alignment query is in the base query names set
        if (query_names.find(alignment.query) != query_names.end()) {
            filtered_alignments.push_back(alignment);
        }
    }

    return filtered_alignments;
}
