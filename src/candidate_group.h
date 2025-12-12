//
// Created by xinwei on 5/19/24.
//

#ifndef CANDIDATE_GROUP_H
#define CANDIDATE_GROUP_H


#include "coord_pair.h"
#include "psl_functions.h"
#include "alignments_chosen.h"


#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <optional>
#include <unordered_map>
#include <algorithm>
#include <utility>

std::vector<alignment_t> single_alignments(const std::vector<alignment_t>& alignments);

std::vector<std::pair<alignment_t, alignment_t>> pair_alignments(const std::vector<alignment_t>& alignments);

std::vector<alignment_t> multiple_alignments(const std::vector<alignment_t>& alignments);

void group_alignments(const std::vector<alignment_t>& alignments);


enum AlignmentCategory {
    CONTAINS_SAME_STRAND,
    CONTAINS_DIFFERENT_STRAND,
    OVERLAPS_SAME_STRAND,
    OVERLAPS_DIFFERENT_STRAND,
    GAP_SAME_STRAND,
    GAP_DIFFERENT_STRAND
};


std::unordered_map<AlignmentCategory, std::vector<std::pair<alignment_t, alignment_t>>> classify_alignments(
        const std::vector<std::pair<alignment_t, alignment_t>>& paired_alignments, const options_t& options);



std::string AlignmentCategoryToString(AlignmentCategory category);


std::vector<alignment_t> filter_and_combine_same_strand(
        const std::unordered_map<AlignmentCategory, std::vector<std::pair<alignment_t, alignment_t>>>& classified_alignments,
        AlignmentCategory category, int max_value, bool is_gap);
std::vector<alignment_t> filter_and_combine_diff_strand(
        const std::unordered_map<AlignmentCategory, std::vector<std::pair<alignment_t, alignment_t>>>& classified_alignments,
        AlignmentCategory category, int max_value, bool is_gap);

std::unordered_set<std::string> extractBaseQueryNames(const std::vector<alignment_t>& alignments);
std::vector<alignment_t> filterAlignmentsByBaseQuery(const std::vector<alignment_t>& alignments, const std::unordered_set<std::string>& base_query_names);
std::vector<alignment_t> filterAlignmentsByQuery(const std::vector<alignment_t>& alignments, const std::unordered_set<std::string>& query_names);

std::pair<std::string, int> parseQuery(const std::string& query);
bool compareAlignments(const alignment_t& a, const alignment_t& b);
std::vector<alignment_t> removeDiscontinuousChromosomes(const std::vector<alignment_t>& alignments);

std::unordered_set<std::string> extractQueryNames(const std::vector<alignment_t>& alignments);

#endif // CANDIDATE_GROUP_H