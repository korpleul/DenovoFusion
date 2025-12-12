//
// Created by xinwei on 12/3/24.
//

// paf_parser.h - Header file for parsing PAF files
#ifndef PAF_H
#define PAF_H

#include <string>
#include <vector>
#include <unordered_map>

// Define a structure to hold the information of a PAF alignment
struct paf_t {
    std::string query_name;
    int query_length;
    int query_start;
    int query_end;
    char strand;
    std::string target_name;
    int target_length;
    int target_start;
    int target_end;
    int match_length;
    int block_length;
    int mapping_quality;
    std::unordered_map<std::string, std::string> optional_fields;
    std::string cigar;
};

// Function to parse a PAF file
void paf_parse(const std::string& filename, std::vector<paf_t>& alignments);

void parse_line(const std::string& line, paf_t& alignment);

void validate_entry(const paf_t& alignment);

// Function to count occurrences of each query name
std::unordered_map<std::string, int> count_qnames(const std::vector<paf_t>& alignments);


#endif // PAF_H
