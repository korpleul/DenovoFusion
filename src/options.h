//
// Created by xinwei on 5/19/24.
//

#ifndef OPTION_PARSER_H
#define OPTION_PARSER_H

#include <string>
#include <vector>
#include <getopt.h>
#include <filesystem>
#include <climits>
#include <cfloat>

#include "common.h"

const std::string DENOVOFUSION_VERSION = "1.0.1";


struct options_t {

    std::string input_type;
    std::string input_file;
    std::string input_assembly;
    std::string output;
    std::string gtf_path;
    std::string prefix;
    std::vector<std::string> input_fastq1;
    std::vector<std::string> input_fastq2;

    int max_alignment_count;
    float min_identity_fract;
    int min_score_total;
    int min_score_each;
    float coverage_differ;

    int max_pair_combination;
    int threads;
    int max_overlap_size;
    int max_gap_size;
    int read_length;
    int min_edge_length;
    float size_ratio_threshold;
    int edge_unaligned;
    float inclusion_fraction_weight;
    float overlap_fraction_weight;
    float size_weight;

    int long_gap_threshold;
    int short_segment_threshold;
    int min_split_reads;
    int min_span_reads;
    int max_itd_length;
    float min_itd_fraction;



    std::unordered_map<std::string,bool> filters;
};


    options_t option_parser(int argc, char** argv);

    bool output_directory_exists(const std::string& output_file);

    std::string wrap_help(const std::string& option, const std::string& text, const unsigned short int max_line_width = 80);

    std::string wrap_help2(const std::string& text, const unsigned short int max_line_width = 80);


    bool validate_int(const char* optarg, unsigned int& value, const unsigned int min_value = 0, const unsigned int max_value = INT_MAX);
    bool validate_float(const char* optarg, float& value, const float min_value = FLT_MIN, const float max_value = FLT_MAX);

    void parse_fastq_files(const char* optarg, std::vector<std::string>& fastq_files);

#endif //OPTIONS_H
