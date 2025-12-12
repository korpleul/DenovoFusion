#include "paf.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

// Function to parse a PAF file
void paf_parse(const std::string& filename, std::vector<paf_t>& alignments) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("无法打开文件: " + filename);
    }

    alignments.clear();  // clear the vector and fill in a new vector
    std::string line;

    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // discard empty lines and comments
        }

        paf_t alignment;
        parse_line(line, alignment);
        validate_entry(alignment);
        if (alignment.optional_fields.count("tp") > 0 && alignment.optional_fields["tp"] == "A:P") {
            alignments.push_back(alignment);  // store alignment in the container only if it's primary
        }
    }

    file.close();
}


// Function to parse a line of PAF format
void parse_line(const std::string& line, paf_t& alignment) {
    std::istringstream ss(line);

    ss >> alignment.query_name
       >> alignment.query_length
       >> alignment.query_start
       >> alignment.query_end
       >> alignment.strand
       >> alignment.target_name
       >> alignment.target_length
       >> alignment.target_start
       >> alignment.target_end
       >> alignment.match_length
       >> alignment.block_length
       >> alignment.mapping_quality;

    // Parse the optional fields
    std::string field;
    while (ss >> field) {
        size_t colon_pos = field.find(':');
        if (colon_pos != std::string::npos) {
            std::string key = field.substr(0, colon_pos);
            std::string value = field.substr(colon_pos + 1);

            if (key == "cg") {
                alignment.cigar = value.substr(2); // Extract CIGAR without "Z:"
            } else {
                alignment.optional_fields[key] = value; // Store in optional fields
            }

        }
    }
}


// Function to validate an entry (optional)
void validate_entry(const paf_t& alignment) {
    if (alignment.query_start > alignment.query_end || alignment.target_start > alignment.target_end) {
        throw std::runtime_error("Invalid PAF alignment range.");
    }
}

// Function to count occurrences of each query name in a vector of PAFs
std::unordered_map<std::string, int> count_qnames(const std::vector<paf_t>& alignments) {
    std::unordered_map<std::string, int> qname_counts;
    for (const auto& alignment : alignments) {
        qname_counts[alignment.query_name]++;
    }
    return qname_counts;
}
