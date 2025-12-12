//
// Created by xinwei on 5/16/24.
//

#include "sam.h"
#include <sstream>
#include <cctype>
#include <vector>
#include <string>
#include <iostream>


// Parse CIGAR string and fill in num_matches, num_insertions, num_deletions, etc.
void parse_cigar(const std::string &cigar, sam_t &sam) {
    int num = 0; // Used to store the number before each operation

    for (char ch : cigar) {
        if (isdigit(ch)) {
            num = num * 10 + (ch - '0'); // Cumulative figures
        } else {
            // Update the corresponding fields according to CIGAR operations
            switch (ch) {
                case 'M': // （match/mismatch）
                case '=': // Exact match
                    sam.num_matches += num;
                break;
                case 'I': // （insertion）
                    sam.num_insertions += num;
                break;
                case 'D': // （deletion）
                    sam.num_deletions += num;
                break;
                case 'S': // （soft clipping）
                    sam.num_soft_clips += num;
                break;
                case 'H': // （hard clipping）
                    sam.num_hard_clips += num;
                break;
                case 'N': // Skipping (often used to describe introns in RNA-seq)
                    sam.num_skipped += num;
                break;
                case 'X': // Mismatch (explicit mismatch)
                    sam.num_matches += num; // You may need another field to store mismatches explicitly if needed
                break;
                default:
                    // Other operations (e.g. P, B, PU, etc.)
                    // Ignore them for now
                break;
            }
            // Reset the number for the next operation
            num = 0;
        }
    }
}



sam_t::sam_t(const std::string &line)
    : num_matches(0), num_insertions(0), num_deletions(0), num_soft_clips(0), num_hard_clips(0), num_skipped(0), tp_label("") {
    std::istringstream ss(line);
    ss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext >> tlen >> seq >> qual;

    // Parsing optional fields
    std::string optional_field;
    bool has_tp_A = false; // Track if "tp:A" is found

    while (ss >> optional_field) {
        optional.push_back(optional_field);
        // Check if tp:A is found in optional fields
        if (optional_field.find("tp:A") != std::string::npos) {
            has_tp_A = true;
            // Check if it's specifically "tp:A:P"
            if (optional_field == "tp:A:P") {
                tp_label = "P";  // It's a primary alignment
            }
        }
    }

    // If tp:A is present but not "tp:A:P", mark as secondary
    if (has_tp_A && tp_label.empty()) {
        tp_label = "S";  // Secondary alignment
    }

    // Fill in num_matches, num_insertions, num_deletions, etc. based on CIGAR string
    parse_cigar(cigar, *this);
}




std::string sam_t::output() const {
    std::ostringstream oss;
    oss << qname << "\t" << flag << "\t" << rname << "\t" << pos << "\t" << mapq << "\t" << cigar << "\t"
        << rnext << "\t" << pnext << "\t" << tlen << "\t" << seq << "\t" << qual;
    for (const auto &opt : optional) {
        oss << "\t" << opt;
    }
    return oss.str();
}


SamFileCls::SamFileCls(const std::string &path, const std::string &fail_msg)
        : finished_(false) {
    openFile(path, fail_msg);
}

SamFileCls::~SamFileCls() {
    close();
}

void SamFileCls::openFile(const std::string &path, const std::string &fail_msg) {
    file_.open(path);
    if (!file_.is_open()) {
        throw SamError(fail_msg);
    }
}

bool SamFileCls::next(sam_t &sam) {
    if (finished_) {
        return false;
    }

    while (std::getline(file_, curr_line_)) {
        if (curr_line_.empty() || curr_line_[0] == '@') {  // Skip header rows or empty lines
            continue;
        }

        sam = sam_t(curr_line_);

        // Skip secondary alignments (tp:A:S)
        if (sam.tp_label == "S") {
            continue;
        }

        // If CIGAR is "*", skip the line and continue reading the next line
        if (sam.cigar == "*") {
            continue;
        }

        return true;  // If a valid SAM record is found, it returns
    }

    finished_ = true;
    return false;
}





void SamFileCls::reset() {
    file_.clear();
    file_.seekg(0);
    curr_line_.clear();
    finished_ = false;
}

void SamFileCls::close() {
    if (file_.is_open()) {
        file_.close();
    }
}

// Function to count occurrences of each query name in a vector of PAFs
std::unordered_map<std::string, int> count_qnames(const std::vector<sam_t>& alignments) {
    std::unordered_map<std::string, int> qname_counts;
    for (const auto& alignment : alignments) {
        qname_counts[alignment.qname]++;
    }
    return qname_counts;
}


// Function to load a SAM file into a vector of sam_t objects
void loadSamFile(const std::string& filePath, std::vector<sam_t>& samEntries) {
    // Open the SAM file
    SamFileCls samFile(filePath);

    sam_t samEntry;  // Temporarily store each entry
    while (samFile.next(samEntry)) {  // Loop through each entry in the SAM file
        samEntries.push_back(samEntry);  // Store the parsed entries into a vector
    }

    // Check if the SAM file is empty
    if (samEntries.empty()) {
        throw SamError("SAM file is empty：" + filePath);
    }
}