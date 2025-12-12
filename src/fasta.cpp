//
// Created by xinwei on 5/16/24.
//

#include "fasta.h"
#include <cctype>
#include <sstream>

// 序列类实现
SequenceCls::SequenceCls(const std::string &id, const std::string &extra)
        : id(id), extra(extra), sequence("") {}

size_t SequenceCls::length() const {
    return sequence.length();
}

std::string SequenceCls::output(bool maintain_case) const {
    std::ostringstream output_str;
    output_str << ">" << id;
    if (!extra.empty()) {
        output_str << " " << extra;
    }
    output_str << "\n";
    if (maintain_case) {
        output_str << sequence;
    } else {
        for (char c : sequence) {
            output_str << (char)std::toupper(c);
        }
    }
    return output_str.str();
}

// FASTA 文件处理类实现
FastaFileCls::FastaFileCls(const std::string &path, const std::string &fail_msg, const std::string &line_delim, bool maintain_case)
        : line_delim_(line_delim), finished_(false), maintain_case_(maintain_case) {
    openFile(path, fail_msg);
}

FastaFileCls::~FastaFileCls() {
    close();
}

void FastaFileCls::openFile(const std::string &path, const std::string &fail_msg) {
    file_.open(path);
    if (!file_.is_open()) {
        throw FastaError(fail_msg);
    }
}

bool FastaFileCls::next(SequenceCls &seq) {
    if (finished_) {
        return false;
    }

    try {
        if (curr_line_.empty()) {
            if (!std::getline(file_, curr_line_)) {
                finished_ = true;
                return false;
            }
        }

        if (curr_line_[0] != '>') {
            throw FastaError("Improperly formatted fasta file: sequence id line must begin with \">\": \"" + curr_line_ + "\".");
        }

        std::string seq_id, seq_extra;
        std::istringstream iss(curr_line_.substr(1));
        iss >> seq_id;
        std::getline(iss, seq_extra);
        seq = SequenceCls(seq_id, seq_extra);

        while (std::getline(file_, curr_line_)) {
            if (curr_line_[0] == '>') {
                break;
            }
            if (!seq.sequence.empty()) {
                seq.sequence += line_delim_;
            }
            seq.sequence += curr_line_;
        }

        if (!maintain_case_) {
            for (char &c : seq.sequence) {
                c = std::toupper(c);
            }
        }

        if (file_.eof()) {
            finished_ = true;
        }

        return true;
    } catch (const std::exception &e) {
        finished_ = true;
        throw;
    }
}

void FastaFileCls::reset() {
    file_.clear();
    file_.seekg(0);
    curr_line_.clear();
    finished_ = false;
}

void FastaFileCls::close() {
    if (file_.is_open()) {
        file_.close();
    }
}


std::unordered_map<std::string, std::string> load_fasta_sequences(const std::string &fasta_path) {
    FastaFileCls fasta_file(fasta_path);
    SequenceCls seq("", "");

    std::unordered_map<std::string, std::string> fasta_sequences;

    while (fasta_file.next(seq)) {
        fasta_sequences[seq.id] = seq.sequence;
    }

    return fasta_sequences;
}



// Collect and merge sequences based on base query name
std::unordered_map<std::string, std::string> collectAndMergeSequences(const std::vector<alignment_t>& alignments, const std::unordered_map<std::string, std::string>& fasta_sequences) {
    std::unordered_map<std::string, std::vector<std::pair<int, std::string>>> sequencesToMerge;

    for (const auto& alignment : alignments) {
        auto [baseName, number] = parseQuery(alignment.query);
        if (fasta_sequences.find(alignment.query) != fasta_sequences.end()) {
            sequencesToMerge[baseName].emplace_back(number, alignment.query);
        }
    }

    std::unordered_map<std::string, std::string> mergedSequences;
    for (auto& [baseName, parts] : sequencesToMerge) {
        // Sort by number
        std::sort(parts.begin(), parts.end());

        std::string combinedSequence;
        std::string lastPart;
        for (const auto& [_, partName] : parts) {
            if (partName != lastPart) {
                combinedSequence += fasta_sequences.at(partName);
                lastPart = partName;
            }
        }
        mergedSequences[baseName] = combinedSequence;
    }

    return mergedSequences;
}

// Collect and merge sequences based on base query name
std::unordered_map<std::string, std::string> collectSequences(const std::vector<alignment_t>& alignments, const std::unordered_map<std::string, std::string>& fasta_sequences) {
    std::unordered_map<std::string, std::string> mergedSequences;

    // Loop through each alignment
    for (const auto& alignment : alignments) {
        // Extract the base query name directly from the alignment
        const std::string& baseName = alignment.query;

        // Only add the sequence if it hasn't been added yet
        if (fasta_sequences.find(baseName) != fasta_sequences.end()) {
            // Check if the baseName has already been added to mergedSequences
            if (mergedSequences.find(baseName) == mergedSequences.end()) {
                // Add the sequence from fasta_sequences to mergedSequences
                mergedSequences[baseName] = fasta_sequences.at(baseName);
            }
        }
    }

    return mergedSequences;
}



// Output the merged sequence to a file
void outputMergedSequences(const std::unordered_map<std::string, std::string>& mergedSequences, const options_t& options) {
    std::ofstream output_file(options.output + "/" + options.prefix + ".chosen.fasta");
    if (!output_file.is_open()) {
        throw std::runtime_error("Failed to open output file: " + options.prefix + ".chosen.fasta");
    }

    for (const auto& [queryName, sequence] : mergedSequences) {
        output_file << ">" << queryName << "\t" << sequence.length() << "\n" << sequence << "\n";
    }

    output_file.close();
}