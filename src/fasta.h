//
// Created by xinwei on 5/16/24.
//

#ifndef FASTA_H
#define FASTA_H

#include <string>
#include <fstream>
#include <stdexcept>
#include <iterator>
#include <unordered_map>
#include <unordered_set>

#include "alignment.h"
#include "options.h"


class FastaError : public std::runtime_error {
public:
    explicit FastaError(const std::string &message)
            : std::runtime_error("FastaError: " + message) {}
};


class SequenceCls {
public:
    SequenceCls(const std::string &id, const std::string &extra);
    size_t length() const;
    std::string output(bool maintain_case = false) const;

    std::string id;
    std::string extra;
    std::string sequence;
};

// FASTA 文件处理类
class FastaFileCls {
public:
    FastaFileCls(const std::string &path, const std::string &fail_msg = "cannot open fasta file", const std::string &line_delim = "", bool maintain_case = false);
    ~FastaFileCls();

    bool next(SequenceCls &seq);
    void reset();
    void close();

private:
    std::ifstream file_;
    std::string line_delim_;
    std::string curr_line_;
    bool finished_;
    bool maintain_case_;

    void openFile(const std::string &path, const std::string &fail_msg);
};

std::unordered_map<std::string, std::string> load_fasta_sequences(const std::string &fasta_path) ;


std::string mergeSequences(const std::vector<std::string>& fragments, const std::unordered_map<std::string, std::string>& fasta_sequences);

std::string getBaseQueryName(const std::string& query);

void outputMergedContigs(const std::vector<alignment_t>& sorted_alignments, const std::unordered_map<std::string, std::string>& fasta_sequences, const options_t& options);

std::pair<std::string, int> parseQuery(const std::string& query);
std::unordered_map<std::string, std::string> collectAndMergeSequences(const std::vector<alignment_t>& alignments, const std::unordered_map<std::string, std::string>& fasta_sequences);
std::unordered_map<std::string, std::string> collectSequences(const std::vector<alignment_t>& alignments, const std::unordered_map<std::string, std::string>& fasta_sequences);
void outputMergedSequences(const std::unordered_map<std::string, std::string>& mergedSequences, const options_t& options);

#endif // FASTA_H