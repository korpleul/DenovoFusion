//
// Created by xinwei on 5/16/24.
//

#ifndef SAM_H
#define SAM_H

#include <string>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <unordered_map>


class SamError : public std::runtime_error {
public:
    explicit SamError(const std::string &message)
            : std::runtime_error("SamError: " + message) {}
};


class sam_t {
public:
    sam_t() {};
    sam_t(const std::string &line);

    std::string output() const;

    std::string qname;
    int flag;
    std::string rname;
    int pos;
    int mapq;
    std::string cigar;
    std::string rnext;
    int pnext;
    int tlen;
    std::string seq;
    std::string qual;
    std::vector<std::string> optional;
    std::string tp_label;

    // additional fields, computed from CIGAR information
    int num_matches;
    int num_insertions;
    int num_deletions;
    int num_soft_clips;
    int num_hard_clips;
    int num_skipped;

    void parseCigar();
};


class SamFileCls {
public:
    SamFileCls(const std::string &path, const std::string &fail_msg = "cannot open SAM file");
    ~SamFileCls();

    bool next(sam_t &sam);
    void reset();
    void close();

    std::ifstream file_;
    std::string curr_line_;
    bool finished_;

    void openFile(const std::string &path, const std::string &fail_msg);

};

std::unordered_map<std::string, int> count_qnames(const std::vector<sam_t>& alignments);

void loadSamFile(const std::string& filePath, std::vector<sam_t>& samEntries);

#endif // SAM_H