#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <regex>
#include <iterator>
#include <cmath>


#include "psl.h"
#include "paf.h"
#include "sam.h"



// class alignment
class alignment_t {
public:

    std::string method_;
    std::string query;
    std::string target;
    int query_len;
    int target_len;
    char query_strand;
    int qstart;
    int qend;
    int tstart;
    int tend;
    int num_bases_aligned;
    int mismatch;
    int qnuminsert;
    int tnuminsert;
    int matches;
    int repmatch;
    int tbaseinsert;
    int qbaseinsert;
    int blockcount;
    double identity;
    double score;
    std::string model;
    std::string pairwise;
    std::string psl_str;
    std::vector<std::pair<int, int>> blocks;
    std::vector<std::pair<int, int>> query_blocks;
    std::vector<std::string> splice_sites;
    char orient;
    int contig;


    // Added constructors to allow direct initialization of query, qstart, and qend
    alignment_t(const std::string& q, int qs, int qe) : query(q), qstart(qs), qend(qe) {}


    alignment_t(const std::string &method, const psl_t &psls)
            : method_(method),
              query(psls.qName),
              target(psls.tName),
              query_len(psls.qSize),
              target_len(psls.tSize),
              query_strand(psls.strand[0]),
              qstart(psls.qStart),
              qend(psls.qEnd),
              tstart(psls.tStart),
              tend(psls.tEnd),
              num_bases_aligned(psls.matches + psls.misMatches + psls.repMatches),
              mismatch(psls.misMatches),
              qnuminsert(psls.qNumInsert),
              tnuminsert(psls.tNumInsert),
              matches(psls.matches),
              repmatch(psls.repMatches),
              tbaseinsert(psls.tBaseInsert),
              qbaseinsert(psls.qBaseInsert),
              blockcount(psls.blockCount)

              {
                  // score add for each alignment
                  identity = setIdentity(psls.qStart,psls.qEnd, psls.tStart, psls.tEnd, psls.qNumInsert, psls.misMatches, num_bases_aligned);
                  score = setScore(psls.matches, psls.misMatches, psls.qNumInsert, psls.tNumInsert, psls.qSize);
                  // set model as default
                  model = "BLAT";
                  // set pairwise as default
                  pairwise = "";
                  // psl_str
                  psl_str = psl();


                  for (int i = 0; i < psls.blockCount; ++i) {
                      blocks.emplace_back(psls.qStarts[i], psls.qEnds[i]);
                      query_blocks.emplace_back(psls.tStarts[i], psls.tEnds[i]);
        }
    }

    // New constructor from PAF
    alignment_t(const std::string &method, const paf_t &pafs)
        : method_(method),
          query(pafs.query_name),
          target(pafs.target_name),
          query_len(pafs.query_length),
          target_len(pafs.target_length),
          query_strand(pafs.strand),
          qstart(pafs.query_start),
          qend(pafs.query_end),
          tstart(pafs.target_start),
          tend(pafs.target_end),
          num_bases_aligned(pafs.match_length),  //
          mismatch(0),  // Mismatch information is not given in PAF and can be calculated during further processing.
          qnuminsert(0), //
          tnuminsert(0), //
          matches(pafs.match_length),
          repmatch(0),
          tbaseinsert(0),
          qbaseinsert(0),
          blockcount(1) {


        // Extract mismatch, qnuminsert, and tnuminsert from optional fields
        if (pafs.optional_fields.find("NM") != pafs.optional_fields.end()) {
            mismatch = std::stoi(pafs.optional_fields.at("NM").substr(2));
        }
        if (pafs.optional_fields.find("XI") != pafs.optional_fields.end()) {
            qnuminsert = std::stoi(pafs.optional_fields.at("XI").substr(2));  //
        }
        if (pafs.optional_fields.find("XT") != pafs.optional_fields.end()) {
            tnuminsert = std::stoi(pafs.optional_fields.at("XT").substr(2));  //
        }



        identity = setIdentity(qstart, qend, tstart, tend, qnuminsert, mismatch, num_bases_aligned);
        score = setScore(matches, mismatch, qnuminsert, tnuminsert, query_len);
        model = "Minimap2";
        pairwise = "";
    }



// New constructor from SAM
alignment_t(const std::string &method, const sam_t &sam, const std::unordered_map<std::string, std::string>& fasta_sequences)
    : method_(method),
      query(sam.qname),
      target(sam.rname),
      query_len(0),  // To be calculated from fasta_sequences
      target_len(0),  // Target length is not given in SAM
      query_strand((sam.flag & 16) ? '-' : '+'),
      qstart(0),
      qend(0),
      tstart(sam.pos - 1),
      tend(0),
      num_bases_aligned(0),
      mismatch(0),
      qnuminsert(0),
      tnuminsert(0),
      matches(0),
      repmatch(0),
      tbaseinsert(0),
      qbaseinsert(0),
      blockcount(1) {

    // Step 1: Extract contig length
    auto it = fasta_sequences.find(query);
    if (it != fasta_sequences.end()) {
        query_len = it->second.length();
    } else {
        std::cerr << "Warning: Could not find contig length for query: " << query << std::endl;
        query_len = 0;
    }

    // Step 2: Parse optional fields for mismatch (NM:i:)
    for (const auto &field : sam.optional) {
        if (field.find("NM:i:") == 0) {
            int nm = std::stoi(field.substr(5));
            mismatch = nm;  // Adjust below after parsing CIGAR
        }
    }

     // Reset variables
int current_qpos = 0;
int current_tpos = tstart;
int num = 0;
bool qstart_set = false;
bool is_reverse = (sam.flag & 16); // Negative strand (FLAG 16)
bool is_softclip_or_hardclip_in_front = false; // Flag to check if softclip or hardclip is before the match

// Check if S or H comes before M in the CIGAR string
if (!is_reverse) {
    // Forward strand: Check if S or H comes before M
    if (sam.cigar.find('S') < sam.cigar.find('M') || sam.cigar.find('H') < sam.cigar.find('M')) {
        is_softclip_or_hardclip_in_front = true;
    }
} else {
    // Reverse strand: Check if S or H comes before M (we reverse the CIGAR string once)
    std::string reversed_cigar = sam.cigar;
    std::reverse(reversed_cigar.begin(), reversed_cigar.end());  // Reverse the string once

    // Check if S or H comes before M in the reversed CIGAR string
    if (reversed_cigar.find('S') < reversed_cigar.find('M') || reversed_cigar.find('H') < reversed_cigar.find('M')) {
        is_softclip_or_hardclip_in_front = true;
    }
}


        // Parse the CIGAR string
        std::string parsed_cigar = sam.cigar;

        // If the strand is reverse, reverse each number and character in the CIGAR
        if (is_reverse) {
            std::string reversed_cigar;
            std::string current_number = "";  // To store the numeric part
            std::vector<std::pair<std::string, char>> segments; // To store number-operation pairs

            // Traverse the CIGAR string from left to right and extract segments
            for (int i = 0; i < parsed_cigar.size(); ++i) {
                char ch = parsed_cigar[i];

                if (isdigit(ch)) {
                    current_number += ch;  // Accumulate the number part
                } else if (isalpha(ch)) {
                    // Once we encounter a character (M, S, etc.), we have a full segment
                    if (!current_number.empty()) {
                        segments.push_back({current_number, ch});  // Store the number and operation pair
                        current_number.clear();  // Reset the number for the next segment
                    }
                }
            }

            // Check if there's an unpaired number at the end of the string
            if (!current_number.empty()) {
                std::cerr << "Error: Unpaired number in CIGAR string!\n";
                return;
            }

            // Reconstruct the reversed CIGAR string
            for (auto it = segments.rbegin(); it != segments.rend(); ++it) {
                reversed_cigar += it->first + it->second; // Append number and operation
            }

            parsed_cigar = reversed_cigar; // Update the reversed CIGAR string
        }





        // Parse CIGAR string
        for (char ch : parsed_cigar) {
            if (isdigit(ch)) {
                num = num * 10 + (ch - '0');
            } else {
                if (ch == 'M' || ch == '=' || ch == 'X') {
                    if (!qstart_set) {
                        qstart = current_qpos; // Set qstart if not set already
                        qstart_set = true;
                    }
                    current_qpos += num;  // Add match length to query position
                    current_tpos += num;  // Add match length to target position
                    matches += num;  // Assume all M are matches until adjusted by NM:i:
                    num_bases_aligned += num;
                } else if (ch == 'I') {
                    if (!qstart_set) {
                        qstart = current_qpos; // Set qstart if not set already
                        qstart_set = true;
                    }
                    current_qpos += num;  // Add insert length only to query position
                    qnuminsert += num;
                } else if (ch == 'D') {
                    current_tpos += num; // Add delete length to target position
                    tnuminsert += num;
                } else if (ch == 'S' || ch == 'H') {
                    // Softclip or hardclip: move qpos if before M
                    if (is_softclip_or_hardclip_in_front) {
                        current_qpos += num;
                    }
                }
                num = 0; // Reset num after processing each operation
            }
        }

        // Calculate qend
        qend = current_qpos;
        tend = current_tpos;

        // Adjust matches and mismatch based on NM:i: (if no 'X' operator is present)
        if (sam.cigar.find('X') == std::string::npos) {
            mismatch = mismatch - qnuminsert - tnuminsert;  // Remove insertions and deletions from NM count
            matches -= mismatch;  // Subtract mismatches from matches
        }


    // Step 5: Calculate identity and score
    identity = setIdentity(qstart, qend, tstart, tend, qnuminsert, mismatch, matches);
    score = setScore(matches, mismatch, qnuminsert, tnuminsert, query_len);
}






    void setNumBasesAligned();
    std::string details() const;
    std::string gff(const std::string &type) const;
    std::string psl() const;
    std::string exon() const;

    void correctBlocks(const std::unordered_map<std::string, bool> &splice_motifs, const std::string &target_seq, const std::string &query_seq);
    void correctSingleGaps(const std::unordered_map<std::string, bool> &splice_motifs, const std::string &target_seq);
    void correctNeighborGaps(const std::unordered_map<std::string, bool> &splice_motifs, const std::string &target_seq);
    void correctUnaligned(const std::string &query_seq, const std::string &target_seq, int max_diff = 5);


    double setIdentity(int qstart, int qend, int tstart, int tend, int qnuminsert, int mismatch, int num_bases_aligned);

    double calcIdentity(int qstart, int qend, int tstart, int tend, int qnuminsert, int mismatch, int num_bases_aligned);

    int setScore(int match, int mismatch, int qnuminsert, int tnuminsert, int query_len);

    int calcScore(int match, int qnuminsert, int tnuminsert, int query_len);


    int numExonsOverlapped() const;

    // Friend function: output alignment information
    friend std::ostream& operator<<(std::ostream &os, const alignment_t &alignment) {
        os << "Method: " << alignment.method_ << "\n"
           << "Query: " << alignment.query << "\n"
           << "Target: " << alignment.target << "\n"
           << "Query Length: " << alignment.query_len << "\n"
           << "Target Length: " << alignment.target_len << "\n"
           << "Query Strand: " << alignment.query_strand << "\n"
           << "Query Start: " << alignment.qstart << "\n"
           << "Query End: " << alignment.qend << "\n"
           << "Target Start: " << alignment.tstart << "\n"
           << "Target End: " << alignment.tend << "\n"
           << "Number of Bases Aligned: " << alignment.num_bases_aligned << "\n"
           << "Mismatch: " << alignment.mismatch << "\n"
           << "QNumInsert: " << alignment.qnuminsert << "\n"
           << "TNumInsert: " << alignment.tnuminsert << "\n"
           << "Matches: " << alignment.matches << "\n"
           << "Repmatch: " << alignment.repmatch << "\n"
           << "TBaseInsert: " << alignment.tbaseinsert << "\n"
           << "QBaseInsert: " << alignment.qbaseinsert << "\n"
           << "identity" << alignment.identity << "\n"
           << "Blocks: " << alignment.blockcount << "\n";
        for (const auto &block : alignment.blocks) {
            os << "(" << block.first << ", " << block.second << ") ";
        }
        os << "\nQuery Blocks: ";
        for (const auto &q_block : alignment.query_blocks) {
            os << "(" << q_block.first << ", " << q_block.second << ") ";
        }
        os << "\n";
        return os;
    }






    std::string fixSingleGap(std::pair<int, int> &tblock1, std::pair<int, int> &tblock2, std::pair<int, int> &qblock1, std::pair<int, int> &qblock2, const std::unordered_map<std::string, bool> &splice_motifs, const std::string &target_seq, char query_strand);
    std::string fixNeighborGaps(std::pair<int, int> &tblock1, std::pair<int, int> &tblock2, std::pair<int, int> &tblock3, std::pair<int, int> &qblock1, std::pair<int, int> &qblock2, std::pair<int, int> &qblock3, const std::unordered_map<std::string, bool> &splice_motifs, const std::string &target_seq, char query_strand);
};


std::unordered_map<std::string, int> cleanFilters(const std::unordered_map<std::string, int> &filters);

std::vector<int> splitIndices(const std::string& indices);

std::string ReverseComplement(const std::string &sequence);



// Calculate the end positions of alignment blocks within a PslEntry based on the start positions and sizes.
void calculate_ends(psl_t& entry);

// Index PSL entries by their query name for quicker access and manipulation within data processing routines.
std::unordered_map<std::string, std::vector<alignment_t>> index_by_qname(const std::vector<alignment_t>& entries);


#endif // ALIGNMENT_H
