//
// Created by xinwei on 5/28/24.
//

#ifndef PSL_H
#define PSL_H

#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <algorithm>


// PSL Entry structure
struct psl_t {

public:
    int matches;
    int misMatches;
    int repMatches;
    int nCount;
    int qNumInsert;
    int qBaseInsert;
    int tNumInsert;
    int tBaseInsert;
    std::string strand;
    std::string qName;
    int qSize;
    int qStart;
    int qEnd;
    std::string tName;
    int tSize;
    int tStart;
    int tEnd;
    int blockCount;
    std::vector<int> blockSizes;
    std::vector<int> qStarts;
    std::vector<int> tStarts;
    std::vector<int> qEnds;
    std::vector<int> tEnds;
};


void psl_parse(const std::string &filename, std::vector<psl_t>& psls);

// Parse a single line from a PSL file to extract alignment information and populate a PslEntry structure.
void parse_line(const std::string& line, psl_t& psl);

// Validate a parsed PslEntry to ensure all fields, such as block counts, are correctly formatted and logical.
void validate_entry(const psl_t& psl);

// Extract block sizes and start positions from strings, updating the given PslEntry with these details.
void get_blocks(psl_t& entry, const std::string& blockSizes, const std::string& qStarts, const std::string& tStarts);

//
std::unordered_map<std::string, int> count_qnames(const std::vector<psl_t>& psls);


#endif // PSL_H