//
// Created by xinwei on 5/28/24.
//

#include "psl.h"
#include "error.h"

// Function to handle opening and line-by-line parsing of files
void psl_parse(const std::string &filename, std::vector<psl_t>& psls) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("can't open file: " + filename);
    }

    psls.clear();  // clear the vector and fill in a new vector
    std::string line;
    bool headerSkipped = false;

    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // discard the space row
        }
        // Skip header lines by checking for specific keywords or patterns
        if (!headerSkipped) {
            std::istringstream iss(line);
            std::string firstWord;
            iss >> firstWord;
            if (firstWord == "psLayout" || firstWord == "match") {
                while (getline(file, line) && line.find("-------") == std::string::npos) {
                    // Skip until the dashed line at the end of the header is encountered
                }
                headerSkipped = true;
                continue;
            }
            headerSkipped = true;
        }


        psl_t psl;
        parse_line(line, psl);

        // Adjust qEnd if it equals qSize, because of psl start end point problem itself
        if (psl.qEnd == psl.qSize) {
            psl.qEnd--; // Adjust qEnd
        }

        validate_entry(psl);    // validate the psl file
        psls.push_back(psl);    // put the file into container
    }

    file.close();  // close the file
}






void parse_line(const std::string& line, psl_t& psl) {
    std::istringstream ss(line);

    ss >> psl.matches
       >> psl.misMatches
       >> psl.repMatches
       >> psl.nCount
       >> psl.qNumInsert
       >> psl.qBaseInsert
       >> psl.tNumInsert
       >> psl.tBaseInsert
       >> psl.strand
       >> psl.qName
       >> psl.qSize
       >> psl.qStart
       >> psl.qEnd
       >> psl.tName
       >> psl.tSize
       >> psl.tStart
       >> psl.tEnd
       >> psl.blockCount;

    std::string blockSizes, qStarts, tStarts;
    ss >> blockSizes >> qStarts >> tStarts;

    get_blocks(psl, blockSizes, qStarts, tStarts);

}

void validate_entry(const psl_t &entry) {
    if (entry.qStarts.size() != entry.blockCount || entry.tStarts.size() != entry.blockCount) {
        throw PslParserError("Mismatch in block count and block starts.");
    }
}


void get_blocks(psl_t &entry, const std::string &blockSizes, const std::string &qStarts, const std::string &tStarts) {
    entry.blockSizes.clear();
    entry.qStarts.clear();
    entry.tStarts.clear();

    std::istringstream blockSizesStream(blockSizes);
    std::istringstream qStartsStream(qStarts);
    std::istringstream tStartsStream(tStarts);

    std::string blockSize, qStart, tStart;

    while (std::getline(blockSizesStream, blockSize, ',') &&
           std::getline(qStartsStream, qStart, ',') &&
           std::getline(tStartsStream, tStart, ',')) {
        entry.blockSizes.push_back(std::stoi(blockSize));
        entry.qStarts.push_back(std::stoi(qStart));
        entry.tStarts.push_back(std::stoi(tStart));
    }

    if (entry.blockSizes.size() != entry.blockCount || entry.qStarts.size() != entry.blockCount || entry.tStarts.size() != entry.blockCount) {
        throw PslParserError("Block count does not match the number of block sizes or start positions.");
    }
}


// Function to count the number of occurrences of each qName in a vector of psls
std::unordered_map<std::string, int> count_qnames(const std::vector<psl_t>& psls) {
    std::unordered_map<std::string, int> qname_counts;
    for (const auto& psl : psls) {
        qname_counts[psl.qName]++;
    }
    return qname_counts;
}
