//
// Created by xinwei on 12/5/25.
//

#ifndef SUPPORT_WRITING_H
#define SUPPORT_WRITING_H
#include <fstream>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include "sam.h"
#include "output_fusions.h"

void writeFastqRecord(std::ofstream& out, const sam_t& read, const std::string& contigName);

void writeEvidenceToFastq(
    const std::vector<result_t>& finalResults, // Input: Final fusion list for filtering
    const std::unordered_map<std::string, std::vector<sam_t>>& splitReadsMap,
    const std::unordered_map<std::string, std::vector<sam_t>>& spanReadsMap,
    const std::string& outputPrefix,
    const options_t& options);



#endif //SUPPORT_WRITING_H
