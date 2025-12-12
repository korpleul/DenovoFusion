//
// Created by xinwei on 6/13/24.
//

#ifndef FUSION_DETECTION_2_OUTPUT_FUSIONS_H
#define FUSION_DETECTION_2_OUTPUT_FUSIONS_H


#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <iomanip>

#include "alignment.h"
#include "options.h"
#include "alignments_chosen.h"
#include "common.h"
#include "annotation.h"
#include "coverage.h"


struct result_t {
    std::string contig;
    std::string gene1;
    std::string gene2;
    std::string gene_id1;
    std::string gene_id2;
    unsigned int tstart1;
    unsigned int tend1;
    unsigned int tstart2;
    unsigned int tend2;
    std::string tstrand1;
    std::string tstrand2;
    std::string chromosome1;
    std::string chromosome2;
    unsigned int splitReadsCount;
    unsigned int spanReadsCount;
    float coverage;
    std::string direction1;
    std::string direction2;
    std::string filter_status;
    std::string regiontype1;
    std::string regiontype2;

};


void processAnnotations(const std::vector<annotation_t>& annotations, std::vector<result_t>& results);

void integrateReadCounts(
        std::vector<result_t>& results,
        const std::unordered_map<std::string, int>& splitReadsCount,
        const std::unordered_map<std::string, int>& spanReadsCount);

confidence_t determineConfidence(const result_t& result);

void removeEmptyGenes(std::vector<result_t>& results);

void writeToCSV(const std::vector<result_t>& final_results, const std::string& filename);
void writeToTSV(const std::vector<result_t>& final_results, const std::string& filename);

#endif //FUSION_DETECTION_2_OUTPUT_FUSIONS_H
