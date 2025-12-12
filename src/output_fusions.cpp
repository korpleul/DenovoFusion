//
// Created by xinwei on 6/13/24.
//

#include "output_fusions.h"







void processAnnotations(const std::vector<annotation_t>& annotations, std::vector<result_t>& results) {

    for (size_t i = 0; i < annotations.size(); i += 2) {
        result_t result;
        const auto& ann1 = annotations[i];
        const auto& ann2 = annotations[i+1];

        if (ann1.contigName == ann2.contigName) {  // 同一个 contig
            result.contig = ann1.contigName;
            result.gene1 = ann1.geneName;
            result.gene2 = ann2.geneName;
            result.gene_id1 = ann1.geneId.empty() ? "" : ann1.geneId;
            result.gene_id2 = ann2.geneId.empty() ? "" : ann2.geneId;
            result.tstart1 = ann1.start;
            result.tend1 = ann1.end;
            result.tstart2 = ann2.start;
            result.tend2 = ann2.end;
            result.chromosome1 = ann1.chromosome;
            result.chromosome2 = ann2.chromosome;
            result.tstrand1 = ann1.strand;
            result.tstrand2 = ann2.strand;
            result.direction1 = ann1.direction;
            result.direction2 = ann2.direction;
            result.regiontype1 = ann1.regionType;
            result.regiontype2 = ann2.regionType;
            result.filter_status = "unfiltered";
            results.push_back(result);
        }
    }
}



void integrateReadCounts(
        std::vector<result_t>& results,
        const std::unordered_map<std::string, int>& splitReadsCount,
        const std::unordered_map<std::string, int>& spanReadsCount) {

    for (auto& result : results) {
        if (splitReadsCount.find(result.contig) != splitReadsCount.end()) {
            result.splitReadsCount = splitReadsCount.at(result.contig);
        }
        if (spanReadsCount.find(result.contig) != spanReadsCount.end()) {
            result.spanReadsCount = spanReadsCount.at(result.contig);
        }
    }
}



// Function to filter out results where gene1 or gene2 is empty
void removeEmptyGenes(std::vector<result_t>& results) {
    results.erase(std::remove_if(results.begin(), results.end(),
                                 [](const result_t& result) {
                                     return result.gene_id1.empty() || result.gene_id2.empty() ||
                                           result.gene_id1 == " " || result.gene_id2 == " ";
                                }),
                  results.end());
}

void writeToCSV(const std::vector<result_t>& final_results, const std::string& filename) {
    // Open the file for writing
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    // Write the header
    outFile << "Contig,Gene1,Gene2,Gene_ID1,Gene_ID2,Split_Reads_Count,Span_Reads_Count,"
               "Chromosome1,Chromosome2,Direction1,Direction2,TStart1,TEnd1,TStart2,TEnd2,regiontype1,regiontype2,filter_status\n";

    // Write the data
    for (const auto& result : final_results) {
        outFile << result.contig << ","
                << result.gene1 << ","
                << result.gene2 << ","
                << result.gene_id1 << ","
                << result.gene_id2 << ","
                << result.splitReadsCount << ","
                << result.spanReadsCount << ","
                << result.chromosome1 << ","
                << result.chromosome2 << ","
                << result.direction1 << ","
                << result.direction2 << ","
                << result.tstart1 << ","
                << result.tend1 << ","
                << result.tstart2 << ","
                << result.tend2 << ","
                << result.regiontype1 << ","
                << result.regiontype2<< ","
                << result.filter_status <<  "\n";
    }

    // Close the file
    outFile.close();
}

void writeToTSV(const std::vector<result_t>& final_results, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    //
    outFile << "Contig\tGene1\tGene2\tGene_ID1\tGene_ID2\tSplit_Reads_Count\tSpan_Reads_Count\t"
               "Chromosome1\tChromosome2\tDirection1\tDirection2\tTStart1\tTEnd1\tTStart2\tTEnd2\tregiontype1\tregiontype2\tfilter_status\n";

    //
    for (const auto& result : final_results) {
        outFile << result.contig << "\t"
                << result.gene1 << "\t"
                << result.gene2 << "\t"
                << result.gene_id1 << "\t"
                << result.gene_id2 << "\t"
                << result.splitReadsCount << "\t"
                << result.spanReadsCount << "\t"
                << result.chromosome1 << "\t"
                << result.chromosome2 << "\t"
                << result.direction1 << "\t"
                << result.direction2 << "\t"
                << result.tstart1 << "\t"
                << result.tend1 << "\t"
                << result.tstart2 << "\t"
                << result.tend2 << "\t"
                << result.regiontype1 << "\t"
                << result.regiontype2 << "\t"
                << result.filter_status << "\n";
    }

    outFile.close();
}



