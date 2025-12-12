//
// Created by xinwei on 1/15/25.
//

#ifndef FILTER_HOMOLOGS_H
#define FILTER_HOMOLOGS_H

#include <vector>
#include <string>
#include <unordered_map>
#include "alignment.h"
#include "output_fusions.h"

// FilterHomologs class definition
class FilterHomologs {
public:
    // Constructor
    FilterHomologs(const std::vector<result_t>& final_results,
                   const std::vector<std::pair<alignment_t, alignment_t>>& pairedAlignments,
                   const std::vector<alignment_t>& alignments);

    // Function to apply the homolog filter
    std::vector<result_t> filter_homologs();

private:
    // Helper function to count occurrences of a coordinate pair (qstart, qend) in alignments
    int countCoordinatePairs(int qstart, int qend, const std::vector<alignment_t>& alignments);

    // Member variables
    const std::vector<result_t>& final_results_; // Reference to the results
    const std::vector<std::pair<alignment_t, alignment_t>>& pairedAlignments_; // Reference to paired alignments
    const std::vector<alignment_t>& alignments_; // Reference to full alignments
};


#endif
