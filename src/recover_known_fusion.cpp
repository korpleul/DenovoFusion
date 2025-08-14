//
// Created by xinwei on 10/1/24.
//

#include "recover_known_fusion.h"


// Function to parse known fusions from the TSV file
std::vector<known_fusion_t> load_known_fusions(const std::string& filename) {
    std::vector<known_fusion_t> known_fusions;
    std::ifstream file(filename);

    if (!file.is_open()) {
        return known_fusions;
    }

    std::string line;
    known_fusion_t current_fusion;

    while (std::getline(file, line)) {
        // Check if the line starts with '#', indicating a new gene pair
        if (!line.empty() && line[0] == '#') {
            // Save the previous fusion if it has valid data
            if (!current_fusion.gene1.empty() && !current_fusion.coordinates1.empty() && !current_fusion.coordinates2.empty()) {
                known_fusions.push_back(current_fusion);
            }

            // Start a new known fusion record
            std::istringstream ss(line.substr(1));  // Skip the '#'
            ss >> current_fusion.gene1 >> current_fusion.gene2;

            // Clear previous coordinate data for the new fusion
            current_fusion.coordinates1.clear();
            current_fusion.coordinates2.clear();
            current_fusion.source.clear();
        } else if (!line.empty()) {
            // Parse the coordinates and source in the current line
            std::istringstream ss(line);
            std::string coord1_str, coord2_str, source;

            if (!(ss >> coord1_str >> coord2_str >> source)) {
                std::cerr << "Failed to parse line: " << line << std::endl;
                continue;
            }

            coordinate_t coord1, coord2;

            // Parse first coordinate
            std::istringstream coord1_ss(coord1_str);
            coord1_ss >> coord1.strand >> coord1.chromosome;
            coord1.chromosome = coord1.chromosome.substr(1);  // Remove leading strand symbol
            coord1_ss.ignore(1, ':');
            coord1_ss >> coord1.start;
            coord1_ss.ignore(1, '-');
            coord1_ss >> coord1.end;

            // Parse second coordinate
            std::istringstream coord2_ss(coord2_str);
            coord2_ss >> coord2.strand >> coord2.chromosome;
            coord2.chromosome = coord2.chromosome.substr(1);  // Remove leading strand symbol
            coord2_ss.ignore(1, ':');
            coord2_ss >> coord2.start;
            coord2_ss.ignore(1, '-');
            coord2_ss >> coord2.end;

            // Store the parsed coordinates and source
            current_fusion.coordinates1.push_back(coord1);
            current_fusion.coordinates2.push_back(coord2);
            current_fusion.source = source;
        }
    }

    // Save the last fusion after loop ends
    if (!current_fusion.gene1.empty() && !current_fusion.coordinates1.empty() && !current_fusion.coordinates2.empty()) {
        known_fusions.push_back(current_fusion);
    }

    return known_fusions;
}


// Function to recover known fusions from discarded results
std::vector<result_t> recover_fusions(const std::vector<result_t>& discarded_results, const std::vector<known_fusion_t>& known_fusions) {
    std::vector<result_t> recovered; // To store recovered results

    for (const auto& known_fusion : known_fusions) {
        for (const auto& result : discarded_results) {
            bool recovered_flag = false; // Flag to indicate if a match is found

            // Check each coordinate pair for gene1 and gene2
            for (const auto& coord1 : known_fusion.coordinates1) {
                for (const auto& coord2 : known_fusion.coordinates2) {
                    // Check if genes and positions match (direct case)
                    bool matches_direct =
                        (result.gene1 == known_fusion.gene1 && result.gene2 == known_fusion.gene2 &&
                         result.tstart1 <= coord1.end && result.tend1 >= coord1.start &&
                         result.tstart2 <= coord2.end && result.tend2 >= coord2.start);

                    // Check if genes and positions match (reversed case)
                    bool matches_reversed =
                        (result.gene1 == known_fusion.gene2 && result.gene2 == known_fusion.gene1 &&
                         result.tstart1 <= coord2.end && result.tend1 >= coord2.start &&
                         result.tstart2 <= coord1.end && result.tend2 >= coord1.start);

                    // If either direct or reversed match is found, mark as recovered
                    if (matches_direct || matches_reversed) {
                        result_t recovered_result = result;
                        recovered_result.filter_status = "recovered"; // Update the filter status
                        recovered.push_back(recovered_result);
                        recovered_flag = true;
                        break;  // Stop checking further coordinates for this result
                    }
                }

                if (recovered_flag) {
                    break;  // Stop checking further coordinates for this known fusion
                }
            }
        }
    }

    return recovered; // Return the list of recovered results
}
