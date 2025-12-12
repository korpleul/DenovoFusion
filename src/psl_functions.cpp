//
// Created by xinwei on 6/6/24.
//

#include "psl_functions.h"

void CalcOverlap(int startA, int endA, int startB, int endB) {
    int leftA = std::min(startA, endA);
    int rightA = std::max(startA, endA);
    int leftB = std::min(startB, endB);
    int rightB = std::max(startB, endB);

    int overlap_left = std::max(leftA, leftB);
    int overlap_right = std::min(rightA, rightB);
    int overlap = (overlap_right - overlap_left) + 1;

    if (overlap > 0) {
        int boardA = std::min(leftA, leftB);
        int boardB = std::max(rightA, rightB);
        int total_span = boardB - boardA + 1;
        double overlap_fraction = static_cast<double>(overlap) / total_span;

        // Perform operations or print the results directly
        std::cout << "Overlap: " << overlap << ", Overlap fraction: " << overlap_fraction << std::endl;
    } else {
        std::cout << "No overlap" << std::endl;
    }
}

void CalcGap(int startA, int endA, int startB, int endB) {
    int leftA = std::min(startA, endA);
    int rightA = std::max(startA, endA);
    int leftB = std::min(startB, endB);
    int rightB = std::max(startB, endB);

    int gap = 0;
    if (rightA < leftB) {
        gap = leftB - rightA - 1;
    } else if (rightB < leftA) {
        gap = leftA - rightB - 1;
    }

    if (gap > 0) {
        int total_span = (rightA - leftA) + 1 + (rightB - leftB) + 1;
        double gap_fraction = static_cast<double>(gap) / total_span;

        // Perform operations or print the results directly
        std::cout << "Gap: " << gap << ", Gap fraction: " << gap_fraction << std::endl;
    } else {
        std::cout << "No gap" << std::endl;
    }
}

void CalcAlignOverlap(const CoordPair& coord1, const CoordPair& coord2) {
    CalcOverlap(coord1.start, coord1.end, coord2.start, coord2.end);
}

void CalcAlignGap(const CoordPair& coord1, const CoordPair& coord2) {
    CalcGap(coord1.start, coord1.end, coord2.start, coord2.end);
}


void CompareCoordPairsInGroups(std::unordered_map<std::string, std::vector<CoordPair>>& coord_pairs_start_end) {
    for (auto& group : coord_pairs_start_end) {
        const std::string& group_name = group.first;
        std::vector<CoordPair>& coord_pairs = group.second;

        std::cout << "Processing group: " << group_name << std::endl;
        // Handle groups with only one CoordPair
        if (coord_pairs.size() == 1) {
            std::cout << "This is a single alignment in this group " << std::endl;
            continue;  // Skip to the next group
        }

        for (size_t i = 0; i < coord_pairs.size(); ++i) {
            for (size_t j = i + 1; j < coord_pairs.size(); ++j) {
                CoordPair& coord1 = coord_pairs[i];
                CoordPair& coord2 = coord_pairs[j];

                // Check if the two CoordPairs are adjacent
                if (coord1.Contains(coord2)) {
                    std::cout << coord1.ToString() << " contains " << coord2.ToString() << std::endl;
                } else if (coord2.Contains(coord1)) {
                    std::cout << coord2.ToString() << " contains " << coord1.ToString() << std::endl;
                } else if (coord1.Overlaps(coord2)) {
                    CalcAlignOverlap(coord1, coord2);
                    std::cout << coord1.ToString() << " overlaps with " << coord2.ToString() << std::endl;
                } else if (coord1.Gap(coord2)){
                    CalcAlignGap(coord1, coord2);
                    std::cout << coord1.ToString() << " gaps with " << coord2.ToString() << std::endl;
                }

                else {
                    std::cout << coord1.ToString() << " is separate from " << coord2.ToString() << std::endl;
                }
            }
        }
    }
}

// Constant for the target name column index, adjust according to your specific PSL format
const int PSL_TNAME_COL = 13; // Example index, usually 13 in a PSL format

// A mock function to format chromosome names
std::string FormatChromosomeName(const std::string& name, bool use_chr, bool use_mt) {
    std::string formatted_name = name;
    if (use_chr && name.substr(0, 3) != "chr") {
        formatted_name = "chr" + name;
    }
    if (use_mt && (name == "M" || name == "MT")) {
        formatted_name = "chrM"; // Example conversion
    }
    return formatted_name;
}

std::string FormatChromosomeNameInPsl(const std::string& psl_str, bool use_chr, bool use_mt) {
    std::istringstream iss(psl_str);
    std::vector<std::string> psl_info(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

    if (psl_info.size() > PSL_TNAME_COL) {
        std::string target_name = psl_info[PSL_TNAME_COL];
        target_name = FormatChromosomeName(target_name, use_chr, use_mt);
        psl_info[PSL_TNAME_COL] = target_name;
    }

    std::ostringstream oss;
    std::copy(psl_info.begin(), psl_info.end(), std::ostream_iterator<std::string>(oss, "\t"));
    std::string formatted_psl_str = oss.str();
    if (!formatted_psl_str.empty()) {
        formatted_psl_str.pop_back(); // Remove the last tab added by std::copy
    }

    return formatted_psl_str;
}

AlignBlockCls::AlignBlockCls(const std::pair<int, int>& block_coords_tuple) {
    // Determine the strand based on the coordinates order
    if (block_coords_tuple.first <= block_coords_tuple.second) {
        strand = '+';
    } else {
        strand = '-';
    }

    // Sort the start and end using std::min and std::max
    start = std::min(block_coords_tuple.first, block_coords_tuple.second);
    end = std::max(block_coords_tuple.first, block_coords_tuple.second);
    span = (end - start) + 1;
}


int AlignBlockCls::Gap(const AlignBlockCls& other) const {
    return std::min(other.start - end, start - other.end) - 1;
}

