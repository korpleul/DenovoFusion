//
// Created by xinwei on 5/19/24.
//

#include "coord_pair.h"
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>

#include "error.h"

CoordPair::CoordPair(int start, int end, bool pos_strand, const std::string& name)
        : start(start), end(end), pos_strand(pos_strand), name(name) {
    min = std::min(start, end);
    max = std::max(start, end);
}

CoordPair::CoordPair(const std::string &coord_str, const std::string& name)
    :name(name){
    std::istringstream iss(coord_str);
    char delim;
    iss >> start >> delim >> end;
    if (start <= end) {
        pos_strand = true;
    } else {
        pos_strand = false;
    }
    min = std::min(start, end);
    max = std::max(start, end);
}

CoordPair::CoordPair(const CoordPair &other, const std::string& name)
        : start(other.start), end(other.end), min(other.min), max(other.max), pos_strand(other.pos_strand), name(other.name) {}

int CoordPair::Span() const {
    return (max - min) + 1;
}

void CoordPair::Union(const CoordPair &other) {
    min = std::min(min, other.min);
    max = std::max(max, other.max);
    if (pos_strand) {
        start = min;
        end = max;
    } else {
        start = max;
        end = min;
    }
}

bool CoordPair::Contains(const CoordPair &other) const {
    return (min <= other.min && max >= other.max);
}

bool CoordPair::Overlaps(const CoordPair &other) const {
    return !(other.max < min || max < other.min);
}
bool CoordPair::Gap(const CoordPair &other) const {
    return  other.max < min || max < other.min ;
}

void CoordPair::Intersect(const CoordPair &other) {
    min = std::max(min, other.min);
    max = std::min(max, other.max);
    if (max < min) {
        throw CoordPairError("coord-pairs do not intersect");
    }
    if (pos_strand) {
        start = min;
        end = max;
    } else {
        start = max;
        end = min;
    }
}

std::string CoordPair::ToString() const {
    std::ostringstream oss;
    oss << start << "-" << end;
    return oss.str();
}

CoordPair CoordPair::Copy() const {
    return CoordPair(start, end, pos_strand);
}

void CoordPair::ResortCoords() {
    std::tie(min, max) = std::minmax(start, end);
}

void CoordPair::SetMin(int new_min) {
    min = new_min;
    if (pos_strand) {
        start = new_min;
    } else {
        end = new_min;
    }
}

void CoordPair::SetMax(int new_max) {
    max = new_max;
    if (pos_strand) {
        end = new_max;
    } else {
        start = new_max;
    }
}

void CoordPair::MoveMin(int delta) {
    min += delta;
    if (pos_strand) {
        start = min;
    } else {
        end = min;
    }
}

void CoordPair::MoveMax(int delta) {
    max += delta;
    if (pos_strand) {
        end = max;
    } else {
        start = max;
    }
}

int CoordPair::operator[](int index) const {
    if (index == 0) {
        return start;
    } else if (index == 1) {
        return end;
    } else {
        throw std::out_of_range("Index out of range");
    }
}

std::ostream &operator<<(std::ostream &os, const CoordPair &coord) {
    os << coord.ToString();
    return os;
}

CoordPair BetweenCoords(const CoordPair &coords1, const CoordPair &coords2) {
    int start = std::max(coords1.min, coords2.min) - 1;
    int end = std::min(coords1.max, coords2.max) + 1;
    return CoordPair(start, end);
}

std::vector<CoordPair> CutOrExtendBlocks(const std::vector<CoordPair> &input_blocks, int target_length, bool from_left) {
    std::vector<CoordPair> output_blocks;
    int remaining = target_length;
    auto blocks = input_blocks;
    if (!from_left) {
        std::reverse(blocks.begin(), blocks.end());
    }
    for (auto &inblock : blocks) {
        if (remaining == 0) {
            break;
        }
        if (remaining > inblock.Span()) {
            output_blocks.push_back(inblock);
            remaining -= inblock.Span();
        } else {
            CoordPair part_block;
            if (from_left) {
                part_block = CoordPair(inblock.min, inblock.min + remaining - 1);
            } else {
                part_block = CoordPair(inblock.max - remaining + 1, inblock.max);
            }
            if (part_block.Span() != remaining) {
                throw CoordPairError("improper calculation of partial block");
            }
            output_blocks.push_back(part_block);
            remaining -= part_block.Span();
        }
    }
    if (remaining > 0) {
        CoordPair &last_block = output_blocks.back();
        if (from_left) {
            last_block.MoveMax(remaining);
        } else {
            last_block.MoveMin(-remaining);
        }
    }
    if (!from_left) {
        std::reverse(output_blocks.begin(), output_blocks.end());
    }
    return output_blocks;
}

std::vector<CoordPair> MergeAdjacentBlocks(const std::vector<CoordPair> &input_blocks) {
    std::vector<CoordPair> output_blocks;
    CoordPair prev_block;
    for (auto &inblock : input_blocks) {
        if (prev_block.min > inblock.min) {
            throw CoordPairError("adjacent blocks cannot be merged when the blocks are out of order");
        }
        if (output_blocks.empty() || output_blocks.back().max < inblock.min - 1) {
            output_blocks.push_back(inblock);
        } else if (output_blocks.back().max == inblock.min - 1) {
            output_blocks.back().SetMax(inblock.max);
        } else {
            throw CoordPairError("adjacent blocks cannot be merged when the blocks overlap");
        }
        prev_block = inblock;
    }
    return output_blocks;
}





bool strand_check(const std::string& strand) {
    return strand == "+";
}
std::string strand_string(bool pos_strand) {
    return pos_strand ? "+" : "-";
}
