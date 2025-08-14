//
// Created by xinwei on 5/19/24.
//

#ifndef COORD_PAIR_H
#define COORD_PAIR_H

#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>

// class for coordination pairs
class CoordPair {
public:
    // initialization
    CoordPair() : start(0), end(0), min(0), max(0), pos_strand(true), name("") {}
    // declaration
    CoordPair(int start, int end, bool pos_strand = true, const std::string& name = "");
    CoordPair(const std::string &coord_str, const std::string& name);

    // create a new pair according to the existed pair
    CoordPair(const CoordPair &other, const std::string& name = "");

    // span distance from a pair
    int Span() const;
    // update the current pair and include the other pair
    void Union(const CoordPair &other);
    // check if the current pair include the other pair
    bool Contains(const CoordPair &other) const;
    // check if the current pair overlap the other pair
    bool Overlaps(const CoordPair &other) const;
    bool Gap(const CoordPair &other) const;
    // make the pair to be an intersection
    void Intersect(const CoordPair &other);
    // return the string
    std::string ToString() const;
    //copy
    CoordPair Copy() const;
    // make sure start <= end;
    void ResortCoords();
    // set values
    void SetMin(int new_min);
    void SetMax(int new_max);
    void MoveMin(int delta);
    void MoveMax(int delta);



    // Allows access to coordinate pairs using indices, possibly returning start or end depending on the index.
    int operator[](int index) const;

    //Output stream operator overloading
    friend std::ostream &operator<<(std::ostream &os, const CoordPair &coord);

    int start;
    int end;
    int min;
    int max;
    bool pos_strand;
    std::string name;
};

// Computes and returns the relationship between two coordinate pairs or a coordinate pair.
CoordPair BetweenCoords(const CoordPair &coords1, const CoordPair &coords2);
// Adjust the length of a series of coordinate pairs to meet a target length.
std::vector<CoordPair> CutOrExtendBlocks(const std::vector<CoordPair> &input_blocks, int target_length, bool from_left);
// Merge adjacent coordinate pairs
std::vector<CoordPair> MergeAdjacentBlocks(const std::vector<CoordPair> &input_blocks);
// deal with strand
bool strand_check(const std::string& strand);
std::string strand_string(bool pos_strand);





#endif // COORD_PAIR_H
