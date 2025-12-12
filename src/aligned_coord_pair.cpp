//
// Created by xinwei on 5/19/24.
//

#include "aligned_coord_pair.h"


// Constructor for AlignedCoordPairCls
AlignedCoordPairCls::AlignedCoordPairCls(const std::string& qname, const CoordPair& qcoords, const std::string& tname, const CoordPair& tcoords)
        : qname(qname), qcoords(qcoords), tname(tname), tcoords(tcoords) {
    if (!qcoords.pos_strand) {
        throw AlignedCoordPairError("Cannot create AlignedCoordPairCls with query coordinates: " + qcoords.ToString() + ". Query coordinates should always be considered positive strand.");
    }
    if (qcoords.Span() != tcoords.Span()) {
        throw AlignedCoordPairError("Cannot create AlignedCoordPairCls with query coordinates: " + qcoords.ToString() + " (" + std::to_string(qcoords.Span()) + "), and target coordinates: " + tcoords.ToString() + " (" + std::to_string(tcoords.Span()) + "). Coordinates must have equal spans.");
    }
}

// Determine strand based on target coordinates
bool AlignedCoordPairCls::pos_strand() const {
    return tcoords.pos_strand;
}

// Convert object to string
std::string AlignedCoordPairCls::ToString() const {
    return qname + ":" + qcoords.ToString() + "=" + tname + ":" + tcoords.ToString();
}

// Constructor for AlignedCoordPairListCls
AlignedCoordPairListCls::AlignedCoordPairListCls(const std::string& qname, const std::string& tname)
        : AlignedCoordPairCls(qname, CoordPair(), tname, CoordPair()) {}

// Add coordinate pair to the list
void AlignedCoordPairListCls::add_coord_pair(CoordPair qcoords, CoordPair tcoords) {
    AlignedCoordPairCls new_pair(qname, qcoords, tname, tcoords);

    // Check strand consistency
    if (!coord_pairs.empty() && tcoords.pos_strand != pos_strand()) {
        throw AlignedCoordPairError("Cannot add target coordinates: " + tcoords.ToString() + ". Strand does not match strand of list: " + strand_string(pos_strand()));
    }

    // Check order of coordinates
    if (!coord_pairs.empty()) {
        if (qcoords.end >= qcoords.start) {
            throw AlignedCoordPairError("Cannot add query coordinates: " + qcoords.ToString() + ", since they do not come after the query coordinates already in the list: " + qcoords.ToString());
        }
        if (pos_strand()) {
            if (tcoords.end >= tcoords.start) {
                throw AlignedCoordPairError("Cannot add target coordinates: " + tcoords.ToString() + ", since they do not come after the target coordinates already in the list: " + tcoords.ToString());
            }
        } else {
            if (tcoords.end <= tcoords.start) {
                throw AlignedCoordPairError("Cannot add target coordinates: " + tcoords.ToString() + ", since they do not come before the target coordinates already in the list: " + tcoords.ToString());
            }
        }
    }

    qcoords.Union(qcoords);
    tcoords.Union(tcoords);
    coord_pairs.push_back(new_pair);
}


