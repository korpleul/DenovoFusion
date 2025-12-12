//
// Created by xinwei on 5/19/24.
//

#ifndef ALIGNED_COORD_PAIR_H
#define ALIGNED_COORD_PAIR_H

#include <string>
#include <vector>
#include <stdexcept>
#include "coord_pair.h"
#include "error.h"

class AlignedCoordPairCls {
public:
    // Constructor
    AlignedCoordPairCls(const std::string& qname, const CoordPair& qcoords, const std::string& tname, const CoordPair& tcoords);

    // Public methods
    bool pos_strand() const;
    std::string ToString() const;

    std::string qname;
    CoordPair qcoords;
    std::string tname;
    CoordPair tcoords;
};

// Exception class
class AlignedCoordPairError : public MyError {
public:
    explicit AlignedCoordPairError(const std::string& message) : MyError(message) {}
};

class AlignedCoordPairListCls : public AlignedCoordPairCls {
public:
    // Constructor
    AlignedCoordPairListCls(const std::string& qname, const std::string& tname);

    // Public methods
    void add_coord_pair(CoordPair qcoords, CoordPair tcoords);

private:
    std::vector<AlignedCoordPairCls> coord_pairs;
};

#endif // ALIGNED_COORD_PAIR_H
