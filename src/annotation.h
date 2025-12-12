//
// Created by xinwei on 7/7/24.
//

#ifndef FUSION_DETECTION_2_ANNOTATION_H
#define FUSION_DETECTION_2_ANNOTATION_H
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iostream>

#include "alignment.h"
#include "realign_support.h"


struct annotation_t {
    std::string contigName;
    std::string geneId;
    std::string geneName;
    int start;
    int end;
    std::string chromosome;
    std::string strand;
    std::string direction;
    std::string regionType;// "exon", "intron", or ""



    annotation_t(const std::string& contig = "", const std::string& id = "",
                 const std::string& name = "No valid annotation", int s = -1, int e = -1,
                 const std::string& chr = "", const std::string& str = "", const std::string& dir = "",
                 const std::string& region = "", int exonNum = -1)
            : contigName(contig), geneId(id), geneName(name), start(s), end(e),
              chromosome(chr), strand(str), direction(dir), regionType(region){}
};


class GeneAnnotator {
public:
    explicit GeneAnnotator(const std::string& gtf_path);
    std::vector<annotation_t> annotateAlignments(const std::vector<coordination_t>& coordinations);

private:
    std::string gtf_path_;
    void parseAttributes(const std::string& attributes, std::unordered_map<std::string, std::string>& attr_map);
};

void calculateDirections(std::vector<annotation_t>& annotations);

std::vector<coordination_t> keepOnlyTwoParts(const std::vector<coordination_t>& coords);

void keepOnlyTwoAnnotations(std::vector<annotation_t>& annotations);

#endif //FUSION_DETECTION_2_ANNOTATION_H
