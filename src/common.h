//
// Created by xinwei on 5/31/24.
//

#ifndef FUSION_DETECTION_2_COMMON_H
#define FUSION_DETECTION_2_COMMON_H

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <list>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <vector>

typedef unsigned char filter_t;
class filters_t: public std::vector<std::string> {
public: filter_t define(const std::string& filter_name) { push_back(filter_name); return size() - 1; };
};
static filters_t FILTERS;
const filter_t FILTER_none = FILTERS.define("");
const filter_t FILTER_duplicates = FILTERS.define("duplicates");
const filter_t FILTER_long_gap = FILTERS.define("long_gap");

const filter_t FILTER_homopolymer = FILTERS.define("internal_tandem_duplication");

const filter_t FILTER_same_gene = FILTERS.define("mt");
const filter_t FILTER_small_insert_size = FILTERS.define("small_insert_size");


const filter_t FILTER_small_fragments = FILTERS.define("small_fragments");
const filter_t FILTER_edge_unaligned = FILTERS.define("edge_unaligned");

const filter_t FILTER_min_support = FILTERS.define("min_support");
const filter_t FILTER_homologs = FILTERS.define("homologs");

// when more than 64 filters are added, the size of the filter member of the fusion_t class needs to be enlarged

inline bool str_to_int(const char* s, int& i) {
    char* end_of_parsing;
    long int result = strtol(s, &end_of_parsing, 10);
    i = result;
    return (*s != ' ' && end_of_parsing != s && *end_of_parsing == '\0' && result != LONG_MAX && result != LONG_MIN);
}

inline bool str_to_float(const char* s, float& f) {
    char* end_of_parsing;
    f = strtof(s, &end_of_parsing);
    return (*s != ' ' && end_of_parsing != s && *end_of_parsing == '\0' && f != HUGE_VALF && f != -HUGE_VALF);
}


typedef int position_t;
typedef char strandedness_t;
typedef bool strand_t;
const strand_t FORWARD = true;
const strand_t REVERSE = false;
typedef unsigned char confidence_t;
const confidence_t CONFIDENCE_LOW = 0;
const confidence_t CONFIDENCE_MEDIUM = 1;
const confidence_t CONFIDENCE_HIGH = 2;


typedef short unsigned int contig_t;


struct annotation_record_t {
    contig_t contig;
    position_t start;
    position_t end;
    strand_t strand;
    // function to sort annotation records by coordinate
    inline bool operator < (const annotation_record_t& x) const {
        if (contig != x.contig) return contig < x.contig;
        if (end != x.end) return end < x.end;
        return start < x.start;
    }
    void copy(const annotation_record_t& x) { *this = x; }
    inline unsigned int length() const { return this->end - this->start; }
};
template <class T> class annotation_set_t: public std::vector<T> {
public:
    typename annotation_set_t<T>::iterator insert(const T& value) {
        typename annotation_set_t<T>::iterator existing_element = lower_bound(this->begin(), this->end(), value);
        if (existing_element == this->end() || *existing_element != value)
            return this->insert(upper_bound(this->begin(), this->end(), value), value);
        else
            return existing_element;
    };
    void insert(typename annotation_set_t<T>::const_iterator first, typename annotation_set_t<T>::const_iterator last) {
        this->reserve(this->size() + distance(first, last));
        for (auto annotation_record = first; annotation_record != last; ++annotation_record)
            this->insert(*annotation_record);
    };
    using std::vector<T>::insert;
};

template <class T> class contig_annotation_index_t: public std::map< position_t, annotation_set_t<T> > {};
template <class T> class annotation_index_t: public std::vector< contig_annotation_index_t<T> > {};

struct gene_annotation_record_t: public annotation_record_t {
    unsigned int id; // ID used internally
    std::string gene_id; // ID specified in the GTF file
    std::string name;
    int exonic_length; // sum of the length of all exons in a gene
    bool is_dummy;
    bool is_protein_coding;
};
typedef gene_annotation_record_t* gene_t;


struct fusion_t {
    // the following members are ordered for mininum struct size

    short unsigned int splitReadsCount:15;

    short unsigned int spanReadsCount:15;

    strand_t predicted_strand1:1, predicted_strand2:1;

    confidence_t confidence:2;

    bool predicted_strands_ambiguous:1;
    short unsigned int discordant_mates:15;
    contig_t contig1, contig2;
    float evalue; // expected number of fusions with the given properties by random chance
    position_t breakpoint1, breakpoint2;
    position_t anchor_start1, anchor_start2;
    position_t closest_genomic_breakpoint1, closest_genomic_breakpoint2;
    gene_t gene1, gene2;

    fusion_t(): splitReadsCount(0),  spanReadsCount(0), predicted_strand1(FORWARD), predicted_strand2(FORWARD), confidence(CONFIDENCE_LOW), predicted_strands_ambiguous(true), discordant_mates(0), contig1(USHRT_MAX), contig2(USHRT_MAX), evalue(0), breakpoint1(-1), breakpoint2(-1), anchor_start1(0), anchor_start2(0), closest_genomic_breakpoint1(-1), closest_genomic_breakpoint2(-1), gene1(NULL), gene2(NULL) {};
    inline unsigned int supporting_reads() const { return splitReadsCount + spanReadsCount + discordant_mates; };
    bool breakpoint_overlaps_both_genes(const unsigned int which_breakpoint = 0) const {
        if (which_breakpoint == 1) return breakpoint1 >= gene2->start && breakpoint1 <= gene2->end;
        if (which_breakpoint == 2) return breakpoint2 >= gene1->start && breakpoint2 <= gene1->end;
        return breakpoint_overlaps_both_genes(1) || breakpoint_overlaps_both_genes(2);
    };

};

// convenience function to print an error message and exit if given condition is true
#define crash(condition,message) { if (condition) { std::cerr << std::string("ERROR: ") + message << std::endl; exit(1); }; }

#endif //FUSION_DETECTION_2_COMMON_H
