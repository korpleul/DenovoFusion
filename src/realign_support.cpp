//
// Created by xinwei on 6/19/24.
//

#include <vector>
#include <string>
#include <unordered_set>
#include <fstream>

#include "realign_support.h"
#include "options.h"


bool isReadSupportingOverlap(const sam_t& read, const OverlapResultCls& overlap, const options_t& options) {
    // Calculate the end position of the read sequence on the reference sequence
    int readEnd = read.pos + read.seq.length();
    // Assume that the CIGAR string represents a complete match (in practice, the CIGAR string needs to be parsed to calculate the exact end position)

    // Checks whether the read sequence is at least min_edge_length away from the edge of the overlap region
    bool isStartFarEnough = read.pos <= overlap.start_ - options.min_edge_length;
    bool isEndFarEnough = readEnd >= overlap.end_ + options.min_edge_length;

    // Check if the read sequence covers the overlapping region
    bool isCoveringOverlap = read.pos <= overlap.end_ && readEnd >= overlap.start_;

    return isCoveringOverlap && isStartFarEnough && isEndFarEnough;
}





// Initializing global variables
std::unordered_map<std::string, std::vector<OverlapResultCls>> overlapMap;
std::unordered_map<std::string, std::vector<sam_t>> samMap;

// Filling overlapMap
void fillOverlapMap(const std::vector<OverlapResultCls>& overlaps) {
    for (const auto& overlap : overlaps) {
        overlapMap[overlap.query_id_].push_back(overlap);
    }
}

// Populate samMap, organizing reads by rname (corresponding to query_id)
void fillSamMap(const std::vector<sam_t>& samEntries) {
    for (const auto& sam : samEntries) {
        samMap[sam.rname].push_back(sam);
    }
}



std::unordered_map<std::string, int> countSplitReads(
        const std::unordered_map<std::string, std::vector<OverlapResultCls>>& overlapMap,
        const std::unordered_map<std::string, std::vector<sam_t>>& samMap,
        const options_t& options) {

    std::unordered_map<std::string, int> splitReadsCount;

    for (const auto& overlapEntry : overlapMap) {
        const std::string& queryName = overlapEntry.first;
        const std::vector<OverlapResultCls>& overlaps = overlapEntry.second;

        if (samMap.find(queryName) != samMap.end()) {
            const std::vector<sam_t>& reads = samMap.at(queryName);
            int count = 0;

            for (const auto& overlap : overlaps) {
                for (const auto& read : reads) {
                    if (isReadSupportingOverlap(read, overlap, options)) {
                        ++count;
                    }
                }
            }

            splitReadsCount[queryName] = count;
        }
    }

    return splitReadsCount;
}




bool isSpanningPair(const sam_t& read1, const sam_t& read2, const OverlapResultCls& overlap) {
    int read1End = read1.pos + read1.seq.length();
    int read2End = read2.pos + read2.seq.length();

    // Check if read1 and read2 are within the contig limits
    bool read1WithinContig = read1.pos >= overlap.getContigStart() && read1End <= overlap.getContigEnd();
    bool read2WithinContig = read2.pos >= overlap.getContigStart() && read2End <= overlap.getContigEnd();


    // If either read is outside of contig boundaries, return false
    if (!read1WithinContig || !read2WithinContig) {
        return false;
    }

    // Check if read1 and read2 are on different sides of the overlap region
    bool read1Left = read1.pos < overlap.getStart();
    bool read2Left = read2.pos < overlap.getStart();

    bool read1Right = read1End > overlap.getEnd();
    bool read2Right = read2End > overlap.getEnd();

    // Ensure reads are exactly on either side of the breakpoints area and do not overlap
    if ((read1Left && read2Right) || (read2Left && read1Right)) {
        // Further ensure that the reads do not overlap the breakpoints region
        if ((read1End <= overlap.getStart() && read2.pos >= overlap.getEnd()) ||
            (read2End <= overlap.getStart() && read1.pos >= overlap.getEnd())) {
            return true;
            }
    }

    return false;
}



std::unordered_map<std::string, int> countSpanReadPairs(
        const std::unordered_map<std::string, std::vector<OverlapResultCls>>& overlapMap,
        const std::unordered_map<std::string, std::vector<sam_t>>& samMap) {

    std::unordered_map<std::string, int> spanPairsCount;

    for (const auto& overlapEntry : overlapMap) {
        const std::string& queryName = overlapEntry.first;
        const std::vector<OverlapResultCls>& overlaps = overlapEntry.second;

        if (samMap.find(queryName) != samMap.end()) {
            const std::vector<sam_t>& reads = samMap.at(queryName);
            int count = 0;

            // Use a mapping to group readings by qname
            std::unordered_map<std::string, std::vector<sam_t>> groupedReads;
            for (const auto& read : reads) {
                groupedReads[read.qname].push_back(read);
            }

            // Iterate over each group and compare only reads with the same qname
            for (const auto& group : groupedReads) {
                const std::vector<sam_t>& pairedReads = group.second;
                // At least two readings are required to form a pair
                for (size_t i = 0; i < pairedReads.size(); ++i) {
                    for (size_t j = i + 1; j < pairedReads.size(); ++j) {
                        for (const auto& overlap : overlaps) {
                            if (isSpanningPair(pairedReads[i], pairedReads[j], overlap)) {
                                count++;  // Each time a valid pair is found, the count increases by 1.
                            }
                        }
                    }
                }
            }

            spanPairsCount[queryName] = count;
        }
    }

    return spanPairsCount;
}



std::unordered_set<std::string> filterQueries(
        const std::unordered_map<std::string, int>& splitReadsCount,
        const std::unordered_map<std::string, int>& spanPairsCount) {

    std::unordered_set<std::string> validQueries;

    for (const auto& entry : splitReadsCount) {
        if (entry.second > 0 || spanPairsCount.at(entry.first) > 0) {
            validQueries.insert(entry.first);
        }
    }

    return validQueries;
}

std::vector<alignment_t> filterAlignmentsByValidBaseQueries(
        const std::vector<alignment_t>& alignments,
        const std::unordered_set<std::string>& validBaseQueries) {

    std::vector<alignment_t> filteredAlignments;
    for (const alignment_t& alignment : alignments) {
        std::string baseQuery = baseQueryName(alignment.query);
        if (validBaseQueries.find(baseQuery) != validBaseQueries.end()) {
            filteredAlignments.push_back(alignment);
        }
    }

    return filteredAlignments;
}

std::vector<alignment_t> filterAlignmentsByValidQueries(
        const std::vector<alignment_t>& alignments,
        const std::unordered_set<std::string>& validQueries) {

    std::vector<alignment_t> filteredAlignments;

    for (const alignment_t& alignment : alignments) {
        // Assumes that query has at least one underscore.
        if (validQueries.find(alignment.query) != validQueries.end()) {
            filteredAlignments.push_back(alignment);
        }
    }

    return filteredAlignments;
}


std::string extractBaseQueryName(const std::string& query) {
    size_t pos = query.rfind('_');
    return (pos != std::string::npos) ? query.substr(0, pos) : query;
}



std::vector<coordination_t> extractCoordinations(const std::vector<alignment_t>& alignments) {
    std::vector<coordination_t> coordinations;
    for (const auto& align : alignments) {
        std::string baseQuery = extractBaseQueryName(align.query);
        coordinations.push_back(coordination_t(baseQuery, align.target, align.tstart, align.tend, std::string(1, align.query_strand)));
    }
    return coordinations;
}

std::vector<coordination_t> mergeContinuousSegments(const std::vector<coordination_t>& segments) {
    if (segments.empty()) return {};

    std::vector<coordination_t> merged;
    coordination_t current = segments[0];

    for (int i = 1; i < segments.size(); ++i) {
        // Check if segments belong to the same target
        if (segments[i].target == current.target) {
            // Allow a gap of up to 5 bases for merging
            if (std::abs(segments[i].tstart - current.tend) <= 3 ||
                std::abs(segments[i].tend - current.tstart) <= 3) {

                // Extend the current segment based on the new one
                current.tstart = std::min(current.tstart, segments[i].tstart);
                current.tend = std::max(current.tend, segments[i].tend);

                // Preserve strand consistency, or set to "/" if they differ
                if (current.strand != segments[i].strand) {
                    current.strand = "/";
                }
                } else {
                    // If not close enough, save the current segment and start a new one
                    merged.push_back(current);
                    current = segments[i];
                }
        } else {
            // If targets differ, save the current segment and start a new one
            merged.push_back(current);
            current = segments[i];
        }
    }

    // Add the last segment
    merged.push_back(current);

    return merged;
}


std::vector<coordination_t> extractSimpleCoordinations(const std::vector<alignment_t>& alignments) {
    std::vector<coordination_t> coordinations;

    for (const auto& align : alignments) {
        // Create a coordination_t object directly using the information from the alignment
        coordinations.push_back(coordination_t(align.query, align.target, align.tstart, align.tend, std::string(1, align.query_strand)));
    }

    return coordinations;
}




void writeFastqRecord(std::ofstream& out, const sam_t& read) {
    // FASTQ standard
    // Line 1: @ReadName
    // Line 2: Sequence
    // Line 3: +
    // Line 4: Quality Score

    out << "@" << read.qname << "\n"
        << read.seq << "\n"
        << "+\n";


    if (read.qual.empty() || read.qual == "*") {
        out << std::string(read.seq.length(), 'I') << "\n";
    } else {
        out << read.qual << "\n";
    }
}



std::unordered_map<std::string, std::vector<sam_t>> collectSplitReads(
    const std::unordered_map<std::string, std::vector<OverlapResultCls>>& overlapMap,
    const std::unordered_map<std::string, std::vector<sam_t>>& samMap,
    const options_t& options) {

    std::unordered_map<std::string, std::vector<sam_t>> splitReadsMap;

    for (const auto& overlapEntry : overlapMap) {
        const std::string& queryName = overlapEntry.first;
        const std::vector<OverlapResultCls>& overlaps = overlapEntry.second;

        if (samMap.find(queryName) != samMap.end()) {
            const std::vector<sam_t>& reads = samMap.at(queryName);
            for (const auto& overlap : overlaps) {
                for (const auto& read : reads) {
                    if (isReadSupportingOverlap(read, overlap, options)) {
                        splitReadsMap[queryName].push_back(read);
                    }
                }
            }
        }
    }
    return splitReadsMap;
}

std::unordered_map<std::string, std::vector<sam_t>> collectSpanReads(
    const std::unordered_map<std::string, std::vector<OverlapResultCls>>& overlapMap,
    const std::unordered_map<std::string, std::vector<sam_t>>& samMap,
    const options_t& options) {

    // Output container: QueryName -> List of Spanning Reads
    std::unordered_map<std::string, std::vector<sam_t>> spanReadsMap;

    for (const auto& overlapEntry : overlapMap) {
        const std::string& queryName = overlapEntry.first; // The Contig Name
        const std::vector<OverlapResultCls>& overlaps = overlapEntry.second;

        // Check if this contig has mapped reads
        if (samMap.find(queryName) != samMap.end()) {
            const std::vector<sam_t>& reads = samMap.at(queryName);

            // ---------------------------------------------------------
            // Step 1: Group reads by Read Name (QNAME) to reconstruct pairs
            // ---------------------------------------------------------
            // We use a temporary map to buffer potential pairs
            std::unordered_map<std::string, std::vector<sam_t>> pairBuffer;
            for (const auto& read : reads) {
                pairBuffer[read.qname].push_back(read);
            }

            // ---------------------------------------------------------
            // Step 2: Evaluate pairs against the overlap region
            // ---------------------------------------------------------
            for (const auto& overlap : overlaps) {
                for (const auto& pairEntry : pairBuffer) {
                    const std::vector<sam_t>& pair = pairEntry.second;

                    // We only evaluate valid pairs (exactly 2 reads with the same name)
                    if (pair.size() == 2) {
                        const sam_t& read1 = pair[0];
                        const sam_t& read2 = pair[1];

                        // Use your custom strict logic
                        if (isSpanningPair(read1, read2, overlap)) {
                            // If they are a spanning pair, add BOTH to the output
                            spanReadsMap[queryName].push_back(read1);
                            spanReadsMap[queryName].push_back(read2);
                        }
                    }
                }
            }
        }
    }

    return spanReadsMap;
}