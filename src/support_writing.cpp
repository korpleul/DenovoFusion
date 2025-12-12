//
// Created by xinwei on 12/5/25.
//

#include "support_writing.h"

// Adds "|ContigName" to the read header so you can identify it in IGV
void writeFastqRecord(std::ofstream& out, const sam_t& read, const std::string& contigName) {
    // Format: @ReadName|ContigName
    out << "@" << read.qname << "|" << contigName << "\n"
        << read.seq << "\n"
        << "+\n";

    if (read.qual.empty() || read.qual == "*") {
        out << std::string(read.seq.length(), 'I') << "\n";
    } else {
        out << read.qual << "\n";
    }
}

// Main function: Filter reads based on final results and write to FASTQ
void writeEvidenceToFastq(
    const std::vector<result_t>& finalResults, // Input: Final fusion list for filtering
    const std::unordered_map<std::string, std::vector<sam_t>>& splitReadsMap,
    const std::unordered_map<std::string, std::vector<sam_t>>& spanReadsMap,
    const std::string& outputPrefix,
    const options_t& options) {

    // 1. Build a Whitelist of valid contigs from final results
    std::unordered_set<std::string> validContigs;
    for (const auto& res : finalResults) {
        // You can add logic here to only include "PASS" filters if needed
        // e.g., if (res.filter_status == "PASS")
        validContigs.insert(res.contig);
    }

    std::cout << "[INFO] Filtering evidence reads for " << validContigs.size() << " final fusion candidates.\n";

    // 2. Define output filenames
    std::string filenameR1 = options.output + "/" + outputPrefix + "_R1.fq";
    std::string filenameR2 = options.output + "/" + outputPrefix + "_R2.fq";

    std::ofstream outR1(filenameR1);
    std::ofstream outR2(filenameR2);

    if (!outR1.is_open() || !outR2.is_open()) {
        std::cerr << "Error: Unable to open FASTQ files for writing.\n";
        return;
    }

    // 3. Lambda to process, filter, and write reads
    auto processAndFilterReads = [&](const std::unordered_map<std::string, std::vector<sam_t>>& sourceMap) {
        for (const auto& entry : sourceMap) {
            const std::string& contigName = entry.first;

            // [FILTERING STEP]: Skip this contig if it's not in the final results
            if (validContigs.find(contigName) == validContigs.end()) {
                continue;
            }

            const std::vector<sam_t>& reads = entry.second;
            for (const auto& read : reads) {
                // 0x40 (64): First in pair (R1)
                // 0x80 (128): Second in pair (R2)
                if (read.flag & 64) {
                    writeFastqRecord(outR1, read, contigName);
                } else if (read.flag & 128) {
                    writeFastqRecord(outR2, read, contigName);
                } else {
                    // Default to R1 for single-end or undefined flags
                    writeFastqRecord(outR1, read, contigName);
                }
            }
        }
    };

    // 4. Apply processing to both Split and Span reads maps
    processAndFilterReads(splitReadsMap);
    processAndFilterReads(spanReadsMap);

    outR1.close();
    outR2.close();
    std::cout << "[INFO] Evidence reads extracted to:\n"
              << "  R1: " << filenameR1 << "\n"
              << "  R2: " << filenameR2 << "\n";
}
