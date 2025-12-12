//
// Created by xinwei on 6/6/24.
//
#include "alignment.h"



std::vector<int> splitIndices(const std::string& indices) {
    std::vector<int> idx;
    std::istringstream iss(indices);
    std::string index;
    while (std::getline(iss, index, ' ')) {  // Assuming space-separated indices
        idx.push_back(std::stoi(index));
    }
    return idx;
}

// Alignment class
double alignment_t::setIdentity(int qstart, int qend, int tstart, int tend, int qnuminsert, int mismatch, int num_bases_aligned) {
    return calcIdentity(qstart, qend, tstart, tend, qnuminsert, mismatch, num_bases_aligned) ;
}

int alignment_t::setScore(int match, int mismatch, int qnuminsert, int tnuminsert, int query_len) {
    return calcScore(match, qnuminsert, tnuminsert, query_len) ;
}

void alignment_t::setNumBasesAligned() {
    if (!repmatch) {
        repmatch = 0;
    }
    num_bases_aligned = matches + mismatch + repmatch;
}

std::string alignment_t::details() const {
    std::ostringstream oss;
    if (blocks.size()) {
        oss << query << " " << qstart << "-" << qend << " " << target << " " << tstart << "-" << tend << " " << blocks.size();
    } else {
        oss << query << " " << qstart << "-" << qend << " " << target << " " << tstart << "-" << tend;
    }
    return oss.str();
}

std::string alignment_t::gff(const std::string &type) const {
    if (blocks.empty()) {
        return "";
    }
    std::ostringstream oss;
    std::string chrom = target.find("chr") == 0 ? target : "chr" + target;
    std::string contig_str = query.substr(0, query.find(" "));
    if (contig && contig) {
        contig_str = std::to_string(contig);
    }
    std::string strand = query_strand ? "+" : "-";
    std::string orient_str = "?";
    if (orient == 'f') {
        orient_str = "+";
    } else if (orient == 'b') {
        orient_str = "-";
    }
    std::string gff_group = target + ":" + std::to_string(tstart) + strand + orient_str + ":" + contig_str;
    if (contig) {
        std::string cov = std::to_string(contig);
        if (!cov.empty()) {
            gff_group += ":" + cov;
        }
    }
    if (orient_str != "?") {
        strand = orient_str;
    }
    for (size_t i = 0; i < blocks.size(); ++i) {
        auto [qstart, qend] = query_blocks[i];
        auto [tstart, tend] = blocks[i];
        std::string splice_site;
        if (blocks.size() > 1 && i < blocks.size() - 1) {
            splice_site = splice_sites[i];
        }
        oss << chrom << "\t" << qstart << ":" << qend << ":" << splice_site << "\t" << type << "\t" << tstart << "\t" << tend << "\t.\t" << strand << "\t.\t" << gff_group << "\n";
    }
    return oss.str();
}

std::string alignment_t::psl() const {
    if (!psl_str.empty()) {
        std::vector<std::string> cols;
        std::istringstream iss(psl_str);
        std::string col;
        while (std::getline(iss, col, '\t')) {
            cols.push_back(col);
        }

        if (!std::regex_search(cols[13], std::regex("^(chr|scaffold)", std::regex::icase))) {
            cols[13] = "chr" + cols[13];
        }

        if (orient == 'f') {
            cols[8] = "+";
        } else if (orient == 'b') {
            cols[8] = "-";
        }
        std::ostringstream oss;
        std::copy(cols.begin(), cols.end(), std::ostream_iterator<std::string>(oss, "\t"));
        std::string psl_str = oss.str();
    } else {
        std::vector<std::string> fields;
        fields.push_back(std::to_string(matches));
        fields.push_back(std::to_string(mismatch));
        fields.push_back("0");
        fields.push_back("0");
        fields.push_back(std::to_string(qnuminsert));
        fields.push_back(std::to_string(qbaseinsert));
        fields.push_back(std::to_string(tnuminsert));
        fields.push_back(std::to_string(tbaseinsert));
        std::string target_str = target.substr(0, 3) == "chr" ? target : "chr" + target;
        fields.push_back(std::string(1, query_strand));
        fields.push_back(query);
        fields.push_back(std::to_string(query_len));
        fields.push_back(std::to_string(qstart - 1));
        fields.push_back(std::to_string(qend));
        fields.push_back(target_str);
        fields.push_back(std::to_string(target_len));
        std::vector<int> tcoords = {tstart, tend};
        std::sort(tcoords.begin(), tcoords.end());
        fields.push_back(std::to_string(tcoords[0] - 1));
        fields.push_back(std::to_string(tcoords[1]));
        fields.push_back(std::to_string(blocks.size()));
        std::ostringstream block_sizes_oss;
        for (const auto &block : blocks) {
            block_sizes_oss << block.second - block.first + 1 << ",";
        }
        fields.push_back(block_sizes_oss.str());
        std::ostringstream qstarts_oss;
        for (const auto &block : query_blocks) {
            qstarts_oss << block.first - 1 << ",";
        }
        fields.push_back(qstarts_oss.str());
        std::ostringstream tstarts_oss;
        for (const auto &block : blocks) {
            tstarts_oss << block.first - 1 << ",";
        }
        fields.push_back(tstarts_oss.str());
        std::ostringstream oss;
        std::copy(fields.begin(), fields.end(), std::ostream_iterator<std::string>(oss, "\t"));
        std::string psl_str = oss.str();
    }
    return psl_str;
}


void alignment_t::correctBlocks(const std::unordered_map<std::string, bool> &splice_motifs, const std::string &target_seq, const std::string &query_seq) {
    if (splice_sites.empty() || splice_motifs.empty()) {
        return;
    }
    correctUnaligned(query_seq, target_seq);
    correctSingleGaps(splice_motifs, target_seq);
    correctNeighborGaps(splice_motifs, target_seq);
}


std::string ReverseComplement(const std::string &sequence) {
    std::string rev_comp(sequence.rbegin(), sequence.rend());
    std::transform(rev_comp.begin(), rev_comp.end(), rev_comp.begin(), [](char base) {
        switch (base) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'G': return 'C';
            case 'C': return 'G';
            case 'a': return 't';
            case 't': return 'a';
            case 'g': return 'c';
            case 'c': return 'g';
            default: return base;
        }
    });
    return rev_comp;
}

void alignment_t::correctSingleGaps(const std::unordered_map<std::string, bool> &splice_motifs, const std::string &target_seq) {
    std::unordered_map<int, int> gaps;
    for (size_t i = 0; i < blocks.size() - 1; ++i) {
        const std::string &ss = splice_sites[i];
        if (!ss.empty() && splice_motifs.find(ss) == splice_motifs.end() && splice_motifs.find(ReverseComplement(ss)) == splice_motifs.end()) {
            if (std::abs(query_blocks[i + 1].first - query_blocks[i].second) == 1) {
                gaps[i] = 0;
            }
        }
    }
    if (!gaps.empty()) {
        std::vector<std::pair<int, int>> target_blocks = blocks;
        std::vector<std::pair<int, int>> query_blocks_copy = query_blocks;
        std::vector<std::string> splice_sites_copy = splice_sites;
        std::vector<int> gap_indices;
        for (const auto &gap : gaps) {
            gap_indices.push_back(gap.first);
        }
        std::sort(gap_indices.begin(), gap_indices.end());
        for (int i : gap_indices) {
            std::pair<int, int> &tblock1 = target_blocks[i];
            std::pair<int, int> &tblock2 = target_blocks[i + 1];
            std::pair<int, int> &qblock1 = query_blocks_copy[i];
            std::pair<int, int> &qblock2 = query_blocks_copy[i + 1];
            std::string splice_site = fixSingleGap(tblock1, tblock2, qblock1, qblock2, splice_motifs, target_seq, query_strand);
            if (!splice_site.empty()) {

                std::cerr << query << " changed blocks " << target << " ["
                          << target_blocks[i].first << ", " << target_blocks[i].second << "] to ["
                          << tblock1.first << ", " << tblock1.second << "] to "
                          << splice_site << "\n";


                std::cerr << query << " changed blocks " << target << " ["
                          << target_blocks[i].first << ", " << target_blocks[i + 1].second << "] to ["
                          << tblock1.first << ", " << tblock2.second << "] to "
                          << splice_site << "\n";

                target_blocks[i] = tblock1;
                target_blocks[i + 1] = tblock2;
                query_blocks_copy[i] = qblock1;
                query_blocks_copy[i + 1] = qblock2;
                splice_sites_copy[i] = splice_site;
            }
        }
        if (target_blocks != blocks) {
            blocks = target_blocks;
            query_blocks = query_blocks_copy;
            splice_sites = splice_sites_copy;
            if (!mismatch || mismatch == 0) {
                mismatch = 1;
            }
        }
    }
}

void alignment_t::correctNeighborGaps(const std::unordered_map<std::string, bool> &splice_motifs, const std::string &target_seq) {
    std::unordered_map<int, int> gaps;
    for (size_t i = 0; i < blocks.size() - 1; ++i) {
        const std::string &ss = splice_sites[i];
        if (!ss.empty() && splice_motifs.find(ss) == splice_motifs.end() && splice_motifs.find(ReverseComplement(ss)) == splice_motifs.end()) {
            if (std::abs(query_blocks[i + 1].first - query_blocks[i].second) == 1) {
                gaps[i] = 0;
            }
        }
    }
    if (!gaps.empty()) {
        std::vector<std::pair<int, int>> target_blocks = blocks;
        std::vector<std::pair<int, int>> query_blocks_copy = query_blocks;
        std::vector<std::string> splice_sites_copy = splice_sites;
        std::vector<int> gap_indices;
        for (const auto &gap : gaps) {
            gap_indices.push_back(gap.first);
        }
        std::sort(gap_indices.begin(), gap_indices.end());
        std::unordered_map<std::string, std::tuple<std::pair<int, int>, std::pair<int, int>, std::pair<int, int>, std::pair<int, int>, std::string>> replaced;
        std::vector<std::string> replaced_ordered;
        for (size_t i = 0; i < gap_indices.size() - 1; ++i) {
            int i1 = gap_indices[i];
            int i2 = i1 + 1;
            if (gaps.find(i2) != gaps.end() && reinterpret_cast<const char *>(gaps[i2]) != "fixed" && i2 + 1 < target_blocks.size()) {
                std::pair<int, int> &tblock1 = target_blocks[i1];
                std::pair<int, int> &tblock2 = target_blocks[i2];
                std::pair<int, int> &tblock3 = target_blocks[i2 + 1];
                std::pair<int, int> &qblock1 = query_blocks_copy[i1];
                std::pair<int, int> &qblock2 = query_blocks_copy[i2];
                std::pair<int, int> &qblock3 = query_blocks_copy[i2 + 1];
                std::string splice_site = fixNeighborGaps(tblock1, tblock2, tblock3, qblock1, qblock2, qblock3, splice_motifs, target_seq, query_strand);
                if (!splice_site.empty()) {
                    std::string idx = std::to_string(i1) + " " + std::to_string(i2) + " " + std::to_string(i2 + 1);
                    replaced[idx] = std::make_tuple(tblock1, tblock3, qblock1, qblock3, splice_site);
                    gaps[i1] += 1;
                    gaps[i2] += 1;
                    replaced_ordered.push_back(idx);
                }
            }
        }
        std::reverse(replaced_ordered.begin(), replaced_ordered.end());
        for (const auto &indices : replaced_ordered) {
            auto new_blocks = replaced[indices];
            bool ok = true;
            for (const auto &index : indices) {
                if (gaps.find(index) != gaps.end() && gaps[index] > 1) {
                    ok = false;
                    break;
                }
            }

            auto idx = splitIndices(indices);

            if (ok) {
                std::cerr << query << " changed blocks " << target << " " << blocks[idx[0]].first << ", " << blocks[idx[0]].second << " to " << std::get<0>(new_blocks).first << ", " << std::get<0>(new_blocks).second << "\n";
                std::cerr << query << " changed blocks " << target << " " << blocks[idx[2]].first << ", " << blocks[idx[2]].second << " to " << std::get<1>(new_blocks).first << ", " << std::get<1>(new_blocks).second << "\n";
                std::cerr << query << " removed block " << target << " " << blocks[idx[1]].first << ", " << blocks[idx[1]].second << "\n";


                target_blocks[idx[0]] = std::get<0>(new_blocks);
                target_blocks[idx[2]] = std::get<1>(new_blocks);
                query_blocks_copy[idx[0]] = std::get<2>(new_blocks);
                query_blocks_copy[idx[2]] = std::get<3>(new_blocks);
                splice_sites_copy[idx[0]] = std::get<4>(new_blocks);

                target_blocks.erase(target_blocks.begin() + idx[1]);
                query_blocks_copy.erase(query_blocks_copy.begin() + idx[1]);
                splice_sites_copy.erase(splice_sites_copy.begin() + idx[1]);
            }

        }
        if (target_blocks != blocks) {
            blocks = target_blocks;
            query_blocks = query_blocks_copy;
            splice_sites = splice_sites_copy;
            if (!mismatch || mismatch == 0) {
                mismatch = 1;
            }
        }
    }
}

void alignment_t::correctUnaligned(const std::string &query_seq, const std::string &target_seq, int max_diff) {
    std::unordered_map<int, std::pair<std::string, int>> expand;
    for (size_t i = 0; i < blocks.size() - 1; ++i) {
        if (std::abs(query_blocks[i + 1].first - query_blocks[i].second) != 1 && std::abs(blocks[i + 1].first - blocks[i].second) != 1) {
            std::string qseq, tseq;
            if (query_strand == '+') {
                qseq = query_seq.substr(query_blocks[i].second, query_blocks[i + 1].first - query_blocks[i].second - 1);
            } else {
                qseq = ReverseComplement(query_seq.substr(query_blocks[i + 1].first, query_blocks[i].second - query_blocks[i + 1].first - 1));
            }
            tseq = target_seq.substr(blocks[i].second, blocks[i + 1].first - blocks[i].second - 1);
            std::string gap_seq = tseq.substr(qseq.size());
            std::string ss_start = gap_seq.substr(0, 2) + gap_seq.substr(gap_seq.size() - 2);
            gap_seq = tseq.substr(0, tseq.size() - qseq.size() - 1);
            std::string ss_end = gap_seq.substr(0, 2) + gap_seq.substr(gap_seq.size() - 2);
            if (qseq.size() > tseq.size() || tseq.size() - qseq.size() < 20) {
                continue;
            }
            int start_diff = 0;
            for (size_t j = 0; j < qseq.size(); ++j) {
                if (qseq[j] != tseq[j]) {
                    ++start_diff;
                }
            }
            int end_diff = 0;
            for (size_t j = 0; j < qseq.size(); ++j) {
                if (qseq[j] != tseq[tseq.size() - qseq.size() + j]) {
                    ++end_diff;
                }
            }
            if (start_diff < end_diff && start_diff < max_diff) {
                expand[i] = {"left", static_cast<int>(qseq.size())};
                mismatch = start_diff;
            } else if (end_diff < start_diff && end_diff < max_diff) {
                expand[i] = {"right", static_cast<int>(qseq.size())};
                mismatch = end_diff;
            }
        }
    }
    for (const auto &[idx, where] : expand) {
        if (where.first == "left") {
            std::cerr << query << " changed blocks (" << blocks[idx].first << ", " << blocks[idx].second << ")\n";

            if (query_strand == '+') {
                query_blocks[idx].second += where.second;
            } else {
                query_blocks[idx].second -= where.second;
            }
            blocks[idx].second += where.second;
            splice_sites[idx] = qstart;
        } else {
            std::cerr << query << " changed blocks (" << blocks[idx + 1].first << ", " << blocks[idx + 1].second << ")\n";

            if (query_strand == '+') {
                query_blocks[idx + 1].first -= where.second;
            } else {
                query_blocks[idx + 1].first += where.second;
            }
            blocks[idx + 1].first -= where.second;
            splice_sites[idx] = qend;
        }
    }
}




std::string alignment_t::fixSingleGap(std::pair<int, int> &tblock1, std::pair<int, int> &tblock2, std::pair<int, int> &qblock1, std::pair<int, int> &qblock2, const std::unordered_map<std::string, bool> &splice_motifs, const std::string &target_seq, char query_strand) {
    const int min_size = 10;
    const int max_size = 10000;
    const std::vector<int> shuffle_sizes = {-2, -1, 1, 2};
    const std::pair<int, int> tgap = {tblock1.second + 1, tblock2.first - 1};
    const int tsize = tgap.second - tgap.first + 1;
    if (tsize < min_size || tsize > max_size) {
        return "";
    }
    bool artefact = false;
    bool ambiguous = false;
    std::string splice_site;
    for (int size : shuffle_sizes) {
        if (std::abs(size) >= (tblock1.second - tblock1.first) || std::abs(size) >= (tblock2.second - tblock2.first)) {
            continue;
        }
        const std::pair<int, int> coord = {tgap.first + size, tgap.second + size};
        const std::string gap_seq = target_seq.substr(coord.first - 1, coord.second - coord.first + 1);
        const std::string ss = gap_seq.substr(0, 2) + gap_seq.substr(gap_seq.size() - 2);
        if (splice_motifs.count(ss) || splice_motifs.count(ReverseComplement(ss))) {
            if (!artefact) {
                artefact = true;
            } else {
                ambiguous = true;
            }
            if (!ambiguous) {
                tblock1.second += size;
                tblock2.first += size;
                if (query_strand == '+') {
                    qblock1.second += size;
                    qblock2.first += size;
                } else {
                    qblock1.second -= size;
                    qblock2.first -= size;
                }
                splice_site = ss;
            }
        }
    }
    return artefact && !ambiguous ? splice_site : "";
}

std::string alignment_t::fixNeighborGaps(std::pair<int, int> &tblock1, std::pair<int, int> &tblock2, std::pair<int, int> &tblock3, std::pair<int, int> &qblock1, std::pair<int, int> &qblock2, std::pair<int, int> &qblock3, const std::unordered_map<std::string, bool> &splice_motifs, const std::string &target_seq, char query_strand) {
    const std::pair<int, int> tgap1 = {tblock1.second + 1, tblock2.first - 1};
    const std::pair<int, int> tgap2 = {tblock2.second + 1, tblock3.first - 1};
    const int tgap1_size = tgap1.second - tgap1.first + 1;
    const int tgap2_size = tgap2.second - tgap2.first + 1;
    const int tblock2_size = tblock2.second - tblock2.first + 1;
    const int max_shuffle_size = 10;
    if (tblock2_size > max_shuffle_size) {
        return "";
    }
    const int small_cutoff = 10;
    const int big_cutoff = 20;
    std::string big;
    if (tgap1_size < small_cutoff && tgap2_size > big_cutoff) {
        big = "right";
    } else if (tgap2_size < small_cutoff && tgap1_size > big_cutoff) {
        big = "left";
    }
    if (big.empty()) {
        return "";
    }
    if (big == "right") {
        tblock3.first -= tblock2_size;
        qblock3.first = qblock2.first;
    } else {
        tblock1.second += tblock2_size;
        qblock1.second = qblock2.second;
    }
    const std::pair<int, int> tgap = {tblock1.second + 1, tblock3.first - 1};
    const std::string gap_seq = target_seq.substr(tgap.first - 1, tgap.second - tgap.first + 1);
    const std::string ss = gap_seq.substr(0, 2) + gap_seq.substr(gap_seq.size() - 2);
    std::string splice_site;
    if (splice_motifs.count(ss) || splice_motifs.count(ReverseComplement(ss))) {
        splice_site = ss;
    } else {
        splice_site = fixSingleGap(tblock1, tblock3, qblock1, qblock3, splice_motifs, target_seq, query_strand);
    }
    return splice_site;
}

double alignment_t::calcIdentity(int qstart, int qend, int tstart, int tend, int qnuminsert, int mismatch, int num_bases_aligned) {
    // Calculate the length of the query and target sequences
    int qspan = qend - qstart;
    int tspan = tend - tstart;

    // Calculate the maximum span
    int maxspan = std::min(qspan, tspan);

    // If the maximum span is less than or equal to 0, return 0.
    if (maxspan <= 0) {
        return 0.0;
    }

    // The query span is assumed to be greater than the target span. Otherwise, they are assumed to be equal.
    int sizediff = qspan - tspan;
    if (sizediff < 0) {
        sizediff = 0;
    }

    // Calculate a size difference factor
    int sizediff_factor = static_cast<int>(round(3 * log(1 + sizediff)));

    // Calculate the original
    int raw_numerator = mismatch + qnuminsert + sizediff_factor;
    double millibad = (1000.0 * raw_numerator) / num_bases_aligned;

    return (100.0 - millibad * 0.1) / 100.0;
}

int alignment_t::calcScore(int match, int qnuminsert, int tnuminsert, int query_length) {

    int alignment_length = match + qnuminsert + tnuminsert;


    if (alignment_length == 0 || query_length == 0) {
        return 0;
    }


    double coverage_percentage = static_cast<double>(alignment_length) / query_length * 100;

    int final_score = static_cast<int>(coverage_percentage);

    return final_score > 0 ? final_score : 0;
}



void calculate_ends(psl_t& entry) {
    entry.qEnds = entry.qStarts;
    entry.tEnds = entry.tStarts;
    for (size_t i = 0; i < entry.qStarts.size(); ++i) {
        entry.qEnds[i] += entry.blockSizes[i] - 1;
        entry.tEnds[i] += entry.blockSizes[i] - 1;
    }
}



std::unordered_map<std::string, std::vector<alignment_t>> index_by_qname(const std::vector<alignment_t>& entries) {
    std::unordered_map<std::string, std::vector<alignment_t>> data_by_qname;
    for (const auto& item : entries) {
        data_by_qname[item.query].push_back(item);
    }
    return data_by_qname;
}

void process_alignment(const paf_t& alignment) {
    if (alignment.optional_fields.find("NM") != alignment.optional_fields.end()) {
        int mismatch_count = std::stoi(alignment.optional_fields.at("NM"));
        std::cout << "Mismatch count: " << mismatch_count << std::endl;
    }

    if (alignment.optional_fields.find("cg") != alignment.optional_fields.end()) {
        std::string cigar_string = alignment.optional_fields.at("cg");

        std::cout << "CIGAR: " << cigar_string << std::endl;
    }
}


void parseCigarString(const std::string& cigar, int& qnuminsert, int& tnuminsert) {
    std::string number = "";
    for (char c : cigar) {
        if (isdigit(c)) {
            number += c;
        } else {
            int count = std::stoi(number);
            number = "";

            if (c == 'I') {
                qnuminsert += count;
            } else if (c == 'D') {
                tnuminsert += count;
            }
        }
    }
}
