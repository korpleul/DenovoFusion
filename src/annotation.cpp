//
// Created by xinwei on 7/7/24.
//

#include "annotation.h"


GeneAnnotator::GeneAnnotator(const std::string& gtf_path)
        : gtf_path_(gtf_path) {}

std::vector<annotation_t> GeneAnnotator::annotateAlignments(const std::vector<coordination_t>& coordinations) {
    std::vector<annotation_t> annotations(coordinations.size());
    std::ifstream gtf_file(gtf_path_);

    if (!gtf_file.is_open()) {
        std::cerr << "Failed to open GTF file: " << gtf_path_ << std::endl;
        return annotations;
    }

    std::unordered_map<std::string, std::unordered_map<std::string, int>> gene_transcript_cds_length;
    std::string line;

    while (getline(gtf_file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::vector<std::string> fields;
        std::string field;
        while (getline(iss, field, '\t')) {
            fields.push_back(field);
        }
        if (fields.size() < 9) continue;

        std::string feature = fields[2];
        int start = std::stoi(fields[3]);
        int end = std::stoi(fields[4]);

        std::unordered_map<std::string, std::string> attr_map;
        parseAttributes(fields[8], attr_map);

        if (feature == "exon" && attr_map.count("transcript_id") && attr_map.count("gene_id")) {
            if (attr_map.count("tag") == 0 || attr_map["tag"].find("basic") == std::string::npos)
                continue;

            std::string gene_id = attr_map["gene_id"];
            std::string transcript_id = attr_map["transcript_id"];
            int exon_len = end - start + 1;

            gene_transcript_cds_length[gene_id][transcript_id] += exon_len;
        }
    }

    std::unordered_map<std::string, std::string> gene_to_transcript;
    for (auto &kv : gene_transcript_cds_length) {
        const std::string &gene_id = kv.first;
        const auto &transcripts = kv.second;

        int max_len = -1;
        std::string selected_tx;
        for (auto &tx : transcripts) {
            if (tx.second > max_len) {
                max_len = tx.second;
                selected_tx = tx.first;
            }
        }
        gene_to_transcript[gene_id] = selected_tx;
    }

    gtf_file.clear();
    gtf_file.seekg(0, std::ios::beg);

    std::unordered_map<std::string, std::vector<std::pair<int, int>>> chrom_exons;
    std::unordered_map<std::string, std::vector<int>> chrom_exon_numbers;

    while (getline(gtf_file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::vector<std::string> fields;
        std::string field;
        while (getline(iss, field, '\t')) {
            fields.push_back(field);
        }
        if (fields.size() < 9) continue;

        std::string chrom = fields[0];
        std::string feature = fields[2];
        int start = std::stoi(fields[3]);
        int end = std::stoi(fields[4]);

        std::unordered_map<std::string, std::string> attr_map;
        parseAttributes(fields[8], attr_map);


        if (feature == "exon" && attr_map.count("transcript_id") && attr_map.count("gene_id")) {
            std::string gene_id = attr_map["gene_id"];
            std::string transcript_id = attr_map["transcript_id"];

            if (gene_to_transcript.count(gene_id) && gene_to_transcript[gene_id] == transcript_id) {
                chrom_exons[chrom].emplace_back(start, end);
                chrom_exon_numbers[chrom].push_back(
                    attr_map.count("exon_number") ? std::stoi(attr_map["exon_number"]) : -1
                );
            }
        }

        for (size_t i = 0; i < coordinations.size(); ++i) {
            const auto &coord = coordinations[i];
            if (coord.target == chrom && coord.tstart <= end && coord.tend >= start) {
                if (annotations[i].contigName.empty()) {
                    std::string geneName = attr_map.count("gene_name") ? attr_map["gene_name"] : " ";
                    std::string geneId = attr_map.count("gene_id") ? attr_map["gene_id"] : " ";
                    std::string strand = fields[6];
                    annotations[i] = annotation_t(coord.query, geneId, geneName,
                                                  coord.tstart, coord.tend,
                                                  chrom, coord.strand, "", "", -1);
                }
            }
        }
    }

    std::unordered_map<std::string, std::vector<size_t>> contigToIndices;
    for (size_t i = 0; i < annotations.size(); ++i) {
        contigToIndices[annotations[i].contigName].push_back(i);
    }

    for (auto& [contig, indices] : contigToIndices) {
        for (size_t i = 0; i + 1 < indices.size(); i += 2) {
            auto &a1 = annotations[indices[i]];
            auto &a2 = annotations[indices[i + 1]];

            if (a1.strand == "+" && a2.strand == "+") {
                a1.direction = "UPSTREAM";
                a2.direction = "DOWNSTREAM";
            } else if (a1.strand == "-" && a2.strand == "-") {
                a1.direction = "DOWNSTREAM";
                a2.direction = "UPSTREAM";
            } else {
                a1.direction = "UNDEFINED";
                a2.direction = "UNDEFINED";
            }
        }
    }

    for (size_t i = 0; i < annotations.size(); ++i) {
        const auto &coord = coordinations[i];
        auto &anno = annotations[i];

        auto it_exons = chrom_exons.find(coord.target);
        auto it_exon_nums = chrom_exon_numbers.find(coord.target);
        if (it_exons == chrom_exons.end() || it_exon_nums == chrom_exon_numbers.end()) {
            anno.regionType = "unknown";
            continue;
        }

        const auto &exons = it_exons->second;
        const auto &exon_numbers = it_exon_nums->second;

        int pos = (anno.direction == "UPSTREAM") ? coord.tend : coord.tstart;

        std::vector<std::tuple<int, int, int>> sorted_exons;
        for (size_t j = 0; j < exons.size(); ++j) {
            sorted_exons.emplace_back(exons[j].first, exons[j].second, exon_numbers[j]);
        }
        std::sort(sorted_exons.begin(), sorted_exons.end(),
                  [](auto &a, auto &b) { return std::get<0>(a) < std::get<0>(b); });

        auto annotate_pos = [&](int pos) -> std::pair<std::string, int> {
            for (auto &e : sorted_exons) {
                if (pos >= std::get<0>(e) && pos <= std::get<1>(e))
                    return {"exon", std::get<2>(e)};
            }
            for (size_t j = 0; j + 1 < sorted_exons.size(); ++j) {
                if (pos > std::get<1>(sorted_exons[j]) && pos < std::get<0>(sorted_exons[j + 1])) {
                    return {"intron", std::get<2>(sorted_exons[j]) };
                }
            }
            if (pos < std::get<0>(sorted_exons.front()))
                return {"upstream", std::get<2>(sorted_exons.front())};
            if (pos > std::get<1>(sorted_exons.back()))
                return {"downstream", std::get<2>(sorted_exons.back())};

            return {"unknown", -1};
        };

        auto [region, exon_number] = annotate_pos(pos);
        anno.regionType = region + "@" + std::to_string(exon_number);
    }

    for (size_t i = 0; i < annotations.size(); ++i) {
        if (annotations[i].contigName.empty()) {
            annotations[i] = annotation_t(coordinations[i].query, " ", " ",
                                          coordinations[i].tstart, coordinations[i].tend,
                                          coordinations[i].target, " ", "", "", -1);
        }
    }

    return annotations;
}





void GeneAnnotator::parseAttributes(const std::string& attributes, std::unordered_map<std::string, std::string>& attr_map) {
    std::istringstream attr_stream(attributes);
    std::string attr;

    while (std::getline(attr_stream, attr, ';')) {
        size_t key_start = attr.find_first_not_of(" \t\n\r\f\v");
        if (key_start == std::string::npos) continue;

        size_t key_end = attr.find(' ', key_start);
        if (key_end == std::string::npos) continue;

        std::string key = attr.substr(key_start, key_end - key_start);
        size_t value_start = attr.find('"', key_end);
        size_t value_end = attr.find('"', value_start + 1);

        if (value_start != std::string::npos && value_end != std::string::npos) {
            std::string value = attr.substr(value_start + 1, value_end - value_start - 1);

            if (key == "tag") {
                if (attr_map.count(key) > 0) {
                    attr_map[key] += "," + value;
                } else {
                    attr_map[key] = value;
                }
            } else {
                attr_map[key] = value;
            }
        }
    }
}


std::vector<coordination_t> keepOnlyTwoParts(const std::vector<coordination_t>& coords) {
    std::unordered_map<std::string, std::vector<coordination_t>> grouped;

    // Group by contig
    for (const auto& coord : coords) {
        grouped[coord.query].push_back(coord);
    }

    std::vector<coordination_t> filtered;

    for (auto& [contig, hits] : grouped) {
        // Only take the first and last element
        if (!hits.empty()) {
            filtered.push_back(hits.front());
            if (hits.size() > 1) {
                filtered.push_back(hits.back());
            }
        }
    }

    return filtered;
}

void keepOnlyTwoAnnotations(std::vector<annotation_t>& annotations) {
    std::unordered_map<std::string, std::vector<annotation_t>> grouped;

    // Group by contig
    for (const auto& ann : annotations) {
        grouped[ann.contigName].push_back(ann);
    }

    std::vector<annotation_t> filtered;

    for (auto& [contig, anns] : grouped) {
        // Only take the first and last annotation
        if (!anns.empty()) {
            filtered.push_back(anns.front());
            if (anns.size() > 1) {
                filtered.push_back(anns.back());
            }
        }
    }

    annotations.swap(filtered); // Replace original with filtered
}

