//
// Created by xinwei on 12/4/24.
//

#include "run_minimap2sam.h"
#include <iostream>
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <limits>
#include <sstream>
#include <cstdlib>
#include <string>
#include <sys/resource.h>
#include <unordered_map>
#include <vector>
#include <thread>
#include <future>
#include <queue>

#include "options.h"
#include "sam.h"
#include "coord_pair.h"
#include "psl_functions.h"
#include "alignments_chosen.h"
#include "alignment.h"
#include "candidate_group.h"
#include "fasta.h"
#include "bowtie2.h"
#include "overlap.h"
#include "realign_support.h"
#include "annotation.h"
#include "output_fusions.h"
#include "coverage.h"
#include "aligned_coord_pair.h"
#include "breakpoint.h"
#include "candidate_group.h"
#include "common.h"
#include "error.h"
#include "log.h"
#include "filter_duplicates.h"
#include "filter_edge_unaligned.h"
#include "filter_homologs.h"
#include "filter_long_gap.h"
#include "filter_internal_tandem_duplication.h"
#include "filter_min_support.h"
#include "filter_mt.h"
#include "filter_small_fragments.h"
#include "recover_known_fusion.h"
#include "support_writing.h"
#include "utils.h"


// We apply this function to integrate all the calculation steps into this function, which includes 4 stages, same with the
// description in article,
// Stage1: chose the best alignments which could represent each contig
// Stage2: group the contigs into 4 subsets, single alignments, gap alignments, paired alignments, multiple alignments
// Stage3: realignment from raw reads to chosen contigs, calculate the support split and span reads
// Stage4: filtering, output the final fusions
// Stage5: calculation of the coverage and  prediction of breakpoints

void run_minimap2sam(const options_t& options) {
    // Stage0: general settings for the programm, initialize log file, clear the existing log file, set a time log
    time_t start_time;
    time(&start_time);

    Logger::init(options);
    std::cout << get_time_string() << " Program DenovoFusion start" << std::endl;
    Logger::Info(get_time_string() + " Program DenovoFusion start");

    std::cout << get_time_string() << " Launching fusion gene analysing program version " << DENOVOFUSION_VERSION << "\n" << std::flush;
    Logger::Info(get_time_string() + " DenovoFusion version = " + DENOVOFUSION_VERSION);

    // Stage1: find the representative alignment for each contig, then find the contigs which include the putative fusion gene
    // Load sam alignment file
    std::cout << get_time_string() << " Stage1: Determine representive alignment for each contig " << std::endl;
    Logger::Info(get_time_string() + " Stage1: Determine representive alignment for each contig ");
    std::cout << get_time_string() << " Loading alignments from SAM file:" << " '" << options.input_file << "' " << "\n" << std::flush;
    Logger::Info(get_time_string() + " Loading alignments from SAM file:" + " '" +  options.input_file + "' ");

    // Load the SAM file
    std::vector<sam_t> sams;  // this is the first sam loading, which is distinguished with realignment
    loadSamFile(options.input_file, sams);

    std::cout << get_time_string() << " The SAM file includes in total " << sams.size() << " alignments\n" << std::flush;

    // Load contig file, remember to make a suitable hash container with key and value
    auto fasta_sequences = load_fasta_sequences(options.input_assembly);
    std::cout << get_time_string() << " The original contig include " << fasta_sequences.size() << " sequences "<< std::endl;
    Logger::Info( get_time_string() + " The original contig include " + std::to_string(fasta_sequences.size()) + " sequences ");

    // Process the SAM file alignments
    auto qname_counts = count_qnames(sams);
    std::vector<alignment_t> alignments;

    // Iterate over original sam entries
    for (const auto& sam : sams) {
        alignment_t alignment("minimap2sam", sam, fasta_sequences);
        alignments.push_back(alignment);
    }

    // according to qName, divide the alignment into groups
    auto groups = index_by_qname(alignments);

    // Convert groups to a vector to facilitate task distribution, and filter out groups with more than 100 members
    std::vector<std::unordered_map<std::string, std::vector<alignment_t>>> group_vector;
    // Initialize minimum and maximum row numbers
    int min_line_number = std::numeric_limits<int>::max();
    int max_line_number = std::numeric_limits<int>::min();
    int total_contigs = 0;  // Record the total number of contigs

    for (const auto& group : groups) {
        int group_size = group.second.size();
        // set record in log file
        // Logger::Info(get_time_string() + " contig name: " + group.first + "，alignments number: " + std::to_string(group.second.size()));

        // Record the minimum and maximum row numbers
        if (group_size < min_line_number) {
            min_line_number = group_size;
        }
        if (group_size > max_line_number) {
            max_line_number = group_size;
        }

        // Increase contig count
        ++total_contigs;

        if (qname_counts[group.first] <= options.max_alignment_count) {
            std::unordered_map<std::string, std::vector<alignment_t>> single_group;
            single_group[group.first] = group.second;
            group_vector.push_back(single_group);
        }
    }

    // Record the total number of contigs and line number range at the end
    Logger::Info(get_time_string() + " This sam file includes in total " + std::to_string(total_contigs) + " contigs，alignments number in each contig range from: " +
    std::to_string(min_line_number) + " to " + std::to_string(max_line_number));

    // Assign tasks for each thread
    int num_threads = std::min(options.threads, static_cast<int>(group_vector.size()));
    int chunk_size = group_vector.size() / num_threads;
    int remainder = group_vector.size() % num_threads;

    if (group_vector.size() < num_threads) {
        num_threads = group_vector.size();
    }


    std::vector<std::vector<alignment_t>> results(num_threads);

    // Define a worker function for each thread
    auto worker = [&](int start, int end, int thread_index) {
        for (int i = start; i < end; ++i) {
            auto result = calculate_alignments_score(group_vector[i], options);
            results[thread_index].insert(results[thread_index].end(), result.begin(), result.end());
        }
    };

    // Creating a Thread Pool
    std::vector<std::thread> threads;
    int start = 0;
    for (int i = 0; i < num_threads; ++i) {
        int end = start + chunk_size + (i < remainder ? 1 : 0);
        threads.emplace_back(worker, start, end, i);
        start = end;
    }

    // Wait for all threads to complete
    for (auto& thread : threads) {
        thread.join();
    }

    // Merge the results of all threads
    std::vector<alignment_t> identity_filtered_alignments;
    for (const auto& result : results) {
        for (const auto& alignment : result) {
            if (alignment.identity > options.min_identity_fract) {  // Apply filters
                identity_filtered_alignments.push_back(alignment);
            }
        }
    }

    // Apply score filter to chosen alignments
    std::unordered_map<std::string, std::vector<alignment_t>> contig_groups = index_by_qname(identity_filtered_alignments);

    std::vector<alignment_t> chosen_alignments;


    // Filter out groups with less than 2 alignments
    for (const auto& group : contig_groups) {
        int totalScore = 0;
        bool validGroup = true; // Used to determine whether the group meets the conditions

        for (const auto& alignment : group.second) {
            // Apply score filter
            if (alignment.score <= options.min_score_each) {
                validGroup = false;
                break; // If any alignment has a score less than or equal to 10, the group is invalid
            }
            totalScore += alignment.score;
        }

        // Check if the group meets the conditions
        if (validGroup && totalScore > options.min_score_total) {
            chosen_alignments.insert(chosen_alignments.end(), group.second.begin(), group.second.end());
        }
    }


    //
    Logger::Info(get_time_string() + " Count the chosen alignments which could represent each contig: " + std::to_string(chosen_alignments.size()));
    std::cout << get_time_string() << " Count the chosen alignments which could represent each contig: " << chosen_alignments.size() << std::endl;

    // Stage2: Grouping the alignments
    std::cout << get_time_string() << " Stage2: Grouping the alignments " << std::endl;
    Logger::Info(get_time_string() + " Stage2: Grouping the alignments ");
    // Using functions to get different types of alignment
    auto singles = single_alignments(chosen_alignments);
    auto pairs = pair_alignments(chosen_alignments);
    auto multiples = multiple_alignments(chosen_alignments);

    // Print the number of alignments in each category
    std::cout << get_time_string() << " All the chosen alignments will be divided into different types: single alignments, gaps alignments, overlaps pairs alignments, multiples alignments " << std::endl;
    Logger::Info(get_time_string() + " All the chosen alignments will be divided into different types: single alignments, gaps alignments, overlaps pairs alignments, multiples alignments");
    std::cout << get_time_string() << " single alignments: " << singles.size() << std::endl;
    Logger::Info(get_time_string() + " single alignments: "  + std::to_string(singles.size()));
    std::cout << get_time_string() << " pairs alignments: " << pairs.size() << std::endl;
    Logger::Info(get_time_string() + " pairs alignments: " + std::to_string(pairs.size()));
    std::cout << get_time_string() << " multiples alignments: " << multiples.size() << std::endl;
    Logger::Info(get_time_string() + " multiples alignments: " + std::to_string(multiples.size()));

    group_alignments(chosen_alignments);

    // Classify pairwise alignments using the classify_alignments function
    auto classified_alignments = classify_alignments(pairs,options);



    // Modify the filter_and_combine_same_strand call to include max_value and is_gap parameters
    auto overlaps_same_strand = filter_and_combine_same_strand(classified_alignments, OVERLAPS_SAME_STRAND, options.max_overlap_size, false);
    auto overlaps_diff_strand = filter_and_combine_diff_strand(classified_alignments, OVERLAPS_DIFFERENT_STRAND, options.max_overlap_size, false);
    auto gaps_same_strand = filter_and_combine_same_strand(classified_alignments, GAP_SAME_STRAND, options.max_gap_size, true);
    auto gaps_diff_strand = filter_and_combine_diff_strand(classified_alignments, GAP_DIFFERENT_STRAND, options.max_gap_size, true);


    // Print the number of alignments in each category
    std::cout << get_time_string() << " Number of overlaps same strand alignments: " << overlaps_same_strand.size() << std::endl;
    std::cout << get_time_string() << " Number of overlaps different strand alignments: " << overlaps_diff_strand.size() << std::endl;
    std::cout << get_time_string() << " Number of gaps same strand alignments: " << gaps_same_strand.size() << std::endl;
    std::cout << get_time_string() << " Number of gaps different strand alignments: " << gaps_diff_strand.size() << std::endl;

    // Log detailed alignments for each category
    Logger::logFile << get_time_string() << " Information for the alignments classified Overlaps Same Strand: " << std::endl;
    for (const auto& alignment : overlaps_same_strand) {
        Logger::logFile << get_time_string() << "\tContig:" << alignment.query
                                             << "\tqlength:" << alignment.query_len
                                             << "\tqstart:" << alignment.qstart
                                             << "\tqend:" << alignment.qend
                                             << "\tchr:" << alignment.target
                                             << "\tstrand:" << alignment.query_strand
                                             << "\tidentity:" << alignment.identity
                                             << "\tscore:" << alignment.score
                                             << std::endl;
    }

    Logger::logFile << get_time_string() << " Information for the alignments classified Overlaps Different Strand: " << std::endl;
    for (const auto& alignment : overlaps_diff_strand) {
        Logger::logFile << get_time_string() << "\tContig:" << alignment.query
                                             << "\tqlength:" << alignment.query_len
                                             << "\tqstart:" << alignment.qstart
                                             << "\tqend:" << alignment.qend
                                             << "\tchr:" << alignment.target
                                             << "\tstrand:" << alignment.query_strand
                                             << "\tidentity:" << alignment.identity
                                             << "\tscore:" << alignment.score
                                             << std::endl;
    }

    Logger::logFile << get_time_string() << " Information for the alignments classified Gaps Same Strand: " << std::endl;
    for (const auto& alignment : gaps_same_strand) {
        Logger::logFile << get_time_string() << "\tContig:" << alignment.query
                                              << "\tqlength:" << alignment.query_len
                                              << "\tqstart:" << alignment.qstart
                                              << "\tqend:" << alignment.qend
                                              << "\tchr:" << alignment.target
                                              << "\tstrand:" << alignment.query_strand
                                              << "\tidentity:" << alignment.identity
                                              << "\tscore:" << alignment.score
                                              << std::endl;
    }

    Logger::logFile << get_time_string() << " Information for the alignments classified Gaps Different Strand: " << std::endl;
    for (const auto& alignment : gaps_diff_strand) {
        Logger::logFile << get_time_string() << "\tContig:" << alignment.query
                                             << "\tqlength:" << alignment.query_len
                                             << "\tqstart:" << alignment.qstart
                                             << "\tqend:" << alignment.qend
                                             << "\tchr:" << alignment.target
                                             << "\tstrand:" << alignment.query_strand
                                             << "\tidentity:" << alignment.identity
                                             << "\tscore:" << alignment.score
                                             << std::endl;
    }

    // Extract base query names from overlaps_same_strand
    std::unordered_set<std::string> query_names_overlap = extractBaseQueryNames(overlaps_same_strand);
    std::unordered_set<std::string> query_names_gap = extractBaseQueryNames(gaps_same_strand);

    // Merge the two sets
    std::unordered_set<std::string> query_names = query_names_overlap;
    query_names.insert(query_names_gap.begin(), query_names_gap.end());

    // Filter alignments based on base query name
    std::vector<alignment_t> filtered_alignments = filterAlignmentsByQuery(chosen_alignments, query_names);

    // Collect and merge sequences
    auto mergedSequences = collectSequences(filtered_alignments, fasta_sequences);

    // Output the merged sequence to a file
    outputMergedSequences(mergedSequences, options);
    std::cout << get_time_string() << " Merged sequences " << mergedSequences.size() << " have been written to " << options.output << "/" << options.prefix << ".chosen.fasta" << std::endl;
    Logger::logFile << get_time_string() << " Merged sequences have been written to " << options.output << "/" << options.prefix << ".chosen.fasta" << std::endl;


    // Stage3: realign the contigs by using bowtie2
    // Check if bowtie2 running
    // Build index
    // Get user's home directory from environment variable
    std::cout << get_time_string() << " Stage3: Start realigning reads to contigs with using bowtie2 " << std::endl;
    Logger::Info(get_time_string() + " Stage3: Start realigning reads to contigs with using bowtie2 ");

    const char* homeDir = getenv("HOME");
    if (!homeDir) {
        std::cout << get_time_string() << " Error: HOME directory not found." << std::endl;
        Logger::Error(get_time_string() + " Error: HOME directory not found.");
        return;  // Return error code
    }

    //Set up the index directory
    std::string directorySetupCommand = "mkdir -p " + std::string(options.output) + "/" + options.prefix +"_idx";
    system(directorySetupCommand.c_str());  // Make sure the directory exists

    // Temporarily set the PATH environment variable
    std::string bowtie2_bin = "./bowtie2-2.5.4-linux-x86_64/bowtie2";
    std::cout << get_time_string() << " Updated PATH for bowtie2 binaries." << std::endl;
    Logger::Info(get_time_string() + " Updated PATH for bowtie2 binaries.");

    // Build the index
    build_bowtie2_index(options);
    std::cout << get_time_string() << " Bowtie2 index built." << std::endl;
    Logger::Info(get_time_string() + " Bowtie2 index built.");

    // Run realignment step
    run_bowtie2(options);
    std::cout << get_time_string() << " Finished bowtie2 alignment. " << std::endl;
    Logger::Info(get_time_string() + " Finished bowtie2 alignment. ");


    // Load sam file
    std::string samFilePath = options.output + options.prefix + ".sam";

    // Create an instance of SamFileCls to handle the SAM file
    SamFileCls samFile(samFilePath);
    std::vector<sam_t> samEntries;  // Vector to store all SAM entries

    sam_t samEntry;  // Temporary object to store each entry

    // Read entries from the SAM file
    while (samFile.next(samEntry)) {
        samEntries.push_back(samEntry);  // Store each entry in the vector
    }

    // Reset file reading position, if needed to re-read the file
    samFile.reset();
    // Close the file
    samFile.close();

    // Paired the exact include fusion information alignments
    auto pairedAlignments_overlap = pairAlignments(overlaps_same_strand);
    auto pairedAlignments_gap = pairAlignments(gaps_same_strand);

    // Assuming pairedAlignments_overlap and pairedAlignments_gap are std::vector<std::pair<alignment_t, alignment_t>>
    std::vector<std::pair<alignment_t, alignment_t>> pairedAlignments = pairedAlignments_overlap;
    // Use insert to add all pairs from pairedAlignments_gap
    pairedAlignments.insert(
        pairedAlignments.end(),
        pairedAlignments_gap.begin(),
        pairedAlignments_gap.end()
    );

    // Get the complete best aligned alignments from the same contig
    auto overlapResults = processAlignments(pairedAlignments);

    // Populating Data Structures
    fillOverlapMap(overlapResults);
    fillSamMap(samEntries);

    // Stage3: calculate the number of support reads for each fusion
    // Calculation support reads
    std::unordered_map<std::string, int> splitReadsCount = countSplitReads(overlapMap, samMap,options);

    auto spanReads = collectSpanReads(overlapMap, samMap, options);

    auto splitReads = collectSplitReads(overlapMap, samMap, options);

    std::unordered_map<std::string, int> spanReadsCount = countSpanReadPairs(overlapMap, samMap);

    // Output the number of split reads and spanning read pairs for all query names
    std::cout << get_time_string() << " Calculate the number of split reads and spanning read pairs for all querys " << std::endl;
    Logger::Info(get_time_string() + " Calculate the number of split reads and spanning read pairs for all querys");
    for (const auto& entry : splitReadsCount) {
        Logger::Info(get_time_string() + "\tQuery:" + entry.first + "\tsplit reads: " + std::to_string(entry.second) + "\tspan reads: " + std::to_string(spanReadsCount[entry.first]));
    }

    // Filter queries that meet the criteria
    std::unordered_set<std::string> validQueries = filterQueries(splitReadsCount, spanReadsCount);

    // Get the alignment corresponding to the valid query
    std::vector<alignment_t> relevantAlignments = filterAlignmentsByValidQueries(overlaps_same_strand, validQueries);

    // Optional: Output relevant alignment information
    std::cout << get_time_string() << " Filtered chosen alignments based on support reads and spanning read pairs: " << relevantAlignments.size() << std::endl;
    Logger::logFile << get_time_string() << " Filtered chosen alignments: " << relevantAlignments.size() << std::endl;


    std::vector<coordination_t> finalcoordinations = extractSimpleCoordinations(relevantAlignments);


    // Print the final coordinations
    std::cout << get_time_string() << " Final coordinations after merging " << std::endl;

    // Stage4: annotation of the gene and output the remain results
    // Annotation
    // Creating an Annotator Instance
    GeneAnnotator annotator(options.gtf_path);

    auto filtered_final_coordinations = keepOnlyTwoParts(finalcoordinations);
    auto annotations = annotator.annotateAlignments(filtered_final_coordinations);

    keepOnlyTwoAnnotations(annotations);



    // apply filters to the final results
    std::vector<result_t> final_results;

    processAnnotations(annotations, final_results);
     std::cout << get_time_string() << " Filtering the invalid annotations " << "(remaining=" << final_results.size() << ")" << std::endl;
    Logger::logFile << get_time_string() << " Filtering the invalid annotations " << "(remaining=" << final_results.size() << ")" << std::endl;
    integrateReadCounts(final_results, splitReadsCount, spanReadsCount);


    removeEmptyGenes(final_results);

    std::cout << get_time_string() << " Filtering the empty annotations " << "(remaining=" << final_results.size() << ")" << std::endl;
    Logger::logFile << get_time_string() << " Filtering the empty annotations " << "(remaining=" << final_results.size() << ")" << std::endl;

    int readLength = options.read_length; // Get the read length from the input
    std::vector<std::pair<int, int>> readLengths;

   // Apply filters to final_results
    std::vector<result_t> discarded_duplicates;
    if (options.filters.at("duplicates")) {
        final_results = filter_duplicates(final_results); // Filter duplicates

        for (const auto& result : final_results) {
            if (result.filter_status == "discarded_duplicates") {
                discarded_duplicates.push_back(result);

            }
        }
        final_results.erase(
               std::remove_if(final_results.begin(), final_results.end(),
                              [](const result_t& r) { return r.filter_status != "kept"; }),
               final_results.end());
        std::cout << get_time_string() << " Filtering fusions with duplicates " << "(remaining=" << final_results.size() << ")" << std::endl;
        Logger::logFile << get_time_string() << " Filtering fusions with duplicates " << "(remaining=" << final_results.size() << ")" << std::endl;
    }

    std::vector<result_t> discarded_mt;
    if (options.filters.at("mt")) {
        final_results = filter_mt(final_results); // Filter by minimum support

        for (const auto& result : final_results) {
            if (result.filter_status == "discarded_mt") {
                discarded_mt.push_back(result);
            }
        }

        final_results.erase(
               std::remove_if(final_results.begin(), final_results.end(),
                              [](const result_t& r) { return r.filter_status != "kept"; }),
               final_results.end());
        std::cout << get_time_string() << " Filtering fusions in MT, mitochondrial " << "(remaining=" << final_results.size() << ")" << std::endl;
        Logger::logFile << get_time_string() << " Filtering fusions in MT, mitochondrial " << "(remaining=" << final_results.size() << ")" << std::endl;
    }

    std::vector<result_t> discarded_long_gap;
    if (options.filters.at("long_gap")) {
        final_results = filter_long_gap(final_results, options.long_gap_threshold, options.short_segment_threshold); // Filter long gap

        for (const auto& result : final_results) {
            if (result.filter_status == "discarded_long_gap") {
                discarded_long_gap.push_back(result);
            }
        }
        final_results.erase(
       std::remove_if(final_results.begin(), final_results.end(),
                      [](const result_t& r) { return r.filter_status != "kept"; }),
       final_results.end());
        std::cout << get_time_string() << " Filtering fusions with long gaps " << "(remaining=" << final_results.size() << ")" << std::endl;
        Logger::logFile << get_time_string() << " Filtering fusions with long gaps " << "(remaining=" << final_results.size() << ")" << std::endl;

    }

    std::vector<result_t> discarded_internal_tandem_duplication;
    if (options.filters.at("internal_tandem_duplication")) {
        final_results = filter_internal_tandem_duplication(final_results, options.max_itd_length, options.min_split_reads, options.min_itd_fraction);

        for (const auto& result : final_results) {
            if (result.filter_status == "discarded_internal_tandem_duplication") {
                discarded_internal_tandem_duplication.push_back(result);
            }
        }
        final_results.erase(
       std::remove_if(final_results.begin(), final_results.end(),
                      [](const result_t& r) { return r.filter_status != "kept"; }),
       final_results.end());
        std::cout << get_time_string() << " Filtering internal tandem duplications " << "(remaining=" << final_results.size() << ")" << std::endl;
        Logger::logFile << get_time_string() << " Filtering internal tandem duplications " << "(remaining=" << final_results.size() << ")"  << std::endl;
    }

    std::vector<result_t> discarded_homologs;
    if (options.filters.at("homologs")) {
        // Create FilterHomologs instance and apply it
        FilterHomologs filter_homologs( final_results, pairedAlignments, alignments);  // Pass pairedAlignments and alignments
        final_results = filter_homologs.filter_homologs(); // Filter fusion events based on the logic

        for (const auto& result : final_results) {
            if (result.filter_status == "discarded_homologs") {
                discarded_homologs.push_back(result);
            }
        }
        final_results.erase(
       std::remove_if(final_results.begin(), final_results.end(),
                      [](const result_t& r) { return r.filter_status != "kept"; }),
       final_results.end());
        std::cout << get_time_string() << " Filtering fusions with homologs " << "(remaining=" << final_results.size() << ")" << std::endl;
        Logger::logFile << get_time_string() << " Filtering fusions with homologs " << "(remaining=" << final_results.size() << ")"  << std::endl;

    }

    std::vector<result_t> discarded_small_fragments;
    if (options.filters.at("small_fragments")) {
        final_results = filter_small_fragments(final_results, options); // Adjust threshold as needed

        for (const auto& result : final_results) {
            if (result.filter_status == "discarded_small_fragments") {
                discarded_small_fragments.push_back(result);
            }
        }
        final_results.erase(
       std::remove_if(final_results.begin(), final_results.end(),
                      [](const result_t& r) { return r.filter_status != "kept"; }),
       final_results.end());
        std::cout << get_time_string() << " Filtering fusions with small fragements " << "(remaining=" << final_results.size() << ")" << std::endl;
        Logger::logFile << get_time_string() << " Filtering fusions with small fragements " << "(remaining=" << final_results.size() << ")"  << std::endl;
    }



    std::vector<result_t> discarded_min_support;
    if (options.filters.at("min_support")) {
        final_results = filter_min_support(final_results, options.min_span_reads, options.min_split_reads); // Filter by minimum support

        for (const auto& result : final_results) {
            if (result.filter_status == "discarded_min_support") {
                discarded_min_support.push_back(result);
            }
        }

        final_results.erase(
               std::remove_if(final_results.begin(), final_results.end(),
                              [](const result_t& r) { return r.filter_status != "kept"; }),
               final_results.end());
        std::cout << get_time_string() << " Filtering fusions under minimum support reads " << "(remaining=" << final_results.size() << ")" << std::endl;
        Logger::logFile << get_time_string() << " Filtering fusions under minimum support reads " << "(remaining=" << final_results.size() << ")" << std::endl;
    }

    // Output discarded results
    std::vector<result_t> discarded_results;
    discarded_results.insert(discarded_results.end(), discarded_duplicates.begin(), discarded_duplicates.end());
    discarded_results.insert(discarded_results.end(), discarded_mt.begin(), discarded_mt.end());
    discarded_results.insert(discarded_results.end(), discarded_long_gap.begin(), discarded_long_gap.end());
    discarded_results.insert(discarded_results.end(), discarded_internal_tandem_duplication.begin(), discarded_internal_tandem_duplication.end());
    discarded_results.insert(discarded_results.end(), discarded_homologs.begin(), discarded_homologs.end());
    discarded_results.insert(discarded_results.end(), discarded_small_fragments.begin(), discarded_small_fragments.end());
    discarded_results.insert(discarded_results.end(), discarded_min_support.begin(), discarded_min_support.end());



    std::string filename = "./known_fusions.tsv";
    std::vector<known_fusion_t> known_fusions = load_known_fusions(filename);

    // Filter out known fusions
    std::vector<result_t> recovered_fusions = recover_fusions(discarded_results, known_fusions);

    // Add recovered fusions to kept_results
    final_results.insert(final_results.end(), recovered_fusions.begin(), recovered_fusions.end());

    // Output final results count
    std::cout << get_time_string() << " Final fusion genes count: " << final_results.size() << std::endl;



    // integrate coverage counts into final_results
    // Maps to hold results
    std::unordered_map<std::string, float> averageCoverageMap;
    std::unordered_map<std::string, std::unordered_map<int, int>> perBaseCoverageMap;
    // This structure will hold per-base coverage in a separate data frame
    std::unordered_map<std::string, std::unordered_map<int, int>> newPerBaseCoverageMap;


    // Process coverage data
    processCoverage(samEntries, overlapResults, averageCoverageMap, perBaseCoverageMap);

    // Iterate through final results and add coverage info
    for (auto& result : final_results) {
        const std::string& contig = result.contig; // Assuming result_t has a contig field

        // If contig is found in averageCoverageMap, add the average coverage
        if (averageCoverageMap.find(contig) != averageCoverageMap.end()) {
            result.coverage = averageCoverageMap[contig]; // Assuming result_t has this field
        }

        // Store per-base coverage separately in newPerBaseCoverageMap
        if (perBaseCoverageMap.find(contig) != perBaseCoverageMap.end()) {
            newPerBaseCoverageMap[contig] = perBaseCoverageMap[contig];
        }
    }

    // Output average coverage only for contigs in final_results
    for (const auto& result : final_results) {
        const std::string& contig = result.contig;
    }


    // Save the data to a TSV file
    savePerBaseCoverageToTSV(newPerBaseCoverageMap, options.output + "/" + options.prefix + ".per_base_coverage.tsv");

    // Output the result and fill in the structure
    for (const auto& result : final_results) {
        Logger::logFile << get_time_string() << " Final fusion genes:\t"
                        << result.contig << "\t"
                        << result.gene1 << "\t"
                        << result.gene2 << "\t"
                        << result.gene_id1 << "\t"
                        << result.gene_id2 << "\t"
                        << result.splitReadsCount << "\t"
                        << result.spanReadsCount  << "\t"
                        << result.chromosome1  << "\t"
                        << result.chromosome2  << "\t"
                        << result.direction1 << "\t"
                        << result.direction2 << "\t"
                        << result.tstart1 << "\t"
                        << result.tend1 << "\t"
                        << result.tstart2 << "\t"
                        << result.tend2 << "\t"
                        << result.regiontype1 << "\t"
                        << result.regiontype2 << "\t"
                        << result.filter_status << "\t"
                        << std::endl;
    }

    for (const auto& result : discarded_results) {
        Logger::logFile << get_time_string() << " Discarded fusion genes:\t"
                       << result.contig << "\t"
                       << result.gene1 << "\t"
                       << result.gene2 << "\t"
                       << result.gene_id1 << "\t"
                       << result.gene_id2 << "\t"
                       << result.splitReadsCount << "\t"
                       << result.spanReadsCount  << "\t"
                       << result.chromosome1  << "\t"
                       << result.chromosome2  << "\t"
                       << result.direction1 << "\t"
                       << result.direction2 << "\t"
                       << result.tstart1 << "\t"
                       << result.tend1 << "\t"
                       << result.tstart2 << "\t"
                       << result.tend2 << "\t"
                       << result.regiontype1 << "\t"
                       << result.regiontype2 << "\t"
                       << result.filter_status << "\t"
                       << std::endl;
    }

    // Write supporting reads into fq files
    writeEvidenceToFastq(final_results, splitReads, spanReads, options.prefix + ".evidence", options);

    // Write to TSV file
    writeToTSV(final_results, options.output + "/" + options.prefix + ".fusion_list.tsv");
    std::cout << get_time_string() << " Write fusion list into file:" << " '" << options.output << "/" << options.prefix << ".fusion_list.tsv" << "' " << std::endl;

    writeToTSV(discarded_results, options.output + "/" + options.prefix + ".discarded_fusions.tsv");
    std::cout << get_time_string() << " Write discarded fusion list into file:" << " '" << options.output << "/" << options.prefix << ".discarded_fusions.tsv" << "' " << std::endl;



    // Stage5: calculation of the coverage and  prediction of breakpoints
    //predict the possible breakpoints

    float averageRate = 1.0; // Average rate of coverage per base
    float qualityThreshold = 0.1; // Minimum quality score for breakpoint detection

    // evaluate breakpoints
    auto breakpoints = evaluateBreakpoints(newPerBaseCoverageMap, readLengths, averageRate, qualityThreshold, options);

    // Stage5: calculation of the coverage and  prediction of breakpoints
    // Predict the possible breakpoints
    std::cout << get_time_string() << " Stage5: Calculation of the coverage and  prediction of breakpoints " << std::endl;
    std::cout << get_time_string() << " Saved per-base coverage to " << " '" << options.output << "/" << options.prefix << ".per_base_coverage.tsv" << "' " << std::endl;

    Logger::Info(get_time_string() + " Stage5: Calculation of the coverage and  prediction of breakpoints ");
    Logger::Info(get_time_string() + " Saved per-base coverage to " + " '" + options.output + "/" + options.prefix + ".per_base_coverage.tsv" + "' ");


    // print resource usage stats end exit
    time_t end_time;
    time(&end_time);
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    #ifdef __APPLE__

    #define RU_MAXRSS_UNIT 1024.0*1024*1024
    #else
    #define RU_MAXRSS_UNIT 1024.0*1024
    #endif



    std::cout << get_time_string() << " Done "
         << "(elapsed time=" << get_hhmmss_string(difftime(end_time, start_time)) << ", "
         << "CPU time=" << get_hhmmss_string(usage.ru_utime.tv_sec + usage.ru_stime.tv_sec) << ", "
         << "peak memory=" << std::setprecision(3) << (usage.ru_maxrss/(RU_MAXRSS_UNIT)) << "gb)" << std::endl;

    // Record before program close
    Logger::Info(get_time_string() + " The program ends");

    // Close log file
    Logger::close();

}