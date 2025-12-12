//
// Created by xinwei on 5/19/24.
//
#include <algorithm>
#include <cstring>
#include <iostream>
#include <libgen.h>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <unordered_map>


#include "options.h"
#include "common.h"



std::string wrap_help(const std::string& option, const std::string& text, const unsigned short int max_line_width) {
    std::istringstream iss(text);
    std::string result = "       " + option + ", ";
    std::string indent(result.length(), ' ');
    unsigned short int line_width = result.length();

    std::string word;
    while (iss >> word) {

        if (line_width + word.length() > max_line_width) {
            result += "\n" + indent;
            line_width = indent.length();
        }
        result += (line_width > indent.length() ? " " : "") + word;
        line_width += word.length() + 1;
    }

    return result;
}


std::string wrap_help2(const std::string& text, const unsigned short int max_line_width) {
    std::istringstream iss(text);
    std::string result = "          ";
    std::string indent(result.length(), ' ');
    unsigned short int line_width = result.length();

    std::string word;
    while (iss >> word) {

        if (line_width + word.length() > max_line_width) {
            result += "\n" + indent;
            line_width = indent.length();
        }
        result += (line_width > indent.length() ? " " : "") + word;
        line_width += word.length() + 1;
    }

    return result + "\n";
}


bool output_directory_exists(const std::string& output_dir) {
    if (output_dir.empty())
        return false;

    struct stat file_info;
    return stat(output_dir.c_str(), &file_info) == 0 && S_ISDIR(file_info.st_mode);
}


bool validate_int(const char* optarg, int& value, const int min_value, const int max_value) {
    if (!str_to_int(optarg, value))
        return false;
    return value >= min_value && value <= max_value;
}

bool validate_int(const char* optarg, unsigned int& value, const unsigned int min_value, const unsigned int max_value) {
    int signed_int;
    if (validate_int(optarg, signed_int, min_value, max_value)) {
        value = signed_int;
        return true;
    } else
        return false;
}

bool validate_float(const char* optarg, float& value, const float min_value, const float max_value) {
    if (!str_to_float(optarg, value))
        return false;
    return value >= min_value && value <= max_value;
}


void parse_fastq_files(const char* optarg, std::vector<std::string>& fastq_files) {
    std::istringstream ss(optarg);
    std::string fastq;
    while (ss >> fastq) {
        fastq_files.push_back(fastq);
    }
}



options_t get_default_options() {
    options_t options;

    options.max_alignment_count = 100;
    options.min_identity_fract = 0.95;
    options.min_score_total = 95;
    options.min_score_each = 10;
    options.coverage_differ = 0.85;

    options.max_pair_combination = 3;
    options.threads = 4;
    options.max_overlap_size = 8;
    options.max_gap_size = 2;
    options.min_edge_length = 20;
    options.size_ratio_threshold = 0.1;
    options.edge_unaligned = 50;
    options.inclusion_fraction_weight = 1;
    options.overlap_fraction_weight = 1;
    options.size_weight = 1;

    options.long_gap_threshold = 200000;
    options.short_segment_threshold = 35;
    options.min_split_reads = 3;
    options.min_span_reads = 1;
    options.max_itd_length = 1000;
    options.min_itd_fraction = 0.5;

    for (size_t i = 0; i < FILTERS.size(); ++i)
        if (i != FILTER_none)
            options.filters[FILTERS[i]] = true;

    return options;
}

void print_usage() {

    options_t default_options = get_default_options();
    std::string valid_filters;

    for (auto i = default_options.filters.begin(); i != default_options.filters.end(); ++i) {
        if (i != default_options.filters.begin())
            valid_filters += "";
        valid_filters += i->first;
    }

    std::cout << std::endl
              << "DenovoFusion Version: " << DENOVOFUSION_VERSION << std::endl
              << "------------------------------------------------------------------------" << std::endl
              << wrap_help2("DenovoFusion is a specialized bioinformatics tool designed for the identification of chimeric sequences, "
                 "including gene fusions, from de novo assembled genomic data.") << std::endl
              << "------------------------------------------------------------------------" << std::endl
              << "Usage: " << std::endl
              << "DenovoFusion -m blat -i Aligned.psl \\" << std::endl
              << "             -a assembly.fa -r read_length -g path/to/.gtf \\" << std::endl
              << "             -1 01.fastq.gz -2 02.fastq.gz \\" << std::endl
              << "             -o path/to/results -p prefix \\" << std::endl
              << "             [OPTIONS]" << std::endl
              << "------------------------------------------------------------------------" << std::endl
              << "Help" << std::endl
              << wrap_help("-h","--help") << std::endl
              << wrap_help2("Display this help message and exit.") << std::endl

              << "Mandatory parameters" << std::endl
              << wrap_help("-1","--fastq1") << std::endl
              << wrap_help2("Specifies the FASTQ file for the forward read (e.g.01.fastq.gz). Multiple input "
                                                "should be divided by comma, as in Bowtie2 input format.") << std::endl
              << wrap_help("-2","--fastq2") << std::endl
              << wrap_help2("Specifies the FASTQ file for the reverse read (e.g.02.fastq.gz). Multiple input "
                                                "should be divided by comma, as in Bowtie2 input format.") << std::endl
              << wrap_help("-a","--input_assembly") << std::endl
              << wrap_help2("The contig file from de novo assembly in FASTA format. When using PSL mode, "
                                                "please use the processed (cut) contigs instead of the original contigs. ") << std::endl
              << wrap_help("-g","--gtf_path") << std::endl
              << wrap_help2("GTF annotation file with gene annotation.") << std::endl
              << wrap_help("-i","--input_file") << std::endl
              << wrap_help2("Alignment files in SAM, PAF, or PSL format, generated by alignment tools, "
                                                "This serves as the primary input for DenovoFusion's analysis.") << std::endl
              << wrap_help("-m"," --method") << std::endl
              << wrap_help2("DenovoFusion applies the same alignment method previously used for contig "
                                                "alignment. The supported methods include BLAT, Minimap2sam (SAM output), and Minimap2paf"
                                                " (PAF output). If the input file is in BAM format, it is recommended to convert it "
                                                "to SAM format before proceeding with the analysis.") << std::endl
              << wrap_help("-o","--output") << std::endl
              << wrap_help2("The output file path, which includes separate files for fusion genes, coverage information, and log details.") << std::endl
              << wrap_help("-p","--prefix") << std::endl
              << wrap_help2("Prefix for all temporary and result files; you can use the sample name here.") << std::endl
              << wrap_help("-r","--read_length") << std::endl
              << wrap_help2("Read length from the raw FASTQ file.") << std::endl

              << "Optional parameters" << std::endl
              << wrap_help("-A","--edge-unaligned") << std::endl
              << wrap_help2("Maximum number of unaligned bases at the head or tail of a contig. "
                                                "(Default: 50).") << std::endl
              << wrap_help("-c","--max_alignment_count") << std::endl
              << wrap_help2("Maximum number of alignments considered for each contig."
                                                " Larger values may increase processing time (Default: 100).") << std::endl
              << wrap_help("-d","--min_identity_fract") << std::endl
              << wrap_help2("Minimum alignment identity threshold (Default: 0.95).") << std::endl
              << wrap_help("-e","--min_score_each") << std::endl
              << wrap_help2("Minimum score for each alignment (Default: 20).") << std::endl
              << wrap_help("-E","--min-edge-length") << std::endl
              << wrap_help2("Minimum edge length of split reads in re-alignment; adjust according to"
                                                "read length (Default: 20).") << std::endl
              << wrap_help("-f","--max-gap-size") << std::endl
              << wrap_help2("Maximum number of gap bases allowed in paired alignments (Default: 2).") << std::endl
              << wrap_help("-G","--short-segment-threshold") << std::endl
              << wrap_help2("Minimum segment length requirement for fusion parts (Default: 35).") << std::endl
              << wrap_help("-I","--inclusion-fraction-weight") << std::endl
              << wrap_help2("Weight of the inclusion fraction in alignment score calculation (Default: 1).") << std::endl
              << wrap_help("-l","--max-overlap-size") << std::endl
              << wrap_help2("Maximum number of overlapping bases allowed in paired alignments (Default: 8).") << std::endl
              << wrap_help("-n","--max-pair-combination") << std::endl
              << wrap_help2("Maximum number of alignment combinations per contig. Larger values may "
                                                "increase processing time (Default: 3).") << std::endl
              << wrap_help("-N","--min-span-reads") << std::endl
              << wrap_help2("Minimum number of spanning reads required as support (Default: 1).") << std::endl
              << wrap_help("-O","--overlap-fraction-weight") << std::endl
              << wrap_help2("Weight of the overlap fraction in alignment score calculation (Default: 1).") << std::endl
              << wrap_help("-P","--min-split-reads") << std::endl
              << wrap_help2("Minimum number of split reads required as support (Default: 3).") << std::endl
              << wrap_help("-q","--threads") << std::endl
              << wrap_help2("Number of threads to use (Default: 4).") << std::endl
              << wrap_help("-s","--min_score_total") << std::endl
              << wrap_help2("Minimum total score for combined alignments (Default: 95).") << std::endl
              << wrap_help("-S","--size-weight") << std::endl
              << wrap_help2("Weight of alignment count in alignment score calculation (Default: 1).") << std::endl
              << wrap_help("-T","--long-gap-threshold") << std::endl
              << wrap_help2("Gap threshold for candidate fusions on the same chromosome (Default: 200000).") << std::endl
              << wrap_help("-v","--coverage-differ") << std::endl
              << wrap_help2("Coverage ratio threshold between adjacent bases in breakpoint prediction (Default: 0.85).") << std::endl
              << wrap_help("-z","--size-ratio-threshold") << std::endl
              << wrap_help2("Proportion of fusion part 1 and part 2 must remain within an "
                                                "acceptable range (Default: 0.1).") << std::endl
    ;
}


options_t option_parser(int argc, char **argv) {

    options_t options = get_default_options();
    // throw error when first argument is not prefixed with a dash
    // for some reason getopt does not detect this error and simply skips the argument
    crash(argc > 1 && (std::string(argv[1]).empty() || argv[1][0] != '-'), "cannot interpret the first argument: " + std::string(argv[1]));

    static struct option long_options[] = {
    {"method", required_argument, nullptr, 'm'},          // --method (short option -m)
    {"input", required_argument, nullptr, 'i'},           // --input (short option -i)
    {"assembly", required_argument, nullptr, 'a'},        // --assembly (short option -a)
    {"output", required_argument, nullptr, 'o'},          // --output (short option -o)
    {"gtf", required_argument, nullptr, 'g'},             // --gtf (short option -f)
    {"prefix", required_argument, nullptr, 'p'},          // --prefix (short option -p)
    {"bowtie2", required_argument, nullptr, 'b'},         // --bowtie2 (short option -b)
    {"fastq1", required_argument, nullptr, '1'},          // --fastq1 (short option -1)
    {"fastq2", required_argument, nullptr, '2'},          // --fastq2 (short option -2)
    {"max-alignment-count", required_argument, nullptr, 'c'},  // --max-alignment-count (short option -c)
    {"min-identity-fract", required_argument, nullptr, 'd'},   // --min-identity-fract (short option -d)
    {"min-score-total", required_argument, nullptr, 's'},  // --min-score-total (short option -s)
    {"min-score-each", required_argument, nullptr, 'e'},   // --min-score-each (short option -e)
    {"coverage-differ", required_argument, nullptr, 'v'},  // --coverage-differ (short option -v)
    {"size-ratio-threshold", required_argument, nullptr, 'z'}, // --size-ratio-threshold (short option -z)
    {"edge-unaligned", required_argument, nullptr, 'A'},   // --edge-unaligned (short option -A)
    {"max-pair-combination", required_argument, nullptr, 'n'}, // --max-pair-combination (short option -n)
    {"threads", required_argument, nullptr, 'q'},          // --threads (short option -q)
    {"max-overlap-size", required_argument, nullptr, 'l'}, // --max-overlap-size (short option -l)
    {"max-gap-size", required_argument, nullptr, 'f'},     // --max-gap-size (short option -g)
    {"read-length", required_argument, nullptr, 'r'},      // --read-length (short option -r)
    {"min-edge-length", required_argument, nullptr, 'E'},  // --min-edge-length (short option -E)
    {"inclusion-fraction-weight", required_argument, nullptr, 'I'}, // --inclusion-fraction-weight (short option -I)
    {"overlap-fraction-weight", required_argument, nullptr, 'O'},   // --overlap-fraction-weight (short option -O)
    {"size-weight", required_argument, nullptr, 'S'},        // --size-weight (short option -S)
    {"long-gap-threshold", required_argument, nullptr, 'T'}, // --long-gap-threshold (short option -T)
    {"short-segment-threshold", required_argument, nullptr, 'G'}, // --short-segment-threshold (short option -G)
    {"min-split-reads", required_argument, nullptr, 'P'},  // --min-split-reads (short option -P)
    {"min-span-reads", required_argument, nullptr, 'N'},   // --min-span-reads (short option -N)
    {"help", no_argument, nullptr, 'h'},                   // --help (short option -h)
    {nullptr, 0, nullptr, 0} // Sentinel value
};



    // parse arguments
    opterr = 0;
    int c;
    std::string junction_suffix(".junction");
    std::unordered_map<char,unsigned int> duplicate_arguments;
    const std::string valid_arguments = "1:2:c:x:q:d:g:r:G:o:w:l:O:t:p:a:b:k:s:i:v:f:E:S:m:L:H:D:R:A:M:K:V:F:U:Q:e:T:C:l:z:Z:uXIh";
    // Use getopt_long to handle both short and long options
    while ((c = getopt_long(argc, argv, valid_arguments.c_str(), long_options, nullptr)) != -1) {
        // Throw error if the same argument is specified more than once
        duplicate_arguments[c]++;
        crash(duplicate_arguments[c] > 1, "option -" + std::string(1, (char)c) + " specified too often");


        switch (c) {
            case 'm':
                options.input_type = optarg;
                // Convert input_type to lowercase
                std::transform(options.input_type.begin(), options.input_type.end(), options.input_type.begin(), ::tolower);
                if (!(options.input_type == "blat" || options.input_type == "minimap2paf" || options.input_type == "minimap2sam")) {
                crash(true, "invalid alignment method: " + options.input_type + ", must be BLAT, Minimap2paf, or Minimap2sam");
                }
                break;
            case 'i':
                options.input_file = optarg;
                crash(access(options.input_file.c_str(), R_OK), "file not found/readable: " + options.input_file);
                break;
            case 'a':
                options.input_assembly = optarg;
                crash(access(options.input_assembly.c_str(), R_OK), "file not found/readable: " + options.input_assembly);
                break;
            case 'o':
                options.output = optarg;
                crash(!output_directory_exists(options.output), "parent directory of output file '" + options.output + "' does not exist");
                break;
            case 'g':
                options.gtf_path = optarg;
                crash(access(options.gtf_path.c_str(), R_OK), "file not found/readable: " + options.gtf_path);
                break;
            case 'p':
                options.prefix = optarg;
                break;

            case '1':
                parse_fastq_files(optarg, options.input_fastq1);
                break;
            case '2':
                parse_fastq_files(optarg, options.input_fastq2);
                break;

            case 'c':
                crash(!validate_int(optarg, options.max_alignment_count, 1,200), "invalid argument to -" + ((char) c));
                break;
            case 'd':
                crash(!validate_float(optarg, options.min_identity_fract, 0.6,0.99), "invalid argument to -" + ((char) c));
                break;
            case 's':
                crash(!validate_int(optarg, options.min_score_total, 60,99), "invalid argument to -" + ((char) c));
                break;
            case 'e':
                crash(!validate_int(optarg, options.min_score_each, 0,20), "invalid argument to -" + ((char) c));
                break;
            case 'v':
                crash(!validate_float(optarg, options.coverage_differ, 0.5,0.99), "invalid argument to -" + ((char) c));
                break;
            case 'z':
                crash(!validate_float(optarg, options.size_ratio_threshold, 0.01,0.2), "invalid argument to -" + ((char) c));
                break;
            case 'A':
                crash(!validate_int(optarg, options.edge_unaligned, 1,100), "invalid argument to -" + ((char) c));
                break;
            case 'n':
                crash(!validate_int(optarg, options.max_pair_combination, 1,5), "invalid argument to -" + ((char) c));
                break;
            case 'q':
                crash(!validate_int(optarg, options.threads, 1,64), "invalid argument to -" + ((char) c));
                break;
            case 'l':
                crash(!validate_int(optarg, options.max_overlap_size, 1,100), "invalid argument to -" + ((char) c));
                break;
            case 'f':
                crash(!validate_int(optarg, options.max_gap_size, 1,100), "invalid argument to -" + ((char) c));
                break;
            case 'r':
                crash(!validate_int(optarg, options.read_length, 20,300), "invalid argument to -" + ((char) c));
                break;
            case 'E':
                crash(!validate_int(optarg, options.min_edge_length, 1,100), "invalid argument to -" + ((char) c));
                break;
            case 'I':
                crash(!validate_float(optarg, options.inclusion_fraction_weight, 0,1), "invalid argument to -" + ((char) c));
                break;
            case 'O':
                crash(!validate_float(optarg, options.overlap_fraction_weight, 0,1), "invalid argument to -" + ((char) c));
                break;
            case 'S':
                crash(!validate_float(optarg, options.size_weight, 0,1), "invalid argument to -" + ((char) c));
                break;
            case 'T':
                crash(!validate_int(optarg, options.long_gap_threshold, 100000,1000000), "invalid argument to -" + ((char) c));
                break;
            case 'G':
                crash(!validate_int(optarg, options.short_segment_threshold, 1,100), "invalid argument to -" + ((char) c));
                break;
            case 'P':
                crash(!validate_int(optarg, options.min_split_reads, 1,50), "invalid argument to -" + ((char) c));
                break;
            case 'N':
                crash(!validate_int(optarg, options.min_span_reads, 1,50), "invalid argument to -" + ((char) c));
                break;
            case 'h':
                print_usage();
                exit(0);
                break;
            default:
                crash(valid_arguments.find(std::string(1, (char) optopt) + ":") != std::string::npos, "option -" + ((char) optopt) + " requires an argument");
                crash(true, "unknown option: -" + ((char) optopt));
                break;
        }

        crash(optind < argc && (std::string(argv[optind]).empty() || argv[optind][0] != '-'), "option -" + ((char) c) + " has too many arguments (arguments with blanks must be wrapped in quotes)");

    }

    // check for mandatory arguments
    if (argc == 1) {
        print_usage();
        crash(true, "no arguments given");
    }
    crash(options.input_type.empty(), "missing mandatory option -m");
    crash(options.input_file.empty(), "missing mandatory option -i");
    crash(options.input_assembly.empty(), "missing mandatory option -a");
    crash(options.output.empty(), "missing mandatory option -o");
    crash(options.gtf_path.empty(), "missing mandatory option -g");
    crash(options.input_fastq1.empty(), "missing mandatory option -1");
    crash(options.input_fastq2.empty(), "missing mandatory option -2");
    crash(options.prefix.empty(), "missing mandatory option -p");
    crash(options.read_length == 0, "missing mandatory option -r");


    return options;
}













