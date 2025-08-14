//
// Created by xinwei on 6/20/24.
//

#include "bowtie2.h"

void build_bowtie2_index(const options_t& options) {
    // Get user's home directory from environment variable
    const char* homeDir = getenv("HOME");
    if (!homeDir) {
        std::cerr << "Error: HOME directory not found." << std::endl;
        return;
    }


    std::string bowtie2Path =  "./bowtie2-2.5.4-linux-x86_64/bowtie2-build";
    std::string fastaFilePath = options.output + "/" + options.prefix + ".chosen.fasta";
    std::string indexOutputPath = options.output + "/" + options.prefix + "_idx/" + options.prefix;


    std::string command = bowtie2Path + " " + fastaFilePath + " " + indexOutputPath;
    int result = system(command.c_str());

    if (result != 0) {
        std::cerr << "Bowtie2 indexing failed." << std::endl;
    }
}

std::string generate_bowtie2_command(const options_t& options) {
    std::stringstream command;

    command << "./bowtie2-2.5.4-linux-x86_64/bowtie2" << " -t -p " << options.threads;
    command << " -x " << options.output << "/" << options.prefix << "_idx/" << options.prefix;

    command << " -1";
    for (const auto& fastq : options.input_fastq1) {
        command << " " << fastq;
    }

    command << " -2";
    for (const auto& fastq : options.input_fastq2) {
        command << " " << fastq;
    }

    // Other parameters remain unchanged
    command << " --no-mixed --no-discordant --no-contain --no-overlap --no-head --no-sq -k1 --no-unal";

    // Set the output file path
    command << " -S " << options.output << "/" << options.prefix << ".sam";

    return command.str();
}

void run_bowtie2(const options_t& options) {
    std::string command = generate_bowtie2_command(options);
    std::cout << "Running command: " << command << std::endl;
    int result = system(command.c_str());

    if (result != 0) {
        std::cerr << "Bowtie2 alignment failed." << std::endl;
    }
}
