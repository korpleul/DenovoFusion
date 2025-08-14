//
// Created by xinwei on 6/20/24.
//

#ifndef FUSION_DETECTION_2_BOWTIE2_H
#define FUSION_DETECTION_2_BOWTIE2_H


#include <iostream>
#include <filesystem>
#include <cstdlib>

#include "options.h"


void build_bowtie2_index(const options_t& options);

std::string generate_bowtie2_command(const options_t& options);

void run_bowtie2(const options_t& options);




#endif //FUSION_DETECTION_2_BOWTIE2_H
