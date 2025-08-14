cmake_minimum_required(VERSION 3.10)
project(DenovoFusion)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(DenovoFusion main.cpp
        src/error.h
        src/log.h
        src/log.cpp
        src/options.cpp
        src/psl.cpp
        src/coord_pair.cpp
        src/alignment.cpp
        src/alignment.h
        src/aligned_coord_pair.h
        src/aligned_coord_pair.cpp
        src/psl_functions.cpp
        src/candidate_group.cpp
        src/candidate_group.h
        src/alignments_chosen.cpp
        src/alignments_chosen.h
        src/output_fusions.cpp
        src/output_fusions.h
        src/fasta.cpp
        src/fasta.h
        src/sam.cpp
        src/sam.h
        src/realign_support.cpp
        src/realign_support.h
        src/bowtie2.cpp
        src/bowtie2.h
        src/overlap.cpp
        src/overlap.h
        src/annotation.cpp
        src/annotation.h
        src/coverage.cpp
        src/coverage.h
        src/breakpoint.cpp
        src/breakpoint.h
        src/filter_duplicates.cpp
        src/filter_duplicates.h
        src/filter_long_gap.cpp
        src/filter_long_gap.h
        src/filter_min_support.cpp
        src/filter_min_support.h
        src/recover_known_fusion.cpp
        src/recover_known_fusion.h
        src/paf.cpp
        src/paf.h
        src/run_blat.cpp
        src/run_blat.h
        src/utils.cpp
        src/utils.h
        src/run_minimap2paf.cpp
        src/run_minimap2paf.h
        src/run_minimap2sam.cpp
        src/run_minimap2sam.h
        src/filter_small_fragments.cpp
        src/filter_small_fragments.h
        src/filter_homologs.cpp
        src/filter_homologs.h
        src/filter_edge_unaligned.cpp
        src/filter_edge_unaligned.h
        src/filter_internal_tandem_duplication.cpp
        src/filter_internal_tandem_duplication.h
        src/filter_mt.cpp
        src/filter_mt.h


)
