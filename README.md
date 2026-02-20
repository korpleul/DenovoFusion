# $\textbf{\color{blue}DenovoFusion}$

Cite : A de novo assembly based fusion gene detection concept based on DNA-seq data of 100 Ewing sarcoma cases
Xinwei Zhao , Marc Hotfilder , Eberhard Korsching
Briefings in Bioinformatics, Volume 27, Issue 1, January 2026, [bbag059](https://doi.org/10.1093/bib/bbag059)
Published: 16 February 2026

## Overview

*DenovoFusion* is a C++17-based bioinformatics tool designed for detecting fusion genes from de novo assembled contigs derived from DNA-seq short reads. Unlike traditional alignment-based approaches, *DenovoFusion* leverages the power of whole-genome de novo assembly to reconstruct longer contiguous sequences, allowing for more accurate identification of structural rearrangements.

This tool is capable of detecting both known and novel fusion gene events and provides precise base-level breakpoint resolution, which is essential for understanding the underlying genomic architecture of complex rearrangements. *DenovoFusion* integrates multiple filtering strategies to reduce false positives, including homologous sequence filtering, edge misalignment filtering, and support-based thresholds, ensuring a high level of accuracy and reliability.

Designed with scalability in mind, *DenovoFusion* can process whole-genome datasets efficiently and is compatible with alignment results in PSL, SAM, and PAF formats. It is particularly suitable for large-scale studies, and can serve as a valuable component of comprehensive genomic analysis pipelines.

Even with long read input files the basic concept remains valid.

## Key Features

- **Broad input format support:** &nbsp; Accepts alignment results in PSL, SAM, and PAF formats, enabling seamless integration with a wide range of alignment tools such as BLAT, minimap2, and Bowtie2.

- **Accurate breakpoint clustering and detection:** &nbsp; Implements robust algorithms to cluster split and discordant alignments, ensuring precise identification of base-level breakpoints even in highly repetitive or complex genomic regions.

- **Comprehensive filtering strategies:** &nbsp; Incorporates multiple layers of filtering to minimize false positives, including homologous sequence filtering, edge misalignment detection, minimum support thresholds, and internal tandem duplication filtering.

- **High scalability and performance:** &nbsp; Optimized for large-scale datasets, *DenovoFusion* can process whole-genome assemblies efficiently with multi-threading support, making it suitable for population-level studies or high-throughput analyses.

- **Detection of both known and novel fusion genes:** &nbsp; Capable of recovering well-characterized fusion events from literature as well as discovering novel gene fusions that may not be present in existing databases.

- **Modular and extensible design:** &nbsp; The architecture allows easy integration into existing analysis pipelines and supports future extensions, such as additional filtering modules or visualization tools.

## Installation procedures

### Prerequisites

**Description**

The prerequisites are fully optional and could be adjusted to the users needs.

Nevertheless, for creating a 'psl' alignment file as well as to assemble
the genomic 'contig' sequences for the *DenovoFusion* input, something like an
alignment program and an assembly program are neccessary prerequisites.

In the present workflow 'pblat' and 'MegaHit' are used.
All work below is taking place in a terminal environment.

**Dependencies**

Linux platform (tested on Debian 12), CMake >= 3.10, g++ >= 7.0, git.

**Installation & Usage**

Installation of [pblat](https://github.com/icebert/pblat)  (version: 2.5.1):
```
git clone https://github.com/icebert/pblat
cd pblat
make
```
Installation of [MegaHit](https://github.com/voutcn/megahit) (version: 1.2.9):
```
git clone https://github.com/voutcn/megahit
cd megahit
git submodule update --init
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release  # add -DCMAKE_INSTALL_PREFIX=MY_PREFIX if needed
make -j4
make simple_test  # will test MEGAHIT with a toy dataset
make install  # if needed
```

Example command lines for 'pblat' and 'MegaHit' (see publication for details):
```
pblat -minIdentity=98 -threads=36 genome.fa assembly.fa -ooc=GRCh38.11.ooc output.psl
megahit -1 reads_R1.fastq -2 reads_R2.fastq --k-list 39,59,79,99 -o assembly_output
```

### DenovoFusion

**Description**

Core workflow of the Fusion detection algorithm *DenovoFusion*.

**Dependencies**

Linux platform (tested on Debian 12), CMake >= 3.10, g++ >= 7.0, 
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version: 2.5.4), git.


**Clone repository**
```
git clone https://github.com/korpleul/DenovoFusion.git
```

**Build and test**
```
cd DenovoFusion
mkdir build
cd build
cmake ~/DenovoFusion
make
cd ..
sh ./test.sh
```

If *DenovoFusion* could be started successfully, a lot of lines are shown in the terminal explaining the ongoing work steps
and the successfull termination of the program.
```
$ sh ./test.sh
[2025-08-12 13:04:28] Program DenovoFusion start
[2025-08-12 13:04:28] Launching fusion gene analysing program version 1.0.0
[2025-08-12 13:04:28] Stage1: Determine representive alignment for each contig 
[2025-08-12 13:04:28] Loading alignments from PSL file: 'test/test.psl'
... ...
[2025-08-12 13:04:32] Coverage per base of each contig is calculated and saved
[2025-08-12 13:04:32] Saved per-base coverage to  'test/test.per_base_coverage.tsv' 
[2025-08-12 13:04:32] Write fusion list into file: 'test/test.fusion_list.tsv' 
[2025-08-12 13:04:32] Write discarded fusion list into file: 'test/test.discarded_fusions.tsv' 
[2025-08-12 13:04:32] Stage5: Calculation of the coverage and  prediction of breakpoints 
[2025-08-12 13:04:32] Done (elapsed time=00:00:04, CPU time=00:00:03, peak memory=0.0327gb)
```
If no error occured, the program is ready for production use.

## Input files for DenovoFusion
- **DNA-seq paired-end reads:** &nbsp; Supported in standard **FASTQ** or compressed **FASTQ.GZ** format. 

- **Alignment files:** &nbsp; Alignment results generated by tools such as BLAT, minimap2, or Bowtie2, provided in **PSL**, **SAM**, or **PAF** format. 

- **Contigs from a de novo assembly:** &nbsp; Assembled contig sequences in a **FASTA** file format produced by short-read or hybrid assembly pipelines. 

- **Gene annotation file:** &nbsp; A reference gene annotation file in **GTF** format, required for accurate gene mapping and breakpoint annotation. 

## Description of the program parameters

### Help

```
       -h, --help
          Display this help message and exit.
```

### Mandatory parameters

```
       -1, --fastq1
          Specifies the FASTQ file for the forward read (e.g.01.fastq.gz).
          Multiple input should be divided by comma, as in Bowtie2 input format.

       -2, --fastq2
          Specifies the FASTQ file for the reverse read (e.g.02.fastq.gz).
          Multiple input should be divided by comma, as in Bowtie2 input format.

       -a, --input_assembly
          The contig file from de novo assembly in FASTA format. When using PSL
          mode, please use the processed (cut) contigs instead of the original
          contigs.

       -g, --gtf_path
          GTF annotation file with gene annotation.

       -i, --input_file
          Alignment files in SAM, PAF, or PSL format, generated by alignment
          tools, This serves as the primary input for DenovoFusion's analysis.

       -m, --method
          DenovoFusion applies the same alignment method previously used for
          contig alignment. The supported methods include BLAT, Minimap2sam (SAM
          output), and Minimap2paf (PAF output). If the input file is in BAM
          format, it is recommended to convert it to SAM format before
          proceeding with the analysis.

       -o, --output
          The output file path, which includes separate files for fusion genes,
          coverage information, and log details.

       -p, --prefix
          Prefix for all temporary and result files; you can use the sample name
          here.

       -r, --read_length
          Read length from the raw FASTQ file.
```

### Optional parameter

```
       -A, --edge-unaligned
          Maximum number of unaligned bases at the head or tail of a contig.
          (Default: 50).

       -c, --max_alignment_count
          Maximum number of alignments considered for each contig. Larger values
          may increase processing time (Default: 100).

       -d, --min_identity_fract
          Minimum alignment identity threshold (Default: 0.95).

       -e, --min_score_each
          Minimum score for each alignment (Default: 20).

       -E, --min-edge-length
          Minimum edge length of split reads in re-alignment; adjust according
          toread length (Default: 20).

       -f, --max-gap-size
          Maximum number of gap bases allowed in paired alignments (Default: 2).

       -G, --short-segment-threshold
          Minimum segment length requirement for fusion parts (Default: 35).

       -I, --inclusion-fraction-weight
          Weight of the inclusion fraction in alignment score calculation
          (Default: 1).

       -l, --max-overlap-size
          Maximum number of overlapping bases allowed in paired alignments
          (Default: 8).

       -n, --max-pair-combination
          Maximum number of alignment combinations per contig. Larger values may
          increase processing time (Default: 3).

       -N, --min-span-reads
          Minimum number of spanning reads required as support (Default: 1).

       -O, --overlap-fraction-weight
          Weight of the overlap fraction in alignment score calculation
          (Default: 1).

       -P, --min-split-reads
          Minimum number of split reads required as support (Default: 3).

       -q, --threads
          Number of threads to use (Default: 4).

       -s, --min_score_total
          Minimum total score for combined alignments (Default: 95).

       -S, --size-weight
          Weight of alignment count in alignment score calculation (Default: 1).

       -T, --long-gap-threshold
          Gap threshold for candidate fusions on the same chromosome (Default:
          200000).

       -v, --coverage-differ
          Coverage ratio threshold between adjacent bases in breakpoint
          prediction (Default: 0.85).

       -z, --size-ratio-threshold
          Proportion of fusion part 1 and part 2 must remain within an
          acceptable range (Default: 0.1).
```

## Description of the output files

The standard output folder includes the following files:
```
*.chosen.fasta
   FASTA file containing contigs with putative fusion alignments.
   These contigs serve as the reference for subsequent realignment steps.

*_idx
   Folder with index files generated by Bowtie2 for the realignment process.

*.sam
   SAM file produced after realigning raw reads to chosen.fasta,
   containing all aligned reads.

*.discarded_fusions.tsv
   Tab-delimited list of fusion candidates that were filtered out
   during post-processing based on predefined criteria.

*.fusion_list.csv
   Final list of high-confidence fusion events with annotated breakpoints
   and supporting evidence.

*.per_base_coverage.tsv
   Coverage information at single-base resolution for the overlap regions
   of confirmed fusion genes.

*.log
   Execution log file summarizing the workflow steps, parameter settings,
   and runtime messages of DenovoFusion.
```

## Usage

A minimal test command out of the file *~/DenovoFusion/test.sh* is shown below.
```
DenovoFusion -m blat -q 12 -i input.psl -a assembly.fa -o path/to/result -p prefix -g Homo_sapiens.GRCh38.105.gtf -1 01.fastq.gz -2 02.fastq.gz 
```

---

## License

This project is licensed under the GPL License – see the LICENSE file for details.

## Contact

Developer: Xinwei Zhao

Email: zhaoxi@uni-muenster.de or korschi@uni-muenster.de

---
