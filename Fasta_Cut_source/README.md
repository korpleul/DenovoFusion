# Fasta_Cut

A high-performance, lightweight C++ utility designed to fragment long FASTA contigs into fixed-length segments.

## Introduction
**Fasta_Cut** was developed to address performance bottlenecks in sequence alignment workflows (specifically with tools like BLAT) where processing extremely long contigs can significantly degrade speed. By segmenting long sequences into smaller, manageable chunks, this tool enables parallel processing and faster alignment without data loss. It is particularly designed as a pre-processing step for fusion detection pipelines DenovoFusion, ensuring that downstream alignment steps run efficiently.

## Features

- **Simplicity:** Minimal dependencies and easy to integrate into pipelines.
- **Preservation:** Handles FASTA headers intelligently to ensure segments can be traced back to the original contig.
- **Efficiency:** optimized for large-scale genomic data.

### Requirements
- Linux (Tested on Debian 12 / Ubuntu 20.04+)
- g++ (v7.0 or higher recommended)
- CMake (v3.0+) or Make

### Build
```bash
cd Fasta_Cut-bin
g++ -std=c++11 -O2 -o Fasta_Cut main.cpp
```
# Check if it works
./Fasta_Cut


## Usage
```
Usage: Fasta_Cut <input_file> <output_file> <sequence_length>

*sequence_length
    The length of segment you want to have for cut contigs.
```

---

## License

This project is licensed under the GPL License – see the LICENSE file for details.

## Contact

Developer: Xinwei Zhao

Email: zhaoxi@uni-muenster.de or korschi@uni-muenster.de

---
