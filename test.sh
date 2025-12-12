#!/bin/bash
# Enter DenovoFusion folder
cd "$(dirname "$0")"

# run test command
./DenovoFusion -m blat -q 4 -r 100 -i test/test.psl -a test/test.fa -o test -p test -g test/test.gtf -1 test/test.01.fastq.gz -2 test/test.02.fastq.gz

