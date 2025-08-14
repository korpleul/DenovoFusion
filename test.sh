#!/bin/bash
# 进入 DenovoFusion 目录，确保相对路径正确
cd "$(dirname "$0")"

# 运行程序，参数中用 test/ 目录路径
./DenovoFusion -m blat -q 4 -r 100 -i test/test.psl -a test/test.fa -o test -p test -g test/test.gtf -1 test/test.01.fastq.gz -2 test/test.02.fastq.gz

