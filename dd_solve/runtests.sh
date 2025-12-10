#!/bin/bash

#run all variants
solve/gpusolver 32 32 32 8 16 16 8 1 >> out_32_32_32_8_16_16.txt
solve/gpusolver 128 128 128 8 16 16 8 1 >> out_128_128_128_8_16_16.txt
solve/gpusolver 256 128 128 8 16 16 8 1 >> out_256_128_128_8_16_16.txt

#profile each variant
rocprof -i bsr_test.txt --basenames on -o variant1_128_3.csv solve/gpusolver 128 128 128 8 16 16 1
rocprof --stats --basenames on -o variant1_128_3_stats.csv solve/gpusolver 128 128 128 8 16 16 1

rocprof -i bsr_test.txt --basenames on -o variant2_128_3.csv solve/gpusolver 128 128 128 8 16 16 2
rocprof --stats --basenames on -o variant2_128_3_stats.csv solve/gpusolver 128 128 128 8 16 16 2

rocprof -i bsr_test.txt --basenames on -o variant3_128_3.csv solve/gpusolver 128 128 128 8 16 16 3
rocprof --stats --basenames on -o variant3_128_3_stats.csv solve/gpusolver 128 128 128 8 16 16 3

rocprof -i bsr_test.txt --basenames on -o variant4_128_3.csv solve/gpusolver 128 128 128 8 16 16 4
rocprof --stats --basenames on -o variant4_128_3_stats.csv solve/gpusolver 128 128 128 8 16 16 4

rocprof -i bsr_test.txt --basenames on -o variant5_128_3.csv solve/gpusolver 128 128 128 8 16 16 5
rocprof --stats --basenames on -o variant5_128_3_stats.csv solve/gpusolver 128 128 128 8 16 16 5

rocprof -i bsr_test.txt --basenames on -o variant6_128_3.csv solve/gpusolver 128 128 128 8 16 16 6
rocprof --stats --basenames on -o variant6_128_3_stats.csv solve/gpusolver 128 128 128 8 16 16 6

rocprof -i bsr_test.txt --basenames on -o variant7_128_3.csv solve/gpusolver 128 128 128 8 16 16 7
rocprof --stats --basenames on -o variant7_128_3_stats.csv solve/gpusolver 128 128 128 8 16 16 7


