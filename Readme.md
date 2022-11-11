- *tests* directory contains fastq files. (One can add more fastq files)
- *kmer_couter.py* python program to solve kmercounting.
    - comandline parameters 
      - k_size
      - path to fastq file
    - dependency *pyfastx*
    - This programs generates _fastq.counts file in tests directory which contains kmer and its counts.
    - This *counts will be used by test_kmerind.cpp to check corretness of kmerind codebase.
- *test_kmerind.cpp* cpp program to check if the kmercounting by kmerind is correct.
    - commandline
        - mpic++ -I../kmerind/src/ -I../kmerind/build/ -I../kmerind/ext/mxx/include/ -I../kmerind/ext -I../kmerind/ext/Nadeau -I../kmerind/ext/sparsehash-c11  --std=c++11  test_kmerind.cpp
        - mpirun -np 1 a.out
    - dpendency mpicc++, kmerind codebase, kmerind build
    - parameters.. k value (global variable), file_name(function papramete) as more than file can be tested at time


***Workalbe a better interface for the test_kmerind.cpp inputs***


