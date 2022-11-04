#include <bits/stdc++.h>
#include "../kmerind/src/common/kmer.hpp"
#include "../kmerind/src/common/alphabets.hpp"
#include "../kmerind/src/index/kmer_index.hpp"
#include "mxx/env.hpp"
#include <mpi.h>

template <typename key>
using bal_mapParm = bliss::index::kmer::BimoleculeHashMapParams<key>;

template <typename tupple_type>
using bal_fastqPaser = bliss::io::FASTQParser<tupple_type>;
//using bal_kmerPaser = bliss::index::kmer::KmerCountTupleParser<tupple_type>;

template<typename bal_Iterator, template<typename> class bal_Parser>
using bal_seqIter = bliss::io::SequencesIterator<bal_Iterator, bal_Parser>;

int main(int argc, char** argv)
{
    //OpenMPI intialization
    mxx::env e(argc, argv);
    mxx::comm comm;
    if (comm.rank() == 0) printf("EXECUTING %s\n", argv[0]);
    //comm.barrier();

    
    //Input data
    std::string input_file = "tests/test_small.fastq";


    //Kmer data type intialization
    typedef bliss::common::alphabet::DNA_T<> alph;
    typedef bliss::common::Kmer<4, alph> bal_kmer;


    //KmerIndex creation
    using bal_MapType = dsc::counting_unordered_map<bal_kmer, unsigned int, bal_mapParm>;
    typedef bliss::index::kmer::CountIndex<bal_MapType> kmer_index;
    
    kmer_index first(comm);
    first.build_mmap<bal_fastqPaser, bal_seqIter>(input_file, comm);
    printf("working\n");
    
    return 0;
}