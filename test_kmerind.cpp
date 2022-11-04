#include "../kmerind/src/common/kmer.hpp"
#include "../kmerind/src/common/alphabets.hpp"
#include "../kmerind/src/index/kmer_index.hpp"
#include "mxx/env.hpp"
#include <bits/stdc++.h>
#include <mpi.h>


//Type definations
template <typename key>
//using map_parm = bliss::index::kmer::BimoleculeHashMapParams<key>;
using map_parm = bliss::index::kmer::SingleStrandHashMapParams<key>;

template <typename tupple_type>
using fastq_paser = bliss::io::FASTQParser<tupple_type>;

template<typename Iterator, template<typename> class Parser>
using seq_iter = bliss::io::SequencesIterator<Iterator, Parser>;

typedef bliss::common::alphabet::DNA_T<> alph;

typedef bliss::common::Kmer<2, alph> Kmer_4;

using MapType = dsc::counting_unordered_map<Kmer_4, unsigned int, map_parm>;
typedef bliss::index::kmer::CountIndex<MapType> kmer_index;

typedef std::pair<std::vector<std::string>, std::vector<int>> kmerAndCounts;





//Function declaration
kmerAndCounts read_kmer_counts(std::string file_name);

std::vector<Kmer_4> get_kmerind_kmers(std::vector<std::string> file_string);




int main(int argc, char** argv)
{
    //OpenMPI intialization
    mxx::env e(argc, argv);
    mxx::comm comm;
    if (comm.rank() == 0) printf("EXECUTING %s\n", argv[0]);
    //comm.barrier();

    
    //Input data
    std::string input_file = "tests/test_small.fastq";


    //KmerIndex creation
    kmer_index first(comm);
    first.build_mmap<fastq_paser, seq_iter>(input_file, comm);


    std::string file_name = input_file;
    kmerAndCounts actual_counts = read_kmer_counts(file_name);
    std::vector<Kmer_4> kmerind_kmers = get_kmerind_kmers(actual_counts.first);

    auto counts = first.find(kmerind_kmers);

    std::vector<int> file_counts = actual_counts.second;
    std::vector<Kmer_4> file_kmers = get_kmerind_kmers(actual_counts.first);
    for(int i = 0; i < counts.size(); i++)
    {
        for(int j = 0; j < file_counts.size(); j++)
        {
            if(counts[i].first == file_kmers[j])
            {
                assert(file_counts[j] == counts[i].second);
            }
        }
    }
    

    std::cout << "workdin" << std::endl;
    
    return 0;
}


kmerAndCounts read_kmer_counts(std::string file_name)
{
    std::string in_file = file_name + ".counts";
    std::ifstream file_stream(in_file);

    
    std::vector <std::string> file_kmers;  
    std::vector <int> file_counts;
    std::string line;
    while (getline(file_stream, line))
    {
        //Where to split
        int pos = line.find("@");

        //extracting kmer and number
        std::string kmer = line.substr(0, pos); // store the substring   
        line.erase(0, pos + 1);  

        //adding kmer and number to function
        file_kmers.push_back(kmer);
        file_counts.push_back(stoi(line));
    }
    
    return std::make_pair(file_kmers, file_counts);
}



std::vector<Kmer_4> get_kmerind_kmers(std::vector<std::string> file_string)
{
    std::vector<Kmer_4> kmerind_kmers;
    for(int i = 0; i < file_string.size(); i++)
    {
        kmerind_kmers.push_back(Kmer_4(file_string[i]));
    }
    return kmerind_kmers;
}