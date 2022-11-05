#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "index/kmer_index.hpp"
#include "mxx/env.hpp"
#include <bits/stdc++.h>
#include <mpi.h>


const int K_value = 2;

//Type definations
template <typename key>
using map_parm = bliss::index::kmer::SingleStrandHashMapParams<key>;

template <typename tupple_type>
using fastq_paser = bliss::io::FASTQParser<tupple_type>;

template<typename Iterator, template<typename> class Parser>
using seq_iter = bliss::io::SequencesIterator<Iterator, Parser>;

typedef bliss::common::alphabet::DNA_T<> alph;

typedef bliss::common::Kmer<K_value, alph> Kmer_k;

using MapType = dsc::counting_unordered_map<Kmer_k, long unsigned int, map_parm>;
typedef bliss::index::kmer::CountIndex<MapType> kmer_index;

typedef std::pair<std::vector<std::string>, std::vector<int>> MykmerAndCountsVec;

typedef std::pair<Kmer_k, unsigned long int> MyKmerCountPair;

typedef decltype(::std::declval<MapType>().count(::std::declval<std::vector<Kmer_k> &>())) kmerAndCountsVec;

kmerAndCountsVec just;
typedef decltype(just[0]) KmerCountPair;






//Function declaration
bool test_for_fastq(std::string file_name, mxx::comm& comm);

bool equal_my(kmerAndCountsVec counts, MykmerAndCountsVec file);

MykmerAndCountsVec read_kmer_counts(std::string file_name);

std::vector<Kmer_k> get_kmerind_kmers(std::vector<std::string> file_string);

bool less_than_sort(KmerCountPair i, KmerCountPair j);

bool less_than_lower(KmerCountPair i, MyKmerCountPair j);







int main(int argc, char** argv)
{
    //OpenMPI intialization
    mxx::env e(argc, argv);
    mxx::comm comm;
    if (comm.rank() == 0) printf("EXECUTING %s\n", argv[0]);
    //comm.barrier();

    
    test_for_fastq("tests/test_small.fastq", comm);
    std::cout << "Ran all test cases successfully" << std::endl;
    
    return 0;
}







//Function definations
bool test_for_fastq(std::string file_name, mxx::comm& comm)
{
    MykmerAndCountsVec actual_counts = read_kmer_counts(file_name);
    std::vector<Kmer_k> kmerind_kmers = get_kmerind_kmers(actual_counts.first);

    //KmerIndex creation
    kmer_index first(comm);
    first.build_mmap<fastq_paser, seq_iter>(file_name, comm);
    auto counts = first.find(kmerind_kmers);
    
    bool passed = equal_my(counts, actual_counts);
    return passed;
}


bool equal_my(kmerAndCountsVec counts, MykmerAndCountsVec file_data)
{
    std::sort(counts.begin(), counts.end(), less_than_sort);
    std::vector<int> file_counts = file_data.second;
    std::vector<Kmer_k> file_kmers = get_kmerind_kmers(file_data.first);

    bool equal = true;    
    bool i_correct = true;
    for(int i = 0; i < file_counts.size(); i++)
    {
        MyKmerCountPair val(file_kmers[i], file_counts[i]);
        auto loc = std::lower_bound(counts.begin(), counts.end(), val, less_than_lower);
        if(loc != counts.end())
        {
            KmerCountPair match = *loc;
            if(match.first == file_kmers[i] && match.second == file_counts[i])
                i_correct = true;
            else
                i_correct = false;
        }
        else
            i_correct = false;


        if(i_correct == false)
        {
            equal = false;
            break;
        }
    } 

    return equal;
}


MykmerAndCountsVec read_kmer_counts(std::string file_name)
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



std::vector<Kmer_k> get_kmerind_kmers(std::vector<std::string> file_string)
{
    std::vector<Kmer_k> kmerind_kmers;
    for(int i = 0; i < file_string.size(); i++)
    {
        kmerind_kmers.push_back(Kmer_k(file_string[i]));
    }
    return kmerind_kmers;
}


bool less_than_sort(KmerCountPair i, KmerCountPair j)
{
    return (i.first < j.first);
}


bool less_than_lower(KmerCountPair i, MyKmerCountPair j)
{
    return i.first < j.first;
}