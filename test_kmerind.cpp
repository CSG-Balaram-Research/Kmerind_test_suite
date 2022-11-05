#include "../kmerind/src/common/kmer.hpp"
#include "../kmerind/src/common/alphabets.hpp"
#include "../kmerind/src/index/kmer_index.hpp"
#include "mxx/env.hpp"
#include <bits/stdc++.h>
#include <mpi.h>


//Type definations
template <typename key>
using map_parm = bliss::index::kmer::SingleStrandHashMapParams<key>;

template <typename tupple_type>
using fastq_paser = bliss::io::FASTQParser<tupple_type>;

template<typename Iterator, template<typename> class Parser>
using seq_iter = bliss::io::SequencesIterator<Iterator, Parser>;

typedef bliss::common::alphabet::DNA_T<> alph;

typedef bliss::common::Kmer<2, alph> Kmer_4;

using MapType = dsc::counting_unordered_map<Kmer_4, long unsigned int, map_parm>;
typedef bliss::index::kmer::CountIndex<MapType> kmer_index;

typedef std::pair<std::vector<std::string>, std::vector<int>> kmerAndCounts;

typedef decltype(::std::declval<MapType>().count(::std::declval<std::vector<Kmer_4> &>())) KmerindCountsType;
KmerindCountsType just;
typedef decltype(just[0]) KmerCountPair;

typedef std::pair<Kmer_4, unsigned long int> MyKmerCount;





//Function declaration
kmerAndCounts read_kmer_counts(std::string file_name);

std::vector<Kmer_4> get_kmerind_kmers(std::vector<std::string> file_string);

bool equal_my(KmerindCountsType counts, kmerAndCounts file);

bool less_than_sort(KmerCountPair i, KmerCountPair j);

bool less_than_lower(KmerCountPair i, MyKmerCount j);

bool test_for_fastq(std::string file_name, mxx::comm& comm);





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
    kmerAndCounts actual_counts = read_kmer_counts(file_name);
    std::vector<Kmer_4> kmerind_kmers = get_kmerind_kmers(actual_counts.first);

    //KmerIndex creation
    kmer_index first(comm);
    first.build_mmap<fastq_paser, seq_iter>(file_name, comm);
    auto counts = first.find(kmerind_kmers);
    
    bool passed = equal_my(counts, actual_counts);
    return passed;
}


bool equal_my(KmerindCountsType counts, kmerAndCounts file_data)
{
    std::sort(counts.begin(), counts.end(), less_than_sort);
    std::vector<int> file_counts = file_data.second;
    std::vector<Kmer_4> file_kmers = get_kmerind_kmers(file_data.first);

    bool equal = true;    
    bool i_correct = true;
    for(int i = 0; i < file_counts.size(); i++)
    {
        MyKmerCount val(file_kmers[i], file_counts[i]);
        auto loc = *std::lower_bound(counts.begin(), counts.end(), val, less_than_lower);
        
        if(loc.first == file_kmers[i] && loc.second == file_counts[i])
            i_correct = true;
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


bool less_than_sort(KmerCountPair i, KmerCountPair j)
{
    return (i.first < j.first);
}


bool less_than_lower(KmerCountPair i, MyKmerCount j)
{
    return i.first < j.first;
}