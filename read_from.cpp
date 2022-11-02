#include <bits/stdc++.h>
#include <iostream>
#include <utility>

typedef std::pair<std::vector<std::string>, std::vector<int>> kmerAndCounts;
kmerAndCounts read_kmer_counts(std::string file_name);

int main()
{
    kmerAndCounts actual = read_kmer_counts("tests/test_small.fastq");
    std::vector <std::string> file_kmers = actual.first;  
    std::vector <int> file_counts = actual.second;

    for(int i = 0; i < file_kmers.size(); i++)
        std::cout << file_kmers[i] << "@" << file_counts[i] << std::endl;
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