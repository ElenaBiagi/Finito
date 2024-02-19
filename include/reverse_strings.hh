//
// Created by Biagi, Elena on 19.12.2023.
//
#pragma once

#include <fstream>
#include <string>
#include <tuple>
#include <set>
#include <string.h>
#include <climits>
#include "sbwt/cxxopts.hpp"
#include "sbwt/globals.hh"
#include "sbwt/SBWT.hh"
#include "sbwt/stdlib_printing.hh"
#include "sbwt/SeqIO.hh"
#include "sbwt/buffered_streams.hh"
#include "PackedStrings.hh"
#include "FinimizerIndex.hh"
#include <vector>
#include <utility>
#include <algorithm>
#include "variants.hh"
#include "sdsl/vectors.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/bit_vectors.hpp"
#include "sbwt/commands.hh"
#include "sbwt/Kmer.hh"
#include "sbwt/suffix_group_optimization.hh"
#include "lcs_basic_parallel_algorithm.hpp"
#include "common.hh"
#include "sbwt/SeqIO.hh"

using namespace std;
using namespace sbwt;

void writeFastaFile(const std::vector<std::string>& sequences, const std::string& outputFilePath) {
    // Open the output FASTA file
    std::ofstream outputFile(outputFilePath);
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file: " << outputFilePath << std::endl;
        return;
    }

    // Write each sequence to the FASTA file
    for (const std::string& sequence : sequences) {
        // Assuming here that each sequence has a unique identifier.
        // You may need to adjust this part based on your specific input.
        std::string header = ">Sequence_" + std::to_string(&sequence - &sequences[0]);
        
        // Write header and sequence to the file
        outputFile << header << "\n";
        outputFile << sequence << "\n";
    }

    // Close the file
    outputFile.close();

    std::cout << "FASTA file created successfully." << std::endl;
}

//write the output as a fasta file
template<typename writer_t>
void write_string_fasta(const vector<string>& strings, string& out_prefix){ // writer_t& out
    writer_t out(out_prefix);

    char newline = '\n';
    for (string& S: strings){
        out.write(">", 1);
        out.write(S, S.length());        
        out.write(&newline, 1);

    }
}


template<typename reader_t>
vector<string> rev_file(const string& infile){
    reader_t reader(infile);
    //writer_t writer(out_prefix);
    vector<string> unitigs;

    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) [[unlikely]] break;
        unitigs.push_back(sbwt::get_rc(reader.read_buf));
    }
    return unitigs;
}

void rev_unitigs(const vector<string>& infiles, const string& out_file){
    vector<string> reversed_unitigs;
    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_no_gzip;

    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = SeqIO::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input){
            vector<string> rev_v = rev_file<in_gzip>(infiles[i]);
            reversed_unitigs.insert(reversed_unitigs.end(), rev_v.begin(), rev_v.end()); 
        }
        else {
            vector<string> rev_v = rev_file<in_gzip>(infiles[i]);
            reversed_unitigs.insert(reversed_unitigs.end(), rev_v.begin(), rev_v.end()); 
        }
    }
    writeFastaFile(reversed_unitigs, out_file);

    return;

}

int reverse_strings(int argc, char** argv) {

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Find all Finimizers of all input reads.");

    options.add_options()

        ("i,in-file",
            "The SPSS in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.",
            cxxopts::value<string>())
        ("o,out-file", 
            "Reverse complement files.", cxxopts::value<string>()->default_value("out.fna"))

        ("h,help", "Print usage");

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")) {
        std::cerr << options.help() << std::endl;
        exit(1);
    }
    

    // input files
    string in_file = opts["in-file"].as<string>();
    vector<string> input_files;
    bool multi_file = in_file.size() >= 4 && in_file.substr(in_file.size() - 4) == ".txt";
    if(multi_file){
        input_files = readlines(in_file);
    } else{
        input_files = {in_file};
    }
    for(string file : input_files) check_readable(file);

    string out_file = opts["out-file"].as<string>();
    
    rev_unitigs(input_files, out_file);

    return 0;
}