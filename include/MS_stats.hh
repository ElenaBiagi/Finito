#pragma once

#include <string>
#include <cstring>
#include "sbwt/cxxopts.hpp"
#include "sbwt/globals.hh"
#include "sbwt/SBWT.hh"
#include "sbwt/SubsetWT.hh"
#include "sbwt/stdlib_printing.hh"
#include "sbwt/SeqIO.hh"
#include "sbwt/SubsetMatrixRank.hh"
#include "sbwt/buffered_streams.hh"
#include "sbwt/variants.hh"
#include "sbwt/commands.hh"
#include <filesystem>
#include <cstdio>
#include <optional>

#include "sbwt/throwing_streams.hh"
#include "PackedStrings.hh"
#include "SeqIO.hh"
#include "common.hh"

//#include <sdsl/elias_fano_vector.hpp>

// Defined in build_fmin.cpp. TODO: TO HEADER
extern pair<int64_t,int64_t> drop_first_char_stats(const int64_t  new_len, const pair<int64_t,int64_t>& I, const sdsl::int_vector<>& LCS, const int64_t n_nodes);

using namespace std;

using namespace sbwt;



template<typename reader_t, typename out_stream_t>
int64_t run_MS_streaming(reader_t& reader, out_stream_t& out, const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, vector<int64_t>& left, vector<int64_t>& right){    
    const int64_t k = sbwt.get_k();
    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    int64_t kmers_count = 0;
    int64_t total_positive = 0;
    vector<int64_t> out_buffer, out_buffer_rev;

    int64_t query_seq=0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        int64_t t0 = cur_time_micros();
        vector<tuple<int64_t, int64_t, int64_t>> MS_answer = MatchingStatistics_stats(sbwt, LCS, reader.read_buf, left, right);
        //FinimizerIndex::QueryMS result = index.MS_statistics(reader.read_buf,left, right);

        for(int64_t i = 0; i < MS_answer.size(); i++){
            int64_t len, freq, Istart;
            std::tie(len,freq,Istart) = MS_answer[i];
            if(i > 0) out << ' ';
            out << len << ", ";
            if (len == k) kmers_count++;
            //out << '(' << unitig << ',' << pos << ')';
        }
        out << '\n';

        int64_t tot_kmers = MS_answer.size();
        number_of_queries += tot_kmers;//result.local_offsets.size();
     
        total_micros += cur_time_micros() - t0;
    }
    write_log("k " + to_string(k), LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    write_log("Found kmers: " + to_string(kmers_count), LogLevel::MAJOR);
    return number_of_queries;
}

template<typename reader_t, typename out_stream_t>
int64_t run_file(const string& infile, out_stream_t& out, const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, vector<int64_t>& left, vector<int64_t>& right){
    reader_t reader(infile);
    write_log("Running streaming queries from input file " + infile, LogLevel::MAJOR);
    return run_MS_streaming(reader, out,sbwt, LCS, left, right);
}

// Returns number of queries executed
int64_t run_MS(const vector<string>& infiles, const optional<vector<string>>& outfiles, const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, vector<int64_t>& left, vector<int64_t>& right){

    if(outfiles.has_value()){
        if(infiles.size() != outfiles.value().size()){
            string count1 = to_string(infiles.size());
            string count2 = to_string(outfiles.value().size());
            throw std::runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
        }
    }

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_no_gzip;

    int64_t n_queries_run = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = SeqIO::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input){
            if(outfiles.has_value()){
                ofstream out(outfiles.value()[i]);
                n_queries_run += run_file<in_gzip>(infiles[i], out, sbwt, LCS, left, right);
            } else { // To stdout
                n_queries_run += run_file<in_gzip>(infiles[i], cout, sbwt, LCS, left, right);
            }
        }
        else {
            if(outfiles.has_value()){
                ofstream out(outfiles.value()[i]);
                n_queries_run += run_file<in_no_gzip>(infiles[i], out, sbwt, LCS, left, right);
            } else{ // To stdout
                n_queries_run += run_file<in_no_gzip>(infiles[i], cout, sbwt, LCS, left, right);
            }
        }
    }
    return n_queries_run;
}

int MS_stats(int argc, char** argv){

    int64_t micros_start = cur_time_micros();

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Matching Statistics + LCS scan distribution.");

    options.add_options()
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("i,index-file", "Index filename prefix.", cxxopts::value<string>())
        ("q,query-file", "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())
        ("lcs", "Provide in input the LCS file if available.", cxxopts::value<string>()->default_value(""))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    // Interpret input file
    string queryfile = opts["query-file"].as<string>();
    vector<string> query_files;
    bool multi_file = queryfile.size() >= 4 && queryfile.substr(queryfile.size() - 4) == ".txt";
    if(multi_file){
        query_files = readlines(queryfile);
    } else{
        query_files = {queryfile};
    }
    for(string file : query_files) check_readable(file);

    // Interpret output file
    optional<vector<string>> output_files;
    string outfile = opts["out-file"].as<string>();
    try{
        string outfile = opts["out-file"].as<string>();
        if(multi_file){
            output_files = readlines(outfile);
        } else{
            output_files = {outfile};
        }
        for(string file : output_files.value()) check_writable(file);
    } catch(cxxopts::option_has_no_value_exception& e){
        write_log("No output file given, writing to stdout", LogLevel::MAJOR);
    }

    // sbwt
    string indexfile = opts["index-file"].as<string>();
    check_readable(indexfile);
    throwing_ifstream in(indexfile, ios::binary);
    vector<string> variants = get_available_variants_fmin();
    string variant = load_string(in.stream); // read variant type
    if (std::find(variants.begin(), variants.end(), variant) == variants.end()) {
       cerr << "Error loading index from file: unrecognized variant specified in the file" << endl;
        return 1;
    }

    write_log("Loading the index variant " + variant, LogLevel::MAJOR);
    if (variant == "plain-matrix") {

        plain_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);

        string LCS_file = opts["lcs"].as<string>();
        if (LCS_file.empty()) {
            std::cerr<< "LCS_file empty" << std::endl;
            LCS_file = indexfile + "LCS.sdsl";
            const sdsl::int_vector<> LCS = lcs_basic_parallel_algorithm(sbwt, 8);
            //const sdsl::int_vector<> LCS = lcs_basic_algorithm(sbwt);
            save_v(LCS_file, LCS);
        }
        sdsl::int_vector<> LCS;
        load_v(LCS_file, LCS);
        std::cerr<< "LCS_file loaded" << std::endl;

    vector<int64_t> left, right;
    int64_t number_of_queries = run_MS(query_files, output_files, sbwt, LCS, left, right);
    int64_t new_total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)new_total_micros / number_of_queries), LogLevel::MAJOR);
    write_log("total number of queries: " + to_string(number_of_queries), LogLevel::MAJOR);
    
    write_csv(outfile + "LCS_Lscan.csv", "left", left);
    write_csv(outfile + "LCS_Rscan.csv", "right", right);
    }
    return 0;

}