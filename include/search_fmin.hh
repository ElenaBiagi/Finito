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
extern pair<int64_t,int64_t> drop_first_char(const int64_t  new_len, const pair<int64_t,int64_t>& I, const sdsl::int_vector<>& LCS, const int64_t n_nodes);

using namespace std;

using namespace sbwt;

template<typename reader_t, typename out_stream_t>
int64_t run_fmin_queries_streaming(reader_t& reader, out_stream_t& out, const FinimizerIndex& index, const string& stats_filename){
    const int64_t k = index.sbwt->get_k();
    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    int64_t kmers_count = 0 , kmers_count_rev = 0;
    int64_t total_positive = 0;
    vector<int64_t> out_buffer, out_buffer_rev;

    int64_t query_seq=0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        int64_t t0 = cur_time_micros();
        //pair<vector<int64_t>, int64_t> final_pair = rarest_fmin_streaming_search(DNA_bitvectors, DNA_rs, sbwt, LCS, reader.read_buf, t, fmin_rs, global_offsets, ef_endpoints, Ustart_rs, found_kmers);
        FinimizerIndex::QueryResult result = index.search(reader.read_buf);

        //reverse compl
        const string reverse = sbwt::get_rc(reader.read_buf);
        FinimizerIndex::QueryResult r_result = index.search(reverse);
        int64_t tot_kmers = result.local_offsets.size();
        int64_t str_len = reverse.length(); // the string and its reverse complement have the same length
        for(int64_t i = 0; i < tot_kmers; i++){
            int64_t unitig, pos;
            if (result.local_offsets[i].first==-1) {
                std::tie(unitig,pos) = r_result.local_offsets[str_len-k-i];
            } else{
                std::tie(unitig,pos) = result.local_offsets[i];
            }
            if(unitig != -1) total_positive++;
            if(i > 0) out << ' ';
            out << '(' << unitig << ',' << pos << ')';
        }
        out << '\n';

        kmers_count += result.n_found;
        kmers_count_rev += r_result.n_found;
        number_of_queries += tot_kmers;//result.local_offsets.size();
     
        total_micros += cur_time_micros() - t0;
    }
    write_log("k " + to_string(k), LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    write_log("Found kmers: " + to_string(kmers_count), LogLevel::MAJOR);
    write_log("Found kmers reverse : " + to_string(kmers_count_rev), LogLevel::MAJOR);
    write_log("Total found kmers: " + to_string(total_positive), LogLevel::MAJOR);

    std::ofstream statsfile;
    statsfile.open(stats_filename, std::ios_base::app); // append instead of overwrite
    statsfile << to_string(k) + "," + to_string(kmers_count+kmers_count_rev) + "," + to_string(number_of_queries);
    statsfile.close();
    return number_of_queries;
}

template<typename reader_t, typename out_stream_t>
int64_t run_fmin_file(const string& infile, out_stream_t& out, const string& stats_filename, const FinimizerIndex& index){
    reader_t reader(infile);
    write_log("Running streaming queries from input file " + infile, LogLevel::MAJOR);
    return run_fmin_queries_streaming(reader, out, index, stats_filename);
}

// Returns number of queries executed
int64_t run_fmin_queries(const vector<string>& infiles, const optional<vector<string>>& outfiles, const string& stats_filename, const FinimizerIndex& index){

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
                n_queries_run += run_fmin_file<in_gzip>(infiles[i], out, stats_filename, index);
            } else { // To stdout
                n_queries_run += run_fmin_file<in_gzip>(infiles[i], cout, stats_filename, index);
            }
        }
        else {
            if(outfiles.has_value()){
                ofstream out(outfiles.value()[i]);
                n_queries_run += run_fmin_file<in_no_gzip>(infiles[i], out, stats_filename, index);
            } else{ // To stdout
                n_queries_run += run_fmin_file<in_no_gzip>(infiles[i], cout, stats_filename, index);
            }
        }
    }
    return n_queries_run;
}

int search_fmin(int argc, char** argv){

    int64_t micros_start = cur_time_micros();

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Query all Finimizers of all input reads.");

    vector<string> types = get_available_types();
    string all_types_string;
    for (string type: types) all_types_string += " " + type;


    options.add_options()
        ("o,out-file", "Output filename, or stdout if not given.", cxxopts::value<string>())
        ("i,index-file", "Index filename prefix.", cxxopts::value<string>())
        ("q,query-file", "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())
        ("type", "Decide which streaming search type you prefer. Available types: " + all_types_string,cxxopts::value<string>()->default_value("rarest"))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    // TODO add type, only rarest now
    string type = opts["type"].as<string>();
    if(std::find(types.begin(), types.end(), type) == types.end()){
        cerr << "Error: unknown type: " << type << endl;
        cerr << "Available types are:" << all_types_string << endl;
        return 1;
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

    string index_prefix = opts["index-file"].as<string>();

    int64_t number_of_queries = 0;

    cerr << "Loading index..." << endl;
    FinimizerIndex index;
    index.load(index_prefix);
    cerr << "Index loaded" << endl;

    const int64_t k = index.sbwt->get_k();
    cerr << "k = "<< to_string(k);
    cerr << " SBWT nodes: "<< to_string(index.sbwt->number_of_subsets())<< " kmers: "<< to_string(index.sbwt->number_of_kmers())<< endl;

    number_of_queries += run_fmin_queries(query_files, output_files, index_prefix + ".stats", index);
    int64_t new_total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)new_total_micros / number_of_queries), LogLevel::MAJOR);
    write_log("total number of queries: " + to_string(number_of_queries), LogLevel::MAJOR);
    
    std::ofstream statsfile2;
    statsfile2.open(index_prefix + "stats.txt", std::ios_base::app); // append instead of overwrite
    string results = to_string(number_of_queries);
    statsfile2 << "," + to_string((double)new_total_micros / number_of_queries);
    
    int64_t bytes = index.size_in_bytes();
    write_log("bytes: " + to_string(bytes), LogLevel::MAJOR);
    statsfile2 <<  "," + to_string(bytes);
    statsfile2 << "," + to_string(static_cast<double>(bytes*8)/index.sbwt->number_of_kmers()) + "\n";
    statsfile2 << "," + to_string(index.sbwt->number_of_kmers()) + "\n";
    statsfile2.close();

    int64_t total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)total_micros / number_of_queries), LogLevel::MAJOR);
    
    return 0;

}