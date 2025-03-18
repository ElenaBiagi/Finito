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
int64_t run_fmin_queries_streaming(reader_t& reader, out_stream_t& out, const FinimizerIndex& index, const string& stats_filename, vector<int64_t>& left, vector<int64_t>& right){
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
        FinimizerIndex::QueryResult result = index.search(reader.read_buf, left, right);
        //FinimizerIndex::QueryResult result = index.MS(reader.read_buf, left, right);
        //vector<tuple<int64_t, int64_t, int64_t>> MS_answer = MatchingStatistics_stats(index.sbwt, index.LCS, reader.read_buf, left, right);


        //reverse compl
        //const string reverse = sbwt::get_rc(reader.read_buf);
        //FinimizerIndex::QueryResult r_result = index.search(reverse, left, right);
        int64_t tot_kmers = result.local_offsets.size();
        /* int64_t str_len = reverse.length(); // the string and its reverse complement have the same length
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
        out << '\n'; */

        kmers_count += result.n_found;
        //kmers_count_rev += r_result.n_found;
        number_of_queries += tot_kmers;//result.local_offsets.size();
     
        total_micros += cur_time_micros() - t0;
    }
    write_log("k " + to_string(k), LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    write_log("Found kmers: " + to_string(kmers_count), LogLevel::MAJOR);
    //write_log("Found kmers reverse : " + to_string(kmers_count_rev), LogLevel::MAJOR);
    write_log("Total found kmers: " + to_string(total_positive), LogLevel::MAJOR);

    std::ofstream statsfile;
    statsfile.open(stats_filename, std::ios_base::app); // append instead of overwrite
    //statsfile << to_string(k) + "," + to_string(kmers_count+kmers_count_rev) + "," + to_string(number_of_queries);
    statsfile.close();
    return number_of_queries;
}

template<typename reader_t, typename out_stream_t>
int64_t run_fmin_file(const string& infile, out_stream_t& out, const string& stats_filename, const FinimizerIndex& index, vector<int64_t>& left, vector<int64_t>& right){
    reader_t reader(infile);
    write_log("Running streaming queries from input file " + infile, LogLevel::MAJOR);
    return run_fmin_queries_streaming(reader, out, index, stats_filename, left, right);
}

// Returns number of queries executed
int64_t run_fmin_queries(const vector<string>& infiles, const optional<vector<string>>& outfiles, const string& stats_filename, const FinimizerIndex& index, vector<int64_t>& left, vector<int64_t>& right ){

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
                n_queries_run += run_fmin_file<in_gzip>(infiles[i], out, stats_filename, index, left, right);
            } else { // To stdout
                n_queries_run += run_fmin_file<in_gzip>(infiles[i], cout, stats_filename, index, left, right);
            }
        }
        else {
            if(outfiles.has_value()){
                ofstream out(outfiles.value()[i]);
                n_queries_run += run_fmin_file<in_no_gzip>(infiles[i], out, stats_filename, index, left, right);
            } else{ // To stdout
                n_queries_run += run_fmin_file<in_no_gzip>(infiles[i], cout, stats_filename, index, left, right);
            }
        }
    }
    return n_queries_run;
}

int search_fmin(int argc, char** argv){

    int64_t micros_start = cur_time_micros();

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Query all Finimizers of all input reads.");

    options.add_options()
        ("o,out-file", "Output filename, or stdout if not given.", cxxopts::value<string>())
        ("i,index-file", "Index filename prefix.", cxxopts::value<string>())
        ("q,query-file", "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())
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

    string index_prefix = opts["index-file"].as<string>();

    int64_t number_of_queries = 0;

    cerr << "Loading index..." << endl;
    FinimizerIndex index;
    index.load(index_prefix);
    cerr << "Index loaded" << endl;

    const int64_t k = index.sbwt->get_k();
    cerr << "k = "<< to_string(k);
    cerr << " SBWT nodes: "<< to_string(index.sbwt->number_of_subsets())<< " kmers: "<< to_string(index.sbwt->number_of_kmers())<< endl;

    vector<int64_t> left, right;
    number_of_queries += run_fmin_queries(query_files, output_files, index_prefix + ".stats", index, left, right);

    int64_t new_total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)new_total_micros / number_of_queries), LogLevel::MAJOR);
    write_log("total number of queries: " + to_string(number_of_queries), LogLevel::MAJOR);
    
    write_csv(outfile + "LCS_Lscan.csv", "left", left);
    write_csv(outfile + "LCS_Rscan.csv", "right", right);

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

    std::ofstream statsfile3;
    statsfile3.open(outfile + "LCSstats.csv", std::ios_base::app); // append instead of overwrite
    // left mean
    int64_t left_sum = std::accumulate(left.begin(), left.end(), 0);
    double left_mean = static_cast<double>(left_sum) / left.size();
    // left standard deviation
    double l_sumOfSquares = 0;
    for (const auto& number : left) {
        double diff = static_cast<double>(number) - left_mean;
        l_sumOfSquares += diff * diff;
    }
    // Calculate the variance and take the square root to get the standard deviation
    double l_variance = l_sumOfSquares / (left.size() - 1);
    double l_standardDeviation = std::sqrt(l_variance);

    auto l_max = max_element(std::begin(left), std::end(left));
    auto l_min = min_element(std::begin(left), std::end(left));
    std::string left_max = std::to_string(*l_max);
    std::string left_min = std::to_string(*l_min);
    // right mean
    int64_t right_sum = std::accumulate(right.begin(), right.end(), 0);

    double right_mean = static_cast<double>(right_sum) / right.size();
    //right sd
    double r_sumOfSquares = 0;
    for (const auto& number : right) {
        double diff = static_cast<double>(number) - right_mean;
        r_sumOfSquares += diff * diff;
    }
    // Calculate the variance and take the square root to get the standard deviation
    double r_variance = r_sumOfSquares / (right.size() - 1);
    double r_standardDeviation = std::sqrt(r_variance);

    auto r_max = max_element(std::begin(right), std::end(right));
    auto r_min = min_element(std::begin(left), std::end(left));
    std::string right_max = std::to_string(*r_max);
    std::string right_min = std::to_string(*r_min);

    std::vector<int64_t> result;
    result.reserve(left.size());
    int64_t tot_sum = right_sum + left_sum;
    for (size_t i = 0; i < left.size(); ++i) {
        result.push_back(left[i] + right[i]);
    }
    //double tot_mean = static_cast<double>((tot_sum))/left.size();

    double mean = static_cast<double>((left_mean+right_mean))/2; // left and right have the same length
    auto t_max = max_element(std::begin(result), std::end(result));
    auto t_min = min_element(std::begin(result), std::end(result));
    std::string tot_max = std::to_string(*t_max);
    std::string tot_min = std::to_string(*t_min);


    string max = (*l_max > *r_max)? left_max: right_max;
    statsfile3 << "left mean = " + to_string((double)left_sum / left.size()) + "\n";
    statsfile3 << "left sd = " + to_string(double(l_standardDeviation)) + "\n";
    statsfile3 << "left max " + left_max + "\n";
    statsfile3 << "left min " + left_min + "\n";
    statsfile3 << "right mean= " +to_string((double)right_sum / right.size()) + "\n";
    statsfile3 << "right sd = " + to_string(double(r_standardDeviation)) + "\n";
    statsfile3 << "right max " + right_max + "\n";
    statsfile3 << "right min " + right_min + "\n";
     statsfile3 << "tot max " + tot_max + "\n";
    statsfile3 << "tot min " + tot_min + "\n";
    statsfile3 << "tot mean = " + to_string((double)mean) + "\n";
    
    statsfile3 << to_string(index.sbwt->number_of_subsets()) + ", " + to_string(k) + "," ;
    statsfile3 << to_string((double)left_sum / left.size()) + ", ";
    statsfile3 << to_string((double)right_sum / right.size()) + ", ";
    statsfile3 << to_string((double)mean) + ",";
    statsfile3 << left_max + ",";
    statsfile3 << right_max + ",";
    statsfile3 << max + ",";
    statsfile3 << to_string(left.size()) + "\n";
    
    statsfile3.close();
    
    return 0;

}