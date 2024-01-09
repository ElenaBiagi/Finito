//
// Created by Biagi, Elena on 23.6.2023.
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
//#include "lcs_basic_algorithm.hpp"

using namespace std;
using namespace sbwt;


std::vector<std::string> get_available_variants_fmin(){
    // This currently works only with the plain matrix representation
    return {"plain-matrix"};//, "rrr-matrix", "mef-matrix", "plain-split", "rrr-split", "mef-split", "plain-concat", "mef-concat", "plain-subsetwt", "rrr-subsetwt"};
}
std::vector<std::string> get_available_types(){
    return {"rarest", "shortest", "optimal", "verify"};
}

void save_v(const std::string& filename, const sdsl::int_vector<>& v) {
    std::ofstream out(filename, std::ios::binary);
    sdsl::serialize(v, out);
    out.close();
}

void load_v(const std::string& filename, sdsl::int_vector<>& v) {
    std::ifstream in(filename, std::ios::binary);
    sdsl::load(v, in);
    in.close();
}

// TODO check if it is used
void save_bv(const std::string& filename, const sdsl::bit_vector& v) {
    std::ofstream out(filename, std::ios::binary);
    sdsl::serialize(v, out);
    out.close();
}

// TODO check if it is used
void load_bv(const std::string& filename, sdsl::bit_vector& v) {
    std::ifstream in(filename, std::ios::binary);
    sdsl::load(v, in);
    in.close();
}


//write the output as a fasta file
template<typename writer_t>
void write_fasta(const pair<string,string>& p, writer_t& out){
    char newline = '\n';
    out.write(">", 1);
    out.write(static_cast<const char*>(p.first.c_str()), p.first.length());        
    out.write(&newline, 1);
    out.write(static_cast<const char*>(p.second.c_str()), p.second.length());
    out.write(&newline, 1);
}

// TODO check if it is used
void write_fasta_old(const pair<string,string>& sequencePair, string& filename) {
    std::ofstream outputFile(filename + ".fa", std::ios_base::app);
    outputFile << ">" + sequencePair.first + '\n'; // Write sequence ID
    outputFile << sequencePair.second + '\n'; // Write sequence data
    outputFile.close();
}

// TODO check if it is used

void write_str_v(const vector<string>& myVector, const string& filename) {
    //if (filesystem::exists(filename)) { filesystem::remove(filename);}
    std::ofstream outputFile(filename);
    while (outputFile.is_open()) {
        for (int i=0; i<myVector.size(); i++){//(const auto& str : myVector) {
            outputFile << i << " " << myVector[i] << '\n';
        }
        outputFile.close();
    }
}

// Not used, but usefull for debugging

void print_LCS(const string& v,const string& fname) {
    std::ofstream csv_file(fname);
    // Write the contents of the int_vector to the CSV file
    for (size_t i = 0; i < v.size(); i++) {
        csv_file << v[i];
        if (i < v.size() - 1) {
            csv_file << "";
        }
    }
    csv_file.close();
}

// TODO check if it is used
template<typename writer_t>
void write_csv(const tuple<string, string, string, string,string>& p, writer_t& out) {
    //assert(out.is_open());
    char newline = '\n';
    std::string line = std::get<0>(p) + ',' + std::get<1>(p) + ',' + std::get<2>(p) + ',' + std::get<3>(p) + ',' + std::get<4>(p) + '\n';
    out.write(line.c_str(), line.length());
    out.write(&newline, 1);
}

set<tuple<int64_t,int64_t, int64_t>> verify_shortest_streaming_search(const plain_matrix_sbwt_t& sbwt, const string& input, const char t) {
    const int64_t n_nodes = sbwt.number_of_subsets();
    const int64_t k = sbwt.get_k();
    const vector<int64_t> &C = sbwt.get_C_array();
    int64_t freq;    
    const int64_t str_len = input.size();
    tuple<int64_t,int64_t, int64_t, int64_t> w_fmin = {k + 1, n_nodes, n_nodes, str_len}; // {len, freq, I start, start}
    tuple<int64_t,int64_t, int64_t, int64_t> new_fmin; // {len, freq, I start, start}

    set<tuple<int64_t,int64_t,int64_t>> count_all_w_fmin;

    pair<int64_t, int64_t> I = {0, n_nodes - 1};
    int64_t I_start;
    char c;
    char char_idx;
    // window
    for (int64_t i = 0; i <= str_len - k; i++) {
        w_fmin = {k + 1, n_nodes, n_nodes, str_len}; // {len, freq, I start, start}
        // starting pos for the window
        for (int64_t start = i; start < k + i; start ++){
            I = {0, n_nodes - 1};
            // ending pos for the window
            for (int64_t end = start; end < k + i; end++) {
                c = static_cast<char>(input[end] &~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
                char_idx = get_char_idx(c);
                I = sbwt.update_sbwt_interval(&c, 1, I);
                freq = (I.second - I.first + 1);
                I_start = I.first;
                if (freq <= t) { // We found something
                    new_fmin = {end - start + 1, freq, I_start,end };// end-k +1
                    if (new_fmin < w_fmin) { w_fmin = new_fmin; }
                }
            }
        }
        count_all_w_fmin.insert({get<0>(w_fmin), get<1>(w_fmin), get<2>(w_fmin)});// (length,freq,colex)
        //cerr << input.substr(i, k) << endl << input.substr(get<3>(w_fmin)-get<0>(w_fmin)+1, get<0>(w_fmin)) << endl;
        //write_fasta({input.substr(i, k)+ ' ' + to_string(get<1>(w_fmin)), input.substr(get<3>(w_fmin)-get<0>(w_fmin)+1, get<0>(w_fmin))}, writer);
    }
    return count_all_w_fmin;
}

set<tuple<int64_t,int64_t, int64_t>> build_shortest_streaming_search(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input, const char t, sdsl::bit_vector& fmin_found){ //const sdsl::bit_vector** DNA_bitvectors,
    const int64_t n_nodes = sbwt.number_of_subsets();
    const int64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();
    int64_t freq;
    deque<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin;
    const int64_t str_len = input.size();
    tuple<int64_t, int64_t, int64_t, int64_t> w_fmin = {k+2,n_nodes,n_nodes,str_len}; // {len, freq, I start, start}
    set<tuple<int64_t,int64_t, int64_t>> count_all_w_fmin;
    int64_t kmer = 0;
    int64_t start = 0;
    int64_t end;
    pair<int64_t, int64_t> I = {0, n_nodes - 1};
    int64_t I_start;
    tuple<int64_t, int64_t, int64_t, int64_t> curr_substr;
    char c;
    char char_idx;
    for (end = 0; end < str_len; end++) {
        c = static_cast<char>(input[end] & ~32); // convert to uppercase using a bitwise operation
        char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]]{
           std::cerr << "Error: unknown character: " << c << endl;
           std::cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
            return {};
        } else {
            //update the sbwt INTERVAL
            I = sbwt.update_sbwt_interval(&c, 1, I);
            freq = (I.second - I.first + 1);
            I_start = I.first;
            
            if (freq <= t){ // 1. rarest
                while (freq <= t){ // 2. shortest
                curr_substr = {end - start + 1,freq, I_start, end};
                //all_fmin.insert(curr_substr); // insert every unique substr

                // (2) drop the first char
                // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                start++;
                I = drop_first_char(end - start + 1, I, LCS, n_nodes);
                freq = (I.second - I.first + 1);
                I_start = I.first;
                }
                if (w_fmin > curr_substr) {
                    all_fmin.clear();
                    w_fmin = curr_substr;
                } else{
                    while (all_fmin.back() > curr_substr) {
                        all_fmin.pop_back();
                    }
                }
                all_fmin.push_back(curr_substr);
 
            }
        }
        if (end >= k -1 ){
            if (!fmin_found[get<2>(w_fmin)]){ // if the kmer never been found before
                count_all_w_fmin.insert({get<0>(w_fmin), get<1>(w_fmin), get<2>(w_fmin) });// (length,freq,colex)
                if (get<3>(w_fmin) >= k-1){
                    fmin_found[get<2>(w_fmin)] = 1;
                }
            }
            //cerr << input.substr(kmer,k) << endl << input.substr(get<3>(w_fmin)-get<0>(w_fmin)+1,get<0>(w_fmin)) << endl;
            //write_fasta({input.substr(kmer,k) + ' ' + to_string(get<0>(w_fmin)),input.substr(get<3>(w_fmin)-get<1>(w_fmin)+1,get<1>(w_fmin))},writer);
            kmer++;
            // Check if the current minimizer is still in this window
            while (get<3>(w_fmin)- get<0>(w_fmin)+1 < kmer) { // start
                all_fmin.pop_front();
                if (all_fmin.empty()){
                    w_fmin={k+1,n_nodes,n_nodes,kmer+k};//str_len+1}; // {len, freq, I start, end}
                }
                else{ 
                    w_fmin = all_fmin.front();
                }
            }
        }
    }
    return count_all_w_fmin;
}

template<typename reader_t>
void print_shortest_finimizer_stats(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, reader_t& reader, int64_t t){
    set<tuple<int64_t,int64_t, int64_t>> count_all_w_fmin;
    const int64_t n_nodes = sbwt.number_of_subsets();
    sdsl::bit_vector fmin_found(n_nodes,0);
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) [[unlikely]] break;
        set<tuple<int64_t,int64_t, int64_t>> w_fmin = build_shortest_streaming_search(sbwt, LCS, reader.read_buf, t, fmin_found);
        count_all_w_fmin.insert(w_fmin.begin(), w_fmin.end());
    }
    print_finimizer_stats(count_all_w_fmin, sbwt.number_of_kmers(), sbwt.number_of_subsets(), t);
}

vector<string> remove_ns(const string& unitig, const int64_t k){
    vector<string> new_unitigs;
    const int64_t str_len = unitig.size();
    int64_t start = 0;
    char c;
    char char_idx;
    for (int64_t i = 0; i < str_len;i++){
        c = static_cast<char>(unitig[i] &~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]] {
            if ((i - start + 1) >= k ){
                string new_seq = unitig.substr(start,(i - start + 1));
                new_unitigs.push_back(new_seq);
            }
            start = i + 1;
        }
    }
    if ((str_len - start) >= k ){
        string new_seq = unitig.substr(start,(str_len - start));
        new_unitigs.push_back(new_seq);
    }
    return new_unitigs;
}

template<typename sbwt_t, typename reader_t>
void run_fmin_streaming(reader_t& f_reader, reader_t& r_reader, const string& index_prefix, unique_ptr<sbwt_t> sbwt, unique_ptr<sdsl::int_vector<>> LCS, const char t, const string& type){


    if(type == "rarest"){
        if(t != 1){
            throw std::runtime_error("t != 1 does not make sense with rarest type");
        }
        FinimizerIndexBuilder builder(move(sbwt), move(LCS), f_reader, r_reader);
        unique_ptr<FinimizerIndex> index = builder.get_index();
        index->serialize(index_prefix);
    } /* else if(type == "shortest"){
        // Just print stats because we don't have an index for this yet
        print_shortest_finimizer_stats(*sbwt, *LCS, reader, t);
    } else if(type == "verify"){
        // Print stats on shortest finimizers based on a reference implementation
        set<tuple<int64_t,int64_t, int64_t>> finimizers;
        while(true){
            int64_t len = reader.get_next_read_to_buffer();
            if(len == 0) [[unlikely]] break;
            vector<string> preprocessed_unitigs = remove_ns(reader.read_buf, sbwt->get_k());
            for (string& seq : preprocessed_unitigs){
                set<tuple<int64_t, int64_t, int64_t>> new_search = verify_shortest_streaming_search(*sbwt, seq, t);
                finimizers.insert(new_search.begin(), new_search.end());
            }
        }
        print_finimizer_stats(finimizers, sbwt->number_of_kmers(), sbwt->number_of_subsets(), t);    
    }*/

}

template<typename sbwt_t, typename reader_t>
int64_t run_file_fmin(const string& f_infile, const string& r_infile, const string& index_prefix, unique_ptr<sbwt_t> sbwt, unique_ptr<sdsl::int_vector<>> LCS, const char t, const string& type){ // const vector<int64_t>& C, const int64_t k,const sdsl::bit_vector** DNA_bitvectors,
    reader_t f_reader(f_infile);
    reader_t r_reader(r_infile);
    //Assume sbwt has streaming support
    write_log("Searching Finimizers from input file " + f_infile + " & " + r_infile + " to index prefix " + index_prefix, LogLevel::MAJOR);
    run_fmin_streaming<sbwt_t, reader_t>(f_reader, r_reader, index_prefix, move(sbwt), move(LCS), t, type); // C,k,DNA_bitvectors,
    return 0;
}

template<typename sbwt_t>
int64_t fmin_search(const vector<string>& f_infiles,const vector<string>& r_infiles, const string& out_prefix, unique_ptr<sbwt_t> sbwt, unique_ptr<sdsl::int_vector<>> LCS, const char t,const string& type){//const vector<int64_t>& C, const int64_t k,const sdsl::bit_vector** DNA_bitvectors,

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_no_gzip;

    int64_t n_fmin = 0;
    for(int64_t i = 0; i < f_infiles.size(); i++){
        bool gzip_input = SeqIO::figure_out_file_format(f_infiles[i]).gzipped;
        if(gzip_input){
            n_fmin += run_file_fmin<sbwt_t, in_gzip>(f_infiles[i], r_infiles[i], out_prefix, move(sbwt), move(LCS),t, type); 
        }
        else {
            n_fmin += run_file_fmin<sbwt_t, in_no_gzip>(f_infiles[i], r_infiles[i], out_prefix, move(sbwt), move(LCS),t, type);
        }
    }

    return n_fmin;

}

int build_fmin(int argc, char** argv) {

    int64_t micros_start = cur_time_micros();

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Find all Finimizers of all input reads.");

    vector<string> types = get_available_types();
    string all_types_string;
    for (string type: types) all_types_string += " " + type;

    options.add_options()
        ("o,out-file", "Output index filename prefix.", cxxopts::value<string>())
        ("i,index-file", "SBWT file. This has to be a binary matrix.", cxxopts::value<string>())
        ("f,f-file",
            "The unitigs in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.",
            cxxopts::value<string>())
        ("r,r-file", "reverse complement of f-file", cxxopts::value<string>())

        ("type", "Decide which streaming search type you prefer. Available types: " + all_types_string, cxxopts::value<string>()->default_value("rarest"))
        ("t", "Maximum finimizer frequency", cxxopts::value<int64_t>())
        ("lcs", "Provide in input the LCS file if available.", cxxopts::value<string>()->default_value(""))
        ("h,help", "Print usage");

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")) {
        std::cerr << options.help() << std::endl;
        exit(1);
    }
    char t = opts["t"].as<int64_t>();

    // forward unitigs
    string f_file = opts["f-file"].as<string>();
    vector<string> f_files;
    //TODO remove multifiles
    bool multi_file = f_file.size() >= 4 && f_file.substr(f_file.size() - 4) == ".txt";
    if(multi_file){
        f_files = readlines(f_file);
    } else{
        f_files = {f_file};
    }
    for(string file : f_files) check_readable(file);

    // reverse unitigs
    string r_file = opts["r-file"].as<string>();
    vector<string> r_files;
    //TODO remove multifiles
    multi_file = r_file.size() >= 4 && r_file.substr(r_file.size() - 4) == ".txt";
    if(multi_file){
        r_files = readlines(r_file);
    } else{
        r_files = {r_file};
    }
    for(string file : r_files) check_readable(file);


    string out_prefix = opts["out-file"].as<string>();

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

        unique_ptr<plain_matrix_sbwt_t> sbwt = make_unique<plain_matrix_sbwt_t>();
        sbwt->load(in.stream);

        string type = opts["type"].as<string>();
        if(std::find(types.begin(), types.end(), type) == types.end()){
            cerr << "Error: unknown type: " << type << endl;
            cerr << "Available types are:" << all_types_string << endl;
            return 1;
        }

        string LCS_file = opts["lcs"].as<string>();
        if (LCS_file.empty()) {
            std::cerr<< "LCS_file empty" << std::endl;
            LCS_file = out_prefix + ".LCS.sdsl";
            const sdsl::int_vector<> LCS = lcs_basic_parallel_algorithm(*sbwt, 8);
            //const sdsl::int_vector<> LCS = lcs_basic_algorithm(*sbwt);
            save_v(LCS_file, LCS);
        }
        unique_ptr<sdsl::int_vector<>> LCS = make_unique<sdsl::int_vector<>>();
        load_v(LCS_file, *LCS);
        std::cerr<< "LCS_file loaded" << std::endl;
        fmin_search(f_files, r_files, out_prefix, move(sbwt), move(LCS), t, type);//DNA_bitvectors,

    }
    return 0;
}