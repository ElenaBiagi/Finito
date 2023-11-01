//
// Created by Biagi, Elena on 23.6.2023.
//
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

void save_bv(const std::string& filename, const sdsl::bit_vector& v) {
    std::ofstream out(filename, std::ios::binary);
    sdsl::serialize(v, out);
    out.close();
}

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


void write_fasta_old(const pair<string,string>& sequencePair, string& filename) {
    std::ofstream outputFile(filename + ".fa", std::ios_base::app);
    outputFile << ">" + sequencePair.first + '\n'; // Write sequence ID
    outputFile << sequencePair.second + '\n'; // Write sequence data
    outputFile.close();
}

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

template<typename writer_t>
void write_csv(const tuple<string, string, string, string,string>& p, writer_t& out) {
    //assert(out.is_open());
    char newline = '\n';
    std::string line = std::get<0>(p) + ',' + std::get<1>(p) + ',' + std::get<2>(p) + ',' + std::get<3>(p) + ',' + std::get<4>(p) + '\n';
    out.write(line.c_str(), line.length());
    out.write(&newline, 1);
}

char get_char_idx(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

pair<int64_t,int64_t> update_sbwt_interval(const int64_t C_char, const pair<int64_t,int64_t>& I, const sdsl::rank_support_v5<>& Bit_rs){
    if(I.first == -1) return I;
    pair<int64_t,int64_t> new_I;
    // both start and end are included
    new_I.first = C_char + Bit_rs(I.first);
    new_I.second = C_char + Bit_rs(I.second+1) -1;
    if(new_I.first > new_I.second){
        return {-1,-1}; // Not found
    } 
    return new_I;
}

pair<int64_t,int64_t> drop_first_char(const int64_t  new_len, const pair<int64_t,int64_t>& I, const sdsl::int_vector<>& LCS, const int64_t n_nodes){
    if(I.first == -1) return I;
    pair<int64_t,int64_t> new_I = I;
    //Check top and bottom w the LCS
    while (LCS[new_I.first] >= new_len ){new_I.first --;}
    while(new_I.second < (n_nodes - 1) && LCS[new_I.second + 1] >= new_len ){
        new_I.second ++;
    }
    return {new_I};
}

template<typename writer_t>
set<tuple<int64_t,int64_t, int64_t>> verify_shortest_streaming_search( const sdsl::rank_support_v5<>** DNA_rs, const plain_matrix_sbwt_t& sbwt, const string& input, const char t, writer_t& writer) {
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
                const sdsl::rank_support_v5<> &Bit_rs = *(DNA_rs[char_idx]);
                I = update_sbwt_interval(C[char_idx], I, Bit_rs);
                freq = (I.second - I.first + 1);
                I_start = I.first;
                if (freq <= t) { // We found something
                    new_fmin = {end - start + 1, freq, I_start, end-k +1};
                    if (new_fmin < w_fmin) { w_fmin = new_fmin; }
                }
            }
        }
        count_all_w_fmin.insert({get<0>(w_fmin), get<1>(w_fmin), get<2>(w_fmin)});// (length,freq,colex)
        write_fasta({input.substr(i, k)+ ' ' + to_string(get<1>(w_fmin)), input.substr(get<3>(w_fmin), get<0>(w_fmin))}, writer);
    }
    return count_all_w_fmin;
}

template<typename writer_t>
set<tuple<int64_t,int64_t, int64_t>> build_rarest_streaming_search( const sdsl::rank_support_v5<>** DNA_rs, const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input, const char t, writer_t& writer, sdsl::bit_vector& fmin_bv, vector<uint32_t>& unitigs_k, const int64_t id){ //const sdsl::bit_vector** DNA_bitvectors,
    const int64_t n_nodes = sbwt.number_of_subsets();
    const int64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();
    int64_t freq;
    set<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin;
    const int64_t str_len = input.size();
    tuple<int64_t, int64_t, int64_t, int64_t> w_fmin = {n_nodes,k+1,n_nodes,str_len}; // {freq, len, I start, start}
    set<tuple<int64_t,int64_t, int64_t>> count_all_w_fmin;

    int64_t kmer = 0;
    int64_t start = 0;
    int64_t end;
    pair<int64_t, int64_t> I = {0, n_nodes - 1};
    int64_t I_start;
    tuple<int64_t, int64_t, int64_t, int64_t> curr_substr;
    char c;
    char char_idx;
    //string writer = "rarest_";
    // the idea is to start from the first pos which is i and move until finding something of ok freq
    // then drop the first char keeping track of which char you are starting from
    // Start is always < k as start <= end and end <k
    // if start == end than the frequency higher than t
    for (end = 0; end < str_len; end++) {
        c = static_cast<char>(input[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]]{
           std::cerr << "Error: unknown character: " << c << endl;
           std::cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
            return {};
        } else {
            //const sdsl::bit_vector& Bit_v = *(DNA_bitvectors[char_idx]);
            const sdsl::rank_support_v5<> &Bit_rs = *(DNA_rs[char_idx]);
            //update the sbwt INTERVAL
            I = update_sbwt_interval(C[char_idx], I, Bit_rs);
            freq = (I.second - I.first + 1);
            I_start = I.first;
            if (freq == 1){ // 1. rarest 
                while (freq == 1) {  //2. shortest
                curr_substr = {freq, end - start + 1, I_start, end-k +1};
                // (2) drop the first char
                // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                start++;
                I = drop_first_char(end - start + 1, I, LCS, n_nodes);
                freq = (I.second - I.first + 1);
                I_start = I.first;
                }
                if (w_fmin > curr_substr) {w_fmin = curr_substr;}
                all_fmin.insert(curr_substr);//({start, end - start + 1, freq, static_cast<int64_t>(I.first)});
            }
            
        }
        if (end > k){
            size_t old_fmin_count = count_all_w_fmin.size();
            count_all_w_fmin.insert({get<1>(w_fmin),get<0>(w_fmin), get<2>(w_fmin) });// (length,freq,colex) freq = 1 thus == (freq, length,colex)
            if (old_fmin_count != count_all_w_fmin.size()){
                fmin_bv[get<2>(w_fmin)]=1;
                if ((id + get<3>(w_fmin))> ULLONG_MAX){
                    std::cerr<< "ISSUE: global offset exceedes the allowed bit range." << std::endl;
                }
                unitigs_k[get<2>(w_fmin)]= id + get<3>(w_fmin);
            }                
            //write_fasta({input.substr(kmer,k) + ' ' + to_string(get<0>(w_fmin)),input.substr(get<3>(w_fmin),get<1>(w_fmin))},writer);
            kmer++;
            // Check if the current minimizer is still in this window
            while ((get<3>(w_fmin)+k-get<1>(w_fmin)) < kmer) { // start
                all_fmin.erase(all_fmin.begin());
                if (all_fmin.empty()){
                    w_fmin={n_nodes,k+1,kmer+1,str_len};
                }
                else{ 
                    w_fmin = *all_fmin.begin();
                }
            }
        }
    }
    return count_all_w_fmin;
}

template<typename writer_t>
set<tuple<int64_t,int64_t, int64_t>> build_shortest_streaming_search( const sdsl::rank_support_v5<>** DNA_rs, const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input, const char t, writer_t& writer){ //const sdsl::bit_vector** DNA_bitvectors,
    const int64_t n_nodes = sbwt.number_of_subsets();
    const int64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();
    int64_t freq;
    set<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin;
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
            //const sdsl::bit_vector& Bit_v = *(DNA_bitvectors[char_idx]);
            const sdsl::rank_support_v5<> &Bit_rs = *(DNA_rs[char_idx]);
            //update the sbwt INTERVAL
            I = update_sbwt_interval(C[char_idx], I, Bit_rs);
            freq = (I.second - I.first + 1);
            I_start = I.first;

            if (freq <= t){ // 1. rarest
                while (freq <= t){ // 2. shortest
                curr_substr = {end - start + 1,freq, I_start, end-k +1};
                //all_fmin.insert(curr_substr); // insert every unique substr

                // (2) drop the first char
                // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                start++;
                I = drop_first_char(end - start + 1, I, LCS, n_nodes);
                freq = (I.second - I.first + 1);
                I_start = I.first;
               //cerr << input.substr(start,end-start+1) << " freq=" << to_string(freq) << endl;
                }
                if (w_fmin > curr_substr) {w_fmin = curr_substr;}
                all_fmin.insert(curr_substr);
            }
        }
        if (end > k){
            count_all_w_fmin.insert({get<0>(w_fmin),get<1>(w_fmin), get<2>(w_fmin) });// (length,freq,colex)
            //write_fasta({input.substr(kmer,k) + ' ' + to_string(get<1>(w_fmin)),input.substr(get<3>(w_fmin),get<0>(w_fmin))},writer);
            kmer++;
            // Check if the current minimizer is still in this window
            while ((get<3>(w_fmin)+k-get<1>(w_fmin)) < kmer) {
                all_fmin.erase(all_fmin.begin());
                if (all_fmin.empty()){
                    w_fmin={k+1,n_nodes,kmer+1,str_len};//place holder, will never be selected
                }
                else{ 
                    w_fmin = *all_fmin.begin();
                }
            }
        }
    }
    return count_all_w_fmin;
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

template<typename sbwt_t, typename reader_t, typename writer_t>
sdsl::bit_vector run_fmin_streaming(reader_t& reader, writer_t& writer, const string& indexfile, const sbwt_t& sbwt, const sdsl::rank_support_v5<>** DNA_rs,  const sdsl::int_vector<>& LCS, const char t, const string& type){ // const vector<int64_t>& C, const int64_t k, const sdsl::bit_vector** DNA_bitvectors,
    int64_t new_number_of_fmin = 0;

    pair<PackedStrings, sdsl::bit_vector> unitig_data = permute_unitigs(sbwt, reader, indexfile);
    PackedStrings& unitigs = unitig_data.first;
    sdsl::bit_vector& Ustart = unitig_data.second;

    set<tuple<int64_t, int64_t, int64_t>> new_search;
    set<tuple<int64_t, int64_t, int64_t>>  finimizers;
    const int64_t k = sbwt.get_k();

    int64_t id = 0;
    int64_t n_nodes = sbwt.number_of_subsets();
    sdsl::bit_vector fmin_bv(n_nodes);
    vector<uint32_t> unitigs_k; // global offsets
    unitigs_k.reserve(n_nodes);
    unitigs_k.resize(n_nodes, 0);
    vector<int64_t> endpoints;
    vector<char> read_buf;
    //sdsl::int_vector endpoints;
    if (type == "rarest"){
        for(int64_t unitig_idx = 0; unitig_idx < unitigs.number_of_strings(); unitig_idx++){
            uint32_t len = unitigs.get(unitig_idx, read_buf);
            vector<string> preprocessed_unitigs = remove_ns(read_buf.data(), k);
            for (string& seq : preprocessed_unitigs){
                cout << seq << endl;
                new_search = build_rarest_streaming_search(DNA_rs, sbwt ,LCS,seq, t, writer, fmin_bv, unitigs_k, id);
                new_number_of_fmin += new_search.size();
                finimizers.insert(new_search.begin(), new_search.end());
            }
            id += len; //end point
            endpoints.push_back(id);// first char of the next string 
        }
       //cerr << to_string(id) << endl;
    } else if (type == "shortest") {
        for(int64_t unitig_idx = 0; unitig_idx < unitigs.number_of_strings(); unitig_idx++){
            uint32_t len = unitigs.get(unitig_idx, read_buf);
            vector<string> preprocessed_unitigs = remove_ns(read_buf.data(), k);
            for (string& seq : preprocessed_unitigs){
                new_search = build_shortest_streaming_search(DNA_rs, sbwt ,LCS,seq, t, writer);
                new_number_of_fmin += new_search.size();
                finimizers.insert(new_search.begin(), new_search.end());
            }    
        }
    } else if (type == "verify"){
        for(int64_t unitig_idx = 0; unitig_idx < unitigs.number_of_strings(); unitig_idx++){
            uint32_t len = unitigs.get(unitig_idx, read_buf);
            vector<string> preprocessed_unitigs = remove_ns(read_buf.data(), k);
            for (string& seq : preprocessed_unitigs){
                new_search = verify_shortest_streaming_search(DNA_rs, sbwt ,seq, t, writer);
                new_number_of_fmin += new_search.size();
                finimizers.insert(new_search.begin(), new_search.end());
            }
        }
    }
    // (length,freq,colex)
    new_number_of_fmin = finimizers.size();
    int64_t sum_freq = 0;
    int64_t sum_len = 0;
    for (auto x : finimizers){
        sum_freq += get<1>(x);
        sum_len += get<0>(x);
        //fmin_bv[get<2>(x)] = 1;
    }
    tuple<string, string, string, string, string> stats = {to_string(t),to_string( new_number_of_fmin),to_string(sum_freq), to_string(static_cast<float>(sum_freq)/static_cast<float>(new_number_of_fmin)), to_string(static_cast<float>(sum_len)/static_cast<float>(new_number_of_fmin)) };
    std::ofstream outfile;
    string results = to_string(new_number_of_fmin) + "," + to_string(sum_freq) + "," + to_string(static_cast<float>(sum_freq) / static_cast<float>(new_number_of_fmin)) + "," + to_string(static_cast<float>(sum_len) / static_cast<float>(new_number_of_fmin)) + "," + to_string(sbwt.number_of_kmers());

    try {
        outfile.open(indexfile + "test.txt", std::ios_base::app);
        if (outfile.is_open()) {
            outfile << to_string(t) + "," + results + "\n";
            outfile.close();
        } else {
            std::cerr << "Error opening the file for writing." << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
    // outfile.open(indexfile + "test.txt", std::ios_base::app); // append instead of overwrite
    // string results = to_string(new_number_of_fmin)+ "," + to_string(sum_freq) + "," + to_string(static_cast<float>(sum_freq)/static_cast<float>(new_number_of_fmin)) + "," + to_string(static_cast<float>(sum_len)/static_cast<float>(new_number_of_fmin));
    // outfile << to_string(t) + "," + results + "\n";
    // outfile.close();
    write_log(to_string(t) + "," +results, LogLevel::MAJOR);
    write_log("#SBWT nodes: " + to_string( n_nodes) , LogLevel::MAJOR);
    write_log("#Distinct finimizers: " + to_string( new_number_of_fmin) , LogLevel::MAJOR);
    write_log("Sum of frequencies: " + to_string(sum_freq) , LogLevel::MAJOR);
    write_log("Avg frequency: " + to_string(static_cast<float>(sum_freq)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);
    write_log("Avg length: " + to_string(static_cast<float>(sum_len)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);

    if(type == "rarest"){
        // endpoints
        sdsl::int_vector<> endpoints_v(endpoints.size(), 0, 64-__builtin_clzll(id));
        for (int64_t x = 0; x < endpoints.size(); x++){
            endpoints_v[x]=endpoints[x];
        }

        // global offsets
        sdsl::int_vector<> unitigs_v(new_number_of_fmin, 0, 64 - __builtin_clzll(id));
        int64_t j=0;
        for (int64_t i = 0; i < n_nodes; i++){
            if (fmin_bv[i] == 1){
                unitigs_v[j] = unitigs_k[i];
                j++;
            }
        }
    
        save_v(indexfile + "O.sdsl", unitigs_v);
        save_v(indexfile + "E.sdsl", endpoints_v);
        save_bv(indexfile + "FBV.sdsl", fmin_bv);

        std::ofstream packed_unitigs_out(indexfile + "packed_unitigs.sdsl");
        sdsl::serialize(unitigs.concat, packed_unitigs_out);
        
        std::ofstream unitig_endpoints_out(indexfile + "unitig_endpoints.sdsl");
        sdsl::serialize(unitigs.ends, unitig_endpoints_out);

        std::ofstream Ustart_out(indexfile + "Ustart.sdsl");
        sdsl::serialize(Ustart, Ustart_out);
    }
    
    return fmin_bv;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_file_fmin(const string& infile, const string& outfile, const string& indexfile, const sbwt_t& sbwt, const sdsl::rank_support_v5<>** DNA_rs, const sdsl::int_vector<>& LCS, const char t, const string& type){ // const vector<int64_t>& C, const int64_t k,const sdsl::bit_vector** DNA_bitvectors,
    reader_t reader(infile);
    writer_t writer(outfile);
    //Assume sbwt has streaming support
    write_log("Searching Finimizers from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
    sdsl::bit_vector fmin_bv = run_fmin_streaming<sbwt_t, reader_t, writer_t>(reader, writer, indexfile, sbwt,  DNA_rs, LCS, t, type); // C,k,DNA_bitvectors,
    return 0;
}

template<typename sbwt_t>
int64_t fmin_search(const vector<string>& infiles, const string& outfile, const string& indexfile, const sbwt_t& sbwt,  const sdsl::rank_support_v5<>** DNA_rs, const sdsl::int_vector<>& LCS, const char t,const string& type, bool gzip_output){//const vector<int64_t>& C, const int64_t k,const sdsl::bit_vector** DNA_bitvectors,

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_no_gzip;

    typedef Buffered_ofstream<zstr::ofstream> out_gzip;
    typedef Buffered_ofstream<std::ofstream> out_no_gzip;

    int64_t n_fmin = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = SeqIO::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input && gzip_output){
            n_fmin += run_file_fmin<sbwt_t, in_gzip, out_gzip>(infiles[i], outfile, indexfile, sbwt,  DNA_rs, LCS,t, type); // DNA_bitvectors,
        }
        if(gzip_input && !gzip_output){
            n_fmin += run_file_fmin<sbwt_t, in_gzip, out_no_gzip>(infiles[i], outfile, indexfile, sbwt,DNA_rs, LCS,t, type);
        }
        if(!gzip_input && gzip_output){
            n_fmin += run_file_fmin<sbwt_t, in_no_gzip, out_gzip>(infiles[i], outfile, indexfile, sbwt,DNA_rs, LCS,t, type);
        }
        if(!gzip_input && !gzip_output){
            n_fmin += run_file_fmin<sbwt_t, in_no_gzip, out_no_gzip>(infiles[i], outfile, indexfile, sbwt, DNA_rs, LCS,t, type);
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
            ("o,out-file", "Output filename.", cxxopts::value<string>())
            ("i,index-file", "Index input file.This has to be a binary matrix.", cxxopts::value<string>())
            ("u,in-file",
             "The SPSS in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.",
             cxxopts::value<string>())
            ("z,gzip-output",
             "Writes output in gzipped form. This can shrink the output files by an order of magnitude.",
             cxxopts::value<bool>()->default_value("false"))
            ("type", "Decide which streaming search type you prefer. Available types: " + all_types_string, cxxopts::value<string>()->default_value("rarest"))
            ("t", "Maximum finimizer frequency", cxxopts::value<int64_t>())
            ("lcs", "Provide in input the LCS file if available.", cxxopts::value<string>()->default_value(""))
            ("h,help", "Print usage");

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")) {
        std:://cerr << options.help() << std::endl;
        exit(1);
    }
    char t = opts["t"].as<int64_t>();

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


    // Interpret output file
    string outfile = opts["out-file"].as<string>() + to_string(t)+".fa";
    bool gzip_output = opts["gzip-output"].as<bool>();
    vector<string> output_files;
    check_writable(outfile);

    // sbwt
    string indexfile = opts["index-file"].as<string>();
    check_readable(indexfile);
    throwing_ifstream in(indexfile, ios::binary);
    vector<string> variants = get_available_variants_fmin();
    string variant = load_string(in.stream); // read variant type
    if (std::find(variants.begin(), variants.end(), variant) == variants.end()) {
       //cerr << "Error loading index from file: unrecognized variant specified in the file" << endl;
        return 1;
    }

    write_log("Loading the index variant " + variant, LogLevel::MAJOR);
    if (variant == "plain-matrix") {

        plain_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);

        string type = opts["type"].as<string>();
        if(std::find(types.begin(), types.end(), type) == types.end()){
            cerr << "Error: unknown type: " << type << endl;
            cerr << "Available types are:" << all_types_string << endl;
            return 1;
        }

        const sdsl::bit_vector& A_bits = sbwt.get_subset_rank_structure().A_bits;
        const sdsl::bit_vector& C_bits = sbwt.get_subset_rank_structure().C_bits;
        const sdsl::bit_vector& G_bits = sbwt.get_subset_rank_structure().G_bits;
        const sdsl::bit_vector& T_bits = sbwt.get_subset_rank_structure().T_bits;

        const sdsl::rank_support_v5<> &A_bits_rs = sbwt.get_subset_rank_structure().A_bits_rs;
        const sdsl::rank_support_v5<> &C_bits_rs = sbwt.get_subset_rank_structure().C_bits_rs;
        const sdsl::rank_support_v5<> &G_bits_rs = sbwt.get_subset_rank_structure().G_bits_rs;
        const sdsl::rank_support_v5<> &T_bits_rs = sbwt.get_subset_rank_structure().T_bits_rs;

        const int64_t n_nodes = sbwt.number_of_subsets();
        const int64_t k = sbwt.get_k();


        const vector<int64_t> &C = sbwt.get_C_array();

        //const sdsl::bit_vector *DNA_bitvectors[4] = {&A_bits, &C_bits, &G_bits, &T_bits};
        const sdsl::rank_support_v5<> *DNA_rs[4] = {&A_bits_rs, &C_bits_rs, &G_bits_rs, &T_bits_rs};

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
        fmin_search(input_files, outfile, indexfile, sbwt, DNA_rs, LCS,t, type, gzip_output);//DNA_bitvectors,

    }
    return 0;
}