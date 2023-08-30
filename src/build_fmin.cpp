//
// Created by Biagi, Elena on 23.6.2023.
//
#include <string>
#include <tuple>
#include <set>
#include <string.h>
#include "sbwt/cxxopts.hpp"
#include "sbwt/globals.hh"
#include "sbwt/SBWT.hh"
#include "sbwt/stdlib_printing.hh"
#include "sbwt/SeqIO.hh"
#include "sbwt/buffered_streams.hh"
#include <vector>
#include <utility>
#include <algorithm>
#include "variants.hh"
#include "sdsl/vectors.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/bit_vectors.hpp"
#include "sbwt/commands.hh"
#include "sbwt/suffix_group_optimization.hh"

using namespace std;
using namespace sbwt;

std::vector<std::string> get_available_variants_fmin(){
    return {"plain-matrix"};//, "rrr-matrix", "mef-matrix", "plain-split", "rrr-split", "mef-split", "plain-concat", "mef-concat", "plain-subsetwt", "rrr-subsetwt"};
}
std::vector<std::string> get_available_types(){
    return {"rarest", "shortest", "optimal", "verify"};
}

//lcs
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

void save_intv(const std::string& filename, const std::vector<uint64_t>& v) {
    std::ofstream out(filename, std::ios::binary);
    for (const auto& p : v) {
        out << p << std::endl;
    }
    out.close();
}
void load_intv(const std::string& filename, std::vector<uint64_t>& v) {
    std::ifstream in(filename, std::ios::binary);
    v.clear();
    int64_t p;
    while (in >> p) {
        v.emplace_back(p);
    }
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

sdsl::int_vector<> get_kmer_lcs(const sdsl::bit_vector& A_bits,
                                const sdsl::bit_vector& C_bits,
                                const sdsl::bit_vector& G_bits,
                                const sdsl::bit_vector& T_bits,
                                int64_t k){
    int64_t n_nodes = A_bits.size();
    vector<int64_t> C_array(4);
    vector<char> last; // last[i] = incoming character to node i
    last.push_back('$');

    C_array[0] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(A_bits[i]) last.push_back('A');
    C_array[1] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(C_bits[i]) last.push_back('C');
    C_array[2] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(G_bits[i]) last.push_back('G');
    C_array[3] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(T_bits[i]) last.push_back('T');
    if(last.size() != n_nodes){
        cerr << "BUG " << last.size() << " " << n_nodes << endl;
        exit(1);
    }

    sdsl::bit_vector mismatch_found_marks(n_nodes, 0);
    sdsl::int_vector<> lcs(n_nodes, 0, 64 - __builtin_clzll(k-1)); // Enough bits per element to store values from 0 to k-1

    for(int64_t round = 0; round < k; round++){
        for(int64_t i = 0; i < n_nodes; i++){
            if(mismatch_found_marks[i] == 0 && (i == 0 || last[i] != last[i-1])){
                mismatch_found_marks[i] = 1;
                lcs[i] = round;
            }
        }

        // Propagate the labels one step forward in the graph
        vector<char> propagated(n_nodes, '$');
        int64_t A_ptr = C_array[0];
        int64_t C_ptr = C_array[1];
        int64_t G_ptr = C_array[2];
        int64_t T_ptr = C_array[3];
        for(int64_t i = 0; i < n_nodes; i++){
            if(A_bits[i]) propagated[A_ptr++] = last[i];
            if(C_bits[i]) propagated[C_ptr++] = last[i];
            if(G_bits[i]) propagated[G_ptr++] = last[i];
            if(T_bits[i]) propagated[T_ptr++] = last[i];
        }
        last = propagated;
    }

    return lcs;
}

int64_t get_char_idx(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

pair<int64_t,int64_t> update_sbwt_interval(char char_idx, const pair<int64_t,int64_t>& I, const sdsl::rank_support_v5<>& Bit_rs, const vector<int64_t>& C){
    //write_log("before updating interval",  LogLevel::MAJOR);
    if(I.first == -1) return I;
    pair<int64_t,int64_t> new_I;
    // both start and end are included
    new_I.first = C[char_idx] + Bit_rs(I.first);
    new_I.second = C[char_idx] + Bit_rs(I.second+1) -1;
    if(new_I.first > new_I.second) return {-1,-1}; // Not found
    return new_I;
}

pair<int64_t,int64_t> drop_first_char(const uint64_t  new_len, const pair<int64_t,int64_t>& I, const sdsl::int_vector<>& LCS, const uint64_t n_nodes){
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
set<tuple<uint64_t,uint64_t, uint64_t>> verify_shortest_streaming_search( const sdsl::rank_support_v5<>** DNA_rs, const plain_matrix_sbwt_t& sbwt, const string& input, const uint64_t t, writer_t& writer) {
    const uint64_t n_nodes = sbwt.number_of_subsets();
    const uint64_t k = sbwt.get_k();
    const vector<int64_t> &C = sbwt.get_C_array();
    uint64_t freq;    
    const uint64_t str_len = input.size();
    tuple<uint64_t,uint64_t, uint64_t, uint64_t> w_fmin = {k + 1, n_nodes, n_nodes, str_len}; // {len, freq, I start, start}
    tuple<uint64_t,uint64_t, uint64_t, uint64_t> new_fmin; // {len, freq, I start, start}

    set<tuple<uint64_t,uint64_t,uint64_t>> count_all_w_fmin;

    pair<int64_t, int64_t> I = {0, n_nodes - 1};
    uint64_t I_start;
    // window
    for (uint64_t i = 0; i <= str_len - k; i++) {
        w_fmin = {k + 1, n_nodes, n_nodes, str_len}; // {len, freq, I start, start}
        // starting pos for the window
        for (uint64_t start = i; start < k + i; start ++){
            I = {0, n_nodes - 1};
            // ending pos for the window
            for (uint64_t end = start; end < k + i; end++) {
                char c = static_cast<char>(input[end] &~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
                int64_t char_idx = get_char_idx(c);
                const sdsl::rank_support_v5<> &Bit_rs = *(DNA_rs[char_idx]);
                I = update_sbwt_interval(char_idx, I, Bit_rs, C);
                freq = (I.second - I.first + 1);
                I_start = I.first;
                if (freq <= t) { // We found something
                    new_fmin = {end - start + 1, freq, I_start, start};
                    if (new_fmin < w_fmin) { w_fmin = new_fmin; }
                }
            }
        }
        count_all_w_fmin.insert({get<0>(w_fmin), get<1>(w_fmin), get<2>(w_fmin)});// (length,freq,colex)
        //count_all_w_fmin.insert({get<2>(w_fmin), get<0>(w_fmin), get<1>(w_fmin)});// Istart, len, freq
        write_fasta({input.substr(i, k)+ ' ' + to_string(get<1>(w_fmin)), input.substr(get<3>(w_fmin), get<0>(w_fmin))}, writer);
    }
    return count_all_w_fmin;
}

template<typename writer_t>
set<tuple<uint64_t,uint64_t, uint64_t>> build_rarest_streaming_search( const sdsl::rank_support_v5<>** DNA_rs, const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input, const uint64_t t, writer_t& writer, sdsl::bit_vector& fmin_bv, vector<uint64_t>& unitigs_k, uint64_t id){ //const sdsl::bit_vector** DNA_bitvectors,
    const uint64_t n_nodes = sbwt.number_of_subsets();
    const uint64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();
    uint64_t freq;
    set<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> all_fmin;
    const uint64_t str_len = input.size();
    tuple<uint64_t, uint64_t, uint64_t, uint64_t> w_fmin = {n_nodes,k+1,n_nodes,str_len}; // {freq, len, I start, start}
    //vector<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> all_w_fmin;
    set<tuple<uint64_t,uint64_t, uint64_t>> count_all_w_fmin;

    uint64_t kmer = 0;
    uint64_t start = 0;
    uint64_t end;
    pair<int64_t, int64_t> I = {0, n_nodes - 1};
    uint64_t I_start;
    tuple<uint64_t, uint64_t, uint64_t, uint64_t> curr_substr;
    //string writer = "rarest_";
    // the idea is to start from the first pos which is i and move until finding something of ok freq
    // then drop the first char keeping track of which char you are starting from
    // Start is always < k as start <= end and end <k
    // if start == end than the frequency higher than t
    for (end = 0; end < str_len; end++) {
        char c = static_cast<char>(input[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        int64_t char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]]{
            cerr << "Error: unknown character: " << c << endl;
            cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
            return {};
        } else {
            //const sdsl::bit_vector& Bit_v = *(DNA_bitvectors[char_idx]);
            const sdsl::rank_support_v5<> &Bit_rs = *(DNA_rs[char_idx]);
            //update the sbwt INTERVAL
            I = update_sbwt_interval(char_idx, I, Bit_rs, C);
            freq = (I.second - I.first + 1);
            I_start = I.first;
            while (freq == 1) { // We found something
                curr_substr = {freq, end - start + 1, I_start, start};
                all_fmin.insert(curr_substr);//({start, end - start + 1, freq, static_cast<uint64_t>(I.first)});

                // Check window fmin
                // 1. rarest (freq=1), 2. shortest
                if (w_fmin > curr_substr) {w_fmin = curr_substr;}

                // (2) drop the first char
                // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                start++;
                // todo check if it is necessary to ensure that start is always <= end or if it is so by construction as if start==end than the freq is >>1
                if (start > end)[[unlikely]] { //todo check if by saying start++ here start increases
                    break;
                }
                I = drop_first_char(end - start + 1, I, LCS, n_nodes);
                freq = (I.second - I.first + 1);
                I_start = I.first;
            }
        }
        // TODO check if this is making it so slow
        if (end >= k-1){
                // Check if the current minimizer is still in this window
                while (get<3>(w_fmin) < kmer) {
                all_fmin.erase(all_fmin.begin());
                w_fmin = *all_fmin.begin();
                }
                size_t old_fmin_count = count_all_w_fmin.size();
                // 
                count_all_w_fmin.insert({get<1>(w_fmin),get<0>(w_fmin), get<2>(w_fmin) });// (length,freq,colex) freq = 1 thus == (freq, length,colex)
                //count_all_w_fmin.insert({get<2>(w_fmin),get<1>(w_fmin), get<0>(w_fmin)});// Istart, len, freq
                if (old_fmin_count != count_all_w_fmin.size()){
                    fmin_bv[get<2>(w_fmin)]=1;
                    if (id >= (1ULL << 36) || start >= (1ULL << 28)) {
                    std::cerr << "ISSUE: One or both numbers exceed the allowed bit range." << std::endl;
                    }
                    unitigs_k[get<2>(w_fmin)]= (id << 28) | get<3>(w_fmin);
                }
                //cerr << input.substr(kmer,k) <<' ' << to_string(get<0>(w_fmin)) << "fmin = " << input.substr(get<3>(w_fmin),get<1>(w_fmin));
                write_fasta({input.substr(kmer,k) + ' ' + to_string(get<0>(w_fmin)),input.substr(get<3>(w_fmin),get<1>(w_fmin))},writer);
                kmer++;
        }
    }
    return count_all_w_fmin;
}

template<typename writer_t>
set<tuple<uint64_t,uint64_t,uint64_t>> build_shortest_streaming_search( const sdsl::rank_support_v5<>** DNA_rs, const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input, const uint64_t t, writer_t& writer){ //const sdsl::bit_vector** DNA_bitvectors,
    //string writer = "shortest_";
    const uint64_t n_nodes = sbwt.number_of_subsets();
    const uint64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();
    uint64_t freq;
    set<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> all_fmin; // later check if [1] and [3] are the same
    const uint64_t str_len = input.size();
    tuple<uint64_t, uint64_t, uint64_t, uint64_t> w_fmin = {k+2,n_nodes,n_nodes,str_len}; // {len, freq, I start, start}
    set<tuple<uint64_t, uint64_t, uint64_t>> count_all_w_fmin;
    uint64_t kmer = 0;

    // (1) Calculate all possible fmin in the FIRST window
    uint64_t start = 0;
    uint64_t end;
    pair<int64_t, int64_t> I = {0, n_nodes - 1};
    uint64_t I_start;
    tuple<uint64_t, uint64_t, uint64_t, uint64_t> curr_substr;
    // the idea is to start from the first pos which is i and move until finding something of ok freq
    // then drop the first char keeping track of which char you are starting from
    // Start is always < k as start <= end and end <k
    // if start == end than the frequency higher than t
    for (end = 0; end <= str_len; end++) {
        char c = static_cast<char>(input[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        int64_t char_idx = get_char_idx(c);
            const sdsl::rank_support_v5<> &Bit_rs = *(DNA_rs[char_idx]);
            //update the sbwt INTERVAL
            I = update_sbwt_interval(char_idx, I, Bit_rs, C);
            freq = (I.second - I.first + 1);
            I_start = I.first;
            while (freq <= t) { // We found something
                curr_substr = {end - start + 1, freq, I_start, start};
                all_fmin.insert(curr_substr);
                // 1. shortest, 2. rarest (freq <= t)
                if (curr_substr < w_fmin){
                    w_fmin = curr_substr;
                }
                // 2. drop the first char
                // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                start++;
                I = drop_first_char(end - start + 1, I, LCS, n_nodes);
                freq = (I.second - I.first + 1);
                I_start = I.first;
            }
            if (end >= k-1){
                while (get<3>(w_fmin) < kmer) { //else keep the current w_fmin
                    all_fmin.erase(all_fmin.begin());
                    w_fmin = *all_fmin.begin();
                }
                count_all_w_fmin.insert({get<0>(w_fmin), get<1>(w_fmin),get<2>(w_fmin),});// (length,freq,colex)
                //count_all_w_fmin.insert({get<2>(w_fmin),get<0>(w_fmin), get<1>(w_fmin)});// Istart, len, freq 
                // add the last interval found to avoid computing it again
                write_fasta({input.substr(kmer,k) + ' ' + to_string(get<1>(w_fmin)),input.substr(get<3>(w_fmin),get<0>(w_fmin))},writer);
                kmer++;
            }

    }
    return count_all_w_fmin;
}

vector<string> remove_ns(const string& unitig, const uint64_t k){
    vector<string> new_unitigs;
    const uint64_t str_len = unitig.size();
    uint64_t start = 0;
    for (uint64_t i = 0; i < str_len;i++){
        char c = static_cast<char>(unitig[i] &~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        int64_t char_idx = get_char_idx(c);
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
sdsl::bit_vector run_fmin_streaming(reader_t& reader, writer_t& writer, const sbwt_t& sbwt, const sdsl::rank_support_v5<>** DNA_rs,  const sdsl::int_vector<>& LCS, const uint64_t t, const string& type){ // const vector<int64_t>& C, const int64_t k, const sdsl::bit_vector** DNA_bitvectors,
    int64_t new_number_of_fmin = 0;
    
    set<tuple<uint64_t, uint64_t, uint64_t>> new_search;
    set<tuple<uint64_t, uint64_t, uint64_t>>  finimizers;
    const uint64_t k = sbwt.get_k();

    int64_t id = 0;
    uint64_t n_nodes = sbwt.number_of_subsets();
    sdsl::bit_vector fmin_bv(n_nodes);
    vector<uint64_t> unitigs_k; // 36 and 28 bits 
    unitigs_k.reserve(n_nodes);
    unitigs_k.resize(n_nodes, 0);
    if (type == "rarest"){
        while(true) {
            int64_t len = reader.get_next_read_to_buffer();
            if(len == 0) [[unlikely]] break;
            id++;
            vector<string> preprocessed_unitigs = remove_ns(reader.read_buf, k);
            for (string& seq : preprocessed_unitigs){
                new_search = build_rarest_streaming_search(DNA_rs, sbwt ,LCS,seq, t, writer, fmin_bv, unitigs_k, id);
                new_number_of_fmin += new_search.size();
                finimizers.insert(new_search.begin(), new_search.end());
            }
        }
    } else if (type == "shortest") {
        while(true) {
            int64_t len = reader.get_next_read_to_buffer();
            if(len == 0) [[unlikely]] break;
            vector<string> preprocessed_unitigs = remove_ns(reader.read_buf, k);
            for (string& seq : preprocessed_unitigs){
                new_search = build_shortest_streaming_search(DNA_rs, sbwt ,LCS,seq, t, writer);
                new_number_of_fmin += new_search.size();
                finimizers.insert(new_search.begin(), new_search.end());
            }
        }
    } else if (type == "verify"){
        while(true) {
            int64_t len = reader.get_next_read_to_buffer();
            if(len == 0) [[unlikely]] break;
            vector<string> preprocessed_unitigs = remove_ns(reader.read_buf, k);
            for (string& seq : preprocessed_unitigs){
                new_search = verify_shortest_streaming_search(DNA_rs, sbwt ,seq, t, writer);
                new_number_of_fmin += new_search.size();
                finimizers.insert(new_search.begin(), new_search.end());
            }
        }
    }

    // (length,freq,colex)
    new_number_of_fmin = finimizers.size();
    uint64_t sum_freq = 0;
    uint64_t sum_len = 0;
    for (auto x : finimizers){
        sum_freq += get<1>(x);
        sum_len += get<0>(x);
        //fmin_bv[get<2>(x)] = 1;
    }
    tuple<string, string, string, string, string> stats = {to_string(t),to_string( new_number_of_fmin),to_string(sum_freq), to_string(static_cast<float>(sum_freq)/static_cast<float>(new_number_of_fmin)), to_string(static_cast<float>(sum_len)/static_cast<float>(new_number_of_fmin)) };
    writer_t writer_stats("stats.csv");
    write_csv(stats, writer_stats);
    write_log("#Distinct finimizers: " + to_string( new_number_of_fmin) , LogLevel::MAJOR);
    write_log("Sum of frequencies: " + to_string(sum_freq) , LogLevel::MAJOR);
    write_log("Avg frequency: " + to_string(static_cast<float>(sum_freq)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);
    write_log("Avg length: " + to_string(static_cast<float>(sum_len)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);

    // Resize the vector for unitigs (id, offset)
    vector<uint64_t> unitigs_v; // 36 and 28 bits 
    unitigs_v.reserve(new_number_of_fmin);
    unitigs_v.resize(new_number_of_fmin);
    uint64_t j=0;
    for (uint64_t i = 0; i < n_nodes; i++){
        if (fmin_bv[i] == 1){
            unitigs_v[j] = unitigs_k[i];
            j++;
        }
    }

    save_intv("fmin_unitigs", unitigs_v);
    save_bv("fmin_bv", fmin_bv);
    return fmin_bv;
}


// find the kmer and look for it streaming all unitigs until you found it
// OR 
// when finding finimizers I know unitig and offset butt I do not know the order
// write them in a larger vector (Istart position) and remove 0 values or neg values
// scan fmin_bv and keep only values of unitig_v where  fmin_bv[i]==1

/* 
void find_rarest_streaming_search( const sdsl::rank_support_v5<>** DNA_rs, const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input, const uint64_t t, const sdsl::bit_vector& fmin_bv, const sdsl::rank_support_v5<>& fmin_rs, vector<uint64_t> unitigs_v, uint64_t id){ //const sdsl::int_vector<>& LCS
    const uint64_t n_nodes = sbwt.number_of_subsets();
    const uint64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();
    pair<int64_t, int64_t> I = {0, n_nodes - 1};
    const uint64_t str_len = input.size();

    int64_t char_idx;
    uint64_t freq;
    uint64_t start = 0;
    uint64_t unitig;
    for (uint64_t end = 0; end < str_len; end++) {
        char c = static_cast<char>(input[end] &~32);
        char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]] {
            cerr << "Error: unknown character: " << c << endl;
            cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
            return;
        } else {
            const sdsl::rank_support_v5<> &Bit_rs = *(DNA_rs[char_idx]);
            //update the sbwt INTERVAL
            I = update_sbwt_interval(char_idx, I, Bit_rs, C);
            freq = (I.second - I.first + 1);
            while (freq == 1 && fmin_bv[I.first]) { // We found something
                if (id >= (1ULL << 36) || start >= (1ULL << 28)) {
                std::cerr << "One or both numbers exceed the allowed bit range." << std::endl;
                return;
                }
                unitig = (id << 28) | start;
                //cerr << to_string(unitig) << " = "<< to_string(id) <<", " << to_string(start) << endl;
                unitigs_v[fmin_rs(I.first)] = unitig;
                start++;
                // todo check if it is necessary to ensure that start is always <= end or if it is so by construction as if start==end than the freq is >>1
                if (start > end)[[unlikely]] {break;}
                I = drop_first_char(end - start + 1, I, LCS, n_nodes);
                freq = (I.second - I.first + 1);
            }
        }
    } 
}
 */

// template<typename sbwt_t, typename reader_t, typename writer_t>
// void find_fmin_streaming(reader_t& reader, const sbwt_t& sbwt, const sdsl::rank_support_v5<>** DNA_rs,  const sdsl::int_vector<>& LCS, const uint64_t t, const string& type, const sdsl::rank_support_v5<>& fmin_rs, const sdsl::bit_vector& fmin_bv){ // const vector<int64_t>& C, const int64_t k, const sdsl::bit_vector** DNA_bitvectors,
//     int64_t number_of_fmin = fmin_rs(sbwt.number_of_subsets()+1);
//     const uint64_t k = sbwt.get_k();
//     int64_t id = 0;
//     vector<uint64_t> unitigs_v; // 36 and 28 bits
//     unitigs_v.reserve(number_of_fmin);
//     unitigs_v.resize(number_of_fmin);
//     //if (type == "rarest"){
//     while(true) {
//             int64_t len = reader.get_next_read_to_buffer();
//             if(len == 0) [[unlikely]] {break;}
//             id++;
//             vector<string> preprocessed_unitigs = remove_ns(reader.read_buf, k);
//             for (string& seq : preprocessed_unitigs){
//                 find_rarest_streaming_search(DNA_rs, sbwt, LCS, seq, t, fmin_bv, fmin_rs, unitigs_v, id);
//             }
//     }
//     //}
//     // ORDERED BASED ON THE ORDER OF APPEARANCE IN THE SBWT -> RANK THEN FILL IN AN ARRAY OF THE SIZE OF SUM OF 1S IN FMIN_BV
//     // TODO modify also SHORTEST AND VERIFY     
//     /* } else if (type == "shortest") {
//         while(true) {
//             int64_t len = reader.get_next_read_to_buffer();
//             if(len == 0) [[unlikely]] break;
//             vector<string> preprocessed_unitigs = remove_ns(reader.read_buf, k);
//             for (string& seq : preprocessed_unitigs){
//                 new_search = build_shortest_streaming_search(DNA_rs, sbwt ,LCS,seq, t, writer);
//                 new_number_of_fmin += new_search.size();
//                 finimizers.insert(new_search.begin(), new_search.end());
//             }
//         }
//     } else if (type == "verify"){
//         while(true) {
//             int64_t len = reader.get_next_read_to_buffer();
//             if(len == 0) [[unlikely]] break;
//             vector<string> preprocessed_unitigs = remove_ns(reader.read_buf, k);
//             for (string& seq : preprocessed_unitigs){
//                 new_search = verify_shortest_streaming_search(DNA_rs, sbwt ,seq, t, writer);
//                 new_number_of_fmin += new_search.size();
//                 finimizers.insert(new_search.begin(), new_search.end());
//             }
//         }
//     }
//     */
//    save_intv("fmin_unitigs", unitigs_v);
// }


template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_file_fmin(const string& infile, const string& outfile, const sbwt_t& sbwt, const sdsl::rank_support_v5<>** DNA_rs, const sdsl::int_vector<>& LCS, const uint64_t t, const string& type){ // const vector<int64_t>& C, const int64_t k,const sdsl::bit_vector** DNA_bitvectors,
    reader_t reader(infile);
    writer_t writer(outfile);
    //Assume sbwt has streaming support
    write_log("Searching Finimizers from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
    sdsl::bit_vector fmin_bv = run_fmin_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt,  DNA_rs, LCS, t, type); // C,k,DNA_bitvectors,
    //sdsl::rank_support_v5<> fmin_rs(&fmin_bv);
    //find_fmin_streaming<sbwt_t, reader_t, writer_t>(reader2, sbwt,  DNA_rs, LCS, t, type, fmin_rs, fmin_bv); // C,k,DNA_bitvectors,
    return 0;
}

template<typename sbwt_t>
int64_t fmin_search(const vector<string>& infiles, const string& outfile, const sbwt_t& sbwt,  const sdsl::rank_support_v5<>** DNA_rs, const sdsl::int_vector<>& LCS, const uint64_t t,const string& type, bool gzip_output){//const vector<int64_t>& C, const int64_t k,const sdsl::bit_vector** DNA_bitvectors,

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_no_gzip;

    typedef Buffered_ofstream<zstr::ofstream> out_gzip;
    typedef Buffered_ofstream<std::ofstream> out_no_gzip;

    int64_t n_fmin = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = SeqIO::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input && gzip_output){
            n_fmin += run_file_fmin<sbwt_t, in_gzip, out_gzip>(infiles[i], outfile, sbwt,  DNA_rs, LCS,t, type); // DNA_bitvectors,
        }
        if(gzip_input && !gzip_output){
            n_fmin += run_file_fmin<sbwt_t, in_gzip, out_no_gzip>(infiles[i], outfile, sbwt,DNA_rs, LCS,t, type);
        }
        if(!gzip_input && gzip_output){
            n_fmin += run_file_fmin<sbwt_t, in_no_gzip, out_gzip>(infiles[i], outfile, sbwt,DNA_rs, LCS,t, type);
        }
        if(!gzip_input && !gzip_output){
            n_fmin += run_file_fmin<sbwt_t, in_no_gzip, out_no_gzip>(infiles[i], outfile, sbwt, DNA_rs, LCS,t, type);
        }
    }

    return n_fmin;

}

std::vector<std::string> dump_all_kmers(const sdsl::bit_vector& A_bits,
                                        const sdsl::bit_vector& C_bits,
                                        const sdsl::bit_vector& G_bits,
                                        const sdsl::bit_vector& T_bits,
                                        int64_t k){
    cerr << "Dumping k-mers"<< endl;
    int64_t n_nodes = A_bits.size();
    vector<int64_t> C_array(4);

    vector<char> last; // last[i] = incoming character to node i
    last.push_back('$');

    C_array[0] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(A_bits[i]) last.push_back('A');

    C_array[1] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(C_bits[i]) last.push_back('C');

    C_array[2] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(G_bits[i]) last.push_back('G');

    C_array[3] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(T_bits[i]) last.push_back('T');

    if(last.size() != n_nodes){
        cerr << "BUG " << last.size() << " " << n_nodes << endl;
        exit(1);
    }

    vector<string> all_kmers(n_nodes);
    for(string& S : all_kmers) S.resize(k);

    for(int64_t round = 0; round < k; round++){
        for(int64_t i = 0; i < n_nodes; i++){
            all_kmers[i][k-1-round] = last[i];
        }

        // Propagate the labels one step forward in the graph
        vector<char> propagated(n_nodes, '$');
        int64_t A_ptr = C_array[0];
        int64_t C_ptr = C_array[1];
        int64_t G_ptr = C_array[2];
        int64_t T_ptr = C_array[3];
        for(int64_t i = 0; i < n_nodes; i++){
            if(A_bits[i]) propagated[A_ptr++] = last[i];
            if(C_bits[i]) propagated[C_ptr++] = last[i];
            if(G_bits[i]) propagated[G_ptr++] = last[i];
            if(T_bits[i]) propagated[T_ptr++] = last[i];
        }
        last = propagated;
    }
    return all_kmers;
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
            ("i,index-file", "Index input file.", cxxopts::value<string>())
            ("u,in-file",
             "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.",
             cxxopts::value<string>())
            ("add-reverse-complements", "Also add the reverse complement of every k-mer to the index. Warning: this creates a temporary reverse-complemented duplicate of each input file before construction. Make sure that the directory at --temp-dir can handle this amount of data. If the input is gzipped, the duplicate will also be compressed, which might take a while.", cxxopts::value<bool>()->default_value("false"))
            ("z,gzip-output",
             "Writes output in gzipped form. This can shrink the output files by an order of magnitude.",
             cxxopts::value<bool>()->default_value("false"))
            ("type", "Decide which streaming search type you prefer. Available types: " + all_types_string,
             cxxopts::value<string>()->default_value("rarest"))
            ("t", "Maximum finimizer frequency", cxxopts::value<uint64_t>())
            ("lcs", "Provide in input the LCS file if available.", cxxopts::value<string>()->default_value(""))
            ("rmq", "For the option --lcs new-rmq provide in input the rmqLCS file if available.",
             cxxopts::value<string>()->default_value(""))
            ("h,help", "Print usage");

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")) {
        std::cerr << options.help() << std::endl;
        exit(1);
    }
    uint64_t t = opts["t"].as<uint64_t>();
    bool revcomps = opts["add-reverse-complements"].as<bool>();

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
    
    //TODO check if the outfile exists

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

        //int64_t new_k = sbwt.get_k();
        //vector<string> all_kmers = dump_all_kmers(A_bits, C_bits, G_bits, T_bits, new_k);
        //write_str_v( all_kmers, "all_kmers.txt");

        //const sdsl::bit_vector *DNA_bitvectors[4] = {&A_bits, &C_bits, &G_bits, &T_bits};
        const sdsl::rank_support_v5<> *DNA_rs[4] = {&A_bits_rs, &C_bits_rs, &G_bits_rs, &T_bits_rs};

        string LCS_file = opts["lcs"].as<string>();
        if (LCS_file.empty()) {
            std::cerr << "LCS_file empty" << std::endl;
            LCS_file = indexfile + "LCS.sdsl";
            const sdsl::int_vector<> LCS = get_kmer_lcs(A_bits, C_bits, G_bits, T_bits, k);
            save_v(LCS_file, LCS);
        }
        sdsl::int_vector<> LCS;
        load_v(LCS_file, LCS);
        std::cerr << "LCS_file loaded" << std::endl;
        fmin_search(input_files, outfile, sbwt, DNA_rs, LCS,t, type, gzip_output);//DNA_bitvectors,

        return 0;
    }
}
