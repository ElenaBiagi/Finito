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

class FinimizerIndex{
public:
    unique_ptr<plain_matrix_sbwt_t> sbwt; // These are smart pointers because they are passed in to the constructor
    unique_ptr<sdsl::int_vector<>> LCS; // These are smart pointers because they are passed in to the constructor
    PackedStrings unitigs;
    sdsl::bit_vector fmin;
    sdsl::rank_support_v5<> fmin_rs;
    sdsl::int_vector<> global_offsets;
    sdsl::bit_vector Ustart;
    sdsl::rank_support_v5<> Ustart_rs;

    FinimizerIndex() {}

    void search(const std::string& query) {
        // TODO
    }

    void serialize(const string& index_prefix) {
        std::ofstream global_offsets_out(index_prefix + ".O.sdsl");
        sdsl::serialize(global_offsets, global_offsets_out);

        std::ofstream fmin_out(index_prefix + ".FBV.sdsl");
        sdsl::serialize(fmin, fmin_out);

        std::ofstream packed_unitigs_out(index_prefix + ".packed_unitigs.sdsl");
        sdsl::serialize(unitigs.concat, packed_unitigs_out);
        
        std::ofstream unitig_endpoints_out(index_prefix + ".unitig_endpoints.sdsl");
        sdsl::serialize(unitigs.ends, unitig_endpoints_out);

        std::ofstream Ustart_out(index_prefix + ".Ustart.sdsl");
        sdsl::serialize(Ustart, Ustart_out);

        std::ofstream LCS_out(index_prefix + ".LCS.sdsl");
        sdsl::serialize(*LCS, LCS_out);

        sbwt->serialize(index_prefix + ".sbwt");
    }

    void load(const string& index_prefix) {

        LCS = make_unique<sdsl::int_vector<>>();
        ifstream LCS_in(index_prefix + ".LCS.sdsl");
        sdsl::load(*LCS, LCS_in);
        std::cerr<< "LCS_file loaded"<<std::endl;

        ifstream fmin_bv_in(index_prefix + ".FBV.sdsl");
        sdsl::load(fmin, fmin_bv_in);
        std::cerr<< "fmin_bv_file loaded"<<std::endl;
        sdsl::rank_support_v5<> fmin_rs(&fmin);

        ifstream global_offsets_in(index_prefix + ".O.sdsl");
        sdsl::load(global_offsets, global_offsets_in);
        std::cerr<< "offsets loaded"<<std::endl;

        std::ifstream packed_unitigs_in(index_prefix + ".packed_unitigs.sdsl");
        sdsl::load(unitigs.concat, packed_unitigs_in);
        std::cerr << "unitigs loaded" << std::endl;

        std::ifstream unitig_endpoints_in(index_prefix + ".unitig_endpoints.sdsl");
        sdsl::load(unitigs.ends, unitig_endpoints_in);
        std::cerr << "unitig endpoints loaded" << std::endl;

        std::ifstream Ustart_in(index_prefix + ".Ustart.sdsl");
        sdsl::load(Ustart, Ustart_in);
        sdsl::rank_support_v5<> Ustart_rs(&Ustart);
        std::cerr << "Ustart loaded" << std::endl;

        sbwt = make_unique<plain_matrix_sbwt_t>();
        sbwt->load(index_prefix + ".sbwt");

    }

    void size_in_bytes(){
        // TODO
    }
};


class FinimizerIndexBuilder{
public:

    unique_ptr<plain_matrix_sbwt_t> sbwt;
    unique_ptr<sdsl::int_vector<>> LCS;

    unique_ptr<FinimizerIndex> index;

    void print_stats(set<tuple<int64_t, int64_t, int64_t>> finimizers, int64_t n_kmers, int64_t n_nodes, int64_t t){
        int64_t new_number_of_fmin = finimizers.size();
        int64_t sum_freq = 0;
        int64_t sum_len = 0;
        for (auto x : finimizers){
            sum_freq += get<1>(x);
            sum_len += get<0>(x);
        }

        string results = to_string(new_number_of_fmin) + "," + to_string(sum_freq) + "," + to_string(static_cast<float>(sum_freq) / static_cast<float>(new_number_of_fmin)) + "," + to_string(static_cast<float>(sum_len) / static_cast<float>(new_number_of_fmin)) + "," + to_string(n_kmers);

        write_log(to_string(t) + "," + results, LogLevel::MAJOR);
        write_log("#SBWT nodes: " + to_string(n_nodes) , LogLevel::MAJOR);
        write_log("#Distinct finimizers: " + to_string(new_number_of_fmin) , LogLevel::MAJOR);
        write_log("Sum of frequencies: " + to_string(sum_freq) , LogLevel::MAJOR);
        write_log("Avg frequency: " + to_string(static_cast<float>(sum_freq)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);
        write_log("Avg length: " + to_string(static_cast<float>(sum_len)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);
    }

    // Takes ownership of sbwt and LCS
    template<typename reader_t, typename writer_t>
    FinimizerIndexBuilder(unique_ptr<plain_matrix_sbwt_t> sbwt, unique_ptr<sdsl::int_vector<>> LCS, reader_t& reader, writer_t& writer) {
        index = make_unique<FinimizerIndex>();
        this->sbwt = move(sbwt); // Take ownership
        this->LCS = move(LCS); // Take ownership

        int64_t n_nodes = this->sbwt->number_of_subsets();
        sdsl::bit_vector fmin_bv(n_nodes, 0); // Finimizer marks
        sdsl::bit_vector fmin_found(n_nodes, 0);
        vector<uint32_t> global_offsets; // TODO: -> 64 bits
        global_offsets.reserve(n_nodes);
        global_offsets.resize(n_nodes, 0);
        
        pair<PackedStrings, sdsl::bit_vector> unitig_data = permute_unitigs(*(this->sbwt), reader, "unused_parameter"); // TODO: remove the unused parameter 
        PackedStrings& unitigs = unitig_data.first;
        sdsl::bit_vector& Ustart = unitig_data.second;

        set<tuple<int64_t, int64_t, int64_t>>  finimizers;
        int64_t total_len = 0;
        vector<char> unitig_buf;
        for(int64_t i = 0; i < unitigs.number_of_strings(); i++){
            int64_t len = unitigs.get(i, unitig_buf);
            set<tuple<int64_t, int64_t, int64_t>> new_search = add_sequence(unitig_buf.data(), fmin_bv, fmin_found, global_offsets, total_len, writer);
            total_len += len; 
            finimizers.insert(new_search.begin(), new_search.end());
        }

        sdsl::int_vector<> packed_global_offsets(global_offsets.size(), 64 - __builtin_clzll(total_len));
        for(int64_t i = 0; i < global_offsets.size(); i++){
            packed_global_offsets[i] = global_offsets[i];
        }
    
        print_stats(finimizers, this->sbwt->number_of_kmers(), this->sbwt->number_of_subsets(), 1);

        index->sbwt = std::move(this->sbwt); // Transfer ownership
        index->LCS = std::move(this->LCS); // Transfer ownership 
        index->unitigs = std::move(unitigs); // Transfer ownership
        index->fmin = std::move(fmin_bv); // Transfer ownership
        index->fmin_rs = sdsl::rank_support_v5<>(&(index->fmin));
        index->global_offsets = std::move(packed_global_offsets); // Transfer ownership
        index->Ustart = std::move(Ustart); // Transfer ownership
        index->Ustart_rs = sdsl::rank_support_v5<>(&(index->Ustart)); 

        
    }

    // TODO: 32 bits for global offsets might not be enough
    template<typename writer_t>
    set<tuple<int64_t, int64_t, int64_t>> add_sequence(const std::string& seq, sdsl::bit_vector& fmin_bv, sdsl::bit_vector& fmin_found, vector<uint32_t>& global_offsets, const int64_t unitig_start, writer_t& writer) {
        const int64_t n_nodes = sbwt->number_of_subsets();
        const int64_t k = sbwt->get_k();
        const vector<int64_t>& C = sbwt->get_C_array();
        int64_t freq;
        set<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin;
        const int64_t str_len = seq.size();
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
        // the idea is to start from the first pos which is i and move until finding something of ok freq
        // then drop the first char keeping track of which char you are starting from
        // Start is always < k as start <= end and end <k
        // if start == end than the frequency higher than t
        for (end = 0; end < str_len; end++) {
            c = static_cast<char>(seq[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
            char_idx = get_char_idx(c);
            if (char_idx == -1) [[unlikely]]{
               std::cerr << "Error: unknown character: " << c << endl;
               std::cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
                return {};
            } else {
                //const sdsl::rank_support_v5<> &Bit_rs = *(DNA_rs[char_idx]);
                //update the sbwt INTERVAL
                I = this->sbwt->update_sbwt_interval(&c, 1, I);
                freq = (I.second - I.first + 1);
                I_start = I.first;
                if (freq == 1){ // 1. rarest 
                    while (freq == 1) {  //2. shortest
                    curr_substr = {freq, end - start + 1, I_start, end};
                    // (2) drop the first char
                    // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                    start++;
                    I = drop_first_char(end - start + 1, I, *(this->LCS), n_nodes);
                    freq = (I.second - I.first + 1);
                    I_start = I.first;
                    }
                    if (w_fmin > curr_substr) {w_fmin = curr_substr;}
                    all_fmin.insert(curr_substr);//({start, end - start + 1, freq, static_cast<int64_t>(I.first)});
                }
            }
            if (end >= k -1 ){
                count_all_w_fmin.insert({get<1>(w_fmin),get<0>(w_fmin), get<2>(w_fmin) });// (length,freq,colex) freq = 1 thus == (freq, length,colex)
                fmin_bv[get<2>(w_fmin)]=1;
                if (!fmin_found[get<2>(w_fmin)]){ // if the kmer never been found before
                    fmin_bv[get<2>(w_fmin)]=1;
                    if ((unitig_start + get<3>(w_fmin))> UINT32_MAX){
                        std::cerr<< "ISSUE: global offset exceedes the allowed bit range." << std::endl;
                    }
                    global_offsets[get<2>(w_fmin)]= unitig_start + get<3>(w_fmin);
                }
                if (get<3>(w_fmin) >= k-1){
                    fmin_found[get<2>(w_fmin)] = 1;
                }
                // write_fasta({input.substr(kmer,k) + ' ' + to_string(get<0>(w_fmin)),input.substr(get<3>(w_fmin)-get<1>(w_fmin)+1,get<1>(w_fmin))},writer);
                kmer++;
                // Check if the current minimizer is still in this window
                while (get<3>(w_fmin)- get<1>(w_fmin)+1 < kmer) { // start
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

    // Transfer ownership of the index out of the builder
    unique_ptr<FinimizerIndex> get_index(){
        return std::move(this->index);
    }
};