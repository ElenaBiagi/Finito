#pragma once

#include <string>
#include <cstring>
#include <algorithm>
#include <bitset>
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

    struct QueryResult{
        vector<pair<int64_t, int64_t>> local_offsets; // Unitig id, distance from the start of the unitig
        int64_t n_found = 0;
    };

private:
    // Forbid copying because we have pointers to our internal data structures
    FinimizerIndex(const FinimizerIndex& other) = delete;
    FinimizerIndex& operator=(const FinimizerIndex& other) = delete;

    void add_to_query_result(int64_t global_kmer_end, QueryResult& answer) const{
        int64_t global_kmer_start = global_kmer_end - sbwt->get_k() + 1;
        pair<int64_t, int64_t> local_start = unitigs.global_offset_to_local_offset(global_kmer_start);
        answer.local_offsets.push_back(local_start);
        answer.n_found++;
    }

    void walk_in_unitigs(const std::string& query, const PackedStrings& unitigs, int64_t global_kmer_end, QueryResult& answer, int64_t& kmer_end, const int64_t k) const{
        //cout << "take a nice unitig walk" << endl;
        int64_t unitig_id = answer.local_offsets.back().first;
        int64_t u_end = unitigs.ends[unitig_id]; // exlusive end
        int64_t max_match = std::min(u_end - global_kmer_end-1, (int64_t)(query.length()-kmer_end-1));

        if (global_kmer_end > u_end or max_match <= 0){return;}
        kmer_end++;
        //convert part of the query to an int_vector<2>
        sdsl::int_vector<2> query_v(max_match);
        
        for (int64_t i = 0; i < max_match;) {
            char c = query[kmer_end+i];
            switch(c){
                case 'A': query_v[i++] = 0; break;
                case 'C': query_v[i++] = 1; break;
                case 'G': query_v[i++] = 2; break;
                case 'T': query_v[i++] = 3; break;
                default: throw std::runtime_error("Invalid character: " + c);
            }
        }

        int64_t word_start = 0;
        int64_t global_kmer_end_copy = global_kmer_end;

        while(max_match > 0){
            int word_len = std::min({(int64_t)32, max_match});
            uint64_t query_word = query_v.get_int(word_start, (word_len)*2); // word start = least significant bit (word_len+1)*2-1
            uint64_t unitig_word = unitigs.concat.get_int(word_start + (global_kmer_end_copy+1)*2, (word_len)*2);
            //cerr << std::bitset<8 * sizeof(int64_t)>(query_word) << endl;
            //cerr << std::bitset<8 * sizeof(int64_t)>(unitig_word) << endl;


            int64_t result = query_word ^ unitig_word;
            //cerr << std::bitset<8 * sizeof(int64_t)>(result) << endl;

            if (result){
                int trailing_zeros = __builtin_ctzll(result);
                for (int i = 0; i < trailing_zeros/2; i++){
                    global_kmer_end++;
                    add_to_query_result(global_kmer_end, answer);
                    kmer_end++; // keep track of how many kmers were found    
                }
                break; // No need to check further. A mismatch has been found.
            }
            int trailing_zeros = (word_len)*2;
            for (int i = 0; i < trailing_zeros/2; i++){
                global_kmer_end++;
                add_to_query_result(global_kmer_end, answer);
                kmer_end++; // keep track of how many kmers were found    
            }
            max_match -= word_len;
            word_start += word_len*2;
        }
        kmer_end--; // will be updated later
    }


public:

    // Note: if you add members, update size_in_bytes(), serialize(), and load()
    unique_ptr<plain_matrix_sbwt_t> sbwt; // These are smart pointers because they are passed in to the constructor
    unique_ptr<sdsl::int_vector<>> LCS; // These are smart pointers because they are passed in to the constructor
    PackedStrings unitigs;
    sdsl::bit_vector fmin;
    sdsl::rank_support_v5<> fmin_rs;
    sdsl::int_vector<> global_offsets;
    sdsl::bit_vector Ustart;
    sdsl::rank_support_v5<> Ustart_rs;

    FinimizerIndex() {}

    QueryResult search(const std::string& query, vector<int64_t>& left, vector<int64_t>& right) const {
        // For each k-mer S that is known to be in the SBWT
        //   - Find the finimizer x.
        //   - Walk forward in the SBWT from the colex rank of x (singleton interval), to the end of S,
        //     recording the rightmost branch point, and the distance from the end of x to this branch
        //     point.
        //   - If a branch point was found, the k-mer containing x is at the unitig after the rightmost
        //     branch. Otherwise, the k-mer containing x is at the unitig that x points to.
        //   - If there was a branch, the k- mer endpoint in that unitig is k + the number of steps taken after the last branch.
        //     

        const plain_matrix_sbwt_t& sbwt = *(this->sbwt.get());
        const int64_t n_nodes = sbwt.number_of_subsets();
        const int64_t k = sbwt.get_k();
        const vector<int64_t>& C = sbwt.get_C_array();
        const int64_t query_len = query.length();

        QueryResult answer{{}, 0};

        if(query.size() < k) answer; 


        //Find kmers and Finimizers together
        tuple<vector<optional<int64_t>>,vector<optional< pair<int64_t, int64_t> > >, vector<optional< pair<int64_t, int64_t> >>>colex_finimizers = rarest_fmin_streaming_search(sbwt, *LCS, query, Ustart);
        vector<optional<int64_t>> kmer_colex_ranks = get<0>(colex_finimizers);
        vector<optional< pair<int64_t, int64_t>>> finimizers_ends_colex = get<1>(colex_finimizers);
        vector<optional< pair<int64_t, int64_t>>> rightmost_Ustart = get<2>(colex_finimizers);

        for(int64_t kmer_end = k-1; kmer_end < query_len; kmer_end++) {

            if(finimizers_ends_colex[kmer_end].has_value()){//kmer_colex_ranks[kmer_end].has_value()){
                // kmer exists
                int64_t global_kmer_end;
                int64_t finimizer_end = finimizers_ends_colex[kmer_end].value().first;
                
                //cout << "Finimizer ends at: " << finimizer_end << endl; //<< " with len " << shortest_unique_lengths[finimizer_end].value() << endl;
                //optional<pair<int64_t, int64_t>> rightmost_branch_end = get_rightmost_Ustart(query, kmer_end, finimizer_end, finimizers_ends_colex, sbwt, Ustart);
                optional<pair<int64_t, int64_t>> rightmost_branch_end = rightmost_Ustart[kmer_end];
                if(rightmost_branch_end.has_value()) {
                    // Look up from the branch dictionary
                    int64_t p = rightmost_branch_end.value().first;
                    int64_t colex = rightmost_branch_end.value().second; 
                    // Get the global off set of the end of the k-mer
                    global_kmer_end = lookup_from_branch_dictionary(colex, k, Ustart_rs, unitigs);
                    global_kmer_end += kmer_end - p; // Shift to the right place in the unitig
                } else {
                    //Look up from Finimizer dictionary
                    int64_t p = finimizer_end;
                    //int64_t colex = shortest_unique_colex_ranks[p].value();
                    int64_t colex = finimizers_ends_colex[kmer_end].value().second;

                    global_kmer_end = lookup_from_finimizer_dictionary(colex, fmin_rs, global_offsets);
                    global_kmer_end += kmer_end - p; // Shift to the right place in the unitig
                }
                // A kmer has been found 
                add_to_query_result(global_kmer_end, answer);
                if ((kmer_end + 1)< query.size()){
                    walk_in_unitigs(query, unitigs, global_kmer_end, answer, kmer_end, k);
                }
                
            } else { 
                answer.local_offsets.push_back({-1, -1});
            } 
        }
        return answer;
    }
 /* 
    QueryResult MS(const string& input, vector<int64_t>& left, vector<int64_t>& right)const{
        const plain_matrix_sbwt_t& sbwt = *(this->sbwt.get());
        const int64_t n_nodes = sbwt.number_of_subsets();
        const int64_t k = sbwt.get_k();
        const vector<int64_t>& C = sbwt.get_C_array();    
        pair<int64_t, int64_t> I = {0, n_nodes - 1};
        pair<int64_t, int64_t> I_new;
        QueryResult answer{{}, 0};

        const int64_t str_len = input.size();
        char c;

        int64_t d = 0; //lenght of the current match
        vector<tuple<int64_t, int64_t, int64_t>> MS;
        for (int64_t end = 0; end < str_len; end++) {
            c = static_cast<char>(input[end] & ~32); // convert to uppercase using a bitwise operation  
            I_new = sbwt.update_sbwt_interval(&c, 1, I);
            while(d > 0 && I_new.first == -1){
                d--;
                // Contract left
                I = drop_first_char_stats(d, I, *LCS, n_nodes, left, right);
                I_new = sbwt.update_sbwt_interval(&c, 1, I);
            }
            if (I_new.first != -1){
                I = I_new;
                d = min(k, d+1);
            }
            MS.push_back({d,(I.second - I.first + 1),I.first });
        }
        return answer;
    }
    */
    void serialize(const string& index_prefix) const {
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
        sdsl::util::init_support(fmin_rs, &fmin);

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
        sdsl::util::init_support(Ustart_rs, &Ustart);
        std::cerr << "Ustart loaded" << std::endl;

        sbwt = make_unique<plain_matrix_sbwt_t>();
        sbwt->load(index_prefix + ".sbwt");

    }

    // This also includes the rank structures which are not serialized
    int64_t size_in_bytes() const{
        int64_t total = 0;
        total += sdsl::size_in_bytes(*LCS);
        total += sdsl::size_in_bytes(fmin);
        total += sdsl::size_in_bytes(fmin_rs);
        total += sdsl::size_in_bytes(global_offsets);
        total += sdsl::size_in_bytes(unitigs.concat);
        total += sdsl::size_in_bytes(unitigs.ends);
        total += sdsl::size_in_bytes(Ustart);
        total += sdsl::size_in_bytes(Ustart_rs);

        sbwt::SeqIO::NullStream ns;
        total += sbwt->serialize(ns);
        return total;
    }
};


class FinimizerIndexBuilder{
public:

    unique_ptr<plain_matrix_sbwt_t> sbwt;
    unique_ptr<sdsl::int_vector<>> LCS;

    unique_ptr<FinimizerIndex> index;


    // Takes ownership of sbwt and LCS
    template<typename reader_t>
    FinimizerIndexBuilder(unique_ptr<plain_matrix_sbwt_t> sbwt, unique_ptr<sdsl::int_vector<>> LCS, reader_t& reader) {
        index = make_unique<FinimizerIndex>();
        this->sbwt = move(sbwt); // Take ownership
        this->LCS = move(LCS); // Take ownership

        int64_t n_nodes = this->sbwt->number_of_subsets();
        sdsl::bit_vector fmin_bv(n_nodes, 0); // Finimizer marks
        //sdsl::bit_vector fmin_found(n_nodes, 0);
        sdsl::int_vector fmin_found(n_nodes, 0);
        
        vector<uint64_t> global_offsets;
        global_offsets.reserve(n_nodes);
        global_offsets.resize(n_nodes, 0);
        
        pair<PackedStrings, sdsl::bit_vector> unitig_data = permute_unitigs(*(this->sbwt), reader);
        PackedStrings& unitigs = unitig_data.first;
        sdsl::bit_vector& Ustart = unitig_data.second;

        set<tuple<int64_t, int64_t, int64_t>>  finimizers;
        int64_t total_len = 0;
        vector<char> unitig_buf;
        for(int64_t i = 0; i < unitigs.number_of_strings(); i++){
            int64_t len = unitigs.get(i, unitig_buf);
            set<tuple<int64_t, int64_t, int64_t>> new_search = add_sequence(unitig_buf.data(), fmin_bv, fmin_found, global_offsets, total_len);
            total_len += len; 
            finimizers.insert(new_search.begin(), new_search.end());
        }

        sdsl::int_vector<> packed_global_offsets(finimizers.size(), 0, 64 - __builtin_clzll(*std::max_element(global_offsets.begin(), global_offsets.end())));
        
        int64_t global_offsets_idx = 0;
        for(int64_t i = 0; i < global_offsets.size(); i++){
            if(fmin_bv[i]) packed_global_offsets[global_offsets_idx++] = global_offsets[i];
        }
    
        print_finimizer_stats(finimizers, this->sbwt->number_of_kmers(), this->sbwt->number_of_subsets(), 1);

        index->sbwt = std::move(this->sbwt); // Transfer ownership
        index->LCS = std::move(this->LCS); // Transfer ownership 
        index->unitigs = std::move(unitigs); // Transfer ownership
        index->fmin = std::move(fmin_bv); // Transfer ownership
        index->fmin_rs = sdsl::rank_support_v5<>(&(index->fmin));
        index->global_offsets = std::move(packed_global_offsets); // Transfer ownership
        index->Ustart = std::move(Ustart); // Transfer ownership
        index->Ustart_rs = sdsl::rank_support_v5<>(&(index->Ustart)); 

    }

    set<tuple<int64_t, int64_t, int64_t>> add_sequence(const std::string& seq, sdsl::bit_vector& fmin_bv, sdsl::int_vector<>& fmin_found, vector<uint64_t>& global_offsets, const int64_t unitig_start) {
        const int64_t n_nodes = sbwt->number_of_subsets();
        const int64_t k = sbwt->get_k();
        const vector<int64_t>& C = sbwt->get_C_array();
        int64_t freq;
        BoundedDeque<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin(seq.size());
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
                if (w_fmin > curr_substr) {
                    all_fmin.clear();
                    w_fmin = curr_substr;
                } else{
                    while (all_fmin.back() > curr_substr) {all_fmin.pop_back();}
                }
                all_fmin.push_back(curr_substr);
            }
            if (end >= k -1 ){
                count_all_w_fmin.insert({get<1>(w_fmin),get<0>(w_fmin), get<2>(w_fmin) });// (length,freq,colex) freq = 1 thus == (freq, length,colex)

                if (fmin_found[get<2>(w_fmin)] == 0 or fmin_found[get<2>(w_fmin)]< get<3>(w_fmin)){ // if the finimizer has been found before in a full kmer
                    fmin_bv[get<2>(w_fmin)]=1;
                    fmin_found[get<2>(w_fmin)] = get<3>(w_fmin);

                    if ((unitig_start + get<3>(w_fmin))> UINT64_MAX){
                        std::cerr<< "ISSUE: global offset exceedes the allowed bit range." << std::endl;
                    }
                    global_offsets[get<2>(w_fmin)]= unitig_start + get<3>(w_fmin);
                }
                // write_fasta({input.substr(kmer,k) + ' ' + to_string(get<0>(w_fmin)),input.substr(get<3>(w_fmin)-get<1>(w_fmin)+1,get<1>(w_fmin))},writer);
                kmer++;
                // Check if the current minimizer is still in this window
                while (get<3>(w_fmin)- get<1>(w_fmin)+1 < kmer) { // start
                    all_fmin.pop_front();
                    w_fmin = (all_fmin.size()==0) ? tuple<int64_t, int64_t, int64_t, int64_t> {n_nodes,k+1,kmer+1,kmer+k} : all_fmin.front();
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