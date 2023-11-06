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

char get_char_idx(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

bool is_branching(const plain_matrix_sbwt_t& sbwt, int64_t colex){
    // Rewind to the start of the suffix group
    while(colex > 0 && sbwt.get_streaming_support()[colex] == 0)
        colex--;

    return sbwt.get_subset_rank_structure().A_bits[colex]
         + sbwt.get_subset_rank_structure().C_bits[colex]
         + sbwt.get_subset_rank_structure().G_bits[colex]
         + sbwt.get_subset_rank_structure().T_bits[colex]
         > 1;
}

// Inclusive ends. Retuns (end, colex of end)
// This function assumes that the k-mer we are looking for exists in the sbwt
optional<pair<int64_t, int64_t>> get_rightmost_branch_end(const std::string& query, int64_t kmer_end, int64_t k, int64_t finimizer_end, const vector<optional<int64_t>>& finimizer_end_colex, const plain_matrix_sbwt_t& sbwt){

    if(!finimizer_end_colex[finimizer_end].has_value()){
        throw std::runtime_error("BUG: get_rightmost_branch_end");
    }

    int64_t colex = finimizer_end_colex[finimizer_end].value();
    optional<pair<int64_t, int64_t>> best = nullopt;
    for(int64_t p = finimizer_end; p < kmer_end - 1; p++){
        if(is_branching(sbwt, colex)){
            colex = sbwt.forward(colex, query[p+1]);
            best = {p+1, colex};
        } else{
            colex = sbwt.forward(colex, query[p+1]);
        }
    }
    return best;

    //for(int64_t p = kmer_end - 1; p >= finimizer_end; p--){
    //    if(finimizer_end_colex[p].has_value() && is_branching(sbwt, finimizer_end_colex[p].value())){
    //        char c = query[p+1];
    //        int64_t colex = sbwt.forward(finimizer_end_colex[p].value(), c); // there should be a branch with input[p+1] because if we are in this function, the k-mer we are looking for should exist
    //        if(colex == -1) throw std::runtime_error("BUG: could not find edge with label " + c);
    //        return optional<pair<int64_t, int64_t>>({p+1, colex}); // First k-mer end after branch
    //    }
    //}
    //return nullopt;
}

// Returns the end point (inclusice) of the first k-mer in the concatenation of the unitigs
int64_t lookup_from_branch_dictionary(int64_t kmer_colex, int64_t k, const sdsl::rank_support_v5<>& Ustart_rs, const PackedStrings& unitigs){
    int64_t unitig_rank = Ustart_rs.rank(kmer_colex);
    assert(unitig_rank < unitigs.ends.size());
    int64_t global_unitig_start = 0;
    if(unitig_rank > 0) global_unitig_start = unitigs.ends[unitig_rank-1];
    return global_unitig_start + k - 1;
}

int64_t lookup_from_finimizer_dictionary(int64_t finimizer_colex, const sdsl::rank_support_v5<>& fmin_rs, const sdsl::int_vector<>& global_offsets){
    int64_t finimizer_id = fmin_rs.rank(finimizer_colex);
    return global_offsets[finimizer_id];
}

vector<optional<int64_t>> get_kmer_colex_ranks(const plain_matrix_sbwt_t& sbwt, const string& query){
    vector<int64_t> colex_ranks = sbwt.streaming_search(query);
    vector<optional<int64_t>> answers;
    for(int64_t i = 0; i < sbwt.get_k()-1; i++) // First k-1 are always null because the k-mer is not full
        answers.push_back(optional<int64_t>());

    for(int64_t x : colex_ranks)
        answers.push_back(x == -1 ? optional<int64_t>() : optional<int64_t>(x));
    
    return answers;
}

// Returns for each endpoint in the query the length of the shortest unique match (if exists)
// ending there, and the colex of that match
pair<vector<optional<int64_t>>, vector<optional<int64_t>>> get_shortest_unique_lengths_and_colex_ranks(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& query){
    const int64_t n_nodes = sbwt.number_of_subsets();
    const vector<int64_t>& C = sbwt.get_C_array();
    int64_t freq;
    const int64_t str_len = query.size();
    
    int64_t start = 0;
    int64_t kmer_start = 0;
    pair<int64_t, int64_t> I = {0, n_nodes - 1}, I_kmer = {0, n_nodes - 1};
    pair<int64_t, int64_t> I_new;
    int64_t I_start;
    pair<int64_t, int64_t> curr_substr;

    auto is_singleton = [](const pair<int64_t, int64_t>& interval){ // Helper to make code more readable
        return interval.second == interval.first;
    };

    vector<optional<int64_t>> shortest_unique_lengths(query.size());
    vector<optional<int64_t>> shortest_unique_colex_ranks(query.size());
    
    for (int64_t end = 0; end < str_len; end++) {
        //cout << "I: " << I.first << " " << I.second << endl;
        char c = static_cast<char>(query[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        char char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]]{
            cerr << "Error: unknown character: " << c << endl;
            cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
            return {};
        } else {
            I_new = sbwt.update_sbwt_interval(&c, 1, I);
            //cout << "I_new: " << I_new.first << " " << I_new.second << endl;


            // (1) Finimizer(subseq) NOT found
            while(I_new.first == -1){
                shortest_unique_lengths[end]= optional<int64_t>{};
                shortest_unique_colex_ranks[end]= optional<int64_t>{};
                I = drop_first_char(end - start, I, LCS, n_nodes); // The result (substr(start++,end)) cannot have freq == 1 as substring(start,end) has freq >1
                I_new = sbwt.update_sbwt_interval(&c, 1, I);
            }
            I = I_new;
            freq = (I.second - I.first + 1);
            I_start = I.first;
            // (2) Finimizer(subseq) freq > 0
          
            // (2b) Finimizer found
            if (freq ==1){ // 1. rarest
                while (freq == 1) { // 2. shortest
                    curr_substr = {end - start + 1, I_start};
                    // 2. drop the first char
                    // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                    start ++;
                    I = drop_first_char(end - start + 1, I, LCS, n_nodes);
                    freq = (I.second - I.first + 1);
                    I_start = I.first;
                }

                
                shortest_unique_lengths[end]= curr_substr.first;
                shortest_unique_colex_ranks[end]= curr_substr.second;
            }else{
                shortest_unique_lengths[end]= optional<int64_t>{};
                shortest_unique_colex_ranks[end]= optional<int64_t>{};
            }
        }
    }
    //for(auto x : shortest_unique_lengths) cout << (x.has_value() ? to_string(x.value()) : "Null") << " "; cout << endl;

    return {shortest_unique_lengths, shortest_unique_colex_ranks};
}

int64_t pick_finimizer(int64_t kmer_end, int64_t k, const vector<optional<int64_t>>& shortest_unique_lengths, const vector<optional<int64_t>>& shortest_unique_colex_ranks){
    if(kmer_end < k-1) throw std::runtime_error("bug: end < k-1");

    const int64_t INFINITE = 1e18;
    int64_t best_end = INFINITE;
    int64_t best_length = INFINITE;
    int64_t colex_tiebreaker = INFINITE;

    int64_t kmer_start = kmer_end - k + 1;

    for(int64_t f_end = kmer_start; f_end <= kmer_end; f_end++){
        auto x = shortest_unique_lengths[f_end].has_value() ? to_string(shortest_unique_lengths[f_end].value()) : "Null";
        if(shortest_unique_lengths[f_end].has_value()){
            int64_t len = shortest_unique_lengths[f_end].value();
            int64_t colex = shortest_unique_colex_ranks[f_end].value();
            if(f_end - len + 1 >= kmer_start){
                // This finimizer is in window
                if(len < best_length || (len == best_length && colex < colex_tiebreaker)){
                    best_end = f_end;
                    best_length = len;
                    colex_tiebreaker = colex;
                }
            }
        }
    }

    if(best_length == INFINITE){
        cerr << "Error: no finimizer found of a k-mer" << endl;
    }

    return best_end;
}