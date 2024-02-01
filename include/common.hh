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
#include <deque>

#include "sbwt/throwing_streams.hh"
#include "PackedStrings.hh"
#include "SeqIO.hh"
#include "BoundedDeque.hh"


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
    if (new_len<=0){return {0, n_nodes - 1};}
    pair<int64_t,int64_t> new_I = I;
    //Check top and bottom w the LCS
    while (new_I.first > 0 && LCS[new_I.first] >= new_len ){new_I.first --;}
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

// Returns the end point (inclusice) of the first k-mer in the concatenation of the unitigs
int64_t lookup_from_branch_dictionary(int64_t kmer_colex, int64_t k, const sdsl::rank_support_v5<>& Ustart_rs, const PackedStrings& unitigs){
    int64_t unitig_rank = Ustart_rs.rank(kmer_colex);
    assert(unitig_rank < unitigs.ends.size());
    int64_t global_unitig_start = (unitig_rank > 0) ? unitigs.ends[unitig_rank-1] :  0;    
    return global_unitig_start + k - 1;
}

int64_t lookup_from_finimizer_dictionary(int64_t finimizer_colex, const sdsl::rank_support_v5<>& fmin_rs, const sdsl::int_vector<>& global_offsets){
    int64_t finimizer_id = fmin_rs.rank(finimizer_colex);
    return global_offsets[finimizer_id];
}


// If a kmer exists, it returns:
// kmers colex rank
// finimizer end, finimizer colex rank
tuple<vector<optional<int64_t>>,vector<optional< pair<int64_t, int64_t> > >, vector<optional< pair<int64_t, int64_t> > >> rarest_fmin_streaming_search(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input, const sdsl::bit_vector& Ustart){ 
    const int64_t n_nodes = sbwt.number_of_subsets();
    const int64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();
    BoundedDeque<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin(input.size());
    const int64_t str_len = input.size();
    tuple<int64_t, int64_t, int64_t, int64_t> w_fmin = {n_nodes,k+1,n_nodes,str_len+1}; // {freq, len, I start, end}
    vector<optional<int64_t>> colex_ranks(str_len, optional<int64_t>());
    vector<optional< pair<int64_t, int64_t>>> finimizers(str_len, optional<pair<int64_t, int64_t>>());

    int64_t freq;
    int64_t count = 0;
    int64_t start = 0;
    int64_t end;
    int64_t kmer_start = 0;
    pair<int64_t, int64_t> I = {0, n_nodes - 1}, I_kmer = {0, n_nodes - 1};
    pair<int64_t, int64_t> I_new, I_kmer_new;
    int64_t I_start;
    tuple<int64_t, int64_t, int64_t, int64_t> curr_substr;
    pair<int64_t, int64_t> best_Ustart = {-1,-1};

    vector<optional<pair<int64_t, int64_t>>> best(str_len, optional<pair<int64_t, int64_t>>());
    
    // the idea is to start from the first pos which is i and move until finding something of ok freq
    // then drop the first char keeping track of which char you are starting from
    // Start is always < k as start <= end and end <k
    // if start == end than the frequency higher than t
    for (end = 0; end < str_len; end++) {
        char c = static_cast<char>(input[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        // 1) fmin interval
        I_new = sbwt.update_sbwt_interval(&c, 1, I); //I_new = update_sbwt_interval(C[char_idx], I, Bit_rs);
        // (1) Finimizer(subseq) NOT found
        // TODO We already know that no kmer will be found
            while(I_new.first == -1){
                kmer_start = ++start;
                if (start>end)[[unlikely]]{
                    I_new = {0, n_nodes - 1};
                    I_kmer = I_new;
                    break;
                }
                I = drop_first_char(end - start, I, LCS, n_nodes); // The result (substr(start++,end)) cannot have freq == 1 as substring(start,end) has freq >1
                I_new = sbwt.update_sbwt_interval(&c, 1, I);// I_new = update_sbwt_interval(C[char_idx], I, Bit_rs);
                I_kmer = I_new;
            }
            I = I_new;
            freq = (I.second - I.first + 1);
            I_start = I.first;
            // (2) Finimizer(subseq) freq > 0
            // Check if the Kmer interval has to be updated
            if ( start != kmer_start){
                I_kmer_new = sbwt.update_sbwt_interval(&c, 1, I_kmer); //I_kmer_new = update_sbwt_interval(C[char_idx], I_kmer, Bit_rs);
                while(I_kmer_new.first == -1){
                    // kmer NOT found
                    kmer_start++;
                    I_kmer = drop_first_char(end - kmer_start, I_kmer, LCS, n_nodes);
                    I_kmer_new = sbwt.update_sbwt_interval(&c, 1, I_kmer);//I_kmer_new = update_sbwt_interval(C[char_idx], I_kmer, Bit_rs);
                } 
                I_kmer = I_kmer_new;
            } else { 
                I_kmer = I;
            }
            // (2b) Finimizer found
            if (freq ==1){ // 1. rarest
                while (freq == 1) { // 2. shortest
                    I_start = I.first;
                    curr_substr = {freq, end - start + 1, I_start, end};
                    // 2. drop the first char
                    // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                    start ++;
                    I = drop_first_char(end - start + 1, I, LCS, n_nodes);
                    freq = (I.second - I.first + 1);
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

            // Ustart
            if (I_kmer.first == I_kmer.second && Ustart[I_kmer.first]==1){ best_Ustart = {end, I_kmer.first};}

            // Check if the kmer is found
            if (end - kmer_start + 1 == k){
            
                count++;
                while ((get<3>(w_fmin)-get<1>(w_fmin) +1) < kmer_start) {
                    all_fmin.pop_front();
                    w_fmin = all_fmin.front();
                }
                colex_ranks[kmer_start+k-1] = optional<int64_t>(I_kmer.first);
                finimizers[kmer_start+k-1] = optional<pair<int64_t, int64_t>>({get<3>(w_fmin), get<2>(w_fmin)});// end , I-start(colex)
                if (best_Ustart.first >= get<3>(w_fmin)){ best[kmer_start+k-1] = optional<pair<int64_t, int64_t>>(best_Ustart);}
                kmer_start++;
                I_kmer = drop_first_char(end - kmer_start + 1, I_kmer, LCS, n_nodes);
            }
        //}
    }
    return tie(colex_ranks, finimizers, best);
}

void print_finimizer_stats(const set<tuple<int64_t, int64_t, int64_t>>& finimizers, int64_t n_kmers, int64_t n_nodes, int64_t t){
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
