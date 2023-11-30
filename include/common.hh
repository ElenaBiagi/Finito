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

//TODO remove
// Inclusive ends. Retuns (end, colex of end)
// This function assumes that the k-mer we are looking for exists in the sbwt
optional<pair<int64_t, int64_t>> get_rightmost_Ustart_old(const std::string& query, int64_t kmer_end, int64_t finimizer_end, const vector<optional<int64_t>>& finimizer_end_colex, const plain_matrix_sbwt_t& sbwt, const sdsl::bit_vector& Ustart){

    if(!finimizer_end_colex[finimizer_end].has_value()){
        throw std::runtime_error("BUG: get_rightmost_branch_end");
    }

    int64_t colex = finimizer_end_colex[finimizer_end].value();
    optional<pair<int64_t, int64_t>> best = nullopt;
    int64_t p;
    for(p = finimizer_end; p < kmer_end; p++){
        // We're not checking the last position because we would extend
        // after the k-mer end.
        if (Ustart[colex]){
            best = {p, colex};  
        } 
        colex = sbwt.forward(colex, query[p+1]);
    }
    // Check also the last position without extending further
    if (Ustart[colex]){
            best = {p, colex};
            
        }
    return best;
}

//TODO remove
// Inclusive ends. Retuns (end, colex of end)
// This function assumes that the k-mer we are looking for exists in the sbwt
optional<pair<int64_t, int64_t>> get_rightmost_Ustart(const std::string& query, int64_t kmer_end, int64_t finimizer_end, const vector<optional<pair<int64_t, int64_t>>>& finimizer_end_colex, const plain_matrix_sbwt_t& sbwt, const sdsl::bit_vector& Ustart){

    if(!finimizer_end_colex[kmer_end].has_value()){
        throw std::runtime_error("BUG: get_rightmost_branch_end");
    }

    int64_t colex = finimizer_end_colex[kmer_end].value().second;
    optional<pair<int64_t, int64_t>> best = nullopt;
    int64_t p;
    for(p = finimizer_end; p < kmer_end; p++){
        // We're not checking the last position here because we would extend
        // after the k-mer end.
        if (Ustart[colex]){
            best = {p, colex};  
        } 
        colex = sbwt.forward(colex, query[p+1]);
    }
    // Check also the last position without extending further
    if (Ustart[colex]){
            best = {p, colex};
            
        }
    return best;
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


// If a kmer exists, it returns:
// kmers colex rank
// finimizer end, finimizer colex rank
tuple<vector<optional<int64_t>>,vector<optional< pair<int64_t, int64_t> > >, vector<optional< pair<int64_t, int64_t> > >> rarest_fmin_streaming_search(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input, const sdsl::bit_vector& Ustart){ 
    const int64_t n_nodes = sbwt.number_of_subsets();
    const int64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();
    deque<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin;
    const int64_t str_len = input.size();
    tuple<int64_t, int64_t, int64_t, int64_t> w_fmin = {n_nodes,k+1,n_nodes,str_len+1}; // {freq, len, I start, end}

    //vector<int64_t> found_kmers(str_len, -1);
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
        int64_t char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]]{
            cerr << "Error: unknown character: " << c << endl;
            cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
            return {};
        } else {
            // 1) fmin interval
            I_new = sbwt.update_sbwt_interval(&c, 1, I); //I_new = update_sbwt_interval(C[char_idx], I, Bit_rs);
            // (1) Finimizer(subseq) NOT found
            // TODO We already know that no kmer will be found
            while(I_new.first == -1){
                kmer_start = ++start;
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
            //cout << I_start << endl;
            //for (auto x:Ustart){ cout << x << ",";}
            //cout << endl;
            if (Ustart[I_kmer.first]){
                best_Ustart = {end, I_kmer.first};  
            }

            // Check if the kmer is found
            if (end - kmer_start + 1 == k){
                count++;
                while ((get<3>(w_fmin)-get<1>(w_fmin) +1) < kmer_start) {
                    all_fmin.pop_front();
                    w_fmin = all_fmin.front();
                }
                colex_ranks[kmer_start+k-1] = optional<int64_t>(I_kmer.first);
                finimizers[kmer_start+k-1] = optional<pair<int64_t, int64_t>>({get<3>(w_fmin), get<2>(w_fmin)});
                if (best_Ustart.first >= get<3>(w_fmin)){
                    best[kmer_start+k-1] = optional<pair<int64_t, int64_t>>(best_Ustart);
                }
                kmer_start++;
                I_kmer = drop_first_char(end - kmer_start + 1, I_kmer, LCS, n_nodes);
            }
        }
    }
    return tie(colex_ranks, finimizers, best);
}

// TODO remove
vector<int64_t> kmer_LCS_streaming_search(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input){
    const uint64_t n_nodes = sbwt.number_of_subsets();
    const uint64_t k = sbwt.get_k();
    const uint64_t str_len = input.size();
    vector<int64_t> colex(str_len -k +1, -1);
    uint64_t start = 0;
    uint64_t end;
    pair<int64_t, int64_t> I = {0, n_nodes - 1};
    pair<int64_t, int64_t> I_new;

    for (end = 0; end < str_len; end++) {

        char c = static_cast<char>(input[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        int64_t char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]]{
            cerr << "Error: unknown character: " << c << endl;
            cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
            return {};
        } else {
            I_new = sbwt.update_sbwt_interval(&c, 1, I);
            // (1) kmer (or subseq) NOT found
            // We already know that no kmer will be found thus we update start
            while(I_new.first == -1){
                if (start == end){
                    start++;
                    I_new = {0, n_nodes - 1};
                    break;
                }
                start++;
                I = drop_first_char(end - start, I, LCS, n_nodes); // The result (substr(start++,end)) cannot have freq == 1 as substring(start,end) has freq >1
                I_new = sbwt.update_sbwt_interval(&c, 1, I);
                
            }
            I = I_new;
        
            // Check if the kmer is found
            if (end - start + 1 == k){
                colex[start]= I.first;
                start++;
                I = drop_first_char(end - start + 1, I, LCS, n_nodes);
            }   
        }
    }
    return colex;
}

// TODO remove
vector<optional<int64_t>> get_kmer_colex_ranks(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& query){
    //vector<int64_t> colex_ranks = sbwt.streaming_search(query);
    vector<int64_t> colex_ranks = kmer_LCS_streaming_search(sbwt, LCS, query);
    vector<optional<int64_t>> answers(sbwt.get_k()-1,optional<int64_t>());

    for(int64_t x : colex_ranks){
        answers.push_back(x == -1 ? optional<int64_t>() : optional<int64_t>(x));
    }
    return answers;
}

// TODO remove
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

    vector<optional<int64_t>> shortest_unique_lengths(query.size());
    vector<optional<int64_t>> shortest_unique_colex_ranks(query.size());
    
    for (int64_t end = 0; end < str_len; end++) {
        char c = static_cast<char>(query[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        char char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]]{
            cerr << "Error: unknown character: " << c << endl;
            cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
            return {};
        } else {
            I_new = sbwt.update_sbwt_interval(&c, 1, I);

            // (1) Finimizer(subseq) NOT found
            while(I_new.first == -1){
                shortest_unique_lengths[end]= optional<int64_t>{};
                shortest_unique_colex_ranks[end]= optional<int64_t>{};
                if (start == end){
                    start++;
                    I_new = {0, n_nodes - 1};
                    break;
                }
                start++;
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

// TODO remove
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
