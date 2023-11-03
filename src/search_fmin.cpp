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

size_t size_in_bytes(const sdsl::int_vector<>& LCS, const sdsl::bit_vector& fmin_bv, const sdsl::rank_support_v5<>& fmin_rs, const sdsl::int_vector<>& unitigs_v,  const sdsl::sd_vector<>& ef_endpoints, const plain_matrix_sbwt_t& sbwt){
            size_t sz = 0;
            // SBWT
            const sdsl::bit_vector& A_bits = sbwt.get_subset_rank_structure().A_bits;                
            const sdsl::bit_vector& C_bits = sbwt.get_subset_rank_structure().C_bits;
            const sdsl::bit_vector& G_bits = sbwt.get_subset_rank_structure().G_bits;
            const sdsl::bit_vector& T_bits = sbwt.get_subset_rank_structure().T_bits;

            const sdsl::rank_support_v5<> &A_bits_rs = sbwt.get_subset_rank_structure().A_bits_rs;
            const sdsl::rank_support_v5<> &C_bits_rs = sbwt.get_subset_rank_structure().C_bits_rs;
            const sdsl::rank_support_v5<> &G_bits_rs = sbwt.get_subset_rank_structure().G_bits_rs;
            const sdsl::rank_support_v5<> &T_bits_rs = sbwt.get_subset_rank_structure().T_bits_rs;
            sz += sdsl::size_in_bytes(A_bits);
            sz += sdsl::size_in_bytes(A_bits_rs);
            sz += sdsl::size_in_bytes(C_bits);
            sz += sdsl::size_in_bytes(C_bits_rs);
            sz += sdsl::size_in_bytes(G_bits);
            sz += sdsl::size_in_bytes(G_bits_rs);
            sz += sdsl::size_in_bytes(T_bits);
            sz += sdsl::size_in_bytes(T_bits_rs);
            cerr << "SBWT size = " << to_string(sz) << endl;

            // LCS
            sz += sdsl::size_in_bytes(LCS);
            cerr << "LCS size = " << to_string(sdsl::size_in_bytes(LCS)) << endl;
            // marksRound 3/3

            // ids
            sz += sdsl::size_in_bytes(unitigs_v);
            cerr << "offsets size = " << to_string(sdsl::size_in_bytes(unitigs_v)) << endl;

            // endpoints
            sz += sdsl::size_in_bytes(ef_endpoints);
            
            cerr << "endpoints size = " << to_string(sdsl::size_in_bytes(ef_endpoints)) << endl;

            cerr << "Total size in bytes = " << to_string(sz) << endl;


            return sz;
        }

// Assumes values of v are -1 or larger
template <typename writer_t>
inline void print_vector(const vector<int64_t>& v, writer_t& out){
    // Fast manual integer-to-string conversion
    char buffer[32];
    char newline = '\n';
    for(int64_t x : v){
        int64_t i = 0;
        if(x == -1){
            buffer[0] = '1';
            buffer[1] = '-';
            i = 2;
        } else{
            while(x > 0){
                buffer[i++] = '0' + (x % 10);
                x /= 10;
            }
        }
        std::reverse(buffer, buffer + i);
        buffer[i] = ' ';
        out.write(buffer, i+1);
    }
    out.write(&newline, 1);
}


void shortest_unique_search_jarno_rewrite(const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input, const char t, const sdsl::rank_support_v5<>& fmin_rs, const sdsl::int_vector<>& global_offsets, const PackedStrings& unitigs, const sdsl::rank_support_v5<>& Ustart_rs){ // writer_t& writer

    // For each k-mer S that is known to be in the SBWT
    //   - Find the finimizer x.
    //   - Walk forward in the SBWT from the colex rank of x (singleton interval), to the end of S,
    //     recording the rightmost branch point, and the distance from the end of x to this branch
    //     point.
    //   - If a branch point was found, the k-mer containing x is at the unitig after the rightmost
    //     branch. Otherwise, the k-mer containing x is at the unitig that x points to.
    //   - If there was a branch, the k- mer endpoint in that unitig is k + the number of steps taken after the last branch.
    //     

    const int64_t n_nodes = sbwt.number_of_subsets();
    const int64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();

    if(input.size() < k) return;

    // For each position of the input, the SBWT colex rank of the k-mer that ends there, if exists.
    // A k-mer does not exist if the ending position is too close to the start of the input, or the SBWT does
    // not contain that k-mer.
    vector<optional<int64_t>> kmer_colex_ranks = get_kmer_colex_ranks(sbwt, input);

    // For each position of the input, the length of the shortest unique substring that ends there, if exists
    // Unique substrings may not exist near the start of the input.
    vector<optional<int64_t>> shortest_unique_lengths;

    // For each position of the input, the colex rank ofthe shortest unique substring that ends there, if exists
    vector<optional<int64_t>> shortest_unique_colex_ranks;

    std::tie(shortest_unique_lengths, shortest_unique_colex_ranks) = get_shortest_unique_lengths_and_colex_ranks(DNA_rs, sbwt, LCS, input);

    for(int64_t kmer_end = k-1; kmer_end < input.size(); kmer_end++) {
        if(kmer_colex_ranks[kmer_end].has_value()){
            // kmer exists

            int64_t finimizer_end = pick_finimizer(kmer_end, k, shortest_unique_lengths, shortest_unique_colex_ranks);
            optional<pair<int64_t, int64_t>> rightmost_branch_end = get_rightmost_branch_end(input, kmer_end, k, finimizer_end, shortest_unique_colex_ranks, sbwt);
            if(rightmost_branch_end.has_value()) {
                // Look up from the branch dictionary
                int64_t p = rightmost_branch_end.value().first;
                int64_t colex = rightmost_branch_end.value().second; 

                // Get the global off set of the end of the k-mer
                int64_t global_kmer_end = lookup_from_branch_dictionary(colex, k, Ustart_rs, unitigs);
                global_kmer_end += kmer_end - p + 1;
            } else {
                // Look up from Finimizer dictionary
                int64_t p = finimizer_end;
                int64_t colex = shortest_unique_colex_ranks[p].value();
                int64_t global_kmer_end = lookup_from_finimizer_dictionary(colex, fmin_rs, global_offsets);
                global_kmer_end += kmer_end - p + 1;
            }
        }        
    }

    // TODO: report results to the caller

    
}




// Here you are noT sure to find the interval as when building fmin
//template<typename writer_t>

pair<vector<int64_t>, int64_t> rarest_fmin_streaming_search(const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input, const char t, const sdsl::rank_support_v5<>& fmin_rs, const  sdsl::int_vector<>& unitigs_v, const sdsl::sd_vector<>& ef_endpoints, const sdsl::rank_support_v5<>& Ustart_rs, vector<int64_t>& found_kmers){ // writer_t& writer
    const int64_t n_nodes = sbwt.number_of_subsets();
    const int64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();
    int64_t freq;
    set<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin;
    const int64_t str_len = input.size();
    tuple<int64_t, int64_t, int64_t, int64_t> w_fmin = {n_nodes,k+1,n_nodes,str_len}; // {freq, len, I start, start}
    char c;
    char char_idx;
    int64_t count = 0;
    int64_t start = 0;
    int64_t end;
    int64_t kmer_start = 0;
    pair<int64_t, int64_t> I = {0, n_nodes - 1}, I_kmer = {0, n_nodes - 1};
    pair<int64_t, int64_t> I_new, I_kmer_new;
    int64_t I_start = -1;
    tuple<int64_t, int64_t, int64_t, int64_t> curr_substr;
    set<pair<int64_t, int64_t>> last_branch; // {index in the query, I_start,}
    int64_t unitig_id;
    int64_t unitig_start;
    
    // the idea is to start from the first pos which is i and move until finding something of ok freq
    // then drop the first char keeping track of which char you are starting from
    // Start is always < k as start <= end and end <k
    // if start == end than the frequency higher than t
    for (end = 0; end < str_len; end++) {
        c = static_cast<char>(input[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]]{
            cerr << "Error: unknown character: " << c << endl;
            cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
            return {};
        } else {
            // WRONG WE NEVER EXTEND RIGHT IF THE FREQ IS 1
            /* // Check the previous branch before searching the next char
            if (freq==1){
                char branch = 0;
                vector<char> bases = {0,1,2,3};
                bases.erase(bases.begin()+char_idx);
                // TODO: remove char_idx
                for (char base :bases){
                const sdsl::bit_vector &Bit_char = *(DNA_bitvectors[base]);
                branch= branch | Bit_char[I_start];
            }
                    if (branch){
                        last_branch = {I_start, end};   
                    }  
            } */

            const sdsl::rank_support_v5<> &Bit_rs = *(DNA_rs[char_idx]);
            // 1) fmin interval
            I_new = update_sbwt_interval(C[char_idx], I, Bit_rs);
            // (1) Finimizer(subseq) NOT found
            // TODO We already know that no kmer will be found
            while(I_new.first == -1){
                kmer_start = ++start;
                I = drop_first_char(end - start, I, LCS, n_nodes); // The result (substr(start++,end)) cannot have freq == 1 as substring(start,end) has freq >1
                I_new = update_sbwt_interval(C[char_idx], I, Bit_rs);
                I_kmer = I_new;
            }
            I = I_new;
            freq = (I.second - I.first + 1);
            I_start = I.first;
            // (2) Finimizer(subseq) freq > 0
            // Check if the Kmer interval has to be updated
            if ( start != kmer_start){
                I_kmer_new = update_sbwt_interval(C[char_idx], I_kmer, Bit_rs);
                while(I_kmer_new.first == -1){
                    // kmer NOT found
                    kmer_start++;
                    I_kmer = drop_first_char(end - kmer_start, I_kmer, LCS, n_nodes);
                    I_kmer_new = update_sbwt_interval(C[char_idx], I_kmer, Bit_rs);
                } 
                I_kmer = I_kmer_new;
            } else { 
                I_kmer = I;
            }
            if (I_kmer.first == I_kmer.second){
                // Check the previous branch before searching the next char
                char branch = 0;
                vector<char> bases = {0,1,2,3}; // still checking all 4 possible char here
                // TODO: remove char_idx without going out of bound!!!
                // c = static_cast<char>(input[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
                // char_idx = get_char_idx(c);
                // bases.erase(bases.begin()+char_idx);
                for (char base :bases){
                    const sdsl::bit_vector &Bit_char = *(DNA_bitvectors[base]);
                    branch = branch + Bit_char[I_start];
                    if (branch > 1) {
                        last_branch.insert({end+1,I_start}); // insert a new branch 
                        break;}
                }
            }
            // (2b) Finimizer found
            if (freq ==1){ // 1. rarest
                while (freq == 1) { // 2. shortest
                    curr_substr = {freq, end - start + 1, I_start, end - k + 1};
                    // 2. drop the first char
                    // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                    start ++;
                    I = drop_first_char(end - start + 1, I, LCS, n_nodes);
                    freq = (I.second - I.first + 1);
                    I_start = I.first;
                }
                if (w_fmin > curr_substr) {w_fmin = curr_substr;}
                all_fmin.insert(curr_substr);
            }
            // Check if the kmer is found
            if (end - kmer_start + 1 == k){
                count++;
                while ((get<3>(w_fmin)+k-get<1>(w_fmin)) < kmer_start) {
                    all_fmin.erase(all_fmin.begin());
                    w_fmin = *all_fmin.begin();
                }
                // only now I know the correct finimizer and the kmer
                if (!last_branch.empty()){
                    // identify the first branch after the finimizer
                    pair<int64_t, int64_t> f_branch = *last_branch.begin();
                    while( f_branch.first <= (get<3>(w_fmin)+k) ){ // remove all branches before the finimizer end
                        last_branch.erase(last_branch.begin());
                        f_branch = *last_branch.begin();
                    }

                    if (!last_branch.empty()){
                        // Check when the branch occured
                        set<pair<int64_t, int64_t>>::reverse_iterator itr;
                        pair<int64_t, int64_t> l_branch;
                        for (itr = last_branch.rbegin(); itr != last_branch.rend(); itr++){
                            l_branch = *itr;
                            if (l_branch.first <= end){break;}
                        }
                        //pair<int64_t, int64_t> l_branch = last_branch.rbegin()[0]; // last element //*last_branch.end();
                        //if (l_branch.first > end){
                        //    pair<int64_t, int64_t> l_branch = last_branch.rbegin()[1];
                        //}
                        if (l_branch.first >= get<3>(w_fmin)+k){ // after the finimizer
                            // we are pointing to a different unitig!!
                            // we have the I_start and the pos of the char preceding the branch in the query
                            c = static_cast<char>(input[l_branch.first] & ~32); // branching char
                            char_idx = get_char_idx(c);
                            pair<int64_t, int64_t> I_branch = update_sbwt_interval(C[char_idx], {l_branch.second,l_branch.second}, Bit_rs);
                            unitig_id= Ustart_rs(unitig_start);
                            ef_endpoints[unitig_id];// correct finimizer
                            //I have the index in the previous finimizer and know where it ends
                            // get<3>(w_fmin) start of the kmer ending with the finimizer
                            // get<3>(w_fmin)+k-1 last pos of the finimizer (- get<1>(w_fmin) +1 = start pos)
                            // start pos of the kmer ending with the finimizer in the wrong unitig unitigs_v[fmin_rs(get<2>(w_fmin))]
                            int64_t f_start = k - (f_branch.first - (get<3>(w_fmin)+k-get<1>(w_fmin))); // starting pos of the finimizer in the last kmer of the WRONG unitig      //(get<3>(w_fmin)+k-get<1>(w_fmin)) - (l_branch.first - 1 - k); 
                            
                            char branches = l_branch.first - f_branch.first + 1; // all the char between the first and the last branch following unitig included!!! 
                            f_start += unitig_id - branches; // starting pos of the finimizer in the new unitig
                            found_kmers[kmer_start]= f_start - ((get<3>(w_fmin)+k-get<1>(w_fmin)) - kmer_start);

                        } else { 
                        // no relevant branches in the kmer after the finimizer 
                        // TODO check if the branch occured before right at the start of the finimizer!!! $$$$xxx(finimizer)
                        found_kmers[kmer_start] = unitigs_v[fmin_rs(get<2>(w_fmin))]+(get<3>(w_fmin) - kmer_start);
                        }
                    }else{
                        found_kmers[kmer_start] = unitigs_v[fmin_rs(get<2>(w_fmin))]+(get<3>(w_fmin) - kmer_start);
                    }
                    
                
                }else{
                    found_kmers[kmer_start] = unitigs_v[fmin_rs(get<2>(w_fmin))]+(get<3>(w_fmin) - kmer_start);

                }
                
                kmer_start++;
                I_kmer = drop_first_char(end - kmer_start + 1, I_kmer, LCS, n_nodes);
            }
        }
    }
    return {found_kmers, count};
}

pair<vector<int64_t>, int64_t> rarest_fmin_streaming_search_r( const sdsl::rank_support_v5<>** DNA_rs, const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input, const char t, const sdsl::rank_support_v5<>& fmin_rs, const  sdsl::int_vector<>& unitigs_v, const sdsl::sd_vector<>& ef_endpoints, vector<int64_t>& found_kmers){ //const sdsl::bit_vector** DNA_bitvectors, writer_t& writer
    char c;
    char char_idx;
    const int64_t n_nodes = sbwt.number_of_subsets();
    const int64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();
    int64_t freq;
    set<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin;
    const int64_t str_len = input.size();
    tuple<int64_t, int64_t, int64_t, int64_t> w_fmin = {n_nodes,k+1,n_nodes,str_len}; // {freq, len, I start, start}
    
    int64_t count = 0;
    int64_t start = 0;
    int64_t end;
    int64_t kmer_start = 0, rev_start;
    pair<int64_t, int64_t> I = {0, n_nodes - 1}, I_kmer = {0, n_nodes - 1};
    pair<int64_t, int64_t> I_new, I_kmer_new;
    int64_t I_start;
    tuple<int64_t, int64_t, int64_t, int64_t> curr_substr;
    
    // the idea is to start from the first pos which is i and move until finding something of ok freq
    // then drop the first char keeping track of which char you are starting from
    // Start is always < k as start <= end and end <k
    // if start == end than the frequency higher than t
    for (end = 0; end < str_len; end++) {
        c = static_cast<char>(input[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]]{
            cerr << "Error: unknown character: " << c << endl;
            cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
            return {};
        } else {
            const sdsl::rank_support_v5<> &Bit_rs = *(DNA_rs[char_idx]);
            // 1) fmin interval
            I_new = update_sbwt_interval(C[char_idx], I, Bit_rs);
            // (1) Finimizer(subseq) NOT found
            // TODO We already know that no kmer will be found
            while(I_new.first == -1){
                kmer_start = ++start;
                I = drop_first_char(end - start, I, LCS, n_nodes); // The result (substr(start++,end)) cannot have freq == 1 as substring(start,end) has freq >1
                I_new = update_sbwt_interval(C[char_idx], I, Bit_rs);
                I_kmer = I_new;
            }
            I = I_new;
            freq = (I.second - I.first + 1);
            I_start = I.first;
            // (2) Finimizer(subseq) freq > 0
            // Check if the Kmer interval has to be updated
            if ( start != kmer_start){
                I_kmer_new = update_sbwt_interval(C[char_idx], I_kmer, Bit_rs);
                while(I_kmer_new.first == -1){
                    // kmer NOT found
                    kmer_start++;
                    I_kmer = drop_first_char(end - kmer_start, I_kmer, LCS, n_nodes);
                    I_kmer_new = update_sbwt_interval(C[char_idx], I_kmer, Bit_rs);
                } 
                I_kmer = I_kmer_new;
            } else { 
                I_kmer = I;
            }
            // (2b) Finimizer found
            if (freq == 1){ // 1. rarest 
                while (freq == 1) { // 2. shortest
                    curr_substr = {freq, end - start + 1, I_start, end - k + 1};
                    // 2. drop the first char
                    // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                    start ++;
                    I = drop_first_char(end - start + 1, I, LCS, n_nodes);
                    freq = (I.second - I.first + 1);
                    I_start = I.first;
                }
                if (w_fmin > curr_substr) {w_fmin = curr_substr;}
                all_fmin.insert(curr_substr);
            }
            
            // Check if the kmer is found
            if (end - kmer_start + 1 == k){
                count++;
                while ((get<3>(w_fmin)+k-get<1>(w_fmin))< kmer_start) {
                    all_fmin.erase(all_fmin.begin());
                    w_fmin = *all_fmin.begin();
                }
                rev_start = str_len - (kmer_start + k);
                found_kmers[rev_start] = unitigs_v[fmin_rs(get<2>(w_fmin))]+(get<3>(w_fmin) - rev_start);
                kmer_start++;
                I_kmer = drop_first_char(end - kmer_start + 1, I_kmer, LCS, n_nodes);
            }
            
        }
        
    }
    return {found_kmers, count};
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_fmin_queries_streaming(reader_t& reader, writer_t& writer, const string& indexfile, const sbwt_t& sbwt, const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const sdsl::int_vector<>& LCS, const sdsl::rank_support_v5<>& fmin_rs, const  sdsl::int_vector<>& global_offsets, const PackedStrings& unitigs, const sdsl::rank_support_v5<>& Ustart_rs, const char t){
    const int64_t k = sbwt.get_k();
    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    int64_t kmers_count = 0 , count, count_rev, kmers_count_rev = 0;
    vector<int64_t> out_buffer, out_buffer_rev;
    //found_kmers.reserve(str_len - k + 1);
    //found_kmers.resize(str_len - k + 1); 

    int64_t query_seq=0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        int64_t t0 = cur_time_micros();
        vector<int64_t> found_kmers(len - k + 1,-1);
        //pair<vector<int64_t>, int64_t> final_pair = rarest_fmin_streaming_search(DNA_bitvectors, DNA_rs, sbwt, LCS, reader.read_buf, t, fmin_rs, global_offsets, ef_endpoints, Ustart_rs, found_kmers);
        shortest_unique_search_jarno_rewrite(DNA_bitvectors, DNA_rs, sbwt, LCS, reader.read_buf, 1, fmin_rs, global_offsets, unitigs, Ustart_rs);

        /* TODO Jarno
        out_buffer = final_pair.first;
        count = final_pair.second;
        number_of_queries += out_buffer.size();
        kmers_count += count;
    
        print_vector(found_kmers, writer);*/
     
       total_micros += cur_time_micros() - t0;
    }
    write_log("k " + to_string(k), LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    write_log("Found kmers: " + to_string(kmers_count), LogLevel::MAJOR);
    //write_log("Found kmers reverse : " + to_string(kmers_count_rev), LogLevel::MAJOR);
    write_log("Total found kmers: " + to_string(kmers_count+kmers_count_rev), LogLevel::MAJOR);

    std::ofstream statsfile;
    statsfile.open(indexfile + "stats.txt", std::ios_base::app); // append instead of overwrite
    statsfile << to_string(k) + "," + to_string(kmers_count+kmers_count_rev) + "," + to_string(number_of_queries);
    statsfile.close();
    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_not_streaming(reader_t& reader, writer_t& writer, const sbwt_t& sbwt){

    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    int64_t k = sbwt.get_k();
    vector<int64_t> out_buffer;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;

        for(int64_t i = 0; i < len - k + 1; i++){
            int64_t t0 = cur_time_micros();
            int64_t ans = sbwt.search(reader.read_buf + i);
            total_micros += cur_time_micros() - t0;
            number_of_queries++;
            out_buffer.push_back(ans);
        }

        print_vector(out_buffer, writer);
        out_buffer.clear();
    }
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_fmin_file(const string& infile, const string& outfile, const string& indexfile, const sbwt_t& sbwt, const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const sdsl::int_vector<>& LCS, const sdsl::rank_support_v5<>& fmin_rs, const sdsl::int_vector<>& unitigs_v, const PackedStrings& unitigs, const sdsl::rank_support_v5<>& Ustart_rs, const char t){
    reader_t reader(infile);
    writer_t writer(outfile);
    if(sbwt.has_streaming_query_support()){
        write_log("Running streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_fmin_queries_streaming<sbwt_t, reader_t, writer_t>(reader, writer, indexfile, sbwt, DNA_bitvectors, DNA_rs, LCS, fmin_rs, unitigs_v, unitigs, Ustart_rs, t);
    }
    else{
        write_log("Running non-streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_not_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
}

// Returns number of queries executed
template<typename sbwt_t>
int64_t run_fmin_queries(const vector<string>& infiles, const vector<string>& outfiles, const string& indexfile, const sbwt_t& sbwt, bool gzip_output, const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const sdsl::int_vector<>& LCS, const sdsl::rank_support_v5<>& fmin_rs, sdsl::int_vector<>& unitigs_v, const PackedStrings& unitigs, const sdsl::rank_support_v5<>& Ustart_rs, const char t){

    if(infiles.size() != outfiles.size()){
        string count1 = to_string(infiles.size());
        string count2 = to_string(outfiles.size());
        //throw std::runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
    }

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_no_gzip;

    typedef Buffered_ofstream<zstr::ofstream> out_gzip;
    typedef Buffered_ofstream<std::ofstream> out_no_gzip;

    int64_t n_queries_run = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = SeqIO::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input && gzip_output){
            n_queries_run += run_fmin_file<sbwt_t, in_gzip, out_gzip>(infiles[i], outfiles[i], indexfile, sbwt, DNA_bitvectors, DNA_rs, LCS, fmin_rs, unitigs_v, unitigs, Ustart_rs, t);
        }
        if(gzip_input && !gzip_output){
            n_queries_run += run_fmin_file<sbwt_t, in_gzip, out_no_gzip>(infiles[i], outfiles[i], indexfile, sbwt, DNA_bitvectors, DNA_rs, LCS, fmin_rs, unitigs_v, unitigs, Ustart_rs, t);
        }
        if(!gzip_input && gzip_output){
            n_queries_run += run_fmin_file<sbwt_t, in_no_gzip, out_gzip>(infiles[i], outfiles[i], indexfile, sbwt, DNA_bitvectors, DNA_rs, LCS, fmin_rs, unitigs_v, unitigs, Ustart_rs, t);
        }
        if(!gzip_input && !gzip_output){
            n_queries_run += run_fmin_file<sbwt_t, in_no_gzip, out_no_gzip>(infiles[i], outfiles[i], indexfile, sbwt, DNA_bitvectors, DNA_rs, LCS, fmin_rs, unitigs_v, unitigs, Ustart_rs, t);
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
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("i,index-file", "Index input file.This has to be a binary matrix.", cxxopts::value<string>())
        ("q,query-file", "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())
        ("z,gzip-output", "Writes output in gzipped form. This can shrink the output files by an order of magnitude.", cxxopts::value<bool>()->default_value("false"))
        ("t", "Maximum finimizer frequency", cxxopts::value<int64_t>())
        ("type", "Decide which streaming search type you prefer. Available types: " + all_types_string,cxxopts::value<string>()->default_value("rarest"))
        ("lcs", "Provide in input the LCS file if available.", cxxopts::value<string>()->default_value(""))
        ("f, fmin_bv", "Provide in input the finimizers binary kmers vector.", cxxopts::value<string>()->default_value(""))
        ("e, endpoints", "Provide in input the endpoints of the concatenated unitigs.", cxxopts::value<string>()->default_value(""))
        ("g, global-offsets", "Provide in input the global offsets of finimizers in the concatenated unitigs.", cxxopts::value<string>()->default_value(""))
        ("s, u-start", "Provide in input the bitvector marking the start kmer of each unitig.", cxxopts::value<string>()->default_value(""))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    char t = opts["t"].as<int64_t>();

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
    string outfile = opts["out-file"].as<string>();
    bool gzip_output = opts["gzip-output"].as<bool>();
    vector<string> output_files;
    if(multi_file){
        output_files = readlines(outfile);
    } else{
        output_files = {outfile};
    }
    for(string file : output_files) check_writable(file);


    // sbwt
    string indexfile = opts["index-file"].as<string>();
    check_readable(indexfile);
    throwing_ifstream in(indexfile, ios::binary);
    vector<string> variants = get_available_variants_fmin();
    string variant = load_string(in.stream); // read variant type
    if(std::find(variants.begin(), variants.end(), variant) == variants.end()){
        cerr << "Error loading index from file: unrecognized variant specified in the file" << endl;
        return 1;
    }

    write_log("Loading the index variant " + variant, LogLevel::MAJOR);
    int64_t number_of_queries = 0;

    if (variant == "plain-matrix"){
        plain_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        string type = opts["type"].as<string>();
        if(std::find(types.begin(), types.end(), type) == types.end()){
            cerr << "Error: unknown type: " << type << endl;
            cerr << "Available types are:" << all_types_string << endl;
            return 1;
        }

        const int64_t k = sbwt.get_k();
        cerr << "k = "<< to_string(k);
        cerr << " SBWT nodes: "<< to_string(sbwt.number_of_subsets())<< " kmers: "<< to_string(sbwt.number_of_kmers())<< endl;

        const sdsl::bit_vector& A_bits = sbwt.get_subset_rank_structure().A_bits;                
        const sdsl::bit_vector& C_bits = sbwt.get_subset_rank_structure().C_bits;
        const sdsl::bit_vector& G_bits = sbwt.get_subset_rank_structure().G_bits;
        const sdsl::bit_vector& T_bits = sbwt.get_subset_rank_structure().T_bits;
        const sdsl::bit_vector* DNA_bitvectors[4] = {&A_bits, &C_bits, &G_bits, &T_bits};


        const sdsl::rank_support_v5<> &A_bits_rs = sbwt.get_subset_rank_structure().A_bits_rs;
        const sdsl::rank_support_v5<> &C_bits_rs = sbwt.get_subset_rank_structure().C_bits_rs;
        const sdsl::rank_support_v5<> &G_bits_rs = sbwt.get_subset_rank_structure().G_bits_rs;
        const sdsl::rank_support_v5<> &T_bits_rs = sbwt.get_subset_rank_structure().T_bits_rs;
        const sdsl::rank_support_v5<> *DNA_rs[4] = {&A_bits_rs, &C_bits_rs, &G_bits_rs, &T_bits_rs};


        // LCS file
        string LCS_file = opts["lcs"].as<string>();
        if (LCS_file.empty()) {
            std::cout<< "LCS_file empty"<<std::endl;
            LCS_file = indexfile + "LCS.sdsl";
            const sdsl::int_vector<> LCS = lcs_basic_parallel_algorithm(sbwt, 8);
            //const sdsl::int_vector<> LCS = get_kmer_lcs(A_bits, C_bits, G_bits, T_bits, k);
           save_v(LCS_file, LCS);
        }
        sdsl::int_vector<> LCS;
        load_v(LCS_file,LCS);
        std::cerr<< "LCS_file loaded"<<std::endl;

        string fmin_bv_file = opts["fmin_bv"].as<string>();
        sdsl::bit_vector fmin_bv;
        load_bv(fmin_bv_file,fmin_bv);
        std::cerr<< "fmin_bv_file loaded"<<std::endl;
        sdsl::rank_support_v5<> fmin_rs(&fmin_bv);

        string unitigs_v_file = opts["global-offsets"].as<string>();
        sdsl::int_vector<> unitigs_v;
        load_v(unitigs_v_file,unitigs_v);
        std::cerr<< "offsets loaded"<<std::endl;

        /* 
        string endpoints_file = opts["endpoints"].as<string>();
        //std::vector<uint32_t> endpoints;
        //load_intv32(endpoints_file,endpoints);
        sdsl::int_vector<> endpoints;
        load_v(endpoints_file,endpoints);
        std::cerr<< "endpoints loaded"<<std::endl;
        sdsl::sd_vector<> ef_endpoints(endpoints.begin(),endpoints.end()); // Elias-Fano

        string u_start_file = opts["unitig-start"].as<string>();
        sdsl::bit_vector unitig_start;
        load_bv(u_start_file,fmin_bv);
        std::cerr<< "unitig-start loaded"<<std::endl;
        sdsl::rank_support_v5<> Ustart_rs(&fmin_bv);
        */

        // Load packed_unitigs
        PackedStrings unitigs;
        std::ifstream packed_unitigs_in(indexfile + "packed_unitigs.sdsl");
        sdsl::load(unitigs.concat, packed_unitigs_in);

        // Load unitig_endpoints
        std::ifstream unitig_endpoints_in(indexfile + "unitig_endpoints.sdsl");
        sdsl::load(unitigs.ends, unitig_endpoints_in);

        // Load Ustart
        sdsl::bit_vector Ustart;
        std::ifstream Ustart_in(indexfile + "Ustart.sdsl");
        sdsl::load(Ustart, Ustart_in);
        sdsl::rank_support_v5<> Ustart_rs(&Ustart);

        number_of_queries += run_fmin_queries(query_files, output_files, indexfile, sbwt, gzip_output, DNA_bitvectors, DNA_rs, LCS, fmin_rs, unitigs_v, unitigs, Ustart_rs, t);
        int64_t new_total_micros = cur_time_micros() - micros_start;
        write_log("us/query end-to-end: " + to_string((double)new_total_micros / number_of_queries), LogLevel::MAJOR);
        write_log("total number of queries: " + to_string(number_of_queries), LogLevel::MAJOR);
        
        std::ofstream statsfile2;
        statsfile2.open(indexfile + "stats.txt", std::ios_base::app); // append instead of overwrite
        string results = to_string(number_of_queries);
        statsfile2 << "," + to_string((double)new_total_micros / number_of_queries);
        
        /* TODO Jarno
        size_t bytes = size_in_bytes(LCS, fmin_bv, fmin_rs, unitigs_v,ef_endpoints, sbwt);
        write_log("bytes: " + to_string(bytes), LogLevel::MAJOR);
        statsfile2 <<  "," + to_string(bytes);
        statsfile2 << "," + to_string(static_cast<float>(bytes*8)/sbwt.number_of_kmers()) + "\n";
        statsfile2 << "," + to_string(sbwt.number_of_kmers()) + "\n";
        statsfile2.close();
        */
    }

    int64_t total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)total_micros / number_of_queries), LogLevel::MAJOR);
    
    return 0;

}