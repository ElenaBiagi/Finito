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
        vector<pair<int64_t, int64_t>> local_offsets; // {Unitig id, distance from the start of the unitig}
        int64_t n_found = 0;
    };

private:
    // Forbid copying because we have pointers to our internal data structures
    FinimizerIndex(const FinimizerIndex& other) = delete;
    FinimizerIndex& operator=(const FinimizerIndex& other) = delete;

    void add_to_query_result(const int64_t global_unitig_rank, const int64_t unitig_rank, const bool rc, const int64_t global_kmer_end, QueryResult& answer) const{
        int64_t global_kmer_start = global_kmer_end - sbwt->get_k() + 1;
        pair<int64_t, int64_t> local_start = unitigs.global_offset_to_local_offset(global_unitig_rank, global_kmer_start, unitig_rank, rc, sbwt->get_k());
        answer.local_offsets.push_back(local_start);
        answer.n_found++;
    }
                
    void walk_in_unitigs(const std::string& query, const PackedStrings& unitigs, int64_t global_kmer_end, QueryResult& answer, int64_t& kmer_end, const int64_t k, int64_t unitig_rank, bool rc, int64_t global_unitig_rank, int64_t local_offset) const{
        
        int64_t u_end = unitigs.ends[global_unitig_rank]; // exlusive end
        assert(u_end > global_kmer_end);
       
        int64_t max_match = std::min(u_end - global_kmer_end-1, (int64_t)((query.length()-kmer_end)-1));
        if (max_match <= 0){return;}
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
        int64_t global_kmer_end_copy = global_kmer_end +1;

        while(max_match > 0){
            int word_len = std::min({(int64_t)32, max_match});
            uint64_t query_word = query_v.get_int(word_start, (word_len)*2); // word start = least significant bit (word_len+1)*2-1
            uint64_t unitig_word = unitigs.concat.get_int(word_start + (global_kmer_end_copy)*2, (word_len)*2);
            for (int64_t i = global_kmer_end_copy; i< global_kmer_end_copy + max_match;i++){
        }
            int64_t result = query_word ^ unitig_word;
            
            if (result){
                int trailing_zeros = __builtin_ctzll(result);
                for (int i = 0; i < trailing_zeros/2; i++){
                    global_kmer_end++;
                    (rc)? local_offset--: local_offset++;
                    answer.local_offsets.push_back({unitig_rank, local_offset});
                    answer.n_found++;
                    kmer_end++; // keep track of how many kmers were found   
                }
                break; // No need to check further. A mismatch has been found. // exit the while loop
            }
            int trailing_zeros = (word_len)*2; // the whole word_len matches
            for (int i = 0; i < trailing_zeros/2; i++){
                global_kmer_end++;
                (rc)? local_offset--: local_offset++;
                answer.local_offsets.push_back({unitig_rank, local_offset});
                answer.n_found++;
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
    sdsl::select_support_mcl<> Ustart_ss;

    sdsl::bit_vector Rstart;
    sdsl::rank_support_v5<> Rstart_rs;
    sdsl::int_vector<> Rpermutation;

    FinimizerIndex() {}

    QueryResult search(const std::string& query) const {
        // For each k-mer S that is known to be in the SBWT
        //   - Find the finimizer x.
        //   - Walk forward in the SBWT from the colex rank of x (singleton interval), to the end of S,
        //     recording the rightmost branch point, and the distance from the end of x to this branch
        //     point.
        //   - If a branch point was found, the k-mer containing x is at the unitig after the rightmost
        //     branch. Otherwise, the k-mer containing x is at the unitig that x points to.
        //   - If there was a branch, the k- mer endpoint in that unitig is k + the number of steps taken after the last branch.
        // Branches are marked by the Ustart bitvector

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

            if(kmer_colex_ranks[kmer_end].has_value()){
                // kmer exists
                int64_t global_kmer_end;
                int64_t finimizer_end = finimizers_ends_colex[kmer_end].value().first;
                
                optional<pair<int64_t, int64_t>> rightmost_branch_end = rightmost_Ustart[kmer_end];
                // Get rightmost Ustart
                int64_t colex = 0 ;
                int64_t global_unitig_rank = 0;
                int64_t unitig_rank = 0; // does not include rc unitigs
                int64_t local_offset = 0;
                int64_t global_start = 0;
                bool rc = false;
                if(rightmost_branch_end.has_value()) {
                    // Look up from the Ustart dictionary
                    int64_t p = rightmost_branch_end.value().first;
                    colex = rightmost_branch_end.value().second; 
                    // Get the global offset of the end of the k-mer

                    global_unitig_rank = Ustart_rs.rank(colex);
                    assert(global_unitig_rank < unitigs.ends.size());
                                        
                    // Check if it is found in a rc unitig
                    std::tie(rc,unitig_rank) = unitigs.check_reverse(colex, Rstart, Rstart_rs, Rpermutation, global_unitig_rank);


                    global_start = (global_unitig_rank == 0) ? 0 : unitigs.ends[global_unitig_rank-1];
                    local_offset = kmer_end - p; // Shift to the right place in the unitig
                    global_kmer_end = global_start + local_offset +k -1; // Shift to the right place in the unitig concatenation
                                    
                } else {
                    //Look up from Finimizer dictionary
                    int64_t p = finimizer_end;
                    colex = finimizers_ends_colex[kmer_end].value().second;
                    global_kmer_end = lookup_from_finimizer_dictionary(colex, fmin_rs, global_offsets);
                    global_kmer_end += kmer_end - p; // Shift to the right place in the unitig
                    // Find unitig rank
                    global_unitig_rank = std::upper_bound(unitigs.ends.begin(), unitigs.ends.end(), global_kmer_end-k+1) - unitigs.ends.begin();
                    assert(global_unitig_rank < unitigs.ends.size());
                    colex = Ustart_ss(global_unitig_rank+1);
                    
                    // Check if it is found in a rc unitig
                    std::tie(rc,unitig_rank) = unitigs.check_reverse(colex, Rstart, Rstart_rs, Rpermutation, global_unitig_rank);

                    // add_to_query_result
                    int64_t global_kmer_start = global_kmer_end - k + 1;
                    global_start = (global_unitig_rank == 0) ? 0 : unitigs.ends[global_unitig_rank-1];
                    local_offset = global_kmer_start - global_start;
                }

                if (rc) {
                    int64_t unitig_len = unitigs.ends[global_unitig_rank] - global_start;
                    local_offset = unitig_len - local_offset - k;
                }        
                answer.local_offsets.push_back({unitig_rank, local_offset});
                answer.n_found++;
                
                // walk in the current unitig                
                if ((kmer_end + 1)< query.size()){
                    walk_in_unitigs(query, unitigs, global_kmer_end, answer, kmer_end, k, unitig_rank, rc, global_unitig_rank, local_offset);
                }  
            } else { 
                answer.local_offsets.push_back({-1, -1});
            } 
        }
        return answer;
    }

    void serialize(const string& index_prefix) const {
        sbwt->serialize(index_prefix + ".sbwt");

        std::ofstream LCS_out(index_prefix + ".LCS.sdsl");
        sdsl::serialize(*LCS, LCS_out);

        std::ofstream packed_unitigs_out(index_prefix + ".packed_unitigs.sdsl");
        sdsl::serialize(unitigs.concat, packed_unitigs_out);

        std::ofstream unitig_endpoints_out(index_prefix + ".unitig_endpoints.sdsl");
        sdsl::serialize(unitigs.ends, unitig_endpoints_out);

        std::ofstream fmin_out(index_prefix + ".FBV.sdsl");
        sdsl::serialize(fmin, fmin_out);

        std::ofstream global_offsets_out(index_prefix + ".O.sdsl");
        sdsl::serialize(global_offsets, global_offsets_out);

        std::ofstream Ustart_out(index_prefix + ".Ustart.sdsl");
        sdsl::serialize(Ustart, Ustart_out);

        std::ofstream Rstart_out(index_prefix + ".Rstart.sdsl");
        sdsl::serialize(Rstart, Rstart_out);

        std::ofstream Rpermutation_out(index_prefix + ".R.sdsl");
        sdsl::serialize(Rpermutation, Rpermutation_out);
    }

    void load(const string& index_prefix) {

        sbwt = make_unique<plain_matrix_sbwt_t>();
        sbwt->load(index_prefix + ".sbwt");

        LCS = make_unique<sdsl::int_vector<>>();
        ifstream LCS_in(index_prefix + ".LCS.sdsl");
        sdsl::load(*LCS, LCS_in);
        std::cerr<< "LCS_file loaded"<<std::endl;

        std::ifstream packed_unitigs_in(index_prefix + ".packed_unitigs.sdsl");
        sdsl::load(unitigs.concat, packed_unitigs_in);
        std::cerr << "unitigs loaded" << std::endl;

        std::ifstream unitig_endpoints_in(index_prefix + ".unitig_endpoints.sdsl");
        sdsl::load(unitigs.ends, unitig_endpoints_in);
        std::cerr << "unitig endpoints loaded" << std::endl;

        ifstream fmin_bv_in(index_prefix + ".FBV.sdsl");
        sdsl::load(fmin, fmin_bv_in);
        std::cerr<< "fmin_bv_file loaded"<<std::endl;
        sdsl::util::init_support(fmin_rs, &fmin);

        ifstream global_offsets_in(index_prefix + ".O.sdsl");
        sdsl::load(global_offsets, global_offsets_in);
        std::cerr<< "offsets loaded"<<std::endl;

        std::ifstream Ustart_in(index_prefix + ".Ustart.sdsl");
        sdsl::load(Ustart, Ustart_in);
        sdsl::util::init_support(Ustart_rs, &Ustart);
        sdsl::util::init_support(Ustart_ss, &Ustart);
        std::cerr << "Ustart loaded" << std::endl;

        std::ifstream Rstart_in(index_prefix + ".Rstart.sdsl");
        sdsl::load(Rstart, Rstart_in);
        sdsl::util::init_support(Rstart_rs, &Rstart);
        std::cerr << "Rstart loaded" << std::endl;

        ifstream Rpermutation_in(index_prefix + ".R.sdsl");
        sdsl::load(Rpermutation, Rpermutation_in);
        std::cerr<< "Rpermutation loaded"<<std::endl;

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
        total += sdsl::size_in_bytes(Ustart_ss);
        total += sdsl::size_in_bytes(Rstart);
        total += sdsl::size_in_bytes(Rstart_rs);
        total += sdsl::size_in_bytes(Rpermutation);

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
    FinimizerIndexBuilder(unique_ptr<plain_matrix_sbwt_t> sbwt, unique_ptr<sdsl::int_vector<>> LCS, reader_t& f_reader, reader_t& r_reader) {
        index = make_unique<FinimizerIndex>();
        this->sbwt = move(sbwt); // Take ownership
        this->LCS = move(LCS); // Take ownership

        int64_t n_nodes = this->sbwt->number_of_subsets();
        sdsl::bit_vector fmin_bv(n_nodes, 0); // Finimizer marks
        sdsl::int_vector fmin_found(n_nodes, 0);
        
        vector<uint64_t> global_offsets;
        global_offsets.reserve(n_nodes);
        global_offsets.resize(n_nodes, 0);
        
        tuple<PackedStrings, sdsl::bit_vector, sdsl::bit_vector, sdsl::int_vector<>> unitig_data = permute_unitigs(*(this->sbwt), f_reader, r_reader);
        PackedStrings& unitigs = get<0>(unitig_data);
        sdsl::bit_vector& Ustart = get<1>(unitig_data);
        sdsl::bit_vector& Rstart = get<2>(unitig_data);
        sdsl::int_vector<>& Rpermutation = get<3>(unitig_data);

        set<tuple<int64_t, int64_t, int64_t>>  finimizers;
        int64_t total_len = 0;
        vector<char> unitig_buf;
        for(int64_t i = 0; i < unitigs.number_of_strings(); i++){
            int64_t len = unitigs.get(i, unitig_buf);
            set<tuple<int64_t, int64_t, int64_t>> new_search = add_sequence(unitig_buf.data(), fmin_bv, fmin_found, global_offsets, total_len);
            total_len += len;
            finimizers.insert(new_search.begin(), new_search.end());
        }
        std::cerr << "All unitigs have been scanned"<< endl;

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
        index->Ustart_ss = sdsl::select_support_mcl<>(&(index->Ustart));
        index->Rstart = std::move(Rstart); // Transfer ownership
        index->Rstart_rs = sdsl::rank_support_v5<>(&(index->Rstart));
        index->Rpermutation = std::move(Rpermutation); // Transfer ownership
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

                    assert ((unitig_start + get<3>(w_fmin)) > UINT64_MAX);
                        
                    global_offsets[get<2>(w_fmin)]= unitig_start + get<3>(w_fmin);
                    // TODO We could add a bit to say if the unitig is forward or rc
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
