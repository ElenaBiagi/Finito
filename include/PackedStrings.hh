#pragma once

#include <fstream>
#include <string>
#include <tuple>
#include <set>
#include <string.h>
#include <climits>
#include "sbwt/globals.hh"
#include "sbwt/SBWT.hh"
#include "sbwt/stdlib_printing.hh"
#include "sbwt/SeqIO.hh"
#include <vector>
#include <utility>
#include <algorithm>
#include "variants.hh"
#include "sdsl/vectors.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/bit_vectors.hpp"
#include "sbwt/Kmer.hh"
//#include "lcs_basic_algorithm.hpp"

using namespace std;
using namespace sbwt;

class PackedStrings{
    public:
    sdsl::int_vector<2> concat;
    sdsl::int_vector<> ends; // Exclusive ends

    static constexpr char ACGT[] = "ACGT";

    PackedStrings() {}

    // Concatenates the strings according to the given permutation
    PackedStrings(const vector<string>& strings, const vector<int64_t>& permutation){
        assert(strings.size() == permutation.size());
        int64_t total_length = 0;
        for(const string& S : strings){
            total_length += S.size();
        }

        this->concat = sdsl::int_vector<2>(total_length);
        this->ends = sdsl::int_vector<>(strings.size(), 64 - __builtin_clzll(total_length));
        
        int64_t ends_idx = 0; // Position in ends vector
        int64_t end = 0; // End of current unitig
        int64_t unitigs_count = strings.size()/2;
        int64_t i = 0; // Position in concatenation

        for(int64_t string_idx : permutation){
            const string& S = strings[string_idx];

            for(char c : S){
                switch(c){
                    case 'A': concat[i++] = 0; break;
                    case 'C': concat[i++] = 1; break;
                    case 'G': concat[i++] = 2; break;
                    case 'T': concat[i++] = 3; break;
                    default: throw std::runtime_error("Invalid character: " + c);
                }
            }
            end += S.size();
            ends[ends_idx++] = end;
        }
    }
    
    // Clears the given buffer and stores string with index string_idx into it,
    // including a null-terminator. Returns the length of the stored string,
    // not counting the null.
    int64_t get(int64_t string_idx, vector<char>& buffer) const{
        assert(string_idx < ends.size());
        int64_t start = 0;
        if(string_idx > 0) start = ends[string_idx-1];

        int64_t end = ends[string_idx]; // Exclusive end
        int64_t len = end - start;

        buffer.clear();
        for(int64_t i = 0; i < len; i++){
            buffer.push_back(ACGT[concat[start+i]]);
        }
        buffer.push_back(0);

        return len;
    }

    int64_t number_of_strings() const{
        return ends.size();
    }

    // Returns pair (unitig_id, offset_in_unitig)
    pair<int64_t, int64_t> global_offset_to_local_offset(int64_t global_offset, const int64_t unitig_rank, bool rc, int64_t k) const{
        assert(global_offset >= 0 && global_offset < concat.size());
        
        // Binary search the smallest index in ends that has value larger than the global offset
        // this cannot be 0

        int64_t ends_idx = std::upper_bound(ends.begin(), ends.end(), global_offset) - ends.begin();
        int64_t global_start = (ends_idx == 0) ? 0 : ends[ends_idx-1];
        int64_t local_offset = global_offset - global_start;
        if (rc) {
            int64_t unitig_len = ends[ends_idx] - global_start;
            local_offset = unitig_len - local_offset - k;
        } 
        return {unitig_rank, local_offset};
    }
};

// Returns packed unitigs + the Ustart bit vector + Rstart bit vector
template <typename reader_t>
tuple<PackedStrings, sdsl::bit_vector, sdsl::bit_vector, sdsl::int_vector<> > permute_unitigs(const plain_matrix_sbwt_t& sbwt, reader_t& f_reader, reader_t& r_reader){
    int64_t k = sbwt.get_k();
    vector<pair<Kmer<MAX_KMER_LENGTH>, int64_t>> first_kmers; // pairs (kmer, unitig id)
    vector<string> f_unitigs;

    // forward unitigs
    int64_t unitig_id = 0;
    while(true){
        int64_t len = f_reader.get_next_read_to_buffer();
        if(len == 0) [[unlikely]] break;

        f_unitigs.push_back(f_reader.read_buf);
        first_kmers.push_back({Kmer<MAX_KMER_LENGTH>(f_unitigs.back().c_str(), k), unitig_id});

        unitig_id++;
    }

    // permutation for sorting forward unitigs
    std::sort(first_kmers.begin(), first_kmers.end()); // Sorts by the kmer comparison operator, which is colexicographic
    vector<int64_t> f_permutation;
    vector<pair<Kmer<MAX_KMER_LENGTH>, int64_t>> sorted_first_kmers; // pairs (kmer, unitig id)    

    unitig_id = 0;
    for(auto& P : first_kmers){
        f_permutation.push_back(P.second);
        // create the correct first_kmers based on the sorted f unitigs
        sorted_first_kmers.push_back({P.first,unitig_id++});
    }

    // rc unitigs
    vector<string> r_unitigs;
    while(true){
        int64_t len = r_reader.get_next_read_to_buffer();
        if(len == 0) [[unlikely]] break;
        r_unitigs.push_back(r_reader.read_buf);
        // todo later after they are sorted
        //r_first_kmers.push_back({Kmer<MAX_KMER_LENGTH>(unitigs.back().c_str(), k), unitig_id});
    }

    // Sort f and rc unitgs based on f_permutation
    vector<string> sorted_f_unitigs;
    vector<string> sorted_r_unitigs;
    for(auto& P : f_permutation){
        sorted_f_unitigs.push_back(f_unitigs[P]);
        sorted_r_unitigs.push_back(r_unitigs[P]);
    }
    
    // Sort rc unitigs based on first kmers
    vector<pair<Kmer<MAX_KMER_LENGTH>, int64_t>> r_first_kmers;
    for (int u=0; u < sorted_r_unitigs.size(); u++){
        r_first_kmers.push_back({Kmer<MAX_KMER_LENGTH>(sorted_r_unitigs[u].c_str(), k), u});
        // Add rc first kmers
        sorted_first_kmers.push_back({Kmer<MAX_KMER_LENGTH>(sorted_r_unitigs[u].c_str(), k), u + f_unitigs.size()});
    }

    // Sort rc unitigs first kmers to get r_permutation
    std::sort(r_first_kmers.begin(), r_first_kmers.end()); // Sorts by the kmer comparison operator, which is colexicographic
    sdsl::int_vector<> r_permutation(r_unitigs.size(), 0, 64 - __builtin_clzll(r_unitigs.size()));
    
    for(int64_t r = 0; r < r_first_kmers.size();r++){
        r_permutation[r]=r_first_kmers[r].second;
    }

    // Concatenate f and r unitigs vectors
    sorted_f_unitigs.insert(sorted_f_unitigs.end(), sorted_r_unitigs.begin(), sorted_r_unitigs.end());
    // get a permutation to sort all unitigs together
    
    std::sort(sorted_first_kmers.begin(), sorted_first_kmers.end()); // Sorts by the kmer comparison operator, which is colexicographic
   
    vector<int64_t> permutation;
    vector <string> sorted_unitigs;
    for(auto& P : sorted_first_kmers){
        permutation.push_back(P.second);
        //sorted_unitigs.push_back(sorted_f_unitigs[P.second]);
    }

    sdsl::bit_vector Ustart(sbwt.number_of_subsets(), 0);
    sdsl::bit_vector Rstart(sbwt.number_of_subsets(), 0);
    for(int64_t i = 0; i < sorted_r_unitigs.size(); i++){ // sorted_r_unitigs only contains rc_unitigs, sorted_f_unitigs = sorted_f + sorted_r
        int64_t f_colex = sbwt.search(sorted_f_unitigs[i].substr(0, k));
        int64_t r_colex = sbwt.search(sorted_r_unitigs[i].substr(0, k));
        if(f_colex == -1) cout << "Error: kmer " + sorted_f_unitigs[i].substr(0, k) + " in forward unitigs but not found in SBWT" << endl;
        if(r_colex == -1) cout << "Error: kmer " + sorted_r_unitigs[i].substr(0, k) + " in rc unitigs but not found in SBWT" << endl;
        Ustart[f_colex] = 1;
        Ustart[r_colex] = 1;
        Rstart[r_colex] = 1;
    }

    return {PackedStrings(sorted_f_unitigs, permutation), Ustart, Rstart, r_permutation};


}
