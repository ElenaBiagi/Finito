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

        int64_t i = 0; // Position in concatenation
        int64_t ends_idx = 0; // Position in ends vector
        int64_t end = 0; // End of current unitig
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
    pair<int64_t, int64_t> global_offset_to_local_offset(int64_t global_offset) const{
        assert(global_offset >= 0 && global_offset < concat.size());

        // Binary search the smallest index in ends that has value larger than the global offset
        int64_t ends_idx = std::upper_bound(ends.begin(), ends.end(), global_offset) - ends.begin();
        int64_t global_start = (ends_idx == 0 ? 0 : ends[ends_idx-1]);

        return {ends_idx, global_offset - global_start};

    }
};

// Returns packed unitigs and the Ustart bit vector
template <typename reader_t>
pair<PackedStrings, sdsl::bit_vector> permute_unitigs(const plain_matrix_sbwt_t& sbwt, reader_t& unitig_reader){
    int64_t k = sbwt.get_k();
    vector<pair<Kmer<MAX_KMER_LENGTH>, int64_t>> first_kmers; // pairs (kmer, unitig id)
    vector<string> unitigs;

    int64_t unitig_id = 0;
    while(true){
        int64_t len = unitig_reader.get_next_read_to_buffer();
        if(len == 0) [[unlikely]] break;

        unitigs.push_back(unitig_reader.read_buf);
        first_kmers.push_back({Kmer<MAX_KMER_LENGTH>(unitigs.back().c_str(), k), unitig_id});

        unitig_id++;
    }

    std::sort(first_kmers.begin(), first_kmers.end()); // Sorts by the kmer comparison operator, which is colexicographic
    vector<int64_t> permutation;
    for(auto& P : first_kmers){
        permutation.push_back(P.second);
    }

    sdsl::bit_vector Ustart(sbwt.number_of_subsets(), 0);
    for(int64_t i = 0; i < unitigs.size(); i++){
        int64_t colex = sbwt.search(unitigs[i].substr(0, k));
        if(colex == -1) cout << "Error: kmer " + unitigs[i].substr(0, k) + " in unitigs but not found in SBWT" << endl;
        Ustart[colex] = 1;
    }

    return {PackedStrings(unitigs, permutation), Ustart};


}