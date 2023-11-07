#pragma once
#include <sdsl/int_vector.hpp>
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"
#include <vector>

using namespace std;

// End is one past end
void lcs_propagate_thread(char** output_starts, const char* input, const sdsl::bit_vector** bits, int64_t start, int64_t end){
    for(int64_t i = start; i < end; i++){
        for(char sigma = 0; sigma < 4; sigma++){
            if((*(bits[sigma]))[i] == 1){
                *(output_starts[sigma]) = input[i];
                output_starts[sigma]++;
            }
        }
    }
}

void lcs_update_thread(uint8_t* lcs_bytes, const char* last, int64_t n_nodes, int64_t k, int64_t round, int64_t start, int64_t end){
    for(int64_t i = start; i < end; i++){
        if(lcs_bytes[i] == k && (i == 0 || last[i] != last[i-1])){
            lcs_bytes[i] = round;
        }
    }
}

void populate_first_column(char* column, vector<int64_t> C_array, const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits){

    
    int64_t n_nodes = A_bits.size();
    C_array.push_back(n_nodes);

    column[0] = '$';

    int64_t last_idx = 1;

    string ACGT = "ACGT";
    for(char sigma = 0; sigma < 4; sigma++){
        for(int64_t i = 0; i < C_array[sigma+1] - C_array[sigma]; i++){
            column[last_idx++] = ACGT[sigma];
        }
    }   

    if(last_idx != n_nodes){
        cerr << "BUG " << last_idx << " " << n_nodes << endl;
        exit(1);
    }
}

sdsl::int_vector<> lcs_basic_parallel_algorithm(const sbwt::plain_matrix_sbwt_t& SBWT, int64_t n_threads){

    if(SBWT.get_k() > 255){
        cerr << "Error: maximum k supported is 255" << endl;
        // Because we are using 1 byte per LCS value
    }

    const sdsl::bit_vector& A_bits = SBWT.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = SBWT.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = SBWT.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = SBWT.get_subset_rank_structure().T_bits;
    int64_t k = SBWT.get_k();
    int64_t n_nodes = SBWT.number_of_subsets();

    vector<int64_t> C_array = SBWT.get_C_array();

    // We need to keep two columns in memory at the same time, so we create
    // two buffers and alternate which one is being written to
    char* buf1 = new char[n_nodes]; 
    char* buf2 = new char[n_nodes];
    
    populate_first_column(buf1, C_array, A_bits, C_bits, G_bits, T_bits);
    // buf1[i] = incoming character to node i

    // One byte per LCS value to get atomic writes and reads.
    // Values LCS[i] == k means the value is not yet set
    vector<uint8_t> lcs_bytes(n_nodes, k);

    for(int64_t round = 0; round < k; round++){
        cerr << "Round " << round << "/" << k-1 << endl;

        char* last = (round % 2 == 0) ? buf1 : buf2;
        char* propagated = (round % 2 == 0) ? buf2 : buf1;

        propagated[0] = '$'; // This node has no incoming edge so it's a special case

        #pragma omp parallel for num_threads(n_threads)
        for(int64_t t = 0; t < n_threads; t++){
            int64_t start = n_nodes * t / n_threads;
            int64_t end = n_nodes * (t+1) / n_threads;
            lcs_update_thread(lcs_bytes.data(), last, n_nodes, k, round, start, end);
        }

        // Propagate the labels one step forward in the graph

        #pragma omp parallel for num_threads(n_threads)
        for(int64_t t = 0; t < n_threads; t++){
            int64_t start = n_nodes * t / n_threads;
            int64_t end = n_nodes * (t+1) / n_threads;
            char* A_ptr = propagated + C_array[0] + SBWT.get_subset_rank_structure().A_bits_rs(start);
            char* C_ptr = propagated + C_array[1] + SBWT.get_subset_rank_structure().C_bits_rs(start);
            char* G_ptr = propagated + C_array[2] + SBWT.get_subset_rank_structure().G_bits_rs(start);
            char* T_ptr = propagated + C_array[3] + SBWT.get_subset_rank_structure().T_bits_rs(start);
            char* ptrs[4] = {A_ptr, C_ptr, G_ptr, T_ptr};
            const sdsl::bit_vector* bits[4] = {&A_bits, &C_bits, &G_bits, &T_bits};
            lcs_propagate_thread(ptrs, last, bits, start, end);
        }
    }

    delete[] buf1;
    delete[] buf2;

    // Compress to sdsl::int_vector
    sdsl::int_vector<> lcs(n_nodes, 0, 64 - __builtin_clzll(k-1)); // Enough bits per element to store values from 0 to k-1
    for(int64_t i = 0; i < n_nodes; i++)
        lcs[i] = lcs_bytes[i];

    return lcs;
}