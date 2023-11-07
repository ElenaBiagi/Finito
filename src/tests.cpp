#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "commands.hh"
#include "globals.hh"
#include "build_fmin.hh"
#include "search_fmin.hh"
#include <filesystem>
#include "FinimizerIndex.hh"
#include "lcs_basic_parallel_algorithm.hpp"
#include "backward.hpp"

string temp_dir = "tests_temp";
vector<string> paper_example_unitigs = {"GTAAGTCT", "ACAGGTA", "GTAGGAAA"};

using namespace std;
using namespace std::filesystem;

template<typename T>
void assert_equal(const T& a, const T& b){
    if(a != b){
        backward::StackTrace st; st.load_here(32);
        backward::Printer p; p.print(st);
        cerr << "Assertion failed: " << a << " != " << b << endl;
        exit(1);
    }
}

void write_as_fasta(const vector<string>& seqs, const string& filename){
    ofstream fasta_out(filename);
    for(const string& S : seqs){
        fasta_out << ">\n" << S << "\n" << endl;
    }
}

unique_ptr<FinimizerIndex> build_example_index(){
    string input_filename = temp_dir + "/paper_example.fna";
    write_as_fasta(paper_example_unitigs, input_filename);

    unique_ptr<plain_matrix_sbwt_t> sbwt = make_unique<plain_matrix_sbwt_t>();
    int64_t k = 4;
    NodeBOSSInMemoryConstructor<plain_matrix_sbwt_t> constructor;
    constructor.build(paper_example_unitigs, *sbwt, k, true);

    unique_ptr<sdsl::int_vector<>> LCS = make_unique<sdsl::int_vector<>>(move(lcs_basic_parallel_algorithm(*sbwt, 3))); // 3 threads

    SeqIO::Reader<> reader(input_filename);
    FinimizerIndexBuilder builder(move(sbwt), move(LCS), reader);
    unique_ptr<FinimizerIndex> index = builder.get_index();
    return move(index);
}

void test_shortest_unique_construction(){

    unique_ptr<FinimizerIndex> index = build_example_index();

    // TODO: fill these in
    sdsl::int_vector<> true_LCS = {0,0,1,2,2,1,1,1,0,1,0,2,2,1,3,0,1,2};
    sdsl::util::bit_compress(true_LCS);
    sdsl::int_vector<2> true_unitig_concat = {2,3,0,0,2,3,1,3, 0,1,0,2,2,3,0, 2,3,0,2,2,0,0,0};
    sdsl::int_vector<> true_unitig_ends = {8,15,23};
    sdsl::bit_vector true_fmin = {0,0,0,0,0,1,1,1,0,1,1,0,0,0,0,0,0,1}; // B
    sdsl::int_vector<> true_global_offsets = {10,20,14,6,4,13};
    sdsl::util::bit_compress(true_global_offsets);
    sdsl::bit_vector true_Ustart =  {0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0};

    assert_equal(true_LCS, *index->LCS);
    assert_equal(true_unitig_concat, index->unitigs.concat);
    assert_equal(true_unitig_ends, index->unitigs.ends);
    assert_equal(true_fmin, index->fmin);
    assert_equal(true_global_offsets, index->global_offsets);
    assert_equal(true_Ustart, index->Ustart);
}

void test_shortest_unique_queries(){
    vector<pair<int64_t, int64_t>> true_unique_ptr = {{1,0},{1,1},{1,2},{1,3}, {2,0},{2,1},{2,2},{2,3},{2,4}, {0,0},{0,1},{0,2},{0,3},{0,4}};// {unitig_id , local kmer_start}
    unique_ptr<FinimizerIndex> index = build_example_index();
}

int main(int argc, char** argv){
    // Create test directory if does not exist
    if (!exists(temp_dir)){
        create_directory(temp_dir);
    }
    test_shortest_unique_construction();

}