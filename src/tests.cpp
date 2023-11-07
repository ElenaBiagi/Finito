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
string paper_example_input = ">3\nGTAAGTCT\n>1\nACAGGTA\n>2\nGTAGGAAA\n";

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

// Returns filename
string write_paper_example_to_disk(){
    // Write fasta data to file
    string fasta_filename = temp_dir + "/" + "paper_example.fasta";
    ofstream fasta_out(fasta_filename);
    fasta_out.write(paper_example_input.c_str(), paper_example_input.size());
    return fasta_filename;
}

void test_shortest_unique_construction(){
    string input_filename = write_paper_example_to_disk();

    unique_ptr<plain_matrix_sbwt_t> sbwt = make_unique<plain_matrix_sbwt_t>();
    int64_t k = 4;
    NodeBOSSInMemoryConstructor<plain_matrix_sbwt_t> constructor;
    constructor.build({paper_example_input}, *sbwt, k, true);

    unique_ptr<sdsl::int_vector<>> LCS = make_unique<sdsl::int_vector<>>(move(lcs_basic_parallel_algorithm(*sbwt, 3))); // 3 threads

    SeqIO::Reader<> reader(input_filename);
    FinimizerIndexBuilder builder(move(sbwt), move(LCS), reader);
    unique_ptr<FinimizerIndex> index = builder.get_index();

    // TODO: fill these in
    sdsl::int_vector<> true_LCS;
    sdsl::int_vector<2> true_unitig_concat;
    sdsl::int_vector<> true_unitig_ends;
    sdsl::bit_vector true_fmin;
    sdsl::int_vector<> true_global_offsets;
    sdsl::bit_vector true_Ustart;

    assert_equal(true_LCS, *index->LCS);
    assert_equal(true_unitig_concat, index->unitigs.concat);
    assert_equal(true_unitig_ends, index->unitigs.ends);
    assert_equal(true_fmin, index->fmin);
    assert_equal(true_global_offsets, index->global_offsets);
    assert_equal(true_Ustart, index->Ustart);

}

int main(int argc, char** argv){
    // Create test directory if does not exist
    if (!exists(temp_dir)){
        create_directory(temp_dir);
    }
    test_shortest_unique_construction();

}