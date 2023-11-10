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
vector<string> paper_example_unitigs = {"ACAGGTA", "GTAGGAAA", "GTAAGTCT"};
//                                          1           2           0
vector<string> paper_example_queries = {"ACAGGTA", "GTAGGAAA", "GTAAGTCT", "TAGGATTTTTTAAGTCTA"};

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

// Takes in a spectrum-preserving string set
unique_ptr<FinimizerIndex> build_index(const vector<string>& spss, int64_t k){

    unique_ptr<plain_matrix_sbwt_t> sbwt = make_unique<plain_matrix_sbwt_t>();
    NodeBOSSInMemoryConstructor<plain_matrix_sbwt_t> constructor;
    constructor.build(spss, *sbwt, k, true);

    unique_ptr<sdsl::int_vector<>> LCS = make_unique<sdsl::int_vector<>>(move(lcs_basic_parallel_algorithm(*sbwt, 3))); // 3 threads

    string input_filename = temp_dir + "/spss.fna";
    write_as_fasta(spss, input_filename);
    SeqIO::Reader<> reader(input_filename);
    FinimizerIndexBuilder builder(move(sbwt), move(LCS), reader);
    unique_ptr<FinimizerIndex> index = builder.get_index();
    return move(index);
}

unique_ptr<FinimizerIndex> build_example_index(){
    return build_index(paper_example_unitigs, 4);
}

void test_shortest_unique_construction(){

    unique_ptr<FinimizerIndex> index = build_example_index();

    sdsl::int_vector<> true_LCS = {0,0,1,2,2,1,1,1,0,1,0,2,2,1,3,0,1,2};
    sdsl::util::bit_compress(true_LCS);
    sdsl::int_vector<2> true_unitig_concat = {2,3,0,0,2,3,1,3, 0,1,0,2,2,3,0, 2,3,0,2,2,0,0,0};
    sdsl::int_vector<> true_unitig_ends = {8,15,23};
    sdsl::bit_vector true_fmin = {0,0,0,0,0,1,1,1,0,1,1,0,0,0,0,0,0,1}; // B
    sdsl::int_vector<> true_global_offsets = {10,20,14,6,4,13};
    sdsl::util::bit_compress(true_global_offsets);
    sdsl::bit_vector true_Ustart =  {0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0};
    //sdsl::bit_vector true_is_branch =  {0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0};

    assert_equal(true_LCS, *index->LCS);
    assert_equal(true_unitig_concat, index->unitigs.concat);
    assert_equal(true_unitig_ends, index->unitigs.ends);
    assert_equal(true_fmin, index->fmin);
    assert_equal(true_global_offsets, index->global_offsets);
    assert_equal(true_Ustart, index->Ustart);
    //assert_equal(true_is_branch, index->is_branch);
}

void test_shortest_unique_queries(){
    vector<vector<pair<int64_t, int64_t>>> true_local_offsets = {{{1,0},{1,1},{1,2},{1,3}}, 
                                                                 {{2,0},{2,1},{2,2},{2,3},{2,4}}, 
                                                                 {{0,0},{0,1},{0,2},{0,3},{0,4}},
                                                                 {{2,1}, {2,2}, {-1,-1}, {-1,-1}, {-1,-1}, {-1,-1}, {-1,-1}, {-1,-1}, {-1,-1}, {-1,-1}, {0,1}, {0,2}, {0,3}, {0,4}, {-1,-1}}};// {unitig_id , local kmer_start}

    // TAGGATTTTTTAAGTCTA
    vector<int64_t> true_hit_counts = {4,5,5,6};
    
    unique_ptr<FinimizerIndex> index = build_example_index();

    for(int64_t i = 0; i < paper_example_queries.size(); i++){
        const string& S = paper_example_queries[i];
        FinimizerIndex::QueryResult res = index->search(S);
        for(auto x : res.local_offsets) cout << x << " "; cout << endl;

        assert_equal(res.n_found, true_hit_counts[i]);
        assert_equal(res.local_offsets, true_local_offsets[i]);
    }
}

void test_finimizer_branch(){
    // GG appears only in ACGG as part of a full kmer but it is not found as it is not the finimizer
    // GG is stored for the k-mer $CGT
    int64_t k = 4;
    vector<string> unitigs = {"ACGG", "CGGT", "GCCGT", "CGGC"} ;
    // Permuted order:            2      3        1       0
    string query = "ACGGC";
    /* 
    0 $$$$ 
    1 $$$A *
    2 $$AC
    3 $GCC *
    4 $$GC
    5 CGGC s
    6 $$$G
    7 $ACG
    8 GCCG s
    9 ACGG * s b
    10 CCGT
    11 CGGT s
    */

    sdsl::int_vector<> true_LCS = {0,0,0,1,1,2,0,1,2,1,0,2};
    sdsl::util::bit_compress(true_LCS);
    sdsl::int_vector<> true_unitig_ends = {4,9,13,17};
    sdsl::bit_vector true_fmin = {0,1,0,1,0,0,0,0,0,1,0,0};
    sdsl::int_vector<> true_global_offsets = {9,6,15};
    sdsl::util::bit_compress(true_global_offsets);
    sdsl::bit_vector true_Ustart = {0,0,0,0,0,1,0,0,1,1,0,1};
    sdsl::bit_vector true_is_branch =  {0,0,0,0,0,0,0,0,0,1,0,0};

    unique_ptr<FinimizerIndex> index = build_index(unitigs, 4);

    cout << (int)index->global_offsets.width() << endl;
    cout << (int)true_global_offsets.width() << endl;
    assert_equal(true_LCS, *index->LCS);
    assert_equal(true_unitig_ends, index->unitigs.ends);
    assert_equal(true_fmin, index->fmin);
    assert_equal(true_global_offsets, index->global_offsets);
    assert_equal(true_Ustart, index->Ustart);
    assert_equal(true_is_branch, index->is_branch);


    vector<pair<int64_t, int64_t>> true_local_offsets = {{2,0}, {0,0}};

    FinimizerIndex::QueryResult res = index->search(query);
    assert_equal(res.local_offsets.size(), true_local_offsets.size());
    assert_equal(res.local_offsets, true_local_offsets);
}

void test_reverse_complement_branch(){
    // ACGG has outgoing edges T and C, but the one with C goes to the reverse complemented k-mer GCCG (CGGC)
    int64_t k = 4;
    //vector<string> unitigs = {"CCGG", "CGGT", "GCCGT"}; // These are not unitigs as CCG is at the start of the first and middle of last
    vector<string> unitigs = {"TCGG", "CGGT", "GCCGTC"}; 
    // Permuted order:            1       2        0
    //string query = "CCGGT";
    //vector<pair<int64_t, int64_t>> true_local_offsets = {{1,0}, {2,0}};
    string query = "TCGGTGCCGTC";
    vector<pair<int64_t, int64_t>> true_local_offsets = {{1,0}, {2,0},{-1,-1},{-1,-1}, {-1,-1}, {0,0}, {0,1}, {0,2}};

    unique_ptr<FinimizerIndex> index = build_index(unitigs, 4);
    FinimizerIndex::QueryResult res = index->search(query);
    assert_equal(res.local_offsets.size(), true_local_offsets.size());
    assert_equal(res.local_offsets, true_local_offsets);
}

void test_leftmost(){
    // ACGG has outgoing edges T and C, but the one with C goes to the reverse complemented k-mer GCCG (CGGC)
    int64_t k = 4;
    vector<string> unitigs = {"CGGT", "GGTT", "TACCCGTA"}; // GG is the finimizer of CGGT and GGTC and it is stored in GGTC // reverse complement wit a C
    // Permuted order:           1       2        0
    //string query = "CCGGT";
    //vector<pair<int64_t, int64_t>> true_local_offsets = {{1,0}, {2,0}};
    string query = "CGGTTACCC";
    vector<pair<int64_t, int64_t>> true_local_offsets = {{1,0}, {2,0}, {-1,-1}, {-1,-1}, {0,0}, {0,1}};

    unique_ptr<FinimizerIndex> index = build_index(unitigs, 4);
    FinimizerIndex::QueryResult res = index->search(query);
    assert_equal(res.local_offsets.size(), true_local_offsets.size());
    assert_equal(res.local_offsets, true_local_offsets);
}



void test_finimizer_selection(){
    // ACGG has outgoing edges T and C, but the one with C goes to the reverse complemented k-mer CGGC
    int64_t k = 4;
    vector<string> unitigs = {"ACGG", "CGGT", "GCCGTA"};
    string query = "GCCGTA";
    // Permuted order:            1       2        0

    unique_ptr<FinimizerIndex> index = build_index(unitigs, 4);
    index->search(query);

    sdsl::bit_vector true_fmin = {0,0,1,1,1,0,0,0,0,1,0,0};
    assert_equal(true_fmin, index->fmin);

    /*
    0 $$$$
    1 $$$A
    2 CGTA*
    3 $$AC*
    4 $GCC*
    5 $$GC
    6 $$$G
    7 $ACG
    8 GCCG
    9 ACGG*
    10 CCGT
    11 CGGT
    */

}

// Return map: unitig -> rank
map<string, int64_t> get_unitig_ranks(const vector<string>& unitigs, int64_t k){

    vector<string> vec = unitigs;

    auto unitig_compare = [k](const string& A, const string& B){
        string Ak = A.substr(0, k);
        string Bk = B.substr(0, k);
        string Ak_rev(Ak.rbegin(), Ak.rend());
        string Bk_rev(Bk.rbegin(), Bk.rend());
        return Ak_rev < Bk_rev;
    };

    sort(vec.begin(), vec.end(), unitig_compare);

    map<string, int64_t> M;
   
    for(int64_t i = 0; i < unitigs.size(); i++){
        M[vec[i]] = i;
    };
    return M;

}


void test_incoming_rc_branch(){
    int64_t k = 10;
    vector<string> unitigs = {"AACAAAAAAA",               "ACAAAAAAAA", "CAAAAAAAAA", 
                 sbwt::get_rc("TACAAAAAAA"), sbwt::get_rc("TCAAAAAAAA")};

    unique_ptr<FinimizerIndex> index = build_index(unitigs, k);

    map<string, int64_t> unitig_ranks = get_unitig_ranks(unitigs, k);

    string query = "CAAAAAAAAA";
    FinimizerIndex::QueryResult ans = index->search(query);

    for(string unitig : unitigs) cout << unitig << " " << unitig_ranks[unitig] << endl;

    vector<pair<int64_t, int64_t>> true_local_offsets = {{unitig_ranks[query], 0}};
    assert_equal(true_local_offsets, ans.local_offsets);

}

int main(int argc, char** argv){
    // Create test directory if does not exist
    if (!exists(temp_dir)){
        create_directory(temp_dir);
    }

    cerr << "Testing shortest unique construction..." << endl;
    test_shortest_unique_construction();
    cerr << "...ok" << endl;

    cerr << "Testing shortest unique queries..." << endl;
    test_shortest_unique_queries();
    cerr << "...ok" << endl;

    cerr << "Testing finimizer branch" << endl;
    test_finimizer_branch();
    cerr << "...ok" << endl;

    cerr << "Testing reverse complement branch" << endl;
    test_reverse_complement_branch();
    cerr << "...ok" << endl;

    cerr << "Testing leftmost" << endl;
    test_leftmost();
    cerr << "...ok" << endl;

    cerr << "Testing Finimizer selection" << endl;
    test_finimizer_selection();
    cerr << "...ok" << endl;
 
    cerr << "Testing incoming rc branch" << endl;
    test_incoming_rc_branch();
    cerr << "...ok" << endl;

    cerr << "ALL TESTS PASSED" << endl;

}
