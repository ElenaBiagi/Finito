#include <string>
#include <string.h>
#include "cxxopts.hpp"
#include "globals.hh"
#include "SBWT.hh"
#include "stdlib_printing.hh"
#include "SeqIO.hh"
#include "buffered_streams.hh"
#include <vector>
#include <utility>
#include <algorithm>
#include "variants.hh"
#include "sdsl/vectors.hpp"
#include "sdsl/int_vector.hpp"
#include "commands.hh"


using namespace std;
using namespace sbwt;

std::vector<std::string> get_available_types(){
    return {"rarest", "shortest", "optimal", "verify"};
}

int build_fmin(int argc, char** argv) {

    set_log_level(LogLevel::MINOR);

    vector<string> types = get_available_types();
    string all_types_string;
    for (string type: types) all_types_string += " " + type;

    options.add_options()
            ("o,out-file", "Output filename.", cxxopts::value<string>())
            ("i,index-file", "Index input file.", cxxopts::value<string>())
            ("u,in-file",
             "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.",
             cxxopts::value<string>())
            ("add-reverse-complements", "Also add the reverse complement of every k-mer to the index. Warning: this creates a temporary reverse-complemented duplicate of each input file before construction. Make sure that the directory at --temp-dir can handle this amount of data. If the input is gzipped, the duplicate will also be compressed, which might take a while.", cxxopts::value<bool>()->default_value("false"))
            ("z,gzip-output",
             "Writes output in gzipped form. This can shrink the output files by an order of magnitude.",
             cxxopts::value<bool>()->default_value("false"))
            ("type", "Decide which streaming search type you prefer. Available types: " + all_types_string,
             cxxopts::value<string>()->default_value("rarest"))
            ("t", "Maximum finimizer frequency", cxxopts::value<uint64_t>())
            ("lcs", "Provide in input the LCS file if available.", cxxopts::value<string>()->default_value(""))
            ("h,help", "Print usage");

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")) {
        std::cerr << options.help() << std::endl;
        exit(1);
    }
    uint64_t t = opts["t"].as<uint64_t>();
    bool revcomps = opts["add-reverse-complements"].as<bool>();

    // SBWT index
    string sbwt_index_file = opts["index-file"].as<string>();
    // check_readable(sbwt_index_file);
    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(variant != "plain-matrix"){
        cerr << "Error: input is not a plain-matrix SBWT" << endl;
    }
    sbwt::plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain matrix SBWT with " << sbwt.number_of_subsets() << " subsets" << endl;

    // input files
    string in_file = opts["in-file"].as<string>();
    vector<string> input_files;
    bool multi_file = in_file.size() >= 4 && in_file.substr(in_file.size() - 4) == ".txt";
    if(multi_file){
        input_files = readlines(in_file);
    } else{
        input_files = {in_file};
    }
    for(string file : input_files) check_readable(file);

    // Interpret output file
    string outfile = opts["out-file"].as<string>();
    bool gzip_output = opts["gzip-output"].as<bool>();
    vector<string> output_files;

    check_writable(outfile);
    vector<string> variants = get_available_variants_fmin();

    int64_t n_finimizers = 0;


        string type = opts["type"].as<string>();
        if(std::find(types.begin(), types.end(), type) == types.end()){
            cerr << "Error: unknown type: " << type << endl;
            cerr << "Available types are:" << all_types_string << endl;
            return 1;
        }

        const sdsl::bit_vector& A_bits = sbwt.get_subset_rank_structure().A_bits;
        const sdsl::bit_vector& C_bits = sbwt.get_subset_rank_structure().C_bits;
        const sdsl::bit_vector& G_bits = sbwt.get_subset_rank_structure().G_bits;
        const sdsl::bit_vector& T_bits = sbwt.get_subset_rank_structure().T_bits;

        const sdsl::rank_support_v5<> &A_bits_rs = sbwt.get_subset_rank_structure().A_bits_rs;
        const sdsl::rank_support_v5<> &C_bits_rs = sbwt.get_subset_rank_structure().C_bits_rs;
        const sdsl::rank_support_v5<> &G_bits_rs = sbwt.get_subset_rank_structure().G_bits_rs;
        const sdsl::rank_support_v5<> &T_bits_rs = sbwt.get_subset_rank_structure().T_bits_rs;

        const int64_t n_nodes = sbwt.number_of_subsets();
        const int64_t k = sbwt.get_k();


        const vector<int64_t> &C = sbwt.get_C_array();

        //int64_t new_k = sbwt.get_k();
        //vector<string> all_kmers = dump_all_kmers(A_bits, C_bits, G_bits, T_bits, new_k);
        //write_str_v( all_kmers, "all_kmers.txt");

        vector<pair<int64_t, int64_t> > ans;
        //const sdsl::bit_vector *DNA_bitvectors[4] = {&A_bits, &C_bits, &G_bits, &T_bits};
        const sdsl::rank_support_v5<> *DNA_rs[4] = {&A_bits_rs, &C_bits_rs, &G_bits_rs, &T_bits_rs};

        string LCS_file = opts["lcs"].as<string>();
        if (LCS_file.empty()) {
            std::cout << "LCS_file empty" << std::endl;
            LCS_file = sbwt_index_file + "LCS.sdsl";
            const sdsl::int_vector<> LCS = get_kmer_lcs(A_bits, C_bits, G_bits, T_bits, k);
            save_v(LCS_file, LCS);
        }
        sdsl::int_vector<> LCS;
        load_v(LCS_file, LCS);
        std::cerr << "LCS_file loaded" << std::endl;
        n_finimizers += fmin_search(input_files, outfile, sbwt, DNA_rs, LCS,t, type, gzip_output);//DNA_bitvectors,

        write_log("total number of finimizers: " + to_string(n_finimizers), LogLevel::MAJOR);

        return 0;
}


