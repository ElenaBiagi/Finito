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

using namespace std;

using namespace sbwt;

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

// Here you are nto sure to find the interval as when buildign fmin
// First update Kmer interval then fmin
template<typename writer_t>
vector<int64_t> rarest_fmin_streaming_search( const sdsl::rank_support_v5<>** DNA_rs, const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input, const char t, const sdsl::rank_support_v5<>& fmin_rs, const  std::vector< uint64_t>& unitigs_v, writer_t& writer){ //const sdsl::bit_vector** DNA_bitvectors,
    //writer_t writer_fmin("rarest_fmin_ecoli.csv");
    const uint64_t n_nodes = sbwt.number_of_subsets();
    const uint64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();
    uint64_t freq;
    set<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> all_fmin;
    const uint64_t str_len = input.size();
    tuple<uint64_t, uint64_t, uint64_t, uint64_t> w_fmin = {n_nodes,k+1,n_nodes,str_len}; // {freq, len, I start, start}
    //vector<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> all_w_fmin;
    //set<tuple<uint64_t,uint64_t, uint64_t>> count_all_w_fmin;
    //vector<bool> found_kmers;
    //unitigs_v.reserve(str_len - k + 1);
    //unitigs_v.resize(str_len - k + 1);
    vector<int64_t> found_kmers(str_len - k + 1,-1); 

    // (1) Calculate all possible fmin in the FIRST window
    uint64_t start = 0;
    uint64_t end;
    uint64_t kmer_start = 0;
    bool found = true;
    pair<int64_t, int64_t> I = {0, n_nodes - 1};
    pair<int64_t, int64_t> I_new, I_kmer, I_kmer_new;
    uint64_t I_start;
    tuple<uint64_t, uint64_t, uint64_t, uint64_t> curr_substr;
    //string writer = "rarest_";
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
            //const sdsl::bit_vector& Bit_v = *(DNA_bitvectors[char_idx]);
            const sdsl::rank_support_v5<> &Bit_rs = *(DNA_rs[char_idx]);
            
            // 1) fmin interval
            I_new = update_sbwt_interval(char_idx, I, Bit_rs, C);
            freq = (I.second - I.first + 1);
            I_start = I.first;
            // (1) Finimizer(subseq) NOT found
            if (I_new.first == -1){
                kmer_start = ++start;
                I = drop_first_char(end - start + 1, I, LCS, n_nodes); // The result canno thave freq == 1
                I_kmer = I;
            } else{
                I = I_new;
            }
            // (2) Finimizer(subseq) freq > 0
            // Check if the Kmer interval has to be updated
            if ( start != kmer_start){
                I_kmer_new = update_sbwt_interval(char_idx, I_kmer, Bit_rs, C);
                if(I_kmer_new.first == -1){
                    // kmer NOT found
                    kmer_start++;
                    I_kmer = drop_first_char(end - kmer_start + 1, I_kmer, LCS, n_nodes);
                    //found = false;
                } else{
                    I_kmer = I_kmer_new;
                }
            } else { I_kmer = I; }
            // (2b) Finimizer found
            while (freq == 1) {
                curr_substr = {freq, end - start + 1, I_start, start};
                all_fmin.insert(curr_substr);
                // 1. rarest (freq=1), 2. shortest
                if (w_fmin > curr_substr) {w_fmin = curr_substr;}

                // 2. drop the first char
                // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                start ++;
                // todo check if it is necessary to ensure that start is always <= end or if it is so by construction as if start==end than the freq is >>1
                if (start > end)[[unlikely]] { break; }
                I = drop_first_char(end - start + 1, I, LCS, n_nodes);
                freq = (I.second - I.first + 1);
                I_start = I.first;
            }
        }
        // Check if the kmer is found
        if (end - kmer_start + 1 == k){
            //cerr << "Kmer found" << endl;
            while (get<3>(w_fmin) < kmer_start) { //else keep the current w_fmin
                all_fmin.erase(all_fmin.begin());
                w_fmin = *all_fmin.begin();
            }
            //found_kmers[kmers_start] == true;
            found_kmers[kmer_start] = unitigs_v[fmin_rs(get<2>(w_fmin))];
            kmer_start++;
            I_kmer = drop_first_char(end - kmer_start + 1, I_kmer, LCS, n_nodes);
        }
    }
    return found_kmers;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_fmin_queries_streaming(reader_t& reader, writer_t& writer, const sbwt_t& sbwt, const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const sdsl::int_vector<>& LCS, const sdsl::rank_support_v5<>& fmin_rs, const  std::vector< uint64_t>& unitigs_v, const char t){
    // sbwt is useless
    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    vector<int64_t> out_buffer;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;

        int64_t t0 = cur_time_micros();
        out_buffer = rarest_fmin_streaming_search( DNA_rs, sbwt, LCS, reader.read_buf, t, fmin_rs, unitigs_v, writer);
    
        //vector<int64_t> out_buffer = sbwt.streaming_search(first_kmer, len);
        //pair<std::string, int64_t> eulertig = unitigs_v[fmin_rs(w_fmin.first)];
        number_of_queries += out_buffer.size();
        print_vector(out_buffer, writer);
    }
    
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
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
int64_t run_fmin_file(const string& infile, const string& outfile, const sbwt_t& sbwt, const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const sdsl::int_vector<>& LCS, const sdsl::rank_support_v5<>& fmin_rs, const std::vector<uint64_t>& unitigs_v, const char t){
    reader_t reader(infile);
    writer_t writer(outfile);
    if(sbwt.has_streaming_query_support()){
        write_log("Running streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_fmin_queries_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt, DNA_bitvectors, DNA_rs, LCS, fmin_rs, unitigs_v, t);
    }
    else{
        write_log("Running non-streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_not_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
}

// Returns number of queries executed
template<typename sbwt_t>
int64_t run_fmin_queries(const vector<string>& infiles, const vector<string>& outfiles, const sbwt_t& sbwt, bool gzip_output, const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const sdsl::int_vector<>& LCS, const sdsl::rank_support_v5<>& fmin_rs, std::vector<uint64_t>& unitigs_v, const char t){

    if(infiles.size() != outfiles.size()){
        string count1 = to_string(infiles.size());
        string count2 = to_string(outfiles.size());
        throw std::runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
    }

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_no_gzip;

    typedef Buffered_ofstream<zstr::ofstream> out_gzip;
    typedef Buffered_ofstream<std::ofstream> out_no_gzip;

    int64_t n_queries_run = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = SeqIO::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input && gzip_output){
            n_queries_run += run_fmin_file<sbwt_t, in_gzip, out_gzip>(infiles[i], outfiles[i], sbwt, DNA_bitvectors, DNA_rs, LCS, fmin_rs, unitigs_v,t);
        }
        if(gzip_input && !gzip_output){
            n_queries_run += run_fmin_file<sbwt_t, in_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt, DNA_bitvectors, DNA_rs, LCS, fmin_rs, unitigs_v, t);
        }
        if(!gzip_input && gzip_output){
            n_queries_run += run_fmin_file<sbwt_t, in_no_gzip, out_gzip>(infiles[i], outfiles[i], sbwt, DNA_bitvectors, DNA_rs, LCS, fmin_rs, unitigs_v, t);
        }
        if(!gzip_input && !gzip_output){
            n_queries_run += run_fmin_file<sbwt_t, in_no_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt, DNA_bitvectors, DNA_rs, LCS, fmin_rs, unitigs_v, t);
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
        ("i,index-file", "Index input file.", cxxopts::value<string>())
        ("q,query-file", "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())
        ("z,gzip-output", "Writes output in gzipped form. This can shrink the output files by an order of magnitude.", cxxopts::value<bool>()->default_value("false"))
        ("t", "Maximum finimizer frequency", cxxopts::value<uint64_t>())
        ("type", "Decide which streaming search type you prefer. Available types: " + all_types_string,cxxopts::value<string>()->default_value("rarest"))
        ("lcs", "Provide in input the LCS file if available.", cxxopts::value<string>()->default_value(""))
        ("f, fmin_bv", "Provide in input the finimizers binary kmers vector.", cxxopts::value<string>()->default_value(""))
        ("unitigs-v", "Provide in input the eulertigs headers and offsets of the finimizers.", cxxopts::value<string>()->default_value(""))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    char t = opts["t"].as<uint64_t>();

    // TODO add type, only rarest now
    string type = opts["type"].as<string>();
    if(std::find(types.begin(), types.end(), type) == types.end()){
        cerr << "Error: unknown type: " << type << endl;
        cerr << "Available types are:" << all_types_string << endl;
        return 1;
    }

    // Interpret input file
    string queryfile = opts["query-file"].as<string>();
    vector<string> input_files;
    bool multi_file = queryfile.size() >= 4 && queryfile.substr(queryfile.size() - 4) == ".txt";
    if(multi_file){
        input_files = readlines(queryfile);
    } else{
        input_files = {queryfile};
    }
    for(string file : input_files) check_readable(file);

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
            const sdsl::int_vector<> LCS = get_kmer_lcs(A_bits, C_bits, G_bits, T_bits, k);
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

        string unitigs_v_file = opts["unitigs-v"].as<string>();
        std::vector<uint64_t> unitigs_v;
        load_intv(unitigs_v_file,unitigs_v);
        std::cerr<< "unitigs_v loaded"<<std::endl;

        number_of_queries += run_fmin_queries(input_files, output_files, sbwt, gzip_output, DNA_bitvectors, DNA_rs, LCS, fmin_rs, unitigs_v,t);
        int64_t new_total_micros = cur_time_micros() - micros_start;
        write_log("us/query end-to-end: " + to_string((double)new_total_micros / number_of_queries), LogLevel::MAJOR);
        write_log("total number of queries: " + to_string(number_of_queries), LogLevel::MAJOR);
    }

    int64_t total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)total_micros / number_of_queries), LogLevel::MAJOR);

    return 0;

}
