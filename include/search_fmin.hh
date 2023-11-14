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


// Here you are noT sure to find the interval as when building fmin
//template<typename writer_t>

// Old code below
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
                            //pair<int64_t, int64_t> I_branch = sbwt.update_sbwt_interval(&c, 1, {l_branch.second,l_branch.second});
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

template<typename reader_t, typename out_stream_t>
int64_t run_fmin_queries_streaming(reader_t& reader, out_stream_t& out, const FinimizerIndex& index, const string& stats_filename){
    const int64_t k = index.sbwt->get_k();
    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    int64_t kmers_count = 0 , kmers_count_rev = 0; // todo chek if we need all of these
    vector<int64_t> out_buffer, out_buffer_rev;

    int64_t query_seq=0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        int64_t t0 = cur_time_micros();
        //pair<vector<int64_t>, int64_t> final_pair = rarest_fmin_streaming_search(DNA_bitvectors, DNA_rs, sbwt, LCS, reader.read_buf, t, fmin_rs, global_offsets, ef_endpoints, Ustart_rs, found_kmers);
        FinimizerIndex::QueryResult result = index.search(reader.read_buf);

        //reverse compl
        const string reverse = sbwt::get_rc(reader.read_buf);
        FinimizerIndex::QueryResult r_result = index.search(reverse);
        int64_t tot_kmers = result.local_offsets.size();
        int64_t str_len = reverse.length(); // the string and its reverse complement have the same length
        for(int64_t i = 0; i < tot_kmers; i++){
            int64_t unitig, pos;
            if (result.local_offsets[i].first==-1) {
                std::tie(unitig,pos) = r_result.local_offsets[str_len-k-i];
            } else{
                std::tie(unitig,pos) = result.local_offsets[i];
            }
            if(i > 0) out << ' ';
            out << '(' << unitig << ',' << pos << ')';
        }
        out << '\n';

        number_of_queries += tot_kmers;//result.local_offsets.size();
        kmers_count += result.n_found;
        kmers_count_rev += r_result.n_found;
     
        total_micros += cur_time_micros() - t0;
    }
    write_log("k " + to_string(k), LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    write_log("Found kmers: " + to_string(kmers_count), LogLevel::MAJOR);
    write_log("Found kmers reverse : " + to_string(kmers_count_rev), LogLevel::MAJOR);
    write_log("Total found kmers: " + to_string(kmers_count+kmers_count_rev), LogLevel::MAJOR);

    std::ofstream statsfile;
    statsfile.open(stats_filename, std::ios_base::app); // append instead of overwrite
    statsfile << to_string(k) + "," + to_string(kmers_count+kmers_count_rev) + "," + to_string(number_of_queries);
    statsfile.close();
    return number_of_queries;
}

template<typename reader_t, typename out_stream_t>
int64_t run_fmin_file(const string& infile, out_stream_t& out, const string& stats_filename, const FinimizerIndex& index){
    reader_t reader(infile);
    write_log("Running streaming queries from input file " + infile, LogLevel::MAJOR);
    return run_fmin_queries_streaming(reader, out, index, stats_filename);
}

// Returns number of queries executed
int64_t run_fmin_queries(const vector<string>& infiles, const optional<vector<string>>& outfiles, const string& stats_filename, const FinimizerIndex& index){

    if(outfiles.has_value()){
        if(infiles.size() != outfiles.value().size()){
            string count1 = to_string(infiles.size());
            string count2 = to_string(outfiles.value().size());
            throw std::runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
        }
    }

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_no_gzip;

    int64_t n_queries_run = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = SeqIO::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input){
            if(outfiles.has_value()){
                ofstream out(outfiles.value()[i]);
                n_queries_run += run_fmin_file<in_gzip>(infiles[i], out, stats_filename, index);
            } else { // To stdout
                n_queries_run += run_fmin_file<in_gzip>(infiles[i], cout, stats_filename, index);
            }
        }
        else {
            if(outfiles.has_value()){
                ofstream out(outfiles.value()[i]);
                n_queries_run += run_fmin_file<in_no_gzip>(infiles[i], out, stats_filename, index);
            } else{ // To stdout
                n_queries_run += run_fmin_file<in_no_gzip>(infiles[i], cout, stats_filename, index);
            }
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
        ("o,out-file", "Output filename, or stdout if not given.", cxxopts::value<string>())
        ("i,index-file", "Index filename prefix.", cxxopts::value<string>())
        ("q,query-file", "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())
        ("type", "Decide which streaming search type you prefer. Available types: " + all_types_string,cxxopts::value<string>()->default_value("rarest"))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

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
    optional<vector<string>> output_files;
    try{
        string outfile = opts["out-file"].as<string>();
        if(multi_file){
            output_files = readlines(outfile);
        } else{
            output_files = {outfile};
        }
        for(string file : output_files.value()) check_writable(file);
    } catch(cxxopts::option_has_no_value_exception& e){
        write_log("No output file given, writing to stdout", LogLevel::MAJOR);
    }

    string index_prefix = opts["index-file"].as<string>();

    int64_t number_of_queries = 0;

    cerr << "Loading index..." << endl;
    FinimizerIndex index;
    index.load(index_prefix);
    cerr << "Index loaded" << endl;

    const int64_t k = index.sbwt->get_k();
    cerr << "k = "<< to_string(k);
    cerr << " SBWT nodes: "<< to_string(index.sbwt->number_of_subsets())<< " kmers: "<< to_string(index.sbwt->number_of_kmers())<< endl;

    number_of_queries += run_fmin_queries(query_files, output_files, index_prefix + ".stats", index);
    int64_t new_total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)new_total_micros / number_of_queries), LogLevel::MAJOR);
    write_log("total number of queries: " + to_string(number_of_queries), LogLevel::MAJOR);
    
    std::ofstream statsfile2;
    statsfile2.open(index_prefix + "stats.txt", std::ios_base::app); // append instead of overwrite
    string results = to_string(number_of_queries);
    statsfile2 << "," + to_string((double)new_total_micros / number_of_queries);
    
    int64_t bytes = index.size_in_bytes();
    write_log("bytes: " + to_string(bytes), LogLevel::MAJOR);
    statsfile2 <<  "," + to_string(bytes);
    statsfile2 << "," + to_string(static_cast<double>(bytes*8)/index.sbwt->number_of_kmers()) + "\n";
    statsfile2 << "," + to_string(index.sbwt->number_of_kmers()) + "\n";
    statsfile2.close();

    int64_t total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)total_micros / number_of_queries), LogLevel::MAJOR);
    
    return 0;

}