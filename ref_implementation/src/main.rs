mod minimizer_index;

use std::io::Read;
use std::io::Write;
use std::path::PathBuf;
use jseqio::reverse_complement;
use minimizer_index::{MinimizerIndex, Kmer};

use clap::ArgAction;
use clap::Subcommand;
use jseqio::reader::*;
use jseqio::writer::*;
use jseqio::record::*;
use clap::{Command, Arg};

use log::{info, error};

use jseqio::seq_db::SeqDB;

// The file format magic number is "KMIDXv01"
const MAGIC_NUMBER: [u8; 8] = [b'K', b'M', b'I', b'D', b'X', b'v', b'0', b'1'];

fn colex_compare(a: &[u8], b: &[u8]) -> std::cmp::Ordering{
    a.iter().rev().cmp(b.iter().rev())
}

#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn test_colex_compare(){
        assert_eq!(colex_compare(b"AAAAX", b"AAAAX"), std::cmp::Ordering::Equal);
        assert_eq!(colex_compare(b"AAAAX", b"AAABX"), std::cmp::Ordering::Less);
        assert_eq!(colex_compare(b"AAABX", b"AAAAX"), std::cmp::Ordering::Greater);
        assert_eq!(colex_compare(b"AAAAA", b"AA"), std::cmp::Ordering::Greater);
        assert_eq!(colex_compare(b"AA", b"AAA"), std::cmp::Ordering::Less);
        assert_eq!(colex_compare(b"", b"AAA"), std::cmp::Ordering::Less);
    }
}

fn get_unitig_rank_to_id(db: &SeqDB, k: usize) -> Vec<usize> {

    let mut unitig_rank_to_unitig_id: Vec<usize> = (0..db.sequence_count()).collect();
    // Sort by the colex order of the first k-mer
    unitig_rank_to_unitig_id.sort_by(|a, b|{
        colex_compare(&db.get(*a).seq()[0..k], &db.get(*b).seq()[0..k])
    });

    unitig_rank_to_unitig_id
}

fn get_permuted_unitig_db(db: SeqDB, k: usize) -> SeqDB{
    let unitig_rank_to_unitig_id = get_unitig_rank_to_id(&db, k);
    let mut new_db = SeqDB::new();
    for unitig_id in unitig_rank_to_unitig_id{
        let rec = db.get(unitig_id);
        new_db.push_record(rec);
    }
    new_db
}

fn build_and_store_index<const KmerWidth64BitWords : usize>(unitig_db: std::sync::Arc<SeqDB>, k: usize, m: usize, out_path: &PathBuf){

    let index = MinimizerIndex::<Kmer<KmerWidth64BitWords>>::new(unitig_db.clone(), k, m);

    eprintln!("Saving index to {}", out_path.display());
    let mut out_file = std::fs::File::create(out_path).unwrap();
    out_file.write_all(&MAGIC_NUMBER).unwrap();
    out_file.write_all(&KmerWidth64BitWords.to_le_bytes()).unwrap();
    index.serialize(out_file);
}

fn run_query<const KmerWidth64BitWords : usize>(mut index_in: std::fs::File, queryfile: &PathBuf, rc: bool){
    let index = MinimizerIndex::<Kmer<KmerWidth64BitWords>>::new_from_serialized(index_in);
    let k = index.get_k();
    let mut query_reader = DynamicFastXReader::from_file(&queryfile).unwrap();
    while let Some(query) = query_reader.read_next().unwrap(){
        let mut answers = Vec::<(isize, isize)>::new(); // End points in unitigs. -1 means does not exist
        for kmer in query.seq.windows(k){
            let mut occurrences = index.lookup_kmer(kmer);
            if rc {
                let rc_kmer = reverse_complement(kmer);
                if rc_kmer != kmer { // Don't re-search self-rc kmers
                    let rc_occs = index.lookup_kmer(&rc_kmer);
                    occurrences.extend_from_slice(rc_occs.as_slice());
                }
            }
            if occurrences.len() > 1{
                eprintln!("Error: k-mer {} occurs in {} unitigs", String::from_utf8_lossy(kmer), occurrences.len());
                std::process::exit(1);
            }
            else if occurrences.is_empty(){
                answers.push((-1,-1));
            }
            else { // Single occurrence
                let (unitig_id, unitig_pos) = occurrences.first().unwrap();
                answers.push((*unitig_id as isize, *unitig_pos as isize));
            }
        } 
        // Build the output line
        let out_line = answers.iter().map(|x| format!("({},{})",x.0, x.1)).collect::<Vec<String>>().join(" ");
        println!("{}", out_line);
    }
}

fn extract_index_unitigs<const KmerWidth64BitWords : usize>(indexfile: &PathBuf, outfile: &PathBuf){
    let index = MinimizerIndex::<Kmer::<KmerWidth64BitWords>>::new_from_serialized(std::fs::File::open(indexfile).unwrap());
    let mut writer = DynamicFastXWriter::new_to_file(outfile).unwrap();
    for record in index.get_unitig_db().iter(){
        writer.write(&record).unwrap();
    }
}

// Reads and checks the magic number and returns the k-mer width
fn read_index_header(index_in: &mut std::fs::File) -> usize{
    // Read and check the magic number
    let magic_number_buf = &mut [0u8; 8];
    index_in.read_exact(&mut magic_number_buf[..]).unwrap();
    if magic_number_buf != &MAGIC_NUMBER[..]{
        eprintln!("Error: Index file does not start with the header bytes {:?}", MAGIC_NUMBER);
        std::process::exit(1);
    }

    // Read k-mer width
    let mut kmer_width_64_bit_words_buf = 0_usize.to_le_bytes();
    index_in.read_exact(&mut kmer_width_64_bit_words_buf).unwrap();

    usize::from_le_bytes(kmer_width_64_bit_words_buf)
}


fn main() {
    if std::env::var("RUST_LOG").is_err(){
        std::env::set_var("RUST_LOG", "info");
    }

    env_logger::init();

    let cli = Command::new("kmer-mapper")
        .about("Mapping k-mers to unitigs")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .arg_required_else_help(true)
        .subcommand(Command::new("build")
            .arg_required_else_help(true)
            .arg(Arg::new("unitigs")
                .help("Input FASTA or FASTQ file, possibly gzipped")
                .long("unitigs")
                .short('u')
                .required(true)
                .value_parser(clap::value_parser!(PathBuf))
            ).arg(Arg::new("nthreads")
                .help("Number of threads")
                .long("nthreads")
                .short('t')
                .default_value("1")
                .value_parser(clap::value_parser!(usize))
            ).arg(Arg::new("outfile")
                .help("Output index file")
                .long("outfile")
                .short('o')
                .value_parser(clap::value_parser!(PathBuf))
            ).arg(Arg::new("k")
                .help("k-mer length")
                .short('k')
                .value_parser(clap::value_parser!(usize))
                .required(true)
            ).arg(Arg::new("m")
                .help("minimizer length (default max(1, k-6))")
                .short('m')
                .value_parser(clap::value_parser!(usize))))
        .subcommand(Command::new("query")
            .arg_required_else_help(true)
            .arg(Arg::new("index")
                .help("Index file")
                .long("index")
                .short('i')
                .required(true)
                .value_parser(clap::value_parser!(PathBuf))
            ).arg(Arg::new("query")
                .help("Input FASTA or FASTQ file, possibly gzipped")
                .long("query")
                .short('q')
                .required(true)
                .value_parser(clap::value_parser!(PathBuf))
            ).arg(Arg::new("reverse-complements")
                .help("Whether to also report reverse complement matches")
                .long("reverse-complements")
                .short('r')
                .action(ArgAction::SetTrue)))
        .subcommand(Command::new("extract-index-unitigs")
            .arg_required_else_help(true)
            .arg(Arg::new("index")
                .help("Index file")
                .long("index")
                .short('i')
                .required(true)
                .value_parser(clap::value_parser!(PathBuf))
            ).arg(Arg::new("outfile")
                .help("Output fasta file")
                .long("outfile")
                .short('o')
                .required(true)
                .value_parser(clap::value_parser!(PathBuf))))
    ;

    // Parse subcommand matches
    match cli.get_matches().subcommand(){
        Some(("build", cli_matches)) => {
            let n_threads = *cli_matches.get_one::<usize>("nthreads").unwrap();
            rayon::ThreadPoolBuilder::new().num_threads(n_threads).build_global().unwrap();

            let unitigs_path: &PathBuf = cli_matches.get_one("unitigs").unwrap();
            let out_path: &PathBuf = cli_matches.get_one("outfile").unwrap();
            let k: usize = *cli_matches.get_one("k").unwrap();

            let m: Option<&usize> = cli_matches.get_one("m");
            let m = match m{
                Some(m) => *m,
                None => std::cmp::max(1, k as isize - 6) as usize
            };

            log::info!("k = {}, m = {}", k, m);

            log::info!("Reading unitigs from {}", unitigs_path.display());
            let unitig_db = DynamicFastXReader::from_file(&unitigs_path).unwrap().into_db().unwrap(); 
            log::info!("Read {} unitigs", unitig_db.sequence_count());
            log::info!("Sorting unitigs by first k-mer");
            let unitig_db = get_permuted_unitig_db(unitig_db, k);
            let unitig_db = std::sync::Arc::new(unitig_db);
            log::info!("Building index");
            let kmer_width_64bit_words = (k + 31) / 32; // Two bits per nucleotide

            match kmer_width_64bit_words {
                1 => build_and_store_index::<1>(unitig_db, k, m, out_path),
                2 => build_and_store_index::<2>(unitig_db, k, m, out_path),
                3 => build_and_store_index::<3>(unitig_db, k, m, out_path),
                4 => build_and_store_index::<4>(unitig_db, k, m, out_path),
                5 => build_and_store_index::<5>(unitig_db, k, m, out_path),
                6 => build_and_store_index::<6>(unitig_db, k, m, out_path),
                7 => build_and_store_index::<7>(unitig_db, k, m, out_path),
                8 => build_and_store_index::<8>(unitig_db, k, m, out_path),
                _ => {
                    eprintln!("Error: k-mer length {} is too long", k);
                    std::process::exit(1);
                }
            }

        },
        Some(("query", cli_matches)) => {
            let indexfile: &PathBuf = cli_matches.get_one("index").unwrap();
            let queryfile: &PathBuf = cli_matches.get_one("query").unwrap();
            let rc = cli_matches.get_flag("reverse-complements");

            let mut index_in = std::fs::File::open(indexfile).unwrap();

            let kmer_width_64_bit_words = read_index_header(&mut index_in);

            // Run queries
            match kmer_width_64_bit_words {
                1 => run_query::<1>(index_in, queryfile, rc),
                2 => run_query::<2>(index_in, queryfile, rc),
                3 => run_query::<3>(index_in, queryfile, rc),
                4 => run_query::<4>(index_in, queryfile, rc),
                5 => run_query::<5>(index_in, queryfile, rc),
                6 => run_query::<6>(index_in, queryfile, rc),
                7 => run_query::<7>(index_in, queryfile, rc),
                8 => run_query::<8>(index_in, queryfile, rc),
                _ => {
                    eprintln!("Error: invalid k-mer word width {} in index file {}", kmer_width_64_bit_words, indexfile.display());
                    std::process::exit(1);
                }
            }

        },
        Some(("extract-index-unitigs", cli_matches)) => {
            let indexfile: &PathBuf = cli_matches.get_one("index").unwrap();
            let outfile: &PathBuf = cli_matches.get_one("outfile").unwrap();
            let mut index_in = std::fs::File::open(indexfile).unwrap();
            let kmer_width_64_bit_words = read_index_header(&mut index_in);
            match kmer_width_64_bit_words {
                1 => extract_index_unitigs::<1>(indexfile, outfile),
                2 => extract_index_unitigs::<2>(indexfile, outfile),
                3 => extract_index_unitigs::<3>(indexfile, outfile),
                4 => extract_index_unitigs::<4>(indexfile, outfile),
                5 => extract_index_unitigs::<5>(indexfile, outfile),
                6 => extract_index_unitigs::<6>(indexfile, outfile),
                7 => extract_index_unitigs::<7>(indexfile, outfile),
                8 => extract_index_unitigs::<8>(indexfile, outfile),
                _ => {
                    eprintln!("Error: invalid k-mer word width {} in index file {}", kmer_width_64_bit_words, indexfile.display());
                    std::process::exit(1);
                }
            }

        },
        _ => {
            eprintln!("Error: No subcommand given");
            std::process::exit(1);
        }
    }



}
