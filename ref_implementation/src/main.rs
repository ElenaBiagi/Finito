mod minimizer_index;

use std::path::PathBuf;
use minimizer_index::MinimizerIndex;

use clap::ArgAction;
use clap::Subcommand;
use jseqio::reader::*;
use jseqio::writer::*;
use jseqio::record::*;
use clap::{Command, Arg};

use log::{info, error};

use jseqio::seq_db::SeqDB;

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


fn main() {
    if std::env::var("RUST_LOG").is_err(){
        std::env::set_var("RUST_LOG", "info");
    }

    env_logger::init();
/*
            )
            .arg(Arg::new("query")
                .help("Input FASTA or FASTQ file, possibly gzipped")
                .long("query")
                .short('q')
                .required(true)
                .value_parser(clap::value_parser!(PathBuf))
            )
 */

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
            let index = minimizer_index::MinimizerIndex::new(unitig_db.clone(), k, m);

            eprintln!("Saving index to {}", out_path.display());
            index.serialize(std::fs::File::create(out_path).unwrap());

        },
        Some(("query", cli_matches)) => {
            let indexfile: &PathBuf = cli_matches.get_one("index").unwrap();
            let queryfile: &PathBuf = cli_matches.get_one("query").unwrap();
            let index = MinimizerIndex::new_from_serialized(std::fs::File::open(indexfile).unwrap());
            let k = index.get_k();
            let mut query_reader = DynamicFastXReader::from_file(&queryfile).unwrap();
            while let Some(query) = query_reader.read_next().unwrap(){
                let mut answers = Vec::<(isize, isize)>::new(); // End points in unitigs. -1 means does not exist
                for kmer in query.seq.windows(k){
                    let occurrences = index.lookup_kmer(kmer);
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
            std::process::exit(1);
        },
        _ => {
            eprintln!("Error: No subcommand given");
            std::process::exit(1);
        }
    }



}
