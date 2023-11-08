mod minimizer_index;

use std::path::PathBuf;

use clap::ArgAction;
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

// Map unitig ids and positions to global offsets in the concatenation of unitigs
// in the order of the colex ranks of their first k-mers
struct MapToGlobalOffset{
    unitig_id_to_unitig_rank: Vec<usize>,
    unitig_rank_to_unitig_id: Vec<usize>,
    permuted_unitig_cumulative_lengths: Vec<usize>,
}

impl MapToGlobalOffset{

    fn new(unitig_db: &SeqDB, k: usize) -> Self{

        let mut unitig_rank_to_unitig_id: Vec<usize> = (0..unitig_db.sequence_count()).collect();
        // Sort by the colex order of the first k-mer
        unitig_rank_to_unitig_id.sort_by(|a, b|{
            colex_compare(&unitig_db.get(*a).seq()[0..k], &unitig_db.get(*b).seq()[0..k])
        });

        let mut unitig_id_to_unitig_rank: Vec<usize> = vec![0; unitig_db.sequence_count()];
        for i in 0..unitig_rank_to_unitig_id.len(){
            unitig_id_to_unitig_rank[unitig_rank_to_unitig_id[i]] = i;
        }

        // Permuted cumulative sum
        // Index i contains the length of unitig i in input order, plus the sums of lengths of all unitigs
        // whose first k-mer is colexicographically smaller than the first k-mer of unitig i
        let mut permuted_unitig_cumulative_lengths: Vec<usize> = vec![0; unitig_db.sequence_count()];
        for i in 0..permuted_unitig_cumulative_lengths.len(){
            if i == 0{
                permuted_unitig_cumulative_lengths[unitig_rank_to_unitig_id[i]] = unitig_db.get(unitig_rank_to_unitig_id[i]).seq().len();
            } else{
                permuted_unitig_cumulative_lengths[unitig_rank_to_unitig_id[i]] = permuted_unitig_cumulative_lengths[unitig_rank_to_unitig_id[i-1]] + unitig_db.get(unitig_rank_to_unitig_id[i]).seq().len();
            }
        }

        MapToGlobalOffset{unitig_id_to_unitig_rank, unitig_rank_to_unitig_id, permuted_unitig_cumulative_lengths}

    }

    fn to_global_offset(&self, unitig_id: usize, unitig_pos: usize) -> usize{
        let unitig_rank = self.unitig_id_to_unitig_rank[unitig_id];
        if unitig_rank > 0{
            unitig_pos + self.permuted_unitig_cumulative_lengths[self.unitig_rank_to_unitig_id[unitig_rank-1]]
        }
        else{
            unitig_pos
        }
    }

    fn to_unitig_rank(&self, unitig_id: usize) -> usize{
        self.unitig_id_to_unitig_rank[unitig_id]
    }
        

}

fn main() {
    if std::env::var("RUST_LOG").is_err(){
        std::env::set_var("RUST_LOG", "info");
    }

    env_logger::init();

    let cli = Command::new("kmer-mapper")
        .about("Mapping k-mers to unitigs")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .arg(Arg::new("unitigs")
            .help("Input FASTA or FASTQ file, possibly gzipped")
            .long("unitigs")
            .short('u')
            .required(true)
            .value_parser(clap::value_parser!(PathBuf))
        )
        .arg(Arg::new("query")
            .help("Input FASTA or FASTQ file, possibly gzipped")
            .long("query")
            .short('q')
            .required(true)
            .value_parser(clap::value_parser!(PathBuf))
        ).arg(Arg::new("nthreads")
            .help("Number of threads")
            .long("nthreads")
            .short('t')
            .default_value("1")
            .value_parser(clap::value_parser!(usize))
        ).arg(Arg::new("k")
            .help("k-mer length")
            .short('k')
            .value_parser(clap::value_parser!(usize))
            .required(true)
        ).arg(Arg::new("m")
            .help("minimizer length (default max(1, k-6))")
            .short('m')
            .value_parser(clap::value_parser!(usize))
    );
    let cli_matches = cli.get_matches();

    let n_threads = *cli_matches.get_one::<usize>("nthreads").unwrap();
    rayon::ThreadPoolBuilder::new().num_threads(n_threads).build_global().unwrap();

    let unitigs_path: &PathBuf = cli_matches.get_one("unitigs").unwrap();
    let query_path: &PathBuf = cli_matches.get_one("query").unwrap();
    let k: usize = *cli_matches.get_one("k").unwrap();

    let m: Option<&usize> = cli_matches.get_one("m");
    let m = match m{
        Some(m) => *m,
        None => std::cmp::max(1, k as isize - 6) as usize
    };

    log::info!("k = {}, m = {}", k, m);

    let unitig_db = DynamicFastXReader::from_file(&unitigs_path).unwrap().into_db().unwrap(); 
    let unitig_db = std::sync::Arc::new(unitig_db);
    let index = minimizer_index::MinimizerIndex::new(unitig_db.clone(), k, m);

    let map_to_global = MapToGlobalOffset::new(unitig_db.as_ref(), k);

    let mut query_reader = DynamicFastXReader::from_file(&query_path).unwrap();
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
                answers.push((map_to_global.to_unitig_rank(*unitig_id) as isize, *unitig_pos as isize));
            }
        } 
        // Build the output line
        let out_line = answers.iter().map(|x| format!("({},{})",x.0, x.1)).collect::<Vec<String>>().join(" ");
        println!("{}", out_line);
    }

}
