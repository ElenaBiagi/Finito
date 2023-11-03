use rayon::prelude::*;
use jseqio::seq_db::SeqDB;

pub struct MinimizerIndex<'a>{
    seq_storage: &'a jseqio::seq_db::SeqDB,
    mphf: boomphf::Mphf<Kmer>, // Minimal perfect hash function
    locations: Vec<(u32, u32)>,
    bucket_starts: Vec<usize>,
    k: usize, // k-mer length
    m: usize, // Minimizer length
    n_mmers: usize, // Number of distinct m-mers stored in the mphf
}

// Keeps only parts of length at least k
#[allow(dead_code)]
fn split_at_non_acgt(v: &[u8], k: usize) -> Vec<Vec<u8>>{
    let mut parts = Vec::<Vec::<u8>>::new();
    let mut start = 0_usize;
    for i in 0 .. v.len() + 1 {
        if i == v.len() || !is_dna(v[i]) { // One-past-the-end of a part
            if i - start >= k {
                parts.push(v[start..i].to_owned());
            }
            start = i+1;
        }
    }
    parts
}

// TODO: UPPER-CASE SEQUENCES
fn is_dna(c: u8) -> bool{
    if c == b'a' || c == b'c' || c == b'g' || c == b't' {
        panic!("ERROR: lower case nucleotides found");
    }
    c == b'A' || c == b'C' || c == b'G' || c == b'T'
}

fn get_minimizer_position(kmer: &[u8], m: usize) -> usize{
    let mut minimizer = &kmer[0 .. m];
    let mut min_pos = 0;
    for j in 1 .. (kmer.len() as i64) - (m as i64) + 1 {
        let j = j as usize;
        if kmer[j..j+m] < *minimizer {
            minimizer = &kmer[j..j+m];
            min_pos = j;
        }
    }
    min_pos
}

// The output will be stored to the positions-vector
// Will not return minimizers for k-mers with non-DNA characters
fn get_minimizer_positions(seq: &[u8], positions: &mut Vec<usize>, k: usize, m: usize){

    positions.clear();

    // Store minimizer mappings
    for i in 0 .. (seq.len() as i64) - (k as i64) + 1 {
        let i = i as usize;
        let kmer = &seq[i..i+k];
        if kmer.iter().any(|&c| !is_dna(c)){
            continue;
        }

        let min_pos = i + get_minimizer_position(kmer, m);

        if positions.is_empty() || positions[positions.len()-1] != min_pos {
            positions.push(min_pos)
        }

    }
}

// The output will be stored to the positions-vector
fn get_minimizer_positions_with_return(seq: &[u8], k: usize, m: usize) -> Vec<usize>{

    let mut positions = Vec::<usize>::new();
    get_minimizer_positions(seq, &mut positions, k, m);
    positions
}

#[derive(Copy, Clone, PartialEq, Eq, Ord, PartialOrd, Hash, Debug)]
struct Kmer{
    data: u64
}

#[derive(Debug)]
enum KmerEncodingError{
    InvalidNucleotide(char), // contains the offending char
    TooLong(usize), // Contains the length of the k-mer which was too long
}

impl Kmer{
    fn from_ascii(ascii: &[u8]) -> Result<Self, KmerEncodingError>{
        if ascii.len() > 32{
            return Err(KmerEncodingError::TooLong(ascii.len()));
        }
        let mut data = 0_u64;
        for c in ascii.iter() {
            data <<= 2;
            data |= match *c{
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => {return Err(KmerEncodingError::InvalidNucleotide(*c as char))}
            };
        }
        Ok(Self{data})
    }
}

impl<'a> MinimizerIndex<'a>{

    // Returns all occurrences of the query k-mer
    pub fn lookup_kmer(&self, kmer: &[u8]) -> Vec<(usize, usize)>{
        assert!(kmer.len() == self.k);
        let min_pos = get_minimizer_position(kmer, self.m);
        let minimizer = &kmer[min_pos .. min_pos + self.m];

        let mut ans: Vec<(usize,usize)> = vec![];
        let minmer = match Kmer::from_ascii(minimizer){
            Err(KmerEncodingError::InvalidNucleotide(c)) => {return vec![]}, // No matches
            Err(KmerEncodingError::TooLong(len)) => {panic!("Sequence length {} shorter than k", len)},
            Ok(x) => x
        };

        if let Some(bucket) = self.mphf.try_hash(&minmer){
            let bucket_range = self.bucket_starts[bucket as usize]..self.bucket_starts[bucket as usize + 1];
            for (seq_id, seq_pos) in self.locations[bucket_range].iter(){
                
                // Start of the k-mer that contains this minimizer occurrence:
                let start = *seq_pos as i64 - min_pos as i64;

                // Check if this occurrence is real
                if start >= 0 && start + self.k as i64 <= self.seq_storage.get(*seq_id as usize).seq.len() as i64 {
                    // k-mer is within bounds of the sequence
                    let start = start as usize;
                    let candidate = &self.seq_storage.get(*seq_id as usize).seq[start .. start + self.k];
                    if candidate == kmer{
                        ans.push((*seq_id as usize, start));
                    }
                }
            }
        }

        ans
    }

    // Searches for all k-mer of the query and 
    // returns pairs (i,j) such that seq might be found in sequence i starting from position j
    pub fn get_exact_alignment_candidates(&self, query: &[u8]) -> Vec<(usize,usize)>{
        let mut align_starts = std::collections::HashSet::<(usize,usize)>::new(); // Pairs (seq_id, seq_pos)
        for (query_pos, kmer) in query.windows(self.k).enumerate(){
            for (target_id, target_pos) in self.lookup_kmer(kmer){
                
                let align_start = target_pos as isize - query_pos as isize;
                if align_start >= 0 && (align_start + query.len() as isize) <= self.seq_storage.get(target_id).seq.len() as isize{
                    // Is within bounds
                    align_starts.insert((target_id, align_start as usize));
                }
            }
        }
        return align_starts.iter().map(|&(i,j)| (i,j)).collect();
    }

    pub fn new(db: &'a SeqDB, k: usize, m: usize) -> Self{
        build::build(db, k, m)
    }

}


mod build{

    use super::*;

    // Returns the longest prefix of the slice that consists of elements that are 
    // all the same according to the given equivalence relation predicate
    fn get_prefix_run<T, F: Fn(&T, &T) -> bool>(slice: &[T], is_same: F) -> &[T]{
        let first_different = slice.iter().position(|x| !is_same(x, &slice[0]));
        match first_different{
            Some(i) => &slice[0..i],
            None => slice
        }
    }

    // L must be sorted
    fn get_bucket_sizes(L: &Vec::<(Kmer, u32, u32)>, h: &boomphf::Mphf<Kmer>, n_minimizers: usize) -> Vec<usize>{
        let mut bucket_sizes: Vec::<usize> = vec![0; n_minimizers]; // Bucket sizes in left-to-right order of buckets

        // Find out the bucket sizes
        let mut L_tail = &L[..];
        while !L_tail.is_empty(){
            let prefix = get_prefix_run(L_tail, |(minmer, _, _), (minmer2, _, _)| minmer == minmer2);
            bucket_sizes[h.hash(&prefix[0].0) as usize] = prefix.len(); // Prefix is nonempty because L_tail was nonempty
            L_tail = L_tail.split_at(prefix.len()).1;
        }

        bucket_sizes
    }

    fn get_bucket_starts(bucket_sizes: &Vec<usize>) -> Vec<usize>{
        // Get the starting positions of buckets
        let mut sum = 0_usize;
        let mut bucket_starts = vec![0];
        for b in bucket_sizes.iter(){
            sum += b;
            bucket_starts.push(sum)
        }
        bucket_starts.shrink_to_fit();
        bucket_starts
    }

    // L must be sorted
    fn store_location_pairs(L: Vec::<(Kmer, u32, u32)>, h: &boomphf::Mphf<Kmer>, bucket_starts: &Vec<usize>) -> Vec<(u32, u32)>{
        let mut locations: Vec::<(u32, u32)> = vec![(0,0); *bucket_starts.last().unwrap()]; // Will have an end sentinel

        let mut L_tail = &L[..];
        while !L_tail.is_empty(){
            let prefix = get_prefix_run(L_tail, |(minmer, _, _), (minmer2, _, _)| minmer == minmer2);
            let bucket_id = h.hash(&prefix[0].0); // Prefix is nonempty because L_tail was nonempty
            let start = bucket_starts[bucket_id as usize];
            for i in 0 .. prefix.len(){
                locations[start + i] = (prefix[i].1, prefix[i].2);
            }
            L_tail = L_tail.split_at(prefix.len()).1;
        }
        locations
    }

    // Returns a list of tuples (minmer, seq_id, pos), where 'pos' is the starting position
    // of 'minmer' in sequence with id 'seq_id'.
    fn build_position_list(db: & SeqDB, k: usize, m: usize) -> Vec<(Kmer, u32, u32)>{

        let parts = (0..db.sequence_count()).into_par_iter()
            // Find out minimizer positions in this sequence
            .map(|i|{ 
                (i, get_minimizer_positions_with_return(db.get(i).seq, k, m))
            })
            // Expand the positions into tuples (minmer, seq_id, seq_pos)
            .map(|(i, pos_list)|{
                let mut tuples = Vec::<(Kmer, u32, u32)>::new();
                for p in pos_list.iter(){
                    let minmer = Kmer::from_ascii(&db.get(i).seq[*p..*p+m]).unwrap();
                    tuples.push((minmer, i as u32, *p as u32));
                }
                tuples
            })
            // Concatenate lists together into just a few lists with a parallel fold
            .fold(Vec::new, |mut x, y| {
                x.extend(y); x
            })
            // Collect the lists from the fold
            .collect::<Vec::<Vec::<(Kmer, u32, u32)>>>();

        // Concatenate together the lists from the fold
        parts.into_iter().fold(Vec::new(), |mut x, y| {
            x.extend(y); x
        })   
    }

    fn collect_distinct_minmers_from_sorted_list(position_list: &Vec<(Kmer, u32, u32)>) -> Vec<Kmer>{
        let mut minimizer_list: Vec<Kmer> = vec![];
        for (minmer, _, _) in position_list.iter(){
            match minimizer_list.last(){
                Some(last) => {
                    if minmer != last{
                        minimizer_list.push(*minmer);
                    };
                },
                None => { minimizer_list.push(*minmer); }
            }
        }
        minimizer_list
    }

    fn compress_sorted_position_list(L: Vec::<(Kmer, u32, u32)>, h: &boomphf::Mphf<Kmer>, n_minimizers: usize) -> (Vec<(u32, u32)>, Vec<usize>){

        log::info!("Computing bucket sizes");
        let bucket_sizes = get_bucket_sizes(&L, h, n_minimizers);
        log::info!("Computing bucket starts");
        let bucket_starts = get_bucket_starts(&bucket_sizes);
        log::info!("Storing location pairs to buckets");
        let locations = store_location_pairs(L, h, &bucket_starts);

        (locations, bucket_starts)
    }


    pub fn build(db: &SeqDB, k: usize, m: usize) -> MinimizerIndex{

        log::info!("Extracting minimizers");
        let mut position_list = build_position_list(db, k, m);

        log::info!("Sorting tuples (minimizer, seq_id, seq_pos)");
        position_list.par_sort_unstable();

        log::info!("Collecting distinct minimizers");
        let mut minimizer_list = collect_distinct_minmers_from_sorted_list(&position_list);

        log::info!("Removing duplicate minimizers");
        minimizer_list.dedup();
        minimizer_list.shrink_to_fit();

        log::info!("Found {} distinct minimizers", minimizer_list.len());

        log::info!("Building an MPHF for the minimizers");
        let n_mmers = minimizer_list.len();
        let mphf = boomphf::Mphf::<Kmer>::new_parallel(1.7, minimizer_list.as_slice(), None);
        drop(minimizer_list);
        
        log::info!("Compressing position lists");
        let (locations, bucket_starts) = compress_sorted_position_list(position_list, &mphf, n_mmers);

        log::info!("Stored {} location pairs", locations.len());
        log::info!("Average bucket size: {:.2}", locations.len() as f64 / (bucket_starts.len() as f64 - 1.0)); // -1 because of end sentinel
    
        MinimizerIndex{seq_storage: db, mphf, locations, bucket_starts, k, m, n_mmers}
    }

}


#[cfg(test)]
mod tests{

    use super::*;
    use jseqio::record::*;
    use rand_chacha::{self, rand_core::{SeedableRng, RngCore}};

    fn number_of_kmers(seq_len: usize, k: usize) -> usize { // TODO: use everywhere
        std::cmp::max(0, (seq_len as i64) - (k as i64) + 1) as usize
    }

    fn to_ascii(S: &[u8]) -> String{
        std::str::from_utf8(S).unwrap().to_owned()
    }

    #[test]
    fn test_get_minimizer_positions(){
        let seq = "ATAGCTAGTCGATGCTGATCGTAGGTTCGTAGCTGTATGCTGACCCTGATGTCTGTAGTCGTGACTGACT";
        let k: usize = 31;
        let m: usize = 10;
        let mut positions: Vec<usize> = vec![];
        get_minimizer_positions(seq.as_bytes(), &mut positions, k, m);

        assert_eq!(positions, vec![2,22,30,42])
    }

    fn get_true_kmer_occurrences(seqs: &[Vec<u8>], k: usize) -> std::collections::HashMap<Vec<u8>, Vec<(usize, usize)>>{
        let mut true_kmer_occurrences = std::collections::HashMap::<Vec<u8>, Vec<(usize, usize)>>::new();
        for (seq_id, seq) in seqs.iter().enumerate(){
            for i in 0 .. number_of_kmers(seq.len(), k){
                let key = &seq[i..i+k];
                let new_entry = (seq_id, i);
                if let Some(vec) = true_kmer_occurrences.get_mut(key){ // Existing k-mer
                    vec.push(new_entry);
                } else{ // New k-mer
                    true_kmer_occurrences.insert(key.to_owned(), vec![new_entry]);
                }
            }
        }
        true_kmer_occurrences
    }

    fn test_vs_hash_table(db: &SeqDB, k: usize, m: usize){
        
        // Build index
        let index = MinimizerIndex::new(db, k, m);
        
        // Read sequences
        let seqs = db.iter().map(|rec| rec.seq.to_owned()).collect::<Vec<Vec<u8>>>();

        // Build true k-mer occurrences map
        let true_kmer_occurrences = get_true_kmer_occurrences(&seqs, k);

        // Look up all k-mers in input sequences
        for seq in seqs.iter(){
            for i in 0 .. number_of_kmers(seq.len(), k){
                let kmer = &seq[i..i+k];
                let occs = index.lookup_kmer(kmer);
                eprintln!("{} {:?} {:?}", to_ascii(kmer), &occs, &true_kmer_occurrences[kmer]);
                assert_eq!(occs, true_kmer_occurrences[kmer]);
            }
        }

    }

    #[test]
    fn test_index_lookup_small(){

        let mut db = SeqDB::new();
        db.push_record(RefRecord{head: b"seq1", seq: b"ATAGCTAGTCGATGCTGATCGTAGGTTCGTAGCTGTATGCTGACCCTGATGTCTGTAGTCGTGACTGACT", qual: None});
        db.push_record(RefRecord{head: b"seq2 (substring of seq1)", seq: b"GTCGATGCTGATCGTAGGTTCGTAGCTGTATGCTGACCCTGATGTCTTGACT", qual: None});
        db.push_record(RefRecord{head: b"seq3 (seq2 with a single change in the middle)", seq: b"GTCGATGCTGATCGTAGGTTCGAAGCTGTATGCTGACCCTGATGTCTTGACT", qual: None});

        test_vs_hash_table(&db, 31, 10);
        test_vs_hash_table(&db, 10, 10);
        test_vs_hash_table(&db, 1, 1);

        // Look up a random k-mer (should not be found)
        let index = MinimizerIndex::new(&db, 31, 10);
        let random_kmer = "ATCTTATCTGGGGCTATTGCTAGGGCTTACA".as_bytes();
        assert_eq!(index.lookup_kmer(random_kmer).len(), 0);
    }

    #[test] 
    fn test_index_lookup_large_random(){

        let mut db = SeqDB::new();

        // Create a random number generator with the fixed seed
        let seed = [123; 32];
        let mut rng = rand_chacha::ChaCha20Rng::from_seed(seed);

        // Generate 100 random sequences of length 100
        for i in 0 .. 100 {
            let mut seq = vec![0; 100];
            #[allow(clippy::needless_range_loop)] // Clearer this way IMO
            for j in 0 .. 100 {
                seq[j] = match rng.next_u64() % 4 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => panic!("Impossible"),
                };
            }
            db.push_record(RefRecord{head: format!("seq{}", i).as_bytes(), seq: &seq, qual: None});
        }

        let seqs = db.iter().map(|rec| rec.seq.to_owned()).collect::<Vec<Vec<u8>>>();

        // We now have 10k base pairs.
        // There are 4^6 = 4096 possible 6-mers.
        // We use k = 6 so that many k-mers are found but not all.

        let k = 6;
        let m = 3;

        // For queries, do all possible 6-mers
        let mut queries = vec![vec![0; k]; 4usize.pow(k as u32)];

        #[allow(clippy::needless_range_loop)] // Clearer this way IMO
        for i in 0 .. 4usize.pow(k as u32){
            let mut j = i;
            for l in 0 .. k {
                queries[i][l] = match j % 4 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => panic!("Impossible"),
                };
                j /= 4;
            }

            eprintln!("{} {}", i, to_ascii(&queries[i]));
        }

        let mut true_occurrences = get_true_kmer_occurrences(&seqs, k);

        // Put in empty lists for k-mers that are not present
        for query in queries.iter(){
            if !true_occurrences.contains_key(query){
                true_occurrences.insert(query.to_owned(), vec![]);
            }
        }

        let index = MinimizerIndex::new(&db, k, m);

        // Verify
        for query in queries {
            let occs = index.lookup_kmer(&query);
            eprintln!("{} {:?} {:?}", to_ascii(&query), &occs, &true_occurrences[&query]);
            assert_eq!(occs, true_occurrences[&query]);
        }
        
    }

    // TODO: test handling Ns
}
