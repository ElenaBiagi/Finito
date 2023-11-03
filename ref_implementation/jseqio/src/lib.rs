
//! This crate provides parsers for sequences in FASTA and FASTQ format.
//! 
//! # Libary design
//! In bioinformatics, sequences are usually stored in files in FASTA or FASTQ format,
//! which are often compressed with gzip. This makes a total of four formats: FASTA with
//! and without gzip, and FASTQ with and without gzip. The purpose of the crate is to provide
//! a parser that can automatically detect the format of the file and parse it without
//! the user having to know beforehand which format is being used. **The file format is detected
//! from the first two bytes of the file, and does not depend on the file extension**.
//! 
//! We use dynamic dispatch to hide the details of the file format from the user. This introduces
//! an overhead of one dynamic dispatch per sequence, which is likely negligible unless the sequences
//! are extremely short. This also allows us to support reading from any byte stream, such as the standard
//! input, without having to attach generic parameters onto the parser. The interface is implemented for 
//! the struct [reader::DynamicFastXReader]. There is also [reader::StaticFastXReader] that takes the input
//! stream as a generic parameter.
//! 
//! A sequence is represented with a [record::RefRecord] struct that points to slices in the internal buffers of the reader. 
//! This is to avoid allocating new memory for each sequence. There also exists [record::OwnedRecord] which owns the memory.
//! 
//! Since the readers stream over the data, we can not implement the Rust Iterator trait. The lifetime constraints
//! on Rust Iterators require that all elements are valid until the end of the iteration. To support iterators,
//! we provide the [seq_db::SeqDB] struct that concatenates all sequences, headers and quality values in memory and provides 
//! an iterator over them.
//! 
//! # Examples
//! 
//! ## Streaming all sequences in a file and printing them to the standard output.
//! 
//! ```
//! use jseqio::reader::*;
//! fn main() -> Result<(), Box<dyn std::error::Error>>{
//!     // Reading from a FASTQ file. Also works for FASTA,
//!     // and seamlessly with/without gzip compression.
//!     let mut reader = DynamicFastXReader::from_file(&"tests/data/reads.fastq.gz")?;
//!     while let Some(rec) = reader.read_next().unwrap() {
//!         // Headers do not include the leading '>' in FASTA or '@' in FASTQ.
//!         eprintln!("Header: {}", std::str::from_utf8(rec.head)?);
//!         eprintln!("Sequence: {}", std::str::from_utf8(rec.seq)?);
//!         if let Some(qual) = rec.qual{
//!             // Quality values are present only in fastq files.
//!             eprintln!("Quality values: {}", std::str::from_utf8(qual)?);
//!         }
//!     }
//!     Ok(())
//! }
//! ```
//! 
//! ## Loading sequences into memory and computing the total length using an iterator.
//! 
//! ```
//! use jseqio::reader::DynamicFastXReader;
//! fn main() -> Result<(), Box<dyn std::error::Error>>{
//!     let reader = DynamicFastXReader::from_file(&"tests/data/reads.fna")?;
//!     let db = reader.into_db()?;
//!     let total_length = db.iter().fold(0_usize, |sum, rec| sum + rec.seq.len());
//!     eprintln!("Total sequence length: {}", total_length);
//!     Ok(())
//! }
//! ```
//! 

use std::path::Path;

pub mod reader;
pub mod writer;
pub mod record;
pub mod seq_db;

#[derive(Copy, Clone, Debug)]
pub enum FileType{
    FASTA,
    FASTQ,
}

#[derive(Copy, Clone, Debug)]
pub enum CompressionType{
    Gzip,
    None,
} // This could be just a bool, but in the future we may want to support other compression types as well

// Returns (file type, is_gzipped)
pub fn figure_out_file_format<P: AsRef<Path>>(filepath: P) -> (FileType, bool){
    let filename = filepath.as_ref().as_os_str().to_str().unwrap();

    let is_gzipped = filename.ends_with(".gz");
    let filename = if is_gzipped{
        &filename[0 .. filename.len()-3] // Drop the .gz suffix
    } else {filename};
    let fasta_extensions = vec![".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa"];
    let fastq_extensions = vec![".fastq", ".fq"];
    if fasta_extensions.iter().any(|&suffix| filename.ends_with(suffix)){
        (FileType::FASTA, is_gzipped)
    } else if fastq_extensions.iter().any(|&suffix| filename.ends_with(suffix)){
        (FileType::FASTQ, is_gzipped)
    } else{
        panic!("Unkown file extension: {}", filename);
    }
}

pub fn complement(c: u8) -> u8{
    match c{
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'a' => b't',
        b't' => b'a',
        b'g' => b'c',
        b'c' => b'g',
        other => other,
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8>{
    seq.iter().rev().map(|&c| complement(c)).collect()
}

pub fn reverse_complement_in_place(seq: &mut [u8]){
    for i in 0..seq.len(){
        seq[i] = complement(seq[i]);
    }
    seq.reverse();
}