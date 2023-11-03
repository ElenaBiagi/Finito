use std::fs::File;
use std::io::BufWriter;
// Unit tests
#[cfg(test)]
use std::io::{BufReader};

use jseqio::reader::*;
use jseqio::record::*;
use jseqio::writer::*;
use jseqio::*;

// Utility function
fn get_test_filenames() -> Vec<&'static str>{
    vec!["tests/data/reads.fastq", "tests/data/reads.fna", "tests/data/reads.fastq.gz", "tests/data/reads.fasta.gz"]
}

#[test]
fn fastq() {
    let headers = vec![
        "SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1",
        "SRR403017.2 HWUSI-EAS108E_0007:3:1:10327:976/1",
        "SRR403017.3 HWUSI-EAS108E_0007:3:1:13569:972/1",
    ];
    let seqs = vec![
        "TTGGACCGGCGCAAGACGGACCAGNGCGAAAGCATTTGCCAAGAANNNN",
        "CAACTTTCTATCTGGCATTCCCTGNGGAGGAAATAGAATGCGCGCNNNN",
        "GATCGGAAGAGCACACGTCTGAACNCCAGTCACTTAGGCATCTCGNNNN",
    ];
    let quals = vec![
        "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQ",
        "RSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~####",
        "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@",
    ];

    let n_seqs = headers.len();
    let mut fastq_data: String = "".to_owned();
    for i in 0..n_seqs {
        fastq_data.push_str(format!("@{}\n", headers[i]).as_str());
        fastq_data.push_str(format!("{}\n", seqs[i]).as_str());
        fastq_data.push_str("+\n");
        fastq_data.push_str(format!("{}\n", quals[i]).as_str());
    }

    let input = BufReader::new(fastq_data.as_bytes());
    let mut reader = StaticFastXReader::new(input).unwrap();

    let mut owned_records: Vec<OwnedRecord> = vec![];
    let mut seqs_read = 0;
    while let Some(record) = reader.read_next().unwrap() {
        assert_eq!(record.head, headers[seqs_read].as_bytes());
        assert_eq!(record.seq, seqs[seqs_read].as_bytes());
        assert_eq!(record.qual.unwrap(), quals[seqs_read].as_bytes());
        owned_records.push(record.to_owned());
        seqs_read += 1;
    }
    assert_eq!(seqs_read, n_seqs);

    // Test writer
    let out_buf: Vec<u8> = vec![];
    let mut writer = FastXWriter::<Vec<u8>>::new(out_buf, FileType::FASTQ);

    for rec in owned_records.iter() {
        writer.write(&rec.as_ref_record()).unwrap();
    }

    let written_data = writer.into_inner().unwrap();

    // Read the records back from written data and compare to originals.

    let mut reader2 = StaticFastXReader::new(written_data.as_slice()).unwrap();
    let mut seqs_read2 = 0;
    while let Some(record) = reader2.read_next().unwrap() {
        dbg!(&record);
        assert_eq!(record.head, headers[seqs_read2].as_bytes());
        assert_eq!(record.seq, seqs[seqs_read2].as_bytes());
        assert_eq!(record.qual.unwrap(), quals[seqs_read2].as_bytes());
        seqs_read2 += 1;
    }
    assert_eq!(seqs_read2, n_seqs);
}

#[test]
fn fasta() {
    let headers: Vec<String> = vec![
        "SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1".to_owned(),
        "SRR403017.2 HWUSI-EAS108E_0007:3:1:10327:976/1".to_owned(),
        "SRR403017.3 HWUSI-EAS108E_0007:3:1:13569:972/1".to_owned(),
    ];
    let seqs: Vec<String> = vec![
        "TTGGACCGGCGCAAGACGGACCAGNGCGAAAGCATTTGCCAAGAANNNN".to_owned(),
        "CAACTTTCTATCTGGCATTCCCTGNGGAGGAAATAGAATGCGCGCNNNN".to_owned(),
        "GATCGGAAGAGCACACGTCTGAACNCCAGTCACTTAGGCATCTCGNNNN".to_owned(),
    ];

    fn split_seq_to_lines(seq: &String, line_length: usize) -> Vec<String> {
        let mut i: usize = 0;
        let mut lines = Vec::<String>::new();
        while line_length * i < seq.len() {
            lines.push(
                seq[line_length * i..std::cmp::min(line_length * (i + 1), seq.len())].to_owned(),
            );
            i += 1;
        }
        lines
    }

    let n_seqs = headers.len();
    let mut fasta_data: String = "".to_owned();
    for i in 0..n_seqs {
        fasta_data.push_str(format!(">{}\n", headers[i].as_str()).as_str());
        for line in split_seq_to_lines(&seqs[i], 11) {
            // Line length 11 should make it so that the last line has a different
            // length than the other lines.
            fasta_data.push_str(format!("{}\n", line.as_str()).as_str());
        }
    }

    dbg!(&fasta_data);

    let input = BufReader::new(fasta_data.as_bytes());
    let mut reader = StaticFastXReader::new(input).unwrap();

    let mut owned_records: Vec<OwnedRecord> = vec![];
    let mut seqs_read = 0;
    while let Some(record) = reader.read_next().unwrap() {
        dbg!(&record);
        assert_eq!(record.head, headers[seqs_read].as_bytes());
        assert_eq!(record.seq, seqs[seqs_read].as_bytes());
        assert_eq!(record.qual, None);
        owned_records.push(record.to_owned());
        seqs_read += 1;
    }
    assert_eq!(seqs_read, n_seqs);

    // Test writer
    let out_buf: Vec<u8> = vec![];
    let mut writer = FastXWriter::<Vec<u8>>::new(out_buf, FileType::FASTA);

    for rec in owned_records.iter() {
        writer.write(&rec.as_ref_record()).unwrap();
    }

    let written_data = writer.into_inner().unwrap();

    // This written data may not exactly equal the original data,
    // because the length of FASTA sequence lines is not fixed.
    // Read the records back from written data and compare to originals.

    let mut reader2 = StaticFastXReader::new(written_data.as_slice()).unwrap();
    let mut seqs_read2 = 0;
    while let Some(record) = reader2.read_next().unwrap() {
        dbg!(&record);
        assert_eq!(record.head, headers[seqs_read2].as_bytes());
        assert_eq!(record.seq, seqs[seqs_read2].as_bytes());
        assert_eq!(record.qual, None);
        seqs_read2 += 1;
    }
    assert_eq!(seqs_read2, n_seqs);
}

#[test]
fn test_figure_out_file_format() {
    assert!(matches!(figure_out_file_format("aa.fna"), (FileType::FASTA, false)));
    assert!(matches!(figure_out_file_format("aa.fq"), (FileType::FASTQ, false)));
    assert!(matches!(figure_out_file_format("bbb.fna.gz"), (FileType::FASTA, true)));
    assert!(matches!(figure_out_file_format("cc.fna.gz"), (FileType::FASTA, true)));
    assert!(matches!(figure_out_file_format(".fna.gz"), (FileType::FASTA, true)));
    assert!(matches!(figure_out_file_format(".fasta"), (FileType::FASTA, false)));
    assert!(matches!(figure_out_file_format(".fq"), (FileType::FASTQ, false)));
}

// Utility function
fn get_test_reader(data: &[u8]) -> Result<DynamicFastXReader, Box<dyn std::error::Error>>{
    std::fs::write("/tmp/test.fastq", data).unwrap();
    DynamicFastXReader::from_file(&String::from("/tmp/test.fastq"))
}

// Utility function
fn get_sequences(filename: &str) -> Vec<Vec<u8>>{
    let reader = DynamicFastXReader::from_file(&String::from(filename)).unwrap();
    let db = reader.into_db().unwrap();
    let db_records: Vec<RefRecord> = db.iter().collect();
    let seqs: Vec<Vec<u8>> = db_records.iter().map(|r| r.seq.to_owned()).collect();
    return seqs;
}

#[test]
fn gzip_auto_detection(){
    let true_seqs = get_sequences("tests/data/reads.fastq");

    for filename in get_test_filenames(){
        eprintln!("{}", filename);
        let input = BufReader::new(File::open(filename).unwrap());
        let reader = DynamicFastXReader::new(input).unwrap();
        let db = reader.into_db().unwrap();
        let seqs = db.iter().map(|r| r.seq.to_owned()).collect::<Vec<Vec<u8>>>();
        assert_eq!(true_seqs, seqs);
    }
}

#[test]
fn fastq_errors() {
    // Not enough quality values
    let reader = get_test_reader("@SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1\nACGT\n+\nIII\n".as_bytes()).unwrap();
    assert!(reader.into_db().is_err_and(|e| e.is::<jseqio::reader::ParseError>()));

    // Too many quality values
    let reader = get_test_reader("@SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1\nACGT\n+\nIIIIIII\n".as_bytes()).unwrap();
    assert!(reader.into_db().is_err_and(|e| e.is::<jseqio::reader::ParseError>()));

    // Not quality line
    let reader = get_test_reader("@SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1\nACGT\n+\n".as_bytes()).unwrap();
    assert!(reader.into_db().is_err_and(|e| e.is::<jseqio::reader::ParseError>()));

    // No plus line
    let reader = get_test_reader("@SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1\nACGT\nIIII\n".as_bytes()).unwrap();
    assert!(reader.into_db().is_err_and(|e| e.is::<jseqio::reader::ParseError>()));

    // No sequence line
    let reader = get_test_reader("@SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1\n+\nIIII\n".as_bytes()).unwrap();
    assert!(reader.into_db().is_err_and(|e| e.is::<jseqio::reader::ParseError>()));

    // No header. Should fail already when creating the reader because it
    // can not detect the file format
    let reader = get_test_reader("\nACGT\n+\nIIII\n".as_bytes());
    assert!(reader.is_err_and(|e| e.is::<jseqio::reader::ParseError>()));
}

// TODO: fasta errors

#[test]
fn test_into_db(){
    let reader1 = DynamicFastXReader::from_file(&String::from("tests/data/reads.fastq")).unwrap();
    let db = reader1.into_db().unwrap();
    let db_records: Vec<RefRecord> = db.iter().collect();
    
    let mut reader2 = DynamicFastXReader::from_file(&String::from("tests/data/reads.fastq")).unwrap();

    let mut i: usize = 0;
    while let Some(record) = reader2.read_next().unwrap() {
        assert_eq!(record.head, db_records[i].head);
        assert_eq!(record.seq, db_records[i].seq);
        assert_eq!(record.qual, db_records[i].qual);
        i += 1;
    }

    assert_eq!(i, db_records.len()); // Check that we read all records and no more than that

    // Test also the non-dynamic version
    let db2 = StaticFastXReader::new(
        std::io::BufReader::new(std::fs::File::open("tests/data/reads.fastq").unwrap())).unwrap().into_db().unwrap();
    let db2_records: Vec<RefRecord> = db2.iter().collect();
    assert_eq!(db_records, db2_records);

}

#[test]
fn test_into_db_with_rc(){
    let reader = DynamicFastXReader::from_file(&String::from("tests/data/reads.fastq")).unwrap();
    let (fw_db, rc_db) = reader.into_db_with_revcomp().unwrap();
    assert_eq!(fw_db.sequence_count(), rc_db.sequence_count());
    for i in 0..fw_db.sequence_count(){
        let fw_rec = fw_db.get(i);
        let rc_rec = rc_db.get(i);
        assert_eq!(fw_rec.head, rc_rec.head);
        assert_eq!(reverse_complement(&fw_rec.seq), rc_rec.seq);
        assert_eq!(fw_rec.qual.unwrap().iter().rev().map(|x| *x).collect::<Vec::<u8>>(), rc_rec.qual.unwrap().to_owned());
    }
}

#[test]
fn test_empty_file() {
    let data = b"";

    // Dynamic reader
    let bufreader = BufReader::new(data.as_slice());
    let reader = DynamicFastXReader::new(bufreader).unwrap();
    let db = reader.into_db().unwrap(); // Should not panic
    assert_eq!(db.iter().count(), 0);

    // non-dynamic FASTA reader
    let bufreader = BufReader::new(data.as_slice());
    let reader = StaticFastXReader::new(bufreader).unwrap();
    let db = reader.into_db().unwrap(); // Should not panic
    assert_eq!(db.iter().count(), 0);

    // non-dynamic FASTQ reader
    let bufreader = BufReader::new(data.as_slice());
    let reader = StaticFastXReader::new(bufreader).unwrap();
    let db = reader.into_db().unwrap(); // Should not panic
    assert_eq!(db.iter().count(), 0);
}

#[test]
fn test_flush_on_drop(){
    let bw = BufWriter::new(File::create("/tmp/test.fna").unwrap());
    let mut writer = FastXWriter::<BufWriter<File>>::new(bw, FileType::FASTA);
    let rec = OwnedRecord{head: b"header".to_vec(), seq: b"ACGT".to_vec(), qual: None};
    writer.write(&rec).unwrap();

    drop(writer);

    // Try to read the record back
    let mut reader = StaticFastXReader::new(BufReader::new(File::open("/tmp/test.fna").unwrap())).unwrap();
    let rec = reader.read_next().unwrap().unwrap();
    assert_eq!(rec.head, b"header");
    assert_eq!(rec.seq, b"ACGT");
    assert!(rec.qual.is_none());
 }
