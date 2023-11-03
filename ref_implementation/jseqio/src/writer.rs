
use std::io;
use std::io::BufWriter;
use std::io::Write;
use ex::fs::File;
use flate2::Compression;
use flate2::write::GzEncoder;

use crate::FileType;
use crate::record::{RefRecord,OwnedRecord,Record};
use crate::figure_out_file_format;

pub trait SeqRecordWriter{
    // We can't use the generic Record trait here because then this can not be made into a trait object
    // so we have separate functions for owned and ref records.
    fn write_owned_record(&mut self, rec: &OwnedRecord) -> Result<(), Box<dyn std::error::Error>>;
    fn write_ref_record(&mut self, rec: &RefRecord) -> Result<(), Box<dyn std::error::Error>>;

    fn flush(&mut self) -> Result<(), Box<dyn std::error::Error>>;
}

// A dynamic writer, i.e. one that takes no generics and uses dyn instead
pub struct DynamicFastXWriter {
    stream: Box<dyn SeqRecordWriter + Send>,
}

// Non-dynamic writer, i.e. a writer that takes the internal stream as a generic parameter
pub struct FastXWriter<W: Write>{
    pub filetype: FileType,
    pub output: BufWriter<W>,
}

impl DynamicFastXWriter{

    pub fn write<Rec: Record>(&mut self, rec: &Rec) -> Result<(), Box<dyn std::error::Error>>{
        let r = RefRecord{head: rec.head(), seq: rec.seq(), qual: rec.qual()};
        self.stream.write_ref_record(&r)?;
        Ok(())
    }

    // No need to give a buffered writer. Buffering is handled internally.
    // If a buffered writer is given, then it will be buffered twice.
    pub fn new<W: Write + 'static + Send>(stream: W, filetype: FileType) -> Self{
        let writer = FastXWriter::<W>::new(stream, filetype);
        DynamicFastXWriter {stream: Box::new(writer)}
    }

    // Write to a file
    pub fn new_to_file<P: AsRef<std::path::Path>>(filename: &P) -> Result<Self, Box<dyn std::error::Error>> {
        let output = File::create(filename)?;
        match figure_out_file_format(filename){
            (FileType::FASTQ, true) =>{
                let gzencoder = GzEncoder::<File>::new(output, Compression::fast());
                Ok(Self::new(gzencoder, FileType::FASTQ))
            },
            (FileType::FASTQ, false) => {
                Ok(Self::new(output, FileType::FASTQ))
            },
            (FileType::FASTA, true) => {
                let gzencoder = GzEncoder::<File>::new(output, Compression::fast());
                Ok(Self::new(gzencoder, FileType::FASTA))
            },
            (FileType::FASTA, false) => {
                Ok(Self::new(output, FileType::FASTA))
            },
        }
    }

    pub fn new_to_stdout(filetype: FileType, compression_type: crate::CompressionType) -> Self {
        match compression_type{
            crate::CompressionType::Gzip => Self::new(GzEncoder::new(io::stdout(), Compression::fast()), filetype),
            crate::CompressionType::None => Self::new(io::stdout(), filetype),
        }
    }
}

impl<W: Write> FastXWriter<W>{
    pub fn write<Rec : Record>(&mut self, rec: &Rec) -> Result<(), std::io::Error> {
        match &self.filetype{
            FileType::FASTA => {
                self.output.write_all(b">")?;
                self.output.write_all(rec.head())?;
                self.output.write_all(b"\n")?;
                self.output.write_all(rec.seq())?;
                self.output.write_all(b"\n")?;
            }
            FileType::FASTQ => {
                self.output.write_all(b"@")?;
                self.output.write_all(rec.head())?;
                self.output.write_all(b"\n")?;
                self.output.write_all(rec.seq())?;
                self.output.write_all(b"\n+\n")?;
                self.output.write_all(rec.qual().expect("Quality values missing"))?;
                self.output.write_all(b"\n")?;
            }
        }
        Ok(())
    }

    pub fn new(output: W, filetype: FileType) -> Self{
        Self{
            filetype,
            output: BufWriter::<W>::new(output)
        }
    }

    pub fn flush(&mut self) -> Result<(), std::io::Error>{
        self.output.flush()?;
        Ok(())
    }

    // Returns the internal write-struct.
    // Also flushes the internal buffer.
    pub fn into_inner(mut self) -> Result<W, std::io::Error>{
        self.flush()?;
        Ok(self.output.into_parts().0)
    }
}

impl<W: Write> SeqRecordWriter for FastXWriter<W>{
    fn write_ref_record(&mut self, rec: &RefRecord) -> Result<(), Box<dyn std::error::Error>>{
        self.write(rec)?;
        Ok(())
    }

    fn write_owned_record(&mut self, rec: &OwnedRecord) -> Result<(), Box<dyn std::error::Error>>{
        self.write(rec)?;
        Ok(())
    }

    fn flush(&mut self) -> Result<(), Box<dyn std::error::Error>>{
        self.output.flush()?;
        Ok(())
    }
}

impl SeqRecordWriter for DynamicFastXWriter{
    fn write_ref_record(&mut self, rec: &RefRecord) -> Result<(), Box<dyn std::error::Error>> {
        self.stream.write_ref_record(rec)?;
        Ok(())
    }

    fn write_owned_record(&mut self, rec: &OwnedRecord) -> Result<(), Box<dyn std::error::Error>> {
        self.stream.write_owned_record(rec)?;
        Ok(())
    }

    fn flush(&mut self) -> Result<(), Box<dyn std::error::Error>>{
        self.stream.flush()?;
        Ok(())
    }
}