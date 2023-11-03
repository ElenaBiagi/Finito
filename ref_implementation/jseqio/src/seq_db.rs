use crate::record::{RefRecord};

pub struct SeqDB {
    headbuf: Vec<u8>,
    seqbuf: Vec<u8>,
    qualbuf: Vec<u8>,
    head_starts: Vec<usize>, // Contains end sentinel at the end
    seq_starts: Vec<usize>, // Contains end sentinel at the end
    qual_starts: Vec<usize>, // Contains end sentinel at the end.

    // A mix of records with and without quality values is allowed. Then
    // the quality value slices will have length 0 for records without quality values.
}

impl SeqDB{

    pub fn iter(&self) -> SeqDBIterator {
        SeqDBIterator{seq_db: self, pos: 0}
    }

    pub fn sequence_count(&self) -> usize{
        self.head_starts.len() - 1
        // ^ The -1 is because we have an end sentinel at the end of the head_starts vector
    }

    pub fn get(&self, seq_index: usize) -> RefRecord{
        if seq_index >= self.head_starts.len(){
            panic!("SeqDB: Sequence index {} not found in database containing {} sequences", seq_index, self.sequence_count());
        }

        let head = &self.headbuf[self.head_starts[seq_index]..self.head_starts[seq_index+1]];
        let seq = &self.seqbuf[self.seq_starts[seq_index]..self.seq_starts[seq_index+1]];
        let qual = {
            let start = self.qual_starts[seq_index];   
            let end = self.qual_starts[seq_index+1];
            if start == end {
                None
            }
            else {
                Some(&self.qualbuf[start..end])
            }
        };
        RefRecord{head, seq, qual}
    }

    pub fn new() -> SeqDB{
        let headbuf: Vec<u8> = Vec::new();
        let seqbuf: Vec<u8> = Vec::new();
        let qualbuf: Vec<u8> = Vec::new();

        let head_starts: Vec<usize> = vec![0];
        let seq_starts: Vec<usize> = vec![0];
        let qual_starts: Vec<usize> = vec![0];
        
        SeqDB{headbuf, seqbuf, qualbuf, head_starts, seq_starts, qual_starts}
    }

    pub fn push_record<R: crate::record::Record>(&mut self, rec: R){
        self.headbuf.extend_from_slice(rec.head());
        self.seqbuf.extend_from_slice(rec.seq());
        self.head_starts.push(self.headbuf.len());
        self.seq_starts.push(self.seqbuf.len());

        if let Some(qual) = rec.qual(){
            // Record has quality values
            self.qualbuf.extend_from_slice(qual);
        }
        self.qual_starts.push(self.qualbuf.len());
    }

    // Push a sequence with no quality values or header
    pub fn push_seq(&mut self, seq: &[u8]){
        self.seqbuf.extend_from_slice(seq);
        self.seq_starts.push(self.seqbuf.len());

        self.head_starts.push(self.headbuf.len()); // Empty header
        self.qual_starts.push(self.qualbuf.len()); // Empty quality values
    }

    pub fn shrink_to_fit(&mut self){
        self.headbuf.shrink_to_fit();
        self.seqbuf.shrink_to_fit();
        self.qualbuf.shrink_to_fit();
    }
}

pub struct SeqDBIterator<'a>{
    seq_db: &'a SeqDB,
    pos: usize,
}

impl<'a> Iterator for SeqDBIterator<'a> {
    type Item = RefRecord<'a>;

    fn next(&mut self) -> Option<RefRecord<'a>> {
        match self.pos{
            i if i < self.seq_db.head_starts.len() - 1 => { // Iteration is not finished yet
                self.pos += 1; // Advance pointer to next element for the next round
                Some(self.seq_db.get(i)) // Should never be out of bounds so we unwrap the error.
            }
            _ => None, // End of iteration
        }
    }
}

impl ExactSizeIterator for SeqDBIterator<'_> {
    fn len(&self) -> usize {
        self.seq_db.sequence_count()
    }
}