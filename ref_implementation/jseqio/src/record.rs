use std::fmt;

pub trait Record{
    fn head(&self) -> &[u8];
    fn seq(&self) -> &[u8];
    fn qual(&self) -> Option<&[u8]>;

    fn name(&self) -> &[u8]{
        self.head().split(|c| *c == b' ').next().expect("Could not parse sequence name (first space-separated token)")
    }
}

// rec.head.split(|c| *c == b' ').next().unwrap()

#[derive(Debug, Hash, PartialEq, Eq)]
pub struct OwnedRecord{
    pub head: Vec<u8>,    
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>, // If FASTA, this is None
}

#[derive(Debug, Hash, PartialEq, Eq, Clone, Copy)]
pub struct RefRecord<'a>{
    pub head: &'a [u8],    
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>, // If FASTA, this is None
}

impl<'a> RefRecord<'a>{
    pub fn to_owned(&self) -> OwnedRecord{
        OwnedRecord { 
            head: self.head.to_vec(), 
            seq: self.seq.to_vec(), 
            qual: self.qual.map(|q| q.to_vec()),
        }
    }
}

impl OwnedRecord{
    pub fn as_ref_record(&self) -> RefRecord{
        RefRecord { 
            head: self.head.as_slice(), 
            seq: self.seq.as_slice(), 
            qual: match &self.qual {
                Some(q) => Some(q.as_slice()), 
                None => None
            }
        }
    }

    pub fn reverse_complement(&mut self){
        crate::reverse_complement_in_place(&mut self.seq);
        if let Some(qual) = &mut self.qual{
            qual.reverse();
        }
    }
}


impl<'a> Record for RefRecord<'a>{
    fn head(&self) -> &[u8]{self.head}
    fn seq(&self) -> &[u8]{self.seq}
    fn qual(&self) -> Option<&[u8]>{self.qual}
}

impl Record for OwnedRecord{
    fn head(&self) -> &[u8]{self.head.as_slice()}
    fn seq(&self) -> &[u8]{self.seq.as_slice()}
    fn qual(&self) -> Option<&[u8]>{
        match &self.qual{
            Some(q) => return Some(q.as_slice()),
            None => None,
        }
    }
}

impl<'a> fmt::Display for RefRecord<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
               "SeqRecord{{ \n  Head: {}\n  Seq:  {}\n  Qual: {}\n}}", 
               std::str::from_utf8(self.head).unwrap(),
               std::str::from_utf8(self.seq).unwrap(),
               match self.qual{
                   Some(q) => std::str::from_utf8(q).unwrap(),
                   None => "", // No quality values
               }
               
        )
    }
}