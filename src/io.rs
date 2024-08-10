// This file is part of babbles which is released under the MIT license.
// See file LICENSE or go to https://github.com/HadrienG/babbles for full license details.
use std::process;
use std::path::PathBuf;
use itertools::Itertools;
use log::{debug, error, info};
use std::fs::{File, OpenOptions};


use niffler::get_reader;
use seq_io::fastq::Reader as FastqReader;
use seq_io::fasta::{Reader as FastaReader, Record as FastaRecord};

use bio::alignment::Alignment;
use bio::pattern_matching::myers::Myers;

pub struct Barcode {
    pub index: usize,
    pub name: String,
    pub pool: String,
    pub sequence: Vec<u8>,
    pub pattern: Myers<u64>
}
impl Barcode {
    pub fn seek(&mut self, record: &[u8], distance: u8) -> Vec<(usize, &String, usize, usize, i32)> {
        // use Myers' algorithm to find the barcodes in a read
        // Ref: Myers, G. (1999). A fast bit-vector algorithm for approximate string
        // matching based on dynamic programming. Journal of the ACM (JACM) 46, 395–415.
        let mut hits: Vec<(usize, &String, usize, usize, i32)> = Vec::new();
        let mut aln = Alignment::default();
        let mut matches = self.pattern.find_all_lazy(record, distance);
        let maybe_matches = matches.by_ref().min_set_by_key(|&(_, dist)| dist);
        if maybe_matches.len() > 0 {
            for (best_end, _) in maybe_matches {
                matches.alignment_at(best_end, &mut aln);
                hits.push((self.index, &self.name, aln.ystart, aln.yend, aln.score));
            }
        }
        hits
    }
}

// "only traits defined in the current crate can be implemented for types
// defined outside of the crate. Define and implement a trait or new type instead"
// If we'd ever want to have all fastx IO with niffler + seqio ...
// impl std::io::Seek for Box<(dyn std::io::Read + 'static)> {
//    println!("unimplemented!");
// }


pub fn open_fastq(file_handle: &PathBuf) -> FastqReader<Box<dyn std::io::Read>> {
    let opened_handle = match File::open(file_handle) {
        Ok(file) => file,
        Err(_) => {
            error!("Could not open file {}", &file_handle.display());
            process::exit(1)
        }
    };
    let (reader, _) = match get_reader(Box::new(opened_handle)) {
        Ok((reader, compression)) => {
            debug!("Opened file {} with compression {:?}", &file_handle.display(), &compression);
            (reader, compression)
        },
        Err(_) => {
            error!("Could read reverse file {}", &file_handle.display());
            process::exit(1)
        }
    };
    let fastq = FastqReader::new(reader);
    fastq
}


pub fn open_fasta(file_handle: &PathBuf) -> FastaReader<Box<dyn std::io::Read>> {
    let opened_handle = match File::open(file_handle) {
        Ok(file) => file,
        Err(_) => {
            error!("Could not open file {}", &file_handle.display());
            process::exit(1)
        }
    };
    let (reader, _) = match get_reader(Box::new(opened_handle)) {
        Ok((reader, compression)) => {
            debug!("Opened file {} with compression {:?}", &file_handle.display(), &compression);
            (reader, compression)
        },
        Err(_) => {
            error!("Could read reverse file {}", &file_handle.display());
            process::exit(1)
        }
    };
    let fasta = FastaReader::new(reader);
    fasta
}


pub fn open_fastq_no_box(file_handle: &PathBuf) -> FastqReader<File> {
    let opened = match File::open(file_handle) {
        Ok(file) => file,
        Err(_) => {
            error!("Could not open forward file {}", &file_handle.display());
            process::exit(1)
        }
    };
    let file = FastqReader::new(opened);
    file
}


pub fn read_barcodes(barcode_files: &Vec<PathBuf>) -> Vec<Barcode> {
    let mut barcodes: Vec<Barcode> = Vec::new();
    for barcode_file in barcode_files {
        let mut reader = open_fasta(barcode_file);
        let mut n_barcodes: usize = 0;
        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            let b = Barcode{
                index: n_barcodes,
                name: record.id().unwrap().to_string(),
                pool: barcode_file.file_stem().unwrap().to_str().unwrap().to_string(),
                sequence: record.seq().to_vec(),
                pattern: Myers::<u64>::new(record.seq().to_vec())
            };
            barcodes.push(b);
            n_barcodes += 1;
        }
    };
    // TODO check the edit distance between barcodes
    info!("Found {} barcodes in specified barcode files", barcodes.iter().count()); 
    barcodes
}


pub fn open_buffer_for_writing(path: &PathBuf, append: bool) -> File {
    let buffer = OpenOptions::new()
        .write(true)
        .append(append)
        .create(true)
        .open(&path);
    let buffer = match buffer {
        Ok(buffer) => buffer,
        Err(error) => {
            debug!("{:?}", error);
            error!("Could not create output file {}", &path.display());
            process::exit(1)
        }
    };
    buffer
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_open_buffer_for_writing() {
        let path = PathBuf::from("tests/data/test.txt");

        let maybe_buffer = open_buffer_for_writing(&path, false);
        assert_eq!(maybe_buffer.metadata().unwrap().is_file(), true);

        // cleanup
        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn test_read_barcodes() {
        // read_barcodes() calls open_fasta() which is therefore not tested separately
        let path = PathBuf::from("tests/data/barcodes.fasta");
        let paths = Vec::from([path]);
        let maybe_barcodes = read_barcodes(&paths);

        assert_eq!(maybe_barcodes.len(), 2);
        assert_eq!(maybe_barcodes[0].name, "A_0");
        assert_eq!(maybe_barcodes[1].sequence, b"TTGAGCCG".to_vec());
    }

    #[test]
    fn test_open_fastq_and_seek() {
        use seq_io::fastq::Record;
        let path = PathBuf::from("tests/data/reads.fastq");
        let mut maybe_reader = open_fastq(&path);
        let maybe_id = maybe_reader.next().unwrap().unwrap().to_owned_record();
        assert_eq!(maybe_id.id().unwrap(), "read_1");
    }

    #[test]
    fn test_seek() {
        let sequence = b"CTGCTTGAGCCGAGGGGATTATCTCGTAAGGCAAGCTCGT";

        let mut barcode = Barcode{
            index: 0,
            name: "test".to_string(),
            pool: "A".to_string(),
            sequence: b"TTGAGCCG".to_vec(),
            pattern: Myers::<u64>::new(b"TTGAGCCG".to_vec())
        };
        let hits = barcode.seek(sequence, 1);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].2, 4);  // start
        assert_eq!(hits[0].3, 12);  // end
    }
}