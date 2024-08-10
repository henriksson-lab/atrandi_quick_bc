
use log::{error, debug}; //, info, trace, warn
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::path::PathBuf;
use std::process;
use std::error::Error;
use std::io::{BufWriter, Write};

use seq_io::fastq::Record as FastqRecord;
use seq_io::fastq::Reader as FastqReader;
use niffler::get_reader;
use csv::ReaderBuilder;
use clap::{Parser, Subcommand};
use gzp::{deflate::Gzip, par::compress::{ParCompress, ParCompressBuilder}, ZWriter};
use env_logger::{Builder, Env};


//////////////////////////////////////////
////////////////////////////////////////// Basic whitelist correction
//////////////////////////////////////////

pub struct BarcodeWhitelist {
    list: Vec<String>,    //List for alignment; not sure if worth having separate from set
    set: HashSet<String>, //Dictionary for fast lookup of exact matches
    bc_length: usize
}

impl BarcodeWhitelist {


    /// Compare to each BC, see which fits best --- each base that matches give 1p, other 0p
    fn closest_bc_basewise(&self, bc_to_match: &String) -> String {
        let mut best_bc = &self.list[0];
        let mut best_bc_score = num_similar_elements(bc_to_match.as_bytes(), best_bc.as_bytes());
        for j in 1..self.list.len() {
            let score = num_similar_elements(bc_to_match.as_bytes(), self.list[j].as_bytes());
            if score>best_bc_score {
                best_bc_score = score;
                best_bc = &self.list[j];
            }
        }
        //println!("best bc basewise {}",best_bc.to_string());
        return best_bc.to_string();
    }

    /// Correct barcode using whitelist
    fn correct_to_whitelist(&self, bc_to_match: &String) -> String {  /////////// do Result instead?
        if bc_to_match.len()==0 {
            //Empty barcode
            return "".to_string();
        } else if self.set.contains(bc_to_match) {
            //See if there is a trivial match
            //println!("trivial match");
            return bc_to_match.to_string();
        } else if self.bc_length==bc_to_match.len() {
            //Compare each base if same length
            return self.closest_bc_basewise(bc_to_match);
        } else {
            //Fail
            return "".to_string();
        }
    }

}



/// Count the number of similar elements in two lists of the same size
fn num_similar_elements(a:&[u8], b:&[u8]) -> i32 {
    let mut count = 0;
    for i in 0..a.len() {
        if a[i] == b[i] {
            count = count + 1;
        }
    }
    return count;
}




/// Structure for Atrandi combinatorial barcodes
pub struct AtrandiBarcodes {
    rounds: Vec<BarcodeWhitelist>
}

impl AtrandiBarcodes {

    /// Read dictionary of Atrandi barcodes from file
    fn read_atrandi_barcodes(filename:&str) -> Result<AtrandiBarcodes, Box<dyn Error>> {
        let mut rdr = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(filename)?;
        let mut bcs_for_well = vec![vec![] as Vec<String>; 4];
        let mut bc_length = 666;
        for result in rdr.records() {
            let record = result?;
            let pos=&record[0];
            //let well=&record[1];
            let bc=&record[2];
            bc_length = bc.len();
            let pos_int = pos.parse::<usize>().unwrap() - 1;
            bcs_for_well[pos_int].push(String::from(bc));        
        }

        let whitelists = bcs_for_well.iter().map(|w| BarcodeWhitelist {
            list: w.to_vec(),
            set: HashSet::from_iter(w.to_vec()),
            bc_length: bc_length
        }  ).collect();
        
        Ok(AtrandiBarcodes {rounds: whitelists})
    }


    /// Take read R1, return correct barcode
    fn get_correct_bc_from_read(&self, read_r1:&str) -> Result<(String,String,String,String),()> {

        //Extract each BC
        //let template_bc = br"********AGGA********ACTC********AAGG********T"; 
        //let barcode_tuple = extract_bc_by_alignment(template_bc, read_r1.as_bytes(), false);  

        let barcode_tuple = extract_bc_optimistic_atrandi(read_r1)?;

        //Note swap here of BCs to match logical order in chemistry. Barcode added last is the first one seen in R1
        let corrected_bc1 = self.rounds[0].correct_to_whitelist(&barcode_tuple.3);
        let corrected_bc2 = self.rounds[1].correct_to_whitelist(&barcode_tuple.2);
        let corrected_bc3 = self.rounds[2].correct_to_whitelist(&barcode_tuple.1);
        let corrected_bc4 = self.rounds[3].correct_to_whitelist(&barcode_tuple.0);
    
        return Ok((corrected_bc1,corrected_bc2,corrected_bc3,corrected_bc4));
    }

}


fn extract_bc_optimistic_atrandi(read_r1:&str) -> Result<(String,String,String,String),()> {

    if read_r1.len() > 36+8 {
        let barcode_1 = &read_r1[0..(0+8)];
        let barcode_2 = &read_r1[(12+0)..(12+8)];
        let barcode_3 = &read_r1[(24+0)..(24+8)];
        let barcode_4 = &read_r1[(36+0)..(36+8)];
        return Ok((barcode_1.to_string(),barcode_2.to_string(),barcode_3.to_string(),barcode_4.to_string()))
    } else {
        return Err(());
    }
}
















//////////////////////////////////////////
////////////////////////////////////////// /// Copied from babbles ; fastq reading
//////////////////////////////////////////


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



//////////////////////////////////////////
////////////////////////////////////////// Parse BC to fastq
//////////////////////////////////////////

/* 
fn write_fastq_str(parz: &mut ParCompress<Gzip>, readname:&str, seq:&str, qual:&str) {
    write_fastq(parz, readname.as_bytes(), seq.as_bytes(), qual.as_bytes());
}
*/

fn write_fastq(parz: &mut ParCompress<Gzip>, readname:&[u8], seq:&[u8], qual:&[u8]) {
    parz.write_all(b"@").unwrap();
    parz.write_all(readname).unwrap();
    parz.write_all(b"\n").unwrap();

    parz.write_all(seq).unwrap();
    parz.write_all(b"\n").unwrap();

    parz.write_all(b"+\n").unwrap();

    parz.write_all(qual).unwrap();
    parz.write_all(b"\n").unwrap();
}





fn parse_to_fastq(
    path_in_r1:&PathBuf,
    path_in_r2:&PathBuf,
    path_out_r1:&PathBuf,
    path_out_r2:&PathBuf,
    histogram_file:&PathBuf
) {

    println!("reading whitelist ");
    let atrandi_barcodes = AtrandiBarcodes::read_atrandi_barcodes("bc.csv").expect("Failed to read barcode file");

    /////////// Set up input
    let mut f_r1 = open_fastq(&path_in_r1);
    let mut f_r2 = open_fastq(&path_in_r2);

    /////////// Set up output
    let output_r1 = File::create(path_out_r1).expect("creation of R1 failed");
    let output_r2 = File::create(path_out_r2).expect("creation of R2 failed");

    let mut parz_r1: ParCompress<Gzip> = ParCompressBuilder::new().from_writer(output_r1);
    let mut parz_r2: ParCompress<Gzip> = ParCompressBuilder::new().from_writer(output_r2);


    let mut barcode_per_cell_count = HashMap::new();


    /////////// Handle all reads
    let mut read_count = 0;
    while let Some(record_r1) = f_r1.next() {

        read_count = read_count + 1;
        if read_count%100000 == 0 {
            println!("Processed reads: {}", read_count);
        }


        //let record_r1 = f_r1.next().expect("No r1");
        let record_r2 = f_r2.next().expect("No r2");
        

        let record_r1: seq_io::fastq::RefRecord = record_r1.expect("Error reading record");
        let record_r2: seq_io::fastq::RefRecord = record_r2.expect("Error reading record");
    
        let seq_r2=String::from_utf8_lossy(record_r2.seq());
        let bc = atrandi_barcodes.get_correct_bc_from_read(&seq_r2);

        match bc {
            Ok(bc) => {

                let concat_bc = format!("{}.{}.{}.{}",bc.0,bc.1,bc.2,bc.3);

                //Count barcodes
                match barcode_per_cell_count.get(&concat_bc) {
                    Some(cnt) => {
                        barcode_per_cell_count.insert(concat_bc.clone(), cnt+1);
                    },
                    None => {
                        barcode_per_cell_count.insert(concat_bc.clone(), 1);
                    }
                }

                //Typical FASTQ record
                //@M03699:228:000000000-LCH6K:1:1102:12164:1000 1:N:0:CAGGTT
                //NCAGTTACTTGCAGGAATCTCCACCTGCTCTCCATCGACTACGTCTTTCGACCTCGCCTTAGGTCCCGACTTACC
                //+
                //#8B<CFDGGGFGGFGGFGGGGGGGGGFGCGFFGGGGGDGFDEGGGGGGGGGGGCGCEGGGGGGGGGGGEFGGFGG


                //Read 1 is just passed through; just change the name
                let new_r1_name = format!("{}_{}",&concat_bc, record_r1.id().unwrap());
                write_fastq(&mut parz_r1, 
                    new_r1_name.as_bytes(),
                    record_r1.seq(),
                    record_r1.qual()
                );

                //For Read 2, we will chop off the BC part. Update name
                let new_r2_name = format!("{}_{}",&concat_bc, record_r2.id().unwrap());

                let from: usize = 36+8;
                let to = record_r1.seq().len();
                let from = if from<to {from} else {to}; //to be on the safe side
                let new_r2_seq = &record_r2.seq()[from..to];
                let new_r2_qual = &record_r2.qual()[from..to];

                write_fastq(&mut parz_r2, 
                    new_r2_name.as_bytes(),
                    new_r2_seq,
                    new_r2_qual
                );


            },
            Err(_) => {
                //println!("Cannot tell BC");
            }
        };
    }

    parz_r1.finish().unwrap();
    parz_r2.finish().unwrap();


    ////// Write barcode histogram
    let output_h = File::create(histogram_file).expect("creation of R1 failed");
    let mut writer_h = BufWriter::new(output_h);
    writer_h.write_all("barcode\tcount\n".as_bytes()).expect("Unable to write data");
    for (bc, cnt) in &barcode_per_cell_count {
        let toprint = format!("{}\t{}\n", bc, cnt);
        writer_h.write_all(toprint.as_bytes()).expect("Unable to write data");
    }



    println!("done");

}








/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// CLI parser ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////



#[derive(Parser)]
#[command(author, version, about, long_about = None)]  // reads from Cargo.toml
struct Cli {
    /// print debug info
    #[arg(short, long, default_value_t = false, global = true)]
    debug: bool,
    #[command(subcommand)]
    command: Option<Commands>,
}


#[derive(Subcommand)]
enum Commands {
    /// Identify BC, make fastq
    ToFastq {
        /// forward reads
        #[arg(long)]
        i1: PathBuf,
        /// reverse reads
        #[arg(long)]
        i2: PathBuf,

        /// forward reads
        #[arg(long)]
        o1: PathBuf,
        /// reverse reads
        #[arg(long)]
        o2: PathBuf,

        /// histogram output
        #[arg(long)]
        h: PathBuf

    }    
}


fn main() {

    let cli = Cli::parse();
    let level = if cli.debug { "debug" } else { "info" };
    Builder::from_env(Env::default().default_filter_or(level)).init();

    match &cli.command {
        Some(Commands::ToFastq { i1, i2, o1, o2, h}) => {
            parse_to_fastq(
                &i1, &i2, 
                &o1, &o2,
                &h
            );
        }
        None => {}
    }



}



