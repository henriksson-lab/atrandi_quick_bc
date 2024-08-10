use std::path::PathBuf;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};

use itertools::Itertools;


pub fn store_counttable(
    path_cnt:&PathBuf,
    counts:HashMap<String, HashMap<usize,i32>>,
    name_of_features:Vec<String>
) -> std::io::Result<()> {


    //Create a folder for the counts
    if !path_cnt.exists() {
        fs::create_dir(path_cnt)?;
    }

    //Figure out name of output files
    let path_count_file =  path_cnt.join("matrix.mtx");
    let path_features_file =  path_cnt.join("features.tsv");
    let path_bc_file =  path_cnt.join("barcodes.tsv");
    

    //Figure size of matrix
    //let num_feature = name_of_features.len();
    let num_cell = counts.len();
    let list_cell = counts.keys().map(|x| x).collect_vec();


    //%%MatrixMarket matrix coordinate integer general
    //89083 974 6075361
    
    ////// Write count table
    let output_h = File::create(path_count_file).expect("creation of R1 failed");
    let mut writer_h = BufWriter::new(output_h);
    //writer_h.write_all("%%MatrixMarket matrix coordinate real general\n".as_bytes()).expect("Unable to write data");
    writer_h.write_all("cell\tfeature\tcount\n".as_bytes()).expect("Unable to write data");

    for cellid in 0..num_cell {

        let cellmap = counts.get(list_cell[cellid]).unwrap();
        for (bc,cnt) in cellmap.iter() {
            let line = format!["{}\t{}\t{}\n", cellid+1, bc+1, cnt];
            writer_h.write_all(line.as_bytes()).expect("Unable to write data");
        }
    }

    ////// Write table with BC names
    let mut writer_cells = BufWriter::new(File::create(path_bc_file).expect("creation of cell table failed"));
    for cellid in 0..num_cell {
        let line = format!["{}\n", list_cell[cellid]];
        writer_cells.write_all(line.as_bytes()).expect("Unable to write data");
    }


    ////// Write table with feature names
    let mut writer_cells = BufWriter::new(File::create(path_features_file).expect("creation of feature table failed"));
    for feature in name_of_features {
        let line = format!["{}\n", feature];
        writer_cells.write_all(line.as_bytes()).expect("Unable to write data");
    }

    Ok(())
}
