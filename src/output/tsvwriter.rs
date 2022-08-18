use crate::bamreader::mismatch::Mismatch;
use std::fs::File;
use std::io::Write;
use std::{collections::BTreeMap, path::PathBuf};

pub fn write_mismatches(
    mismatches: &BTreeMap<Mismatch, usize>,
    file: PathBuf,
) -> std::io::Result<()> {
    //open the writer
    let mut file_writer = File::create(file).expect("Could not open file to write mismatches");

    for (mm, count) in mismatches {
        writeln!(file_writer, "{}\t{}", mm, count)?;
    }
    //
    file_writer.flush()?;

    return Ok(());
}
