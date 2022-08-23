use crate::bamreader::mismatch::Mismatch;
use bgzip::BGZFWriter;
use std::io::Write;
use std::{collections::BTreeMap, fs::File, path::PathBuf};

const VCF_HEADER: &[u8] = "##fileformat=VCFv4.2
##source==MisMatchFinder
##FILTER=<ID=PASS,Description=\"Mismatch of high quality\">
##INFO=<ID=MULTI,Number=1,Type=Integer,Description=\"Number of reads supporting the mismatch\">
##INFO=<ID=SOMATIC,Number=1,Type=String,Description=\"This variant is somatic after germline check\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    .as_bytes();

pub fn write_vcf(
    mismatches: &BTreeMap<Mismatch, usize>,
    file: PathBuf,
    include_small_vars: bool,
    somatic: bool,
) -> std::io::Result<()> {
    let mut vcf_fh = File::create(file).expect("Could not open file to write vcf");

    let mut writer = BGZFWriter::new(&mut vcf_fh, flate2::Compression::default());

    //write the header here
    writer.write_all(VCF_HEADER)?;

    for (mm, count) in mismatches {
        let line = match mm.typ {
            crate::bamreader::mismatch::MismatchType::SBS => {
                //we only write this if requested
                if include_small_vars {
                    let mut line = format!(
                        "{}\t{}\t.\t{}\t{}\t.\tPASS\tMULTI={}",
                        mm.chromosome,
                        mm.position,
                        std::str::from_utf8(&mm.reference[1..2]).unwrap(),
                        std::str::from_utf8(&mm.alternative[1..2]).unwrap(),
                        count,
                    );

                    //finish up the line
                    if somatic {
                        line += ";SOMATIC\n";
                    } else {
                        line += "\n";
                    }

                    Some(line)
                } else {
                    None
                }
            }
            crate::bamreader::mismatch::MismatchType::DBS => {
                //we only write this if requested
                if include_small_vars {
                    let mut line = format!(
                        "{}\t{}\t.\t{}\t{}\t.\tPASS\tMULTI={}",
                        mm.chromosome,
                        mm.position,
                        std::str::from_utf8(&mm.reference[..]).unwrap(),
                        std::str::from_utf8(&mm.alternative[..]).unwrap(),
                        count,
                    );

                    //finish up the line
                    if somatic {
                        line += ";SOMATIC\n";
                    } else {
                        line += "\n";
                    }
                    Some(line)
                } else {
                    None
                }
            }
            crate::bamreader::mismatch::MismatchType::INS
            | crate::bamreader::mismatch::MismatchType::DEL => {
                //we always write indels
                let mut line = format!(
                    "{}\t{}\t.\t{}\t{}\t.\tPASS\tMULTI={}",
                    mm.chromosome,
                    mm.position,
                    std::str::from_utf8(&mm.reference[..]).unwrap(),
                    std::str::from_utf8(&mm.alternative[..]).unwrap(),
                    count,
                );

                //finish up the line
                if somatic {
                    line += ";SOMATIC\n";
                } else {
                    line += "\n";
                }
                Some(line)
            }
        };
        if let Some(l) = line {
            writer.write_all(l.as_bytes())?;
        }
    }
    writer.close()?;
    return Ok(());
}
