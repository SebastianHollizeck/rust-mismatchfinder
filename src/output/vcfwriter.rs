use std::{collections::BTreeMap, fs::File, io::Write, path::PathBuf};

use crate::bamreader::mismatch::Mismatch;

const VCF_HEADER: &str = "##fileformat=VCFv4.2
##source==MisMatchFinder
##FILTER=<ID=PASS,Description=\"Mismatch of high quality\">
##INFO=<ID=MULTI,Number=1,Type=Integer,Description=\"Number of reads supporting the mismatch\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

pub fn write_vcf(
    mismatches: &BTreeMap<Mismatch, usize>,
    file: PathBuf,
    include_small_vars: bool,
    somatic: bool,
) -> std::io::Result<()> {
    let mut vcf_fh = File::create(file).expect("Could not open file to write vcf");

    //write the header here
    writeln!(vcf_fh, "{VCF_HEADER}")?;

    for (mm, count) in mismatches {
        match mm.typ {
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

                    if somatic {
                        line += ";SOMATIC";
                    }
                    writeln!(vcf_fh, "{line}")?;
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

                    if somatic {
                        line += ";SOMATIC";
                    }
                    writeln!(vcf_fh, "{line}")?;
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

                if somatic {
                    line += ";SOMATIC";
                }
                writeln!(vcf_fh, "{line}")?;
            }
        }
    }
    return Ok(());
}
