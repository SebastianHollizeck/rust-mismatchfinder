use std::{
    ops::{Range, RangeInclusive},
    path::PathBuf,
};

use clap::{Parser, ValueHint};
use mismatchfinder::bamreader;
use rust_htslib::bam;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Options {
    /// Reference to be used for analysis.
    #[clap(short='T', long="reference", value_hint = ValueHint::FilePath)]
    reference_file: PathBuf,

    /// Bed file for genomic regions to ignore in the analysis. White list bed file regions overwrite black list regions
    #[clap(long="blacklist_bed", value_hint = ValueHint::FilePath)]
    blacklist_file: Option<PathBuf>,

    /// Bed file for genomic regions to include in the analysis. White list bed file regions overwrite black list regions [default: all regions]
    #[clap(long="whitelist_bed", value_hint = ValueHint::FilePath)]
    whitelist_file: Option<PathBuf>,

    /// Mimimum mapping quality for a read to be considered
    #[clap(
        short = 'q',
        long = "mininum_mapping_quality",
        value_parser,
        default_value_t = 20
    )]
    min_mq: u8,

    /// Mimimum base quality for a read to be considered
    #[clap(
        short = 'Q',
        long = "minimum_base_quality",
        value_parser,
        default_value_t = 25
    )]
    min_bq: u16,

    /// Mimimum average base quality in the whole read for a read to be considered
    #[clap(
        long = "minimum_average_basequality",
        value_parser,
        default_value_t = 20
    )]
    min_avg_bq: u16,

    /// Maximum mismatches we allow in the read
    #[clap(
        long = "maximum_mismatches_per_read",
        value_parser,
        default_value_t = 7
    )]
    max_mismatches_per_read: u8,

    /// Mimimum mismatches we require in the read
    #[clap(
        long = "minimum_mismatches_per_read",
        value_parser,
        default_value_t = 1
    )]
    min_mismatches_per_read: u8,

    /// Length of fragments to be considered in the analysis
    #[clap(long= "fragment_length_intervals", value_parser, action = clap::ArgAction::Append, multiple_values=true, default_values(&["74-155", "250-325"]))]
    fragment_length_intervals: Vec<String>,

    /// only use the overlap of the two reads
    #[clap(long, value_parser, default_value_t = false , action = clap::ArgAction::SetTrue)]
    only_overlaps: bool,

    /// Bam to analyse
    #[clap(value_hint = ValueHint::FilePath)]
    bam: PathBuf,
}

fn main() {
    let cli = Options::parse();

    // let gl_filter = bamreader::filter::germline::ZarrStorage::load_zarr_path(
    //     "/Volumes/bioinf/data/reference/dawson_labs/gnomad/3/zarr/chr1",
    // );

    //create white list object
    let bed = match cli.whitelist_file {
        Some(file) => bamreader::filter::region::BedObject::lapper_from_bed(file),
        None => bamreader::filter::region::BedObject::lapper_from_bed(PathBuf::from("/Volumes/bioinf/data/reference/dawson_labs/bed_files/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.100mer.highMappability.bed")),
    };

    let mut bam = bam::Reader::from_path(cli.bam).unwrap();

    let mut fragment_length_intervals: Vec<RangeInclusive<i64>> = Vec::new();

    for ivl in cli.fragment_length_intervals {
        let ivl: Vec<&str> = ivl.split("-").collect();
        if ivl.len() != 2 {
            panic!("We require the fragment length intervals to be in the shape of <x>-<y>")
        } else {
            //try to convert the intervals into numerics
            let start = match ivl[0].parse::<i64>() {
                Ok(v) => v,
                Err(_) => panic!("Could not convert start value to numeric"),
            };
            let end = match ivl[1].parse::<i64>() {
                Ok(v) => v,
                Err(_) => panic!("Could not convert end value to numeric"),
            };
            fragment_length_intervals.push(start..=end);
        }
    }

    mismatchfinder::bamreader::find_mismatches(&mut bam, &bed, &fragment_length_intervals);
}
