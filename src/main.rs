use std::{
    fmt::format,
    fs,
    ops::{Range, RangeInclusive},
    path::PathBuf,
};

use clap::{Parser, ValueHint};
use mismatchfinder::{
    bamreader::{self, mismatch},
    output,
};
use rust_htslib::bam;

/// Reference to be used for analysis. (should only be relevant for crams)
// #[clap(short='T', long="reference", value_hint = ValueHint::FilePath)]
// reference_file: PathBuf,

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Options {
    /// Output folder to write files to
    #[clap(short='o', long="output", value_hint = ValueHint::FilePath)]
    output_folder: PathBuf,

    /// Bed file for genomic regions to ignore in the analysis. White list bed file regions overwrite black list regions
    #[clap(long="blacklist_bed", value_hint = ValueHint::FilePath)]
    blacklist_file: Option<PathBuf>,

    /// Bed file for genomic regions to include in the analysis. White list bed file regions overwrite black list regions [default: all regions]
    #[clap(long="whitelist_bed", value_hint = ValueHint::FilePath)]
    whitelist_file: Option<PathBuf>,

    /// File to read germline information from default is an echtvar file
    #[clap(long="germline_file", value_hint = ValueHint::FilePath)]
    germline_file: Option<PathBuf>,

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
        long = "maximum_edit_distance_per_read",
        value_parser,
        default_value_t = 7
    )]
    max_mismatches_per_read: u8,

    /// Mimimum mismatches we require in the read
    #[clap(
        long = "minimum_edit_distance_per_read",
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

    /// Bams to analyse
    #[clap(value_hint = ValueHint::FilePath, action = clap::ArgAction::Append, multiple_values=true)]
    bams: Vec<PathBuf>,
}

fn main() {
    let cli = Options::parse();

    env_logger::init();

    //check if we can write to the output folder, or if it doesnt exist yet, if we can create it
    match fs::create_dir_all(&cli.output_folder) {
        Ok(_) => { // we are happy and dont worry
        }
        Err(e) => panic!("Couldnt open output folder: {}", e),
    };

    // let gl_filter = bamreader::filter::germline::ZarrStorage::load_zarr_path(
    //     "/Volumes/bioinf/data/reference/dawson_labs/gnomad/3/zarr/chr1",
    // );

    // make this an option with lapper being default
    let lapper = true;
    //create white list object
    let bed;
    if lapper {
        // obviously this heavily depends on the network speed etc, but this was run multiple time to ensure a representative sampling
        // real	11m30.015s, 11m50.762s, 12m14.666s, 11m22.795s
        // user	8m1.104s, 8m3.983s, 7m57.696s, 7m35.625s
        // sys	0m38.480s, 0m42.985s, 0m38.536s, 0m26.196s
        bed = match cli.whitelist_file {
        Some(file) => bamreader::filter::region::BedObject::lapper_from_bed(file),
        None => bamreader::filter::region::BedObject::lapper_from_bed(PathBuf::from("/Volumes/bioinf/data/reference/dawson_labs/bed_files/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.100mer.highMappability.bed")),
    };
    } else {
        // real	12m18.631s, 11m57.539s, 12m48.091s, 13m0.463s, 12m13.922s
        // user	8m17.179s, 7m54.323s, 8m21.796s, 8m15.832s, 8m2.991s
        // sys	0m52.488s, 0m32.826s, 0m42.453s, 0m47.496s, 0m34.951s
        bed = match cli.whitelist_file {
            Some(file) => bamreader::filter::region::BedObject::interval_tree_from_bed(file),
            None => bamreader::filter::region::BedObject::interval_tree_from_bed(PathBuf::from("/Volumes/bioinf/data/reference/dawson_labs/bed_files/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.100mer.highMappability.bed")),
        };
    }

    let mut gnomad = match cli.germline_file {
        Some(file) => bamreader::filter::germline::GermlineResource::load_echtvars_file(file.to_str().unwrap()),
        None => bamreader::filter::germline::GermlineResource::load_echtvars_file("/Volumes/bioinf/data/reference/dawson_labs/gnomad/3/echtvar/gnomad.v3.1.2.echtvar.v2.zip"),
    };

    // let backend = match gnomad.get_backend_mut() {
    //     bamreader::filter::germline::GermlineBackend::Echt(v) => {
    //         v.set_position(1, String::from("chr8"), 75576664)
    //             .expect("Could not change the pointer in germline reference");
    //         v
    //     }
    //     bamreader::filter::germline::GermlineBackend::Zarr(_) => todo!(),
    // };

    // let enc = echtvar_lib::var32::encode(
    //     // we seem to need genomic pos -1 for the echtvar cache to find it
    //     75576663, b"C", b"A", &mut 0,
    // );

    // match backend.var32s.binary_search(&enc) {
    //     Ok(_) => println!("Found"),
    //     Err(_) => println!("not found"),
    // }

    // let test_mismatch_present = bamreader::mismatch::Mismatch {
    //     chromosome: String::from("chr1"),
    //     position: 55051215,
    //     reference: "AGT".as_bytes().to_vec(),
    //     alternative: "AAT".as_bytes().to_vec(),
    //     quality: 37,
    //     typ: bamreader::mismatch::MismatchType::SBS,
    //     rid: 1,
    // };

    // let test_mismatch_absent = bamreader::mismatch::Mismatch {
    //     chromosome: String::from("chr1"),
    //     position: 55051215,
    //     reference: "AGT".as_bytes().to_vec(),
    //     alternative: "ACT".as_bytes().to_vec(),
    //     quality: 37,
    //     typ: bamreader::mismatch::MismatchType::SBS,
    //     rid: 1,
    // };

    // let test_dbs_present1 = bamreader::mismatch::Mismatch {
    //     chromosome: String::from("chr8"),
    //     position: 75576664,
    //     reference: "CCA".as_bytes().to_vec(),
    //     alternative: "ATA".as_bytes().to_vec(),
    //     quality: 74,
    //     typ: bamreader::mismatch::MismatchType::DBS,
    //     rid: 1,
    // };

    // let test_dbs_present2 = bamreader::mismatch::Mismatch {
    //     chromosome: String::from("chr8"),
    //     position: 75576664,
    //     reference: "CCA".as_bytes().to_vec(),
    //     alternative: "GTA".as_bytes().to_vec(),
    //     quality: 74,
    //     typ: bamreader::mismatch::MismatchType::DBS,
    //     rid: 1,
    // };

    // let test_dbs_absent = bamreader::mismatch::Mismatch {
    //     chromosome: String::from("chr8"),
    //     position: 75576664,
    //     reference: "CCA".as_bytes().to_vec(),
    //     alternative: "GGA".as_bytes().to_vec(),
    //     quality: 74,
    //     typ: bamreader::mismatch::MismatchType::DBS,
    //     rid: 1,
    // };

    // let test_indel_short_present = bamreader::mismatch::Mismatch {
    //     chromosome: String::from("chr8"),
    //     position: 1055795,
    //     reference: "CAG".as_bytes().to_vec(),
    //     alternative: "C".as_bytes().to_vec(),
    //     quality: 74,
    //     typ: bamreader::mismatch::MismatchType::DEL,
    //     rid: 8,
    // };

    // let test_indel_short_absent = bamreader::mismatch::Mismatch {
    //     chromosome: String::from("chr8"),
    //     position: 1055795,
    //     reference: "CA".as_bytes().to_vec(),
    //     alternative: "C".as_bytes().to_vec(),
    //     quality: 74,
    //     typ: bamreader::mismatch::MismatchType::DEL,
    //     rid: 8,
    // };

    // let test_indel_long_present = bamreader::mismatch::Mismatch {
    //     chromosome: String::from("chr8"),
    //     position: 1056318,
    //     reference: "GTTCCGTGCATGTGGTGGAGGGTGGAA".as_bytes().to_vec(),
    //     alternative: "G".as_bytes().to_vec(),
    //     quality: 74,
    //     typ: bamreader::mismatch::MismatchType::DEL,
    //     rid: 8,
    // };

    // let test_indel_long_absent = bamreader::mismatch::Mismatch {
    //     chromosome: String::from("chr8"),
    //     position: 1055795,
    //     reference: "GTTCCGTGCATGTGGTGGAGGGTGG".as_bytes().to_vec(),
    //     alternative: "G".as_bytes().to_vec(),
    //     quality: 74,
    //     typ: bamreader::mismatch::MismatchType::DEL,
    //     rid: 8,
    // };

    // println!(
    //     "Variant was found in gnomad {}",
    //     gnomad.is_present(test_mismatch_present)
    // );
    // println!(
    //     "Variant was found in gnomad {}",
    //     gnomad.is_present(test_mismatch_absent)
    // );

    // println!(
    //     "Variant was found in gnomad {}",
    //     gnomad.is_present(&test_dbs_present1)
    // );
    // println!(
    //     "Variant was found in gnomad {}",
    //     gnomad.is_present(&test_dbs_present2)
    // );
    // println!(
    //     "Variant was found in gnomad {}",
    //     gnomad.is_present(&test_dbs_absent)
    // );

    // println!(
    //     "Variant was found in gnomad {}",
    //     gnomad.is_present(&test_dbs_absent)
    // );

    // println!(
    //     "Variant was found in gnomad {}",
    //     gnomad.is_present(&test_indel_short_present)
    // );

    // println!(
    //     "Variant was found in gnomad {}",
    //     gnomad.is_present(&test_indel_short_absent)
    // );

    // println!(
    //     "Variant was found in gnomad {}",
    //     gnomad.is_present(&test_indel_long_present)
    // );

    // println!(
    //     "Variant was found in gnomad {}",
    //     gnomad.is_present(&test_indel_long_absent)
    // );

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

    // println!(
    //     "Generated fragment_length_intervals: {:?} and {:?}",
    //     fragment_length_intervals[0], fragment_length_intervals[1]
    // );

    // println!(
    //     "And both 112 and 260 are in the intervals ({}, {})",
    //     fragment_length_intervals[0].contains(&112),
    //     fragment_length_intervals[1].contains(&260)
    // );

    // let header = bam::Read::header(&bam);

    // // this has a nice "normal" mismatch
    // let record = bam::Record::from_sam(header, "A00130:138:H5KNFDSXY:3:2648:5502:1658:GAATCTACGGA\t99\tchr1\t1099541\t21\t72M3D75M\t=\t1099657\t262\tCTAAGAATTGGAGAGAGAGAGGTTAAAATCTCCGACTATGATTCTGGGTTGTCTATTGATGTTTTTGGTCTATTCTAAGAATTGGAGAGAGAGTGGTTAAAATCTCCGACTATGATTGTGGATTGTCTACTGATGTTTTTGGTCTAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF:FFFFFFFFFFFFFF\tXA:Z:chr1,+1099009,147M,9;chr1,+1099387,70M3D77M,10;chr1,+1099313,147M,10;chr1,+1099464,72M3D75M,11;chr1,+1099245,10S43M1D19M3D75M,10;\tMC:Z:48M2D96M3S\tMD:Z:7A35G28^TTG21A12T21G18\tPG:Z:MarkDuplicates\tRG:Z:MBCB196-pDNA\tNM:i:8\tAS:i:113\tXS:i:105".as_bytes()).unwrap();

    // let mate = bam::Record::from_sam(header, "A00130:138:H5KNFDSXY:3:2648:5502:1658:GAATCTACGGA\t147\tchr1\t1099657\t21\t48M2D96M3S\t=\t1099541\t-262\tGATTGTTGATTGTCTACTGATGTTTTTGGTCTAGTGTTCTAAGAATTGGACACAGAGAGGTTAAAATCTCCGATTATGATTGTGGATTGTCTATTGATGTTTTTGGTCTATTGTTCTAAGAATCAGAGAGAGAAGTTAAAATCTACC\tFFFF,:FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:chr1,-1099511,8S115M2D21M3S,11;chr1,-1099352,33M3I87M2D21M3S,14;chr1,-1099122,49M4I80M14S,13;\tMC:Z:72M3D75M\tMD:Z:6G8G17T14^GA2G1G20C22C27G8G10\tPG:Z:MarkDuplicates\tRG:Z:MBCB196-pDNA\tNM:i:11\tAS:i:91\tXS:i:83".as_bytes()).unwrap();

    // // this has a deletion in the overlap
    // let record = bam::Record::from_sam(header, "A00130:138:H5KNFDSXY:3:2550:19144:15515:TCTTATGAGCC\t163\tchr1\t8276826\t60\t51M2I11M2D19M4I33M1D19M7S\t=\t8276830\t136\tTCTCATTCAAAGAAATACACATACACACACCACACACACACAAACACACACATACACCACACACACACATACACCACACACATACACACACCACACACATACACACACCACACACACACACAAACACACATCACACACATACAAAC\t:FFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFF:FFFFFFFFFF:,:FFF,FFFFFFFFFFFFFF:F,FFFFF:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFF:FFFFFFF:,FFFFFF\tMC:Z:47M2I11M2D19M4I33M1D19M12S\tMD:Z:62^AA36C1C9T3^C2C7C8\tPG:Z:MarkDuplicates\tRG:Z:MBCB196-pDNA\tNM:i:14\tAS:i:75\tXS:i:54".as_bytes()).unwrap();

    // let mate = bam::Record::from_sam(header, "A00130:138:H5KNFDSXY:3:2550:19144:15515:TCTTATGAGCC\t83\tchr1\t8276830\t57\t47M2I11M2D19M4I33M1D19M12S\t=\t8276826\t-136\tATTCAAATAAATACACATACACACACCACACACACACAAACACACACATACACCAAACACACACATACACCACACACATACACACACCACACACATACACACACCACACACACACACAAACACACATCACACACATACAAACACACA\tFFF,:FF,:FF:FF,::FFFF,:FFFFFFFFFFFFFFF:FFFF,:F:FFFF:F:F,F:F:F,FFFFFF:FFFF:,,:FFF:F,F,FFFF,FFFF,FFFF:,:FF::,FF::,:FFFFFF::FF:,FFFFFF,,FFFF:FF,FFFF,,\tXA:Z:chr13,+24682149,32S59M56S,2;\tMC:Z:51M2I11M2D19M4I33M1D19M7S\tMD:Z:7G45C4^AA36C1C9T3^C2C7C8\tPG:Z:MarkDuplicates\tRG:Z:MBCB196-pDNA\tNM:i:16\tAS:i:61\tXS:i:49".as_bytes()).unwrap();

    // let r1 = mismatchfinder::bamreader::cigar::parse_cigar_str(&record, &String::from("chr1"));
    // let r2 = mismatchfinder::bamreader::cigar::parse_cigar_str(&mate, &String::from("chr1"));

    // println!(
    //     "{:?}\n{:?}",
    //     String::from_utf8(r1.get_read_seq().to_vec()).unwrap(),
    //     r1.get_read_qual().to_vec()
    // );
    // println!(
    //     "{:?}\n{:?}",
    //     String::from_utf8(r2.get_read_seq().to_vec()).unwrap(),
    //     r2.get_read_qual().to_vec()
    // );

    // let frag =
    //     bamreader::fragment::Fragment::make_fragment(r1, r2, true, true, &String::from("chr1"))
    //         .unwrap();

    // println!(
    //     "read1: {}\n{:?}",
    //     String::from_utf8((frag.get_read1()).get_read_seq().to_vec()).unwrap(),
    //     frag.get_read1().get_read_qual()
    // );

    // println!(
    //     "read2: {}\n{:?}",
    //     String::from_utf8((frag.get_read2()).get_read_seq().to_vec()).unwrap(),
    //     frag.get_read2().get_read_qual()
    // );

    // let mms = frag.get_mismatches(37);

    // for mm in mms {
    //     println!("{}", mm);
    // }

    for bam in cli.bams {
        //prepare file names
        let base = bam.file_stem().unwrap();
        let mut bam = bam::Reader::from_path(&bam).unwrap();
        let vcf_file = cli
            .output_folder
            .join(format!("{}_bamsites.vcf.gz", base.to_str().unwrap()));
        let tsv_file = cli
            .output_folder
            .join(format!("{}_bamsites.tsv", base.to_str().unwrap()));

        if !vcf_file.exists() {
            let mut mismatches = mismatchfinder::bamreader::find_mismatches(
                &mut bam,
                &bed,
                &fragment_length_intervals,
            );
            println!("Found {} mismatches ", mismatches.len());

            // we could potentially only annotate the germline status, but then we still have to write about a million germline vars
            gnomad.filter_germline_cooccurance(&mut mismatches);

            println!("Found {} somatic mismatches", mismatches.len());

            match output::write_mismatches(&mismatches, tsv_file) {
                Err(_) => println!("Could not write mismatch file"),
                Ok(_) => {}
            }

            match output::write_vcf(&mismatches, vcf_file, true, true) {
                Err(_) => println!("Could not write mismatch file"),
                Ok(_) => {}
            }
        }
    }
}
