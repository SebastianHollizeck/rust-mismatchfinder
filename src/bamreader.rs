use std::{
    cmp::{max, min},
    collections::HashMap,
    hash::Hash,
    ops::RangeInclusive,
    process::exit,
};

use rust_htslib::bam::{self, record::Aux, Read, Record};

use crate::bamreader::{
    cigar::{parse_cigar_str, GappedRead},
    fragment::Fragment,
};

use self::filter::region::BedObject;

mod cigar;
pub mod filter;
mod fragment;
mod mismatch;

pub fn find_mismatches(
    bam: &mut bam::Reader,
    white_list: &BedObject,
    fragment_length_intervals: &Vec<RangeInclusive<i64>>,
) {
    let mut read_cache: HashMap<String, Record> = HashMap::new();

    //get the header for the mapping of tid to chromosome name
    let mut tid_map: HashMap<i32, String> = HashMap::new();
    {
        let header = bam.header();
        let seqnames = header.target_names();

        for (id, chr) in seqnames.iter().enumerate() {
            tid_map.insert(
                id as i32,
                String::from_utf8(chr.to_vec()).expect("We couldnt fonvert the tid to String"),
            );
        }
    }

    let mut counter = 0;
    for r in bam.records() {
        let record = r.unwrap();

        let strict = true;
        let overlap_only = true;

        let min_bq = 37;
        let min_avg_bq = 25.;
        let min_mq = 20;

        let qname = std::str::from_utf8(record.qname()).unwrap().to_owned();

        //println!("Read name: {qname}");

        if read_cache.contains_key(&qname) {
            // get and delete the record from the cache, because we dont want to bloat the storage
            let mate = read_cache.remove(&qname).unwrap();

            println!("r1 {:?} ({counter})", std::str::from_utf8(mate.qname()),);
            //println!("r2 {:?}({counter})", record,);

            if !(record.is_unmapped()
                || record.is_duplicate()
                || record.is_supplementary()
                || record.is_secondary()
                || record.is_quality_check_failed()
                || mate.is_unmapped()
                || record.tid() != record.mtid()
                // we use the average mapping quality of the read as the mapping quality of the fragment
                || (record.mapq() + mate.mapq()) / 2 < min_mq)
            {
                // we check if the reads actually have any changes (edit distance) because otherwise what is the point
                let read1_edit = match record.aux(b"NM") {
                    Ok(rust_htslib::bam::record::Aux::I32(nm_i)) => nm_i,
                    _ => 0,
                };

                let read2_edit = match mate.aux(b"NM") {
                    Ok(rust_htslib::bam::record::Aux::I32(nm_i)) => nm_i,
                    _ => 0,
                };

                if read1_edit > 0 || read2_edit > 0 {
                    //now we check if the fragment has the right size
                    let mut skip = true;
                    for ivl in fragment_length_intervals {
                        if ivl.contains(&(record.insert_size().abs())) {
                            //we are good here
                            skip = false;
                        }
                    }
                    if skip {
                        continue;
                    }

                    //then we check for the average base quality of the reads
                    if (Fragment::average(record.qual()) + Fragment::average(mate.qual())) / 2.
                        < min_avg_bq
                    {
                        continue;
                    }

                    let read1 = parse_cigar_str(&record);
                    let read2 = parse_cigar_str(&mate);

                    // we check if the theoretical fragment overlaps the white list
                    if white_list.has_overlap(
                        tid_map.get(&read1.tid()).unwrap(),
                        min(read1.start() as usize, read2.start() as usize),
                        max(read1.end() as usize, read2.end() as usize),
                    ) {
                        // println!("Found overlap with whitelist");

                        let res;
                        if read1.get_read_pos() <= read2.get_read_pos() {
                            res = Fragment::make_fragment(read1, read2, overlap_only, strict);
                        } else {
                            res = Fragment::make_fragment(read2, read1, overlap_only, strict);
                        }

                        let mismatches = match res {
                            Some(v) => v.get_mismatches(min_bq),
                            None => Vec::new(),
                        };

                        // for mm in mismatches {
                        //     println!("{:?}", mm);
                        // }
                    } else {
                        println!("No overlap");
                    }
                } else {
                    println!("No mismatches")
                }
            }
        } else {
            read_cache.insert(qname, record);
        }
        counter += 1;
        if counter >= 2000000 {
            exit(0);
        }
    }
}
