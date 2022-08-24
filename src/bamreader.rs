use std::{
    cmp::{max, min},
    collections::{BTreeMap, HashMap},
    ops::RangeInclusive,
};

use log::{debug, info, warn};

use rust_htslib::bam::{self, record::Aux, Read, Record};

use crate::bamreader::{
    cigar::{parse_cigar_str},
    fragment::Fragment,
};

use self::{
    filter::{region::BedObject},
    mismatch::Mismatch,
};

pub mod cigar;
pub mod filter;
pub mod fragment;
pub mod mismatch;

pub fn find_mismatches(
    bam: &mut bam::Reader,
    white_list: &BedObject,
    fragment_length_intervals: &Vec<RangeInclusive<i64>>,
) -> BTreeMap<Mismatch, usize> {
    //we create a read cache and set the initial capacity to 10k to reduce the reallocation operations
    let mut read_cache: HashMap<String, Record> = HashMap::with_capacity(10000);

    //get the header for the mapping of tid to chromosome name
    let mut tid_map: HashMap<i32, String> = HashMap::new();
    {
        let header = bam.header();
        let seqnames = header.target_names().into_iter();

        for (id, chr) in seqnames.enumerate() {
            tid_map.insert(
                id as i32,
                String::from_utf8(chr.to_vec()).expect("We couldnt fonvert the tid to String"),
            );
        }
    }

    //store all mismatches found and how often they were found (we also know its going to be a big hash)
    let mut mismatch_store: BTreeMap<Mismatch, usize> = BTreeMap::new();

    let mut counter = 0;
    let mut last_chr = &String::from("*");
    let mut last_pos = -1;
    for r in bam.records() {
        let record = r.unwrap();

        let strict = true;
        let overlap_only = true;

        let min_bq = 65;
        let min_avg_bq = 25.;
        let min_mq = 20;

        let qname = std::str::from_utf8(record.qname()).unwrap().to_owned();

        //println!("Read name: {qname}");
        counter += 1;
        if counter % 500000 == 0 {
            info!(
                    "Read through {counter} reads - last position: {}:{}",
                    last_chr ,
                    last_pos
                );
            
        
        }

        if read_cache.contains_key(&qname) {

            // we have to skip out of this if the read is a secondary (before we get the read from the cache)
            if record.is_secondary() || record.is_supplementary(){
                continue;
            }

            // get and delete the record from the cache, because we dont want to bloat the storage
            let mate = read_cache.remove(&qname).unwrap();

            // println!("r1 {:?} ({counter})", std::str::from_utf8(mate.qname()),);
            //println!("r2 {:?}({counter})", record,);
            debug!("Working on read: {} and {}", std::str::from_utf8(record.qname()).unwrap(), std::str::from_utf8(mate.qname()).unwrap());

            // we check all the quality of the read AND a few for the mate, because we need to be sure we have a proper
            // pair 
            if !(record.is_unmapped() 
                || mate.is_unmapped()
                || record.is_duplicate()
                || mate.is_duplicate()
                || record.is_quality_check_failed()
                || mate.is_quality_check_failed()
                || record.tid() != record.mtid()
                // we use the average mapping quality of the read as the mapping quality of the fragment
                || (record.mapq() + mate.mapq()) / 2 < min_mq)
            {

                // get the chromosome the record is on (because we know they are both on the same)
                let chrom = tid_map.get(&record.tid()).unwrap();
                last_chr = chrom;
                last_pos = record.pos();

                //if the read has any additional mapping locations, we cant really trust the alignment as much
                if let Ok(_) = record.aux(b"XA") {
                    continue;
                }
                if let Ok(_) = mate.aux(b"XA") {
                    continue;
                }

                // we check if the reads actually have any changes (edit distance), this contains both mismatches and
                // insertions or deletions
                let read1_edit = match record.aux(b"NM") {
                    // while technically the thing should only be I8, we accept other integer types as well
                    Ok(Aux::I8(nm_i)) => nm_i,
                    Ok(Aux::I16(nm_i)) => nm_i as i8,
                    Ok(Aux::I32(nm_i)) => nm_i as i8,
                    Ok(Aux::U8(nm_i)) => nm_i as i8,
                    Ok(Aux::U16(nm_i)) => nm_i as i8,
                    Ok(Aux::U32(nm_i)) => nm_i as i8,
                    // we just return the cigar ops necessary to make the read
                    _ => 0,
                };

                // same as above, but with read2
                let read2_edit = match mate.aux(b"NM") {
                    Ok(Aux::I8(nm_i)) => nm_i,
                    Ok(Aux::I16(nm_i)) => nm_i as i8,
                    Ok(Aux::I32(nm_i)) => nm_i as i8,
                    Ok(Aux::U8(nm_i)) => nm_i as i8,
                    Ok(Aux::U16(nm_i)) => nm_i as i8,
                    Ok(Aux::U32(nm_i)) => nm_i as i8,
                    _ => 0,
                };



                if read1_edit > 0 || read2_edit > 0 {
                    //now we check if the fragment has the right size
                    let frag_size = record.insert_size().abs();
                    // println!("Found fragment of size {frag_size}");
                    let mut skip = true;
                    for ivl in fragment_length_intervals {
                        if ivl.contains(&frag_size) {
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

                    // get only the aligned part of the read, without insertions
                    let read1 = parse_cigar_str(&record, chrom);
                    let read2 = parse_cigar_str(&mate, chrom);




                    if white_list.has_overlap(
                        chrom,
                        min(read1.start() as usize, read2.start() as usize),
                        max(read1.end() as usize, read2.end() as usize),
                    ) {
                        // println!("Found overlap in whitelist");

                        let res;
                        if read1.get_read_pos() <= read2.get_read_pos() {
                            res =
                                Fragment::make_fragment(read1, read2, overlap_only, strict, chrom);
                        } else {
                            res =
                                Fragment::make_fragment(read2, read1, overlap_only, strict, chrom);
                        }

                        let mismatches = match res {
                            Some(v) => v.get_mismatches(min_bq),
                            None => Vec::new(),
                        };

                        debug!("Found {} mismatches in fragment", mismatches.len());

                        for mm in mismatches {
                            // add the mismatch to the storage if it wasnt already
                            let val = mismatch_store.entry(mm).or_insert(0);
                            //and add one to the count
                            *val +=1;
                            
                        }
                    } else {
                        // println!("Not in whitelist");
                    }
                } else {
                    // println!("No mismatches");
                }
            } else {
                // println!("Low quality");
            }

            debug!("Done with read {}", &qname);
        // we dont check if the read is mapped here, so we get mapped and unmapped pairs as well
        } else if !(record.is_supplementary()
            || record.is_secondary()) {
            read_cache.insert(qname, record);
        } else {
            // we do discard any additional reads, which are not primary alignment or not mapped to the same chromosome
        }
    }

    // if we have done everything correctly, we should have no reads in the cache anymore
    if read_cache.len() != 0 {
        warn!("Read cache contained unpaired read at the end of the analysis, this shouldnt happen with a well formed bam");
        for (qname, _) in read_cache {
            warn!("{qname}")
        }
    }

    return mismatch_store;
}
