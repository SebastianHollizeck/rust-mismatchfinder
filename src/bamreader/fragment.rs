pub struct Fragment {
    read1: GappedRead,
    read2: GappedRead,
    start: i64,
    end: i64,
    single_end: bool,
}

use std::cmp::{max, min};
use std::ptr::read;

use lazy_static::lazy_static;
use regex::Regex;

use crate::bamreader::mismatch::MismatchType;

use super::cigar::GappedRead;
use super::mismatch::Mismatch;

impl Fragment {
    pub fn make_fragment(
        mut read1: GappedRead,
        mut read2: GappedRead,
        only_overlap: bool,
        strict: bool,
    ) -> Option<Fragment> {
        //calculate the offset and the overlap for both reads
        let read1_start = read1.start();
        let read1_end = read1.end();

        let read2_start = read2.start();
        let read2_end = read2.end();

        let overlap = max(read1_start, read2_start)..min(read1_end, read2_end);
        //if there is no overlap, but the only_overlap flag is set, we can just return None here
        let full_fragment;
        if only_overlap && overlap.is_empty() {
            // we just give the original thing back in the end
            return None;
        } else {
            //otherwise we continue building the full fragment
            full_fragment = min(read1_start, read2_start)..max(read1_end, read2_end);
        }

        for ref_pos in full_fragment {
            let read1_pos = read1.get_read_pos_from_ref(ref_pos);
            let read2_pos = read2.get_read_pos_from_ref(ref_pos);

            // if read1_pos.is_none() && read2_pos.is_none() {
            //     //we skip this, as it seems to be a gap in both reads
            //     continue;
            // }
            let read1_pos = match read1_pos {
                Some(v) => *v as i64,
                None => -1,
            };
            let read2_pos = match read2_pos {
                Some(v) => *v as i64,
                None => -1,
            };

            // println!("{read1_pos}, {read2_pos} : {ref_pos}");

            // if one of them is <0 means we have no overlap
            if read1_pos < 0 || read2_pos < 0 {
                //if we only want to analyse the overlap, we just reduce the quality of these positions to 0
                if only_overlap {
                    if read1_pos >= 0 {
                        read1.get_read_qual_mut()[read1_pos as usize] = 0u8;
                    }
                    if read2_pos >= 0 {
                        read2.get_read_qual_mut()[read2_pos as usize] = 0u8;
                    }
                }
                continue;
            }

            //now both reads are overlapping and we can make a consensus

            let r1_seq = read1.get_read_seq()[read1_pos as usize];
            let r2_seq = read2.get_read_seq()[read2_pos as usize];

            // println!(" = {},{}", r1_seq, r2_seq);
            if r1_seq == r2_seq {
                // we update the quality to reflect the additional support
                read1.get_read_qual_mut()[read1_pos as usize] = read1.get_read_qual()
                    [read1_pos as usize]
                    + read2.get_read_qual()[read2_pos as usize];
            } else if strict {
                //if we have the strict rule and there is an overlap, we set the quality to 0 because we do not want
                // to use thise change
                read1.get_read_qual_mut()[read1_pos as usize] = 0u8;
                read2.get_read_qual_mut()[read2_pos as usize] = 0u8;
            } else {
                // here we need to find a compromise between the two reads.
                let r1_qual = read1.get_read_qual()[read1_pos as usize];
                let r2_qual = read2.get_read_qual()[read2_pos as usize];

                if r1_qual > r2_qual {
                    // here we only want to analyse read1, so we set read2's quality to 0, but also reduce quality in 2
                    // we could think about error profile stuff here like in https://github.com/harvardinformatics/NGmerge/blob/master/qual_profile.txt

                    read1.get_read_qual_mut()[read1_pos as usize] = r1_qual - (r2_qual / 2);
                    read2.get_read_qual_mut()[read2_pos as usize] = 0u8;

                    //we can already put the rest if read1 is reverse in here, because if the first fails we want to evaluate that anyways
                    // this makes it less code duplication
                } else if r2_qual > r1_qual || read1.is_reverse() {
                    //same like on top but other way around
                    read1.get_read_qual_mut()[read1_pos as usize] = 0u8;
                    read2.get_read_qual_mut()[read2_pos as usize] = r2_qual - (r1_qual / 2);
                } else {
                    //we know at this part, the quality is exactly the same, so we trust the forward read, which should be less error prone.
                    // we also know that read1 is not reverse, otherwise we would have evaluated the if above, so here we trust read1
                    read1.get_read_qual_mut()[read1_pos as usize] = r1_qual - (r2_qual / 2);
                    read2.get_read_qual_mut()[read2_pos as usize] = 0u8;
                }
            }
        }

        return Some(Fragment {
            read1,
            read2,
            start: min(read1_end, read2_end),
            end: max(read1_end, read2_end),
            single_end: false,
        });
    }

    pub fn get_mismatches(&self, min_base_qual: u8) -> Vec<Mismatch> {
        let mut ret: Vec<Mismatch> = Vec::new();

        let reads_to_check;
        if self.single_end {
            reads_to_check = Vec::from([&self.read1]);
        } else {
            reads_to_check = Vec::from([&self.read1, &self.read2]);
        }

        for read in reads_to_check {
            let md = read.get_md_str();
            //if we have a full integer as MD string we can skip, because we have no changes
            let md_parse: Result<i64, _> = md.parse();
            match md_parse {
                Ok(_) => return ret,
                Err(v) => v, //we have to work
            };

            lazy_static! {
                static ref RE: Regex =
                    Regex::new(r"(?P<digit>[0-9]+)(\^?)(?P<change>[A-Z]+)?").unwrap();
            }

            let ref_pos = read.get_read_pos();

            let mut ref_index = 0;
            //let mut ref_seq = read.seq().encoded.to_vec();

            //because we cant make the iterator peekable, we store the intermediate result, and only insert at the end;
            let mut prev_mismatch: Option<Mismatch> = None;
            let mut prev_pos: usize = 0;

            let upper_bound = read.len() - 1;

            // println!("Upper bound for read length: {upper_bound}");

            for elem in RE.captures_iter(md) {
                //skip all positions which are just matched
                ref_index += &elem["digit"].parse().unwrap();

                // println!("ref_index: {ref_index}");
                // println!(
                //     "Jumping to position {ref_index}: {}",
                //     read.get_start_pos() + ref_index as i64
                // );

                if ref_index >= upper_bound || ref_index == 0 {
                    //there is no way we can make a tri nucleotide here
                    ref_index += 1;
                    continue;
                }

                //if the quality of the base is low, we stop here
                if read.get_read_qual()[ref_index] < min_base_qual {
                    ref_index += 1;
                    continue;
                }

                //store the change
                // if the change is more than one position we have to insert and shift
                // we need to check if we have a change, because we might also just be at the end
                let change = match elem.get(3) {
                    None => {
                        break; //nothing more to do
                    }
                    Some(v) => v.as_str().as_bytes(),
                };

                //if we have a center group, we have a deletion
                let deletion = match elem.get(2) {
                    None => false,
                    Some(v) => {
                        if v.as_str().len() > 0 {
                            true
                        } else {
                            false
                        }
                    }
                };

                // println!("{:?}; deletion: {}", elem, deletion);
                if deletion {
                    // if we have something stored, we push that as well, as there is no DBS with
                    // deletions
                    if let Some(v) = prev_mismatch {
                        ret.push(v);
                    }

                    // we dont have to worry about making a dbs with a deletion, so we push it right away
                    ret.push(Mismatch {
                        chromosome: read.tid(),
                        position: read.start() + ref_index as i64,
                        reference: change.to_vec(),
                        alternative: Vec::new(),
                        quality: 0u8,
                        typ: MismatchType::DEL,
                    });
                    prev_mismatch = None;
                    prev_pos = 0;
                    ref_index -= 1;
                } else {
                    //otherwise we have single base substitutions
                    let change = change[0];

                    //however, we need a continous stretch of dna on the reference
                    if ref_pos[ref_index - 1] + 1 == ref_pos[ref_index]
                        && ref_pos[ref_index] + 1 == ref_pos[ref_index + 1]
                    {
                        //if we have that, we need to check if there is another base change just behind
                        // so we have a DBS instead so we first store the mismatch
                        if let Some(v) = prev_mismatch {
                            //if the mismatch was just before the current one, we join them together
                            if v.position == read.start() + ref_index as i64 - 1 {
                                //unless they were already joined, in which case we drop all of them
                                if v.typ == MismatchType::DBS {
                                    prev_mismatch = None;
                                    prev_pos = ref_index;
                                } else {
                                    //get the previous reference without the first base
                                    let mut tmp_seq = v.reference[1..].to_vec();
                                    //add the next nucleotide
                                    tmp_seq.push(read.get_read_seq()[ref_index + 1]);
                                    //and adjust according to the md string
                                    tmp_seq[1] = change;

                                    prev_mismatch = Some(Mismatch {
                                        quality: (v.quality + read.get_read_qual()[ref_index]) / 2,
                                        chromosome: v.chromosome,
                                        position: read.start() + ref_index as i64,
                                        reference: tmp_seq,
                                        alternative: read.get_read_seq()
                                            [ref_index - 1..ref_index + 2]
                                            .to_vec(),
                                        typ: MismatchType::DBS,
                                    });

                                    prev_pos = ref_index;
                                }
                            } else {
                                // if it isnt, then we push the old and build a new mismatch
                                ret.push(v);

                                //get the read sequence
                                let mut ref_seq =
                                    read.get_read_seq()[ref_index - 1..ref_index + 2].to_vec();
                                //and change according to MD to get the ref string
                                ref_seq[1] = change;

                                prev_mismatch = Some(Mismatch {
                                    quality: read.get_read_qual()[ref_index],
                                    chromosome: read.tid(),
                                    position: read.start() + ref_index as i64,
                                    reference: ref_seq,
                                    alternative: read.get_read_seq()[ref_index - 1..ref_index + 2]
                                        .to_vec(),
                                    typ: MismatchType::SBS,
                                })
                            }
                        } else {
                            //if we have no previous one, we make one but only if the read position ref_index is not exactly the same,
                            //in which case we skip here as well this can really only happen if there are three or more continous changes
                            if ref_index - 1 == prev_pos {
                                prev_pos = ref_index;
                            } else {
                                //get the read sequence
                                let mut ref_seq =
                                    read.get_read_seq()[ref_index - 1..ref_index + 2].to_vec();
                                //and change according to MD to get the ref string
                                ref_seq[1] = change;

                                prev_mismatch = Some(Mismatch {
                                    quality: read.get_read_qual()[ref_index],
                                    chromosome: read.tid(),
                                    position: read.start() + ref_index as i64,
                                    reference: ref_seq,
                                    alternative: read.get_read_seq()[ref_index - 1..ref_index + 2]
                                        .to_vec(),
                                    typ: MismatchType::SBS,
                                })
                            }
                        }
                    } else {
                        // in this case we can also push the old
                        if let Some(v) = prev_mismatch {
                            ret.push(v);
                        }

                        //but we dont really build a new one
                        prev_mismatch = None;
                    }
                    //jump over the nucleotide
                    ref_index += 1;
                }
            }
            // and if we have something left in the end, we push that as well
            if let Some(v) = prev_mismatch {
                ret.push(v);
            }
        }
        return ret;
    }

    pub fn average(numbers: &[u8]) -> f32 {
        //need to cast to a higher number so that we dont overflow
        numbers.iter().map(|x| *x as u64).sum::<u64>() as f32 / numbers.len() as f32
    }
}
