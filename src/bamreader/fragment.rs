pub struct Fragment {
    read1: GappedRead,
    read2: GappedRead,
    chrom: String,
    start: i64,
    end: i64,
    single_end: bool,
    insertions: Vec<Mismatch>,
    // we stroe this so we know if we only need to look at read1 instead of both
    only_overlap: bool,
}

use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};

use lazy_static::lazy_static;
use regex::Regex;

use crate::bamreader::mismatch::MismatchType;

use super::cigar::GappedRead;
use super::mismatch::Mismatch;

impl Fragment {
    pub fn get_read1(&self) -> &GappedRead {
        &self.read1
    }

    pub fn get_read2(&self) -> &GappedRead {
        &self.read2
    }

    pub fn make_fragment(
        mut read1: GappedRead,
        mut read2: GappedRead,
        only_overlap: bool,
        strict: bool,
        seqname: &String,
    ) -> Option<Fragment> {
        //calculate the offset and the overlap for both reads
        let read1_start = read1.start();
        let read1_end = read1.end();

        let read2_start = read2.start();
        let read2_end = read2.end();

        // println!("r1: {read1_start} - {read1_end}\nr2: {read2_start} - {read2_end}");

        let overlap = max(read1_start, read2_start)..min(read1_end, read2_end);

        // println!("Overlap: {:?}", overlap);
        //if there is no overlap, but the only_overlap flag is set, we can just return None here
        let full_fragment;
        if only_overlap && overlap.is_empty() {
            // we just give the original thing back in the end
            return None;
        } else {
            //otherwise we continue building the full fragment
            full_fragment = min(read1_start, read2_start)..max(read1_end, read2_end);
        }
        // println!("Full fragment: {:?}", full_fragment);

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

            // print!("{read1_pos}, {read2_pos} : {ref_pos}");

            // if one of them is <0 means we have no overlap
            if read1_pos < 0 || read2_pos < 0 {
                //if we only want to analyse the overlap, we just reduce the quality of these positions to 0
                if only_overlap {
                    if read1_pos >= 0 {
                        read1.get_read_qual_mut()[read1_pos as usize] = 0u8;
                        // println!(
                        //     " = {},-",
                        //     std::str::from_utf8(&[read1.get_read_seq()[read1_pos as usize]])
                        //         .unwrap()
                        // );
                    } else if read2_pos >= 0 {
                        read2.get_read_qual_mut()[read2_pos as usize] = 0u8;
                        // println!(
                        //     " = -,{}",
                        //     std::str::from_utf8(&[read2.get_read_seq()[read2_pos as usize]])
                        //         .unwrap()
                        // );
                    } else {
                        //this is a hap situation in the middle of a read
                        // println!(" = -,-");
                    }
                }
                continue;
            }

            //now both reads are overlapping and we can make a consensus

            let r1_seq = read1.get_read_seq()[read1_pos as usize];
            let r2_seq = read2.get_read_seq()[read2_pos as usize];

            // println!(
            //     " = {},{}",
            //     std::str::from_utf8(&[r1_seq]).unwrap(),
            //     std::str::from_utf8(&[r2_seq]).unwrap()
            // );
            if r1_seq == r2_seq {
                // we update the quality to reflect the additional support in read1
                read1.get_read_qual_mut()[read1_pos as usize] = read1.get_read_qual()
                    [read1_pos as usize]
                    + read2.get_read_qual()[read2_pos as usize];
                // we set the quality to 0 in the second read, to not get the mismatch twice
                read2.get_read_qual_mut()[read2_pos as usize] = 0u8;
            } else if strict {
                //if we have the strict rule and there is an overlap, we set the quality to 0 because we do not want
                // to use thise change
                // println!("setting both qualities to 0 ");
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

        //and now we check the inserts for overlaps
        // create final storage
        let mut fragment_ins: HashMap<u32, Mismatch> = HashMap::new();

        // we need this because we cant be sure that we find something that we have in read1 but not in read2 otherwise
        let mut insert_pos: HashSet<u32> = HashSet::new();

        for mm in read1.get_insertions().iter().cloned() {
            // if we only want the overlap, we dont store things that arent in the overlap
            if only_overlap && !overlap.contains(&(mm.position as i64)) {
                continue;
            } else {
                insert_pos.insert(mm.position);
                fragment_ins.insert(mm.position, mm);
            }
        }

        //going through read2 and checking if we already had it
        for mm in read2.get_insertions().iter().cloned() {
            // if we only want the overlap, we dont store things that arent in the overlap
            if only_overlap && !overlap.contains(&(mm.position as i64)) {
                continue;
            } else {
                if fragment_ins.contains_key(&mm.position) {
                    // here we update the mismatch quality to reflect the additional support
                    let entry = fragment_ins.get_mut(&mm.position).unwrap();

                    entry.quality = entry.quality + mm.quality;

                    //remove the insert_pos entry, so we dont delete it afterwards
                    insert_pos.remove(&mm.position);
                } else {
                    if strict {
                        // we only want double support
                        continue;
                    } else {
                        // add the insert in
                        fragment_ins.insert(mm.position, mm);
                    }
                }
            }
        }

        //go through the insert_pos map and delete all remaining insertions
        for pos in insert_pos {
            fragment_ins.remove(&pos);
        }

        return Some(Fragment {
            read1,
            read2,
            chrom: seqname.to_string(),
            start: min(read1_end, read2_end),
            end: max(read1_end, read2_end),
            single_end: false,
            insertions: fragment_ins.into_values().collect(),
            only_overlap: only_overlap,
        });
    }

    pub fn get_mismatches(&self, min_base_qual: u8) -> Vec<Mismatch> {
        let mut ret: Vec<Mismatch> = Vec::new();

        let reads_to_check;
        // if we only trust the overlap, or have a single end read library, we only look at read1
        if self.single_end || self.only_overlap {
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

            //because we cant make the iterator peekable, we store the intermediate result, and only insert at the end;
            let mut prev_mismatch: Option<Mismatch> = None;
            let mut prev_pos: usize = 0;

            let upper_bound = read.len() - 1;

            // println!("Upper bound for read length: {upper_bound}");

            for elem in RE.captures_iter(md) {
                //skip all positions which are just matched
                ref_index += &elem["digit"].parse().unwrap();

                // println!(
                //     "Jumping to position {ref_index}: {}",
                //     read.start() + ref_index as i64
                // );

                if ref_index >= upper_bound || ref_index == 0 {
                    //there is no way we can make a tri nucleotide here
                    ref_index += 1;
                    continue;
                }

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

                //if the quality of the base is low, we stop here
                if !deletion && read.get_read_qual()[ref_index] < min_base_qual {
                    // println!("Low basequal {}", read.get_read_qual()[ref_index]);
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

                // println!("{:?}; deletion: {}", elem, deletion);
                if deletion {
                    // if we have something stored, we push that as well, as there is no DBS with
                    // deletions
                    if let Some(v) = prev_mismatch {
                        ret.push(v);
                    }

                    //we use the quality of the neighbouring bases as an indication of the qual of the deletion
                    let qual = (read.get_read_qual()[ref_index - 1]
                        + read.get_read_qual()[ref_index + change.len()])
                        / 2;
                    if qual > min_base_qual {
                        // we take the base just before the deletion as the alt
                        let alt_seq = read.get_read_seq()[(ref_index - 1)..ref_index].to_vec();
                        // and for the reference, we have to append what the read is missing
                        let mut ref_seq = alt_seq.clone();
                        ref_seq.append(change.to_vec().as_mut());

                        // we dont have to worry about making a dbs with a deletion, so we push it right away
                        ret.push(Mismatch {
                            chromosome: self.chrom.to_string(),
                            rid: read.tid(),
                            position: (read.start() + ref_index as i64) as u32,
                            // we use the read sequence and add on the change
                            reference: ref_seq,
                            alternative: alt_seq,

                            quality: qual,
                            typ: MismatchType::DEL,
                        });
                    }
                    prev_mismatch = None;
                    prev_pos = 0;
                    //step over the deletion to be at the base after
                    ref_index += change.len();
                } else {
                    //otherwise we have single base substitutions
                    let change = change[0];

                    //however, we need a continous stretch of dna on the reference
                    if ref_pos[ref_index - 1] + 1 == ref_pos[ref_index]
                        && ref_pos[ref_index] + 1 == ref_pos[ref_index + 1]
                    {
                        //get the read sequence (which is the alternative)
                        let mut alt_seq =
                            read.get_read_seq()[(ref_index - 1)..(ref_index + 2)].to_vec();

                        // if we have a - on either side it means there is a deletion right next to the SNP
                        // which means we wont trust it anymore
                        if alt_seq[0] == b'-' || alt_seq[2] == b'-' {
                            continue;
                        }

                        let mut ref_seq = alt_seq.clone();
                        //and change according to MD to get the ref string
                        ref_seq[1] = change;

                        //if we have that, we need to check if there is another base change just behind
                        // so we have a DBS instead so we first store the mismatch
                        if let Some(mm) = prev_mismatch {
                            //if the mismatch was just before the current one, we join them together
                            if mm.position == (read.start() + ref_index as i64 - 1) as u32 {
                                //unless they were already joined, in which case we drop all of them
                                if mm.typ == MismatchType::DBS {
                                    prev_mismatch = None;
                                    prev_pos = ref_index;
                                } else {
                                    //alt_seq will be shortened to only 2 bases
                                    alt_seq.pop();
                                    // but we need to update the ref_seq with the center from the previous
                                    ref_seq[0] = mm.reference[1];
                                    // and shorten as well
                                    ref_seq.pop();

                                    prev_mismatch = Some(Mismatch {
                                        quality: (mm.quality + read.get_read_qual()[ref_index]) / 2,
                                        chromosome: mm.chromosome,
                                        rid: read.tid(),
                                        position: (read.start() + ref_index as i64) as u32,
                                        reference: ref_seq,
                                        alternative: alt_seq,
                                        typ: MismatchType::DBS,
                                    });

                                    prev_pos = ref_index;
                                }
                            } else {
                                // if it isnt, then we push the old and build a new mismatch
                                ret.push(mm);

                                prev_mismatch = Some(Mismatch {
                                    quality: read.get_read_qual()[ref_index],
                                    chromosome: self.chrom.to_string(),
                                    rid: read.tid(),
                                    position: (read.start() + ref_index as i64 + 1) as u32,
                                    reference: ref_seq,
                                    alternative: alt_seq,
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
                                    chromosome: self.chrom.to_string(),
                                    rid: read.tid(),
                                    position: (read.start() + ref_index as i64 + 1) as u32,
                                    reference: ref_seq,
                                    alternative: alt_seq,
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

                    // step over the mismatch
                    ref_index += 1;
                }
            }
            // and if we have something left in the end, we push that as well
            if let Some(v) = prev_mismatch {
                ret.push(v);
            }
        }

        // go through all of the insertions we already collected before and add them to the final result if they are of high enough qual
        for mm in self.insertions.iter().cloned() {
            if mm.quality >= min_base_qual {
                ret.push(mm);
            }
        }

        return ret;
    }

    pub fn average(numbers: &[u8]) -> f32 {
        //need to cast to a higher number so that we dont overflow
        numbers.iter().map(|x| *x as u64).sum::<u64>() as f32 / numbers.len() as f32
    }
}
