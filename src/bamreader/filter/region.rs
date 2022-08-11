//filter based on a bed file

use std::{collections::HashMap, hash::Hash, path::PathBuf};

use bio::io::bed;
use rust_lapper::{Interval, Lapper};

enum Backend {
    // the lapper backend uses bool as a value, as it is the smallest and we dont care about the value
    LapperBackend(HashMap<String, Lapper<usize, bool>>),
    IntervalTreeBackend,
}

pub struct BedObject {
    backend: Backend,
}

impl BedObject {
    pub fn has_overlap(&self, chr: &str, query_start: usize, query_end: usize) -> bool {
        match &self.backend {
            // should probably check if we can do seek here, as we should have coordinate sorted reads and beds
            Backend::LapperBackend(core) => {
                let lapper = core.get(chr);

                match lapper {
                    // if we have a lapper, we really just need to check if the iterator is not empty, because we dont care where it overlaps
                    Some(l) => l.find(query_start, query_end).next().is_some(),
                    //if we dont have that chromosome in the storage, we can say we dont overlap
                    None => false,
                }
            }
            Backend::IntervalTreeBackend => todo!("interval tree not supported yet"),
        }
    }

    pub fn lapper_from_bed(bed: PathBuf) -> BedObject {
        let mut bed = bed::Reader::from_file(bed).expect("Could not open bed file");

        // we could assume that the bed file is sorted, but we never know, so we instead store the data per chrsomosome
        // in a hashmap before building the final lapper instead of building it right away

        let mut intervals: HashMap<String, Vec<Interval<usize, bool>>> = HashMap::new();
        for record in bed.records() {
            let record = record.unwrap();

            let mut int_list = intervals
                .entry(record.chrom().to_string())
                .or_insert(Vec::new());
            int_list.push(Interval {
                start: record.start() as usize,
                stop: record.end() as usize,
                val: true,
            })
        }

        //go through all the chromosomes stored and build the lapper representation
        let mut lappers: HashMap<String, Lapper<usize, bool>> = HashMap::new();
        for chrom in intervals.keys() {
            lappers.insert(
                chrom.to_string(),
                Lapper::new(intervals.get(chrom).unwrap().to_vec()),
            );
        }
        //finally just return the full thing
        return BedObject {
            backend: Backend::LapperBackend(lappers),
        };
    }

    pub fn interval_tree_from_bed(bed: PathBuf) -> BedObject {
        todo!("not yet implemented");
    }
}
