use std::collections::BTreeMap;

use echtvar_lib::kmer16;
// function to filter mismatches with a germline resource
// use zarr::{prelude::*, storage::ListableStore};

use crate::bamreader::mismatch::Mismatch;
use crate::bamreader::mismatch::MismatchType;

use echtvar_lib::echtvar;
use echtvar_lib::var32;

pub struct ZarrStorage {}

impl ZarrStorage {
    pub fn filter(&self, _input: Vec<Mismatch>) {}

    pub fn load_zarr_path(root: &str) {
        println!("Trying to open folder {root}");
        todo!("zarr is not su&pported in rust yet");
        // let root = FilesystemHierarchy::open(root).expect("Could not open zarr root");
        // println!("{:?}", root.list_dir("/chr").unwrap());
    }
}

pub enum GermlineBackend {
    Echt(echtvar::EchtVars),
    Zarr(zarr::ArrayMetadata),
}

pub struct GermlineResource {
    backend: GermlineBackend,
}

impl GermlineResource {
    pub fn get_backend_mut(&mut self) -> &mut GermlineBackend {
        &mut self.backend
    }
    pub fn load_echtvars_file(file: &str) -> GermlineResource {
        GermlineResource {
            backend: GermlineBackend::Echt(echtvar::EchtVars::open(file)),
        }
    }

    pub fn is_present(&mut self, mismatch: &Mismatch) -> bool {
        match &mut self.backend {
            GermlineBackend::Echt(backend) => {
                //update where we are in the file
                let _ = backend
                    .set_position(
                        mismatch.rid,
                        mismatch.chromosome.to_string(),
                        mismatch.position,
                    )
                    .expect("Could not change the pointer in germline reference");

                // then we need to check depending on the type of mismatch
                match mismatch.typ {
                    MismatchType::SBS => {
                        // we need to encode the center allele change and search it
                        let enc = var32::encode(
                            // index is also 0 based so we need to substruct 1
                            mismatch.position - 1,
                            &mismatch.reference[1..2],
                            &mismatch.alternative[1..2],
                            &mut 0,
                        );

                        // println!(
                        //     "Checking mismatch: {} with chr: {} start: {} pos: {} and ref: {} to alt: {} with enc: {}",
                        //     mismatch,
                        //     backend.chrom,
                        //     backend.start,
                        //     mismatch.position + 1,
                        //     std::str::from_utf8(&mismatch.reference[1..2]).unwrap(),
                        //     std::str::from_utf8(&mismatch.alternative[1..2]).unwrap(), enc,
                        // );
                        return match backend.var32s.binary_search(&enc) {
                            Ok(_) => true,
                            Err(_) => false,
                        };
                    }
                    MismatchType::DBS => {
                        // here we do the same as for SBS, but for both the first and the second position
                        let mut ret_val = false;
                        for internal_pos in 0usize..2 {
                            let enc = var32::encode(
                                // we need to add  the internal pos but also substract 1 as the first position of the
                                // DBS is already position-1
                                (mismatch.position + internal_pos as u32) - 1,
                                &mismatch.reference[(internal_pos)..(internal_pos + 1)],
                                &mismatch.alternative[internal_pos..(internal_pos + 1)],
                                &mut 0,
                            );

                            // println!(
                            //     "Checking mismatch: {} with chr: {} start: {} pos: {} and ref: {} to alt: {} with enc: {}",
                            //     mismatch,
                            //     backend.chrom,
                            //     backend.start,
                            //     mismatch.position -1 + internal_pos as u32,
                            //     std::str::from_utf8(&mismatch.reference[(internal_pos)..(internal_pos + 1)]).unwrap(),
                            //     std::str::from_utf8(&mismatch.alternative[(internal_pos)..(internal_pos + 1)]).unwrap(), enc,
                            // );

                            match backend.var32s.binary_search(&enc) {
                                Ok(_) => {
                                    ret_val = true;
                                    //we can also leave here
                                    break;
                                }
                                Err(_) => {
                                    //nothing to do here
                                    continue;
                                }
                            };
                        }
                        return ret_val;
                    }
                    MismatchType::DEL | MismatchType::INS => {
                        //in this case, we will have to hope that our deletions are normalised (until we normalise them)
                        // also we need to most likely encode this with the long setting
                        if &mismatch.reference.len() + &mismatch.alternative.len()
                            <= var32::MAX_COMBINED_LEN
                        {
                            let enc = var32::encode(
                                mismatch.position - 1,
                                &mismatch.reference,
                                &&mismatch.alternative,
                                &mut 0,
                            );
                            return match backend.var32s.binary_search(&enc) {
                                Ok(_) => true,
                                Err(_) => false,
                            };
                        } else {
                            let l = var32::LongVariant {
                                position: mismatch.position - 1,
                                sequence: kmer16::encode_var(
                                    &mismatch.reference,
                                    &mismatch.alternative,
                                ),
                                idx: 0,
                            };
                            return match backend.longs.binary_search(&l) {
                                Ok(_) => true,
                                Err(_) => false,
                            };
                        }
                    }
                }
                //then we encode the variant and then we search for it
            }
            GermlineBackend::Zarr(_v) => todo!("Not yet implemented"),
        }
    }

    pub fn filter_germline_cooccurance(&mut self, mismatches: &mut BTreeMap<Mismatch, usize>) {
        mismatches.retain(|key, _| !self.is_present(key));
    }
}
