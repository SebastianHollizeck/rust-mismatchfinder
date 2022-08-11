// function to filter mismatches with a germline resource
use zarr::{prelude::*, storage::ListableStore};

use crate::bamreader::mismatch::Mismatch;

pub struct ZarrStorage {}

impl ZarrStorage {
    pub fn filter(&self, input: Vec<Mismatch>) {}

    pub fn load_zarr_path(root: &str) {
        println!("Trying to open folder {root}");
        todo!("zarr is not supported in rust yet");
        let root = FilesystemHierarchy::open(root).expect("Could not open zarr root");
        println!("{:?}", root.list_dir("/chr").unwrap());
    }
}
