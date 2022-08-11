#[derive(Debug)]
pub enum MismatchType {
    SBS,
    DBS,
    INS,
    DEL,
}

impl PartialEq for MismatchType {
    fn eq(&self, other: &Self) -> bool {
        core::mem::discriminant(self) == core::mem::discriminant(other)
    }

    fn ne(&self, other: &Self) -> bool {
        !self.eq(other)
    }
}

#[derive(Debug)]
pub struct Mismatch {
    pub chromosome: i32,
    pub position: i64,
    pub reference: Vec<u8>,
    pub alternative: Vec<u8>,
    pub quality: u8,
    pub typ: MismatchType,
}

impl PartialEq for Mismatch {
    fn eq(&self, other: &Mismatch) -> bool {
        //we start with the easy things, this will make the method a bit bloaty, but more performant
        if self.typ != other.typ {
            return false;
        }
        // we do the position before the chromosome, as we wont compare mismatches between chromosomes
        if self.position != other.position {
            return false;
        }
        if self.chromosome != other.chromosome {
            return false;
        }
        // each of the things needs to fit
        for (n_self, n_other) in self.alternative.iter().zip(other.alternative.iter()) {
            if n_self != n_other {
                return false;
            }
        }
        for (n_self, n_other) in self.reference.iter().zip(other.reference.iter()) {
            if n_self != n_other {
                return false;
            }
        }

        //we dont really care about the quality being the same, so we ignore that and return true if we got to here
        return true;
    }

    fn ne(&self, other: &Self) -> bool {
        !self.eq(other)
    }
}
