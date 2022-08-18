use std::{cmp::Ordering, fmt};

#[derive(Debug, Hash, Clone, Eq, Ord, PartialOrd)]
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

#[derive(Debug, Hash, Clone)]
pub struct Mismatch {
    pub chromosome: String,
    pub position: u32,
    pub reference: Vec<u8>,
    pub alternative: Vec<u8>,
    pub quality: u8,
    pub typ: MismatchType,
    pub rid: i32,
}

impl fmt::Display for Mismatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{:?}",
            self.chromosome,
            self.position,
            String::from_utf8(self.reference.to_vec()).unwrap(),
            String::from_utf8(self.alternative.to_vec()).unwrap(),
            self.quality,
            self.typ
        )
    }
}

impl PartialOrd for Mismatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Mismatch {
    fn cmp(&self, other: &Self) -> Ordering {
        // check first if the chromsosomes are the same
        if self.chromosome == other.chromosome {
            if self.position < other.position {
                return Ordering::Less;
            } else if self.position > other.position {
                return Ordering::Greater;
            } else {
                if self.reference == other.reference && self.alternative == other.alternative {
                    return Ordering::Equal;
                } else {
                    // if they are the same position, we just compare the length of the variant
                    let len_cmp = (self.reference.len() + self.alternative.len())
                        .cmp(&(other.reference.len() + other.alternative.len()));

                    match len_cmp {
                        Ordering::Equal => self.reference.cmp(&other.reference),
                        Ordering::Greater => Ordering::Greater,
                        Ordering::Less => Ordering::Less,
                    }
                }
            }
        } else {
            self.chromosome.cmp(&other.chromosome)
        }
    }

    fn max(self, other: Self) -> Self {
        // HACK(fee1-dead): go back to using `self.max_by(other, Ord::cmp)`
        // when trait methods are allowed to be used when a const closure is
        // expected.
        match self.cmp(&other) {
            Ordering::Less | Ordering::Equal => other,
            Ordering::Greater => self,
        }
    }

    fn min(self, other: Self) -> Self {
        // HACK(fee1-dead): go back to using `self.min_by(other, Ord::cmp)`
        // when trait methods are allowed to be used when a const closure is
        // expected.
        match self.cmp(&other) {
            Ordering::Less | Ordering::Equal => self,
            Ordering::Greater => other,
        }
    }

    fn clamp(self, min: Self, max: Self) -> Self {
        assert!(min <= max);
        if self < min {
            min
        } else if self > max {
            max
        } else {
            self
        }
    }
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

// our partial eq already satisfies this
impl Eq for Mismatch {}
