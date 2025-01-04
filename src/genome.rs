//! Functionality for referencing genome sequences
//! 
//! Currently, this includes sequence ids (chromosome, scaffold id etc.) and genomic ranges.

use serde::{Deserialize, Serialize};
use std::{cmp::Ordering, fmt, ops::{Range, RangeInclusive}};

use crate::format::{Gff3Row, BedRow};

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("invalid arguments: {0}")]
    InvalidArguments(String),
}

/// Refers to a genomic sequence with an ID e.g. chromosome, scaffold, contig etc. 
/// 
/// Currently, any string is accepted, although strings that can be converted to unsigned integers
/// are ordered numerically before the rest.
#[derive(Debug, Clone, Hash, PartialEq, Eq, Deserialize, Serialize)]
pub struct SeqId(String);

impl SeqId {
    pub fn as_str(&self) -> &str {
        self.0.as_str()
    }
}

impl From<&str> for SeqId {
    fn from(value: &str) -> Self {
        SeqId(value.to_string())
    }
}

impl From<String> for SeqId {
    fn from(value: String) -> Self {
        SeqId(value)
    }
}

impl fmt::Display for SeqId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// [`SeqId`]s are ordered numerically when they can be converted to an unsigned integer
impl Ord for SeqId {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.0.cmp(&other.0) {
            Ordering::Equal => Ordering::Equal,
            ord => {
                // comare numerically if convertable to an unsigned integer
                if let (Ok(seqnum), Ok(other_seqnum)) = (self.0.parse::<u32>(), other.0.parse::<u32>()) {
                    seqnum.cmp(&other_seqnum) 
                } else {
                    ord
                }
            }
        }
    }
}

impl PartialOrd for SeqId {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Stores a genomic range on a specific sequence
///
/// Provides methods for accessing coordinates that require the caller
/// to be explicit about the type of range required. This is to avoid a 
/// common source of errors in genomics where different formats and standards
/// are used to specify ranges.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GenomicRange {
    // these are stored as 0-based, open on the right
    seqid: SeqId,
    start: u64,
    end: u64
}

impl TryFrom<&str> for GenomicRange {
    type Error = Error;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let values: Vec<_> = value.split(':').collect();
        if values.len() == 2 {
            let seqid = values[0];
            let range = values[1];
            let bounds: Vec<_> = range.split('-').collect();
            if bounds.len() == 2 {
                let start = bounds[0].parse::<u64>();
                let end = bounds[1].parse::<u64>();
                if let (Ok(start), Ok(end)) = (start, end) {
                    let location = GenomicRange::from_1closed(seqid, start..=end)?;
                    return Ok(location);
                }
            }
        } 
        Err(Error::InvalidArguments("location should be in the form <seqid>:<start>-<end>".to_string()))
    }
}

impl Ord for GenomicRange {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.seqid.cmp(&other.seqid) {
            Ordering::Equal => {}
            ord => return ord 
        }
        match self.start.cmp(&other.start) {
            core::cmp::Ordering::Equal => {}
            ord => return ord,
        }
        self.end.cmp(&other.end)
    }
}

impl PartialOrd for GenomicRange {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl GenomicRange {
    pub fn combine(&self, other: &GenomicRange) -> Result<GenomicRange, Error> {
        if self.seqid != other.seqid {
            Err(Error::InvalidArguments("can only combine GenomicRanges with the same seqid".to_string()))
        } else {
            let start = self.start.min(other.start);
            let end = self.end.max(other.end);
            Ok(GenomicRange {
                seqid: self.seqid.clone(),
                start,
                end
            })
        }
    }

    pub fn from_0halfopen<T: Into<SeqId>>(seqid: T, range: Range<u64>) -> Result<GenomicRange, Error> {
        let grange = GenomicRange {
            seqid: seqid.into(),
            start: range.start,
            end: range.end,
        };
        Ok(grange)
    }

    pub fn from_1closed<T: Into<SeqId>>(seqid: T, range: RangeInclusive<u64>) -> Result<GenomicRange, Error> {
        if *range.start() == 0 {
            return Err(Error::InvalidArguments("1-based coordinates can't start with 0".to_string()));
        }
        let grange = GenomicRange {
            seqid: seqid.into(),
            start: *range.start() - 1,
            end: *range.end(),
        };
        Ok(grange)
    }
    
    pub fn from_gff_row<T>(row: &Gff3Row<T>) -> GenomicRange {
        GenomicRange {
            seqid: row.seqid.clone(),
            start: row.start - 1,
            end: row.end,
        }
    }

    pub fn from_bed_row<T: BedRow>(row: &T) -> GenomicRange {
        GenomicRange {
            seqid: row.chrom().to_owned(),
            start: row.chrom_start(),
            end: row.chrom_end(),
        }
    }

    pub fn seqid(&self) -> &SeqId {
        &self.seqid
    }

    pub fn range_1closed(&self) -> RangeInclusive<u64> {
        (self.start + 1)..=(self.end)
    }

    pub fn range_0halfopen(&self) -> Range<u64> {
        (self.start)..(self.end)
    }

    pub fn range_0closed(&self) -> RangeInclusive<u64> {
        assert!(self.end > self.start);
        (self.start)..=(self.end - 1) 
    }

}


