//! Supports different genomic formats
//!
//! This module provides basic support for reading and writing GFF3 and BED files.
//! These implementations are not intended to be general and comprehensive.
//! 
use std::fmt;
use indexmap::IndexMap;
use num_traits::NumOps;
use serde::{Deserialize, Deserializer, Serialize};

use crate::genome::{GenomicRange, SeqId};

/// The standard fields of GFF3
///
/// Limitations are:
/// - The text encoding must be UTF-8
/// - Percent encodings are not converted on deserialize and are not used in serialize
/// - Validation is limited to type e.g. String, u64 etc.
///     - [seqid](Gff3Row::seqid), [score](Gff3Row::score), and [phase](Gff3Row::phase) allow any string
#[derive(Debug, Deserialize, Serialize)]
pub struct Gff3Row<T> {
    pub seqid: SeqId,
    pub source: String,
    pub feature_type: T,
    pub start: u64,
    pub end: u64,
    pub score: String,
    pub strand: Strand,
    pub phase: String,
    #[serde(deserialize_with = "deserialize_attributes")]
    pub attributes: IndexMap<String, String>,
}

fn deserialize_attributes<'de, D>(deserializer: D) -> Result<IndexMap<String, String>, D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = Deserialize::deserialize(deserializer)?;
    let mut map = IndexMap::new();

    for kv in s.split(';') {
        let mut iter = kv.splitn(2, '=');
        if let (Some(key), Some(value)) = (iter.next(), iter.next()) {
            map.insert(key.to_string(), value.to_string());
        }
    }

    Ok(map)
}

/// The genome strand the annotation is associated with
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum Strand {
  #[serde(rename = "+")]
  Plus,
  #[serde(rename = "-")]
  Minus,
  #[serde(rename = ".")]
  None,  
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", match self {
            Strand::Plus => "+",
            Strand::Minus => "-",
            Strand::None => ".",
        })
    }
}

/// The standard fields shared by all BED files
pub trait BedRow {
    fn chrom(&self) -> &SeqId;
    fn chrom_start(&self) -> u64;
    fn chrom_end(&self) -> u64;
}

/// The standard fields of BED6
///
/// Validation is by type only e.g. u64 or u16.
#[derive(Deserialize)]
pub struct Bed6Row {
    pub chrom: SeqId,
    pub chrom_start: u64,
    pub chrom_end: u64,
    pub name: String,
    pub score: u16,
    pub strand: Strand,
}

impl BedRow for Bed6Row {
    fn chrom(&self) -> &SeqId {
        &self.chrom
    }

    fn chrom_start(&self) -> u64 {
        self.chrom_start
    }

    fn chrom_end(&self) -> u64 {
        self.chrom_end
    }
}

impl Bed6Row {
    pub fn name(&self) -> &String {
        &self.name
    }

    pub fn score(&self) -> u16 {
        self.score
    }

    pub fn strand(&self) -> Strand {
        self.strand
    }
}

/// A genomic range with zero or more associated data values
pub struct DataInterval<T: NumOps + Copy> {
    range: GenomicRange,
    values: Vec<Option<T>>,
}

impl<T: NumOps + Copy> DataInterval<T> {
    pub fn new(range: GenomicRange, values: Vec<Option<T>>) -> DataInterval<T> {
        DataInterval {
            range,
            values,
        }
    }

    pub fn values(&self) -> &[Option<T>] {
        &self.values
    }

    pub fn range(&self) -> &GenomicRange {
        &self.range
    }
}

/// The standard fields of BedGraph
#[derive(Deserialize, Serialize)]
pub struct BedGraphRow<T: NumOps + Copy> {
    pub chrom: SeqId,
    pub chrom_start: u64,
    pub chrom_end: u64,
    pub data_value: T,
}

impl<T> From<BedGraphRow<T>> for DataInterval<T> where T: NumOps + Copy {

    fn from(row: BedGraphRow<T>) -> Self {
        let range = GenomicRange::from_0halfopen(row.chrom, row.chrom_start..row.chrom_end).unwrap();
        DataInterval {
            range,
            values: vec![Some(row.data_value)]
        }
    }
}

/// BedGraph Extended, supporting zero or more values per row
#[derive(Deserialize, Serialize)]
pub struct BedGraphExtRow<T: NumOps + Copy> {
    pub chrom: SeqId,
    pub chrom_start: u64,
    pub chrom_end: u64,
    #[serde(flatten)]
    pub data_values: Vec<Option<T>>,
}

impl<T> From<BedGraphExtRow<T>> for DataInterval<T> where T: NumOps + Copy {

    fn from(row: BedGraphExtRow<T>) -> Self {
        let range = GenomicRange::from_0halfopen(row.chrom, row.chrom_start..row.chrom_end).unwrap();
        DataInterval {
            range,
            values: row.data_values
        }
    }
}