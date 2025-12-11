use std::io;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FastxErr {
    #[error("I/O error")]
    IOError(#[from] io::Error),

    #[error("invalid FASTA/Q format")]
    InvalidFormat,

    #[error("unequal lengths of sequence ({0}) and quality ({1})")]
    UnequalSeqAndQual(usize, usize),
}
