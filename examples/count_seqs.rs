//! Read FASTA/Q records from stdin and write the number of records, total bases
//!
//! Usage:
//!     cat input.fa          | cargo run --release --example count_seqs
//!     echo -ne ">s\nATCG\n" | cargo run --release --example count_seqs
//!     zcat input.fq.gz      | cargo run --release --example count_seqs

use std::error::Error;

use fastseq::Reader;

fn main() -> Result<(), Box<dyn Error>> {
    // Read FASTA/Q records from stdin and write the number of records, total bases

    let in_file = "-";

    let mut reader = Reader::new(in_file)
        .map_err(|e| format!("failed to parse input file: {}: {}", in_file, e))?;

    let mut total_records: u64 = 0;
    let mut total_bases: u64 = 0;

    while let Some(res) = reader.next() {
        let seq = res?;

        total_records += 1;
        total_bases += seq.seq.len() as u64;
    }

    print!("{}\t{}\n", total_records, total_bases);

    Ok(())
}
