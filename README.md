[![Docs](https://img.shields.io/docsrs/fastseq)](https://docs.rs/fastseq)
[![Crates.io](https://img.shields.io/crates/v/fastseq.svg)](https://crates.io/crates/fastseq)

## fastseq

fastseq is a high-performance Rust crate for parsing FASTA/Q sequences, for learning purposes.

It seemlessly parses FASTA/Q records from either plain files, (gzip, xz, bzip2, zstd, and lz4) compressed files or STDIN.

## Examples

```rust
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
```


## Benchmarking

see [fastx-tools](https://github.com/shenwei356/fastx-tools)


## Reference

- [needletail](https://github.com/onecodex/needletail)

## License

[MIT License](https://github.com/shenwei356/fastseq/blob/master/LICENSE)
