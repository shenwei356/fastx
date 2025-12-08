use bzip2::read::BzDecoder;
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use xz2::read::XzDecoder;
use zstd::stream::read::Decoder as ZstdDecoder;

pub fn xopen(file: &str) -> io::Result<Box<dyn BufRead>> {
    let mut r: Box<dyn BufRead> = if file == "-" {
        Box::new(BufReader::new(io::stdin().lock()))
    } else {
        Box::new(BufReader::new(File::open(file)?))
    };

    // check compression formats
    let buf = r.fill_buf()?; // peek without consuming

    let reader: Box<dyn BufRead> = if buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b {
        // gzip
        Box::new(BufReader::new(GzDecoder::new(r)))
    } else if buf.starts_with(&[0xFD, b'7', b'z', b'X', b'Z', 0x00]) {
        // xz
        Box::new(BufReader::new(XzDecoder::new(r)))
    } else if buf.starts_with(b"BZh") {
        // bzip2
        Box::new(BufReader::new(BzDecoder::new(r)))
    } else if buf.starts_with(&[0x28, 0xB5, 0x2F, 0xFD]) {
        // zstd
        Box::new(BufReader::new(ZstdDecoder::new(r)?))
    } else {
        // no compression
        r
    };

    Ok(reader)
}

// pub(crate) fn trim_cr(line: &[u8]) -> &[u8] {
//     if let Some((&b'\r', remaining)) = line.split_last() {
//         remaining
//     } else {
//         line
//     }
// }

pub(crate) fn trim_crlf(line: &[u8]) -> &[u8] {
    let mut line = line;

    if let Some((&b'\n', remaining)) = line.split_last() {
        line = remaining;
    }

    if let Some((&b'\r', remaining)) = line.split_last() {
        line = remaining;
    }

    line
}
