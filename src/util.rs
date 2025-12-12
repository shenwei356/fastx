use bzip2::read::{BzDecoder, BzEncoder};
use flate2::Compression;
use flate2::read::{GzDecoder, GzEncoder};
use liblzma::read::{XzDecoder, XzEncoder};
use std::fs::File;
use std::io::IsTerminal;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use zstd::stream::read::Decoder as ZstdDecoder;
use zstd::stream::write::Encoder as ZstdEncoder;

pub fn xopen(file: &str, buf_size: usize) -> io::Result<Box<dyn BufRead>> {
    let buf_size = buf_size.max(4096);

    let mut r: Box<dyn BufRead> = if file == "-" {
        if io::stdin().is_terminal() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "no data detected in STDIN",
            ));
        }
        Box::new(BufReader::with_capacity(buf_size, io::stdin().lock()))
    } else {
        Box::new(BufReader::with_capacity(buf_size, File::open(file)?))
    };

    // check compression formats
    let buf = r.fill_buf()?; // peek without consuming

    let reader: Box<dyn BufRead> = if buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b {
        // gzip
        Box::new(BufReader::with_capacity(buf_size, GzDecoder::new(r)))
    } else if buf.starts_with(&[0xFD, b'7', b'z', b'X', b'Z', 0x00]) {
        // xz
        Box::new(BufReader::with_capacity(buf_size, XzDecoder::new(r)))
    } else if buf.starts_with(b"BZh") {
        // bzip2
        Box::new(BufReader::with_capacity(buf_size, BzDecoder::new(r)))
    } else if buf.starts_with(&[0x28, 0xB5, 0x2F, 0xFD]) {
        // zstd
        Box::new(BufReader::with_capacity(buf_size, ZstdDecoder::new(r)?))
    } else {
        // no compression
        r
    };

    Ok(reader)
}

pub fn xwrite(path: &str, buf_size: usize) -> io::Result<Box<dyn Write>> {
    let buf_size = buf_size.max(4096);

    if path == "-" {
        return Ok(Box::new(io::BufWriter::with_capacity(
            buf_size,
            io::stdout(),
        )));
    }

    let file = File::create(path)?;

    let writer: Box<dyn Write> = if path.ends_with(".gz") {
        Box::new(BufWriter::with_capacity(
            buf_size,
            GzEncoder::new(file, Compression::default()),
        ))
    } else if path.ends_with(".xz") {
        Box::new(BufWriter::with_capacity(buf_size, XzEncoder::new(file, 6)))
    } else if path.ends_with(".bz2") {
        Box::new(BufWriter::with_capacity(
            buf_size,
            BzEncoder::new(file, bzip2::Compression::default()),
        ))
    } else if path.ends_with(".zst") || path.ends_with(".zstd") {
        let encoder = ZstdEncoder::new(file, 0)?; // level 0 = default
        Box::new(BufWriter::with_capacity(buf_size, encoder.auto_finish()))
    } else {
        // no compression
        Box::new(BufWriter::with_capacity(buf_size, file))
    };

    Ok(writer)
}

// pub(crate) fn trim_cr(line: &[u8]) -> &[u8] {
//     if let Some((&b'\r', remaining)) = line.split_last() {
//         remaining
//     } else {
//         line
//     }
// }

#[inline(always)]
pub(crate) fn trim_crlf(line: &[u8]) -> &[u8] {
    // let mut line = line;
    // if let Some((&b'\n', remaining)) = line.split_last() {
    //     line = remaining;
    // }
    // if let Some((&b'\r', remaining)) = line.split_last() {
    //     line = remaining;
    // }
    // line

    let mut end = line.len();
    if end > 0 && line[end - 1] == b'\n' {
        end -= 1;
    }
    if end > 0 && line[end - 1] == b'\r' {
        end -= 1;
    }
    &line[..end]
}

// #[inline(always)]
// unsafe fn trim_crlf_mut(v: &mut Vec<u8>) -> &[u8] {
//     let mut len = v.len();
//     if len > 0 && *v.get_unchecked(len - 1) == b'\n' {
//         len -= 1;
//     }
//     if len > 0 && *v.get_unchecked(len - 1) == b'\r' {
//         len -= 1;
//     }
//     v.get_unchecked(0..len)
// }
