use crate::errors::*;
use crate::seq::*;
use crate::util::*;
use crate::xopen::*;
use memchr::{memchr, memchr2};
use std::io::BufRead;

/// A FASTA/Q reader that can read from any BufRead.
/// It supports both FASTA and FASTQ formats,
/// and can automatically detect the format based on the first non-empty line.
pub struct Reader<R: BufRead> {
    reader: R,

    is_fastq: bool,

    record_buf: Vec<u8>,
    line_buf: Vec<u8>,
    lookahead_line: Vec<u8>,
    has_lookahead: bool,

    parse_id: bool,
}

impl Reader<Box<dyn BufRead>> {
    /// Creates a new Reader from a file path.
    /// Stdin is supported by passing "-" as the file path.
    /// The reader will take ownership of the file reader.
    pub fn new(file: &str) -> Result<Self, FastxErr> {
        Self::new_with_buffer_alignment(file, 65536, DEFAULT_IO_BUFFER_ALIGNMENT)
    }

    /// Creates a new Reader from a file path with specified buffer size and alignment.
    pub fn new_with_buffer_alignment(
        file: &str,
        buf_size: usize,
        buf_align: usize,
    ) -> Result<Self, FastxErr> {
        let r: Box<dyn BufRead> =
            xopen_with_alignment(file, buf_size, buf_align).map_err(FastxErr::IOError)?;
        Ok(Self::from_reader(r))
    }
}

impl<R: BufRead> Reader<R> {
    /// Creates a new Reader from any BufRead. The reader will take ownership of the provided BufRead.
    pub fn from_reader(reader: R) -> Self {
        Self {
            reader,
            is_fastq: false,
            record_buf: Vec::with_capacity(1 << 20),
            line_buf: Vec::with_capacity(1024),
            lookahead_line: Vec::with_capacity(1024),
            has_lookahead: false,
            parse_id: true,
        }
    }

    pub fn skip_id_parsing(&mut self) {
        self.parse_id = false
    }

    // read a line into line_buf, return the number of bytes consumed from the reader
    #[inline(always)]
    fn read_line_fill_buf(&mut self) -> std::io::Result<usize> {
        self.line_buf.clear();

        let mut total = 0;
        loop {
            let (consumed, done) = {
                let buf = self.reader.fill_buf()?;
                if buf.is_empty() {
                    // EOF
                    return Ok(total);
                }

                // find the position of the first '\n' in the buffer
                match memchr(b'\n', buf) {
                    Some(pos) => {
                        // found a line ending
                        let end = pos + 1;
                        self.line_buf.extend_from_slice(&buf[..end]); // copy the line (including '\n') into line_buf
                        (end, true)
                    }
                    None => {
                        // no line ending found, consume the entire buffer and continue reading
                        self.line_buf.extend_from_slice(buf); // copy the last part of the line into line_buf
                        (buf.len(), false)
                    }
                }
            };

            self.reader.consume(consumed);
            total += consumed;

            // found a '\n' or reached EOF
            if done {
                return Ok(total);
            }
            // else continue reading until we find a line ending or reach EOF
        }
    }

    // only for the first record
    #[inline(always)]
    fn read_next_nonempty_line(&mut self) -> Result<bool, FastxErr> {
        loop {
            match self.read_line_fill_buf() {
                Ok(0) => return Ok(false),                                 // EOF
                Ok(_) if trim_crlf(&self.line_buf).is_empty() => continue, // skip blank lines
                Ok(_) => return Ok(true), // non-empty line read successfully
                Err(e) => return Err(FastxErr::IOError(e)), // I/O error occurred
            }
        }
    }

    // read the next non-empty line and append it to record_buf,
    // return the outcome of the read operation
    #[inline(always)]
    fn read_next_nonempty_line_into_record_buf(
        &mut self,
        stop_on_fasta_header: bool,
        stop_on_fastq_sep: bool,
    ) -> Result<ReadLineOutcome, FastxErr> {
        loop {
            let fast_path = {
                // read the buffer without consuming it yet
                let buf = self.reader.fill_buf().map_err(FastxErr::IOError)?;
                if buf.is_empty() {
                    return Ok(ReadLineOutcome::Eof);
                }

                // find the position of the first '\n' in the buffer and determine if we can take the fast path
                memchr(b'\n', buf).map(|pos| {
                    let consumed = pos + 1;
                    let line = trim_crlf(&buf[..consumed]);
                    let first_char = if line.is_empty() { None } else { Some(line[0]) };
                    // fasta header line must be copied to lookahead_line to avoid borrowing issues,
                    // while other lines can be directly processed from the buffer without copying
                    let needs_copy = matches!(first_char, Some(b'>') if stop_on_fasta_header);
                    (consumed, line.len(), first_char, needs_copy)
                })
            };

            // if the fast path is available, we can process the line directly from the buffer without copying to line_buf
            if let Some((consumed, line_len, first_char, needs_copy)) = fast_path {
                // skip blank lines
                if line_len == 0 {
                    self.reader.consume(consumed); // do not forget to consume the buffer for blank lines
                    continue;
                }

                // if it's a header line that requires lookahead,
                //we need to copy it to lookahead_line to avoid borrowing issues
                if needs_copy {
                    // we have to read the buffer again because we haven't consumed it yet,
                    // and we need to copy the header line to lookahead_line for lookahead
                    let buf = self.reader.fill_buf().map_err(FastxErr::IOError)?;
                    self.lookahead_line.clear();
                    self.lookahead_line.extend_from_slice(&buf[..consumed]);
                    self.has_lookahead = true;
                    self.reader.consume(consumed);
                    return Ok(ReadLineOutcome::NextHeader);
                }

                // if it's a FASTQ separator line, we can directly return the outcome without copying
                if matches!(first_char, Some(b'+') if stop_on_fastq_sep) {
                    self.reader.consume(consumed);
                    return Ok(ReadLineOutcome::FastqSep);
                }

                // normal line, directly append to record_buf without copying to line_buf
                let buf = self.reader.fill_buf().map_err(FastxErr::IOError)?;
                self.record_buf
                    .extend_from_slice(trim_crlf(&buf[..consumed]));
                self.reader.consume(consumed);
                return Ok(ReadLineOutcome::Appended(line_len));
            }

            // no '\n' is found in the buffer

            // if the first character indicates a potential LONG header or FASTQ separator,
            // we need to read the line using the slow path to properly handle
            // lookahead and avoid borrowing issues.
            // Otherwise, we can directly read the line into record_buf using the
            // slow path without worrying about lookahead.
            let first_char = {
                let buf = self.reader.fill_buf().map_err(FastxErr::IOError)?;
                if buf.is_empty() {
                    return Ok(ReadLineOutcome::Eof);
                }
                buf[0]
            };
            if !matches!(first_char, b'>' if stop_on_fasta_header)
                && !matches!(first_char, b'+' if stop_on_fastq_sep)
            {
                return self.read_long_line_into_record_buf();
            }

            // otherwise, we need to read the line using the slow path to
            // properly handle lookahead and avoid borrowing issues.
            match self.read_line_fill_buf() {
                Ok(0) => return Ok(ReadLineOutcome::Eof),
                Ok(_) => {
                    let line = trim_crlf(&self.line_buf);
                    if line.is_empty() {
                        continue;
                    }

                    match line[0] {
                        b'>' if stop_on_fasta_header => {
                            std::mem::swap(&mut self.line_buf, &mut self.lookahead_line);
                            self.has_lookahead = true;
                            return Ok(ReadLineOutcome::NextHeader);
                        }
                        b'+' if stop_on_fastq_sep => return Ok(ReadLineOutcome::FastqSep),
                        _ => {
                            self.record_buf.extend_from_slice(line);
                            return Ok(ReadLineOutcome::Appended(line.len()));
                        }
                    }
                }
                Err(e) => return Err(FastxErr::IOError(e)),
            }
        }
    }

    #[inline(always)]
    fn read_long_line_into_record_buf(&mut self) -> Result<ReadLineOutcome, FastxErr> {
        let mut line_len = 0;
        let mut pending_cr = false;

        loop {
            let step = {
                let buf = self.reader.fill_buf().map_err(FastxErr::IOError)?;
                if buf.is_empty() {
                    if line_len == 0 && !pending_cr {
                        return Ok(ReadLineOutcome::Eof);
                    }
                    None
                } else if let Some(pos) = memchr(b'\n', buf) {
                    let line = &buf[..pos];
                    let trimmed_len = if line.ends_with(b"\r") {
                        line.len() - 1
                    } else {
                        line.len()
                    };
                    Some((pos + 1, trimmed_len, true, pos == 0, line.ends_with(b"\r")))
                } else {
                    let trimmed_len = if buf.ends_with(b"\r") {
                        buf.len() - 1
                    } else {
                        buf.len()
                    };
                    Some((buf.len(), trimmed_len, false, false, buf.ends_with(b"\r")))
                }
            };

            let Some((consumed, trimmed_len, has_lf, lf_only, ends_with_cr)) = step else {
                return Ok(ReadLineOutcome::Appended(line_len));
            };

            let buf = self.reader.fill_buf().map_err(FastxErr::IOError)?;
            let data = if has_lf {
                &buf[..consumed - 1]
            } else {
                &buf[..consumed]
            };
            let data = &data[..trimmed_len];

            if pending_cr && !(has_lf && lf_only) {
                self.record_buf.push(b'\r');
                line_len += 1;
            }

            if !data.is_empty() {
                self.record_buf.extend_from_slice(data);
                line_len += data.len();
            }

            pending_cr = ends_with_cr && !has_lf;
            self.reader.consume(consumed);

            if has_lf {
                if line_len == 0 {
                    continue;
                }
                return Ok(ReadLineOutcome::Appended(line_len));
            }
        }
    }

    // returns None if EOF is reached, otherwise returns Some(Ok(Seq)) or Some(Err(e))
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Option<Result<Seq<'_>, FastxErr>> {
        self.record_buf.clear();

        // --- Step 1: load or read Header into self.line_buf ---

        if !self.has_lookahead {
            // first record
            match self.read_next_nonempty_line() {
                Ok(false) => return None,
                Ok(true) => {}
                Err(e) => return Some(Err(e)),
            }
            self.is_fastq = match trim_crlf(&self.line_buf)[0] {
                b'>' => false,
                b'@' => true,
                _ => return Some(Err(FastxErr::InvalidFormat)), // not a valid fasta/q record
            };
        } else {
            // not the first record
            std::mem::swap(&mut self.line_buf, &mut self.lookahead_line);
            self.has_lookahead = false;
        }

        // extract header from the header line and store it into record_buf
        let header_line = trim_crlf(&self.line_buf);
        let header: &[u8] = &header_line[1..];
        self.record_buf.extend_from_slice(header); // copy from header_line (excluding the leading '>' or '@') into record_buf
        let header_end = self.record_buf.len();

        // --- Step 2: read Sequence ---

        loop {
            match self.read_next_nonempty_line_into_record_buf(true, self.is_fastq) {
                Ok(
                    ReadLineOutcome::Eof | ReadLineOutcome::NextHeader | ReadLineOutcome::FastqSep,
                ) => break,
                Ok(ReadLineOutcome::Appended(_)) => {}
                Err(e) => return Some(Err(e)),
            }
        }

        let seq_end = self.record_buf.len();

        // --- Step 3: read Quality ---

        if self.is_fastq {
            let seq_len = seq_end - header_end;
            let mut qual_read_len = 0;

            while qual_read_len < seq_len {
                match self.read_next_nonempty_line_into_record_buf(false, false) {
                    Ok(ReadLineOutcome::Eof) => break, // Unexpected EOF in FASTQ
                    Ok(ReadLineOutcome::Appended(len)) => {
                        qual_read_len += len;

                        if qual_read_len > seq_len {
                            return Some(Err(FastxErr::UnequalSeqAndQual(seq_len, qual_read_len)));
                        }
                    }
                    Ok(ReadLineOutcome::NextHeader | ReadLineOutcome::FastqSep) => {
                        return Some(Err(FastxErr::UnequalSeqAndQual(seq_len, qual_read_len)));
                    }
                    Err(e) => return Some(Err(e)),
                }
            }

            if qual_read_len != seq_len {
                return Some(Err(FastxErr::UnequalSeqAndQual(seq_len, qual_read_len)));
            }
        }

        let buf_slice: &Vec<u8> = &self.record_buf;
        let id_slice: &[u8] = &buf_slice[0..header_end];
        let seq_slice: &[u8] = &buf_slice[header_end..seq_end];

        let qual_slice: Option<&[u8]> = if self.is_fastq {
            Some(&buf_slice[seq_end..])
        } else {
            None
        };

        if self.parse_id {
            let (id, desc) = parse_header(id_slice);
            return Some(Ok(Seq {
                id,
                desc,
                seq: seq_slice,
                qual: qual_slice,
            }));
        }
        Some(Ok(Seq {
            id: id_slice,
            desc: &[],
            seq: seq_slice,
            qual: qual_slice,
        }))
    }
}

enum ReadLineOutcome {
    Eof,
    Appended(usize),
    NextHeader,
    FastqSep,
}

#[inline]
fn parse_header(line: &[u8]) -> (&[u8], &[u8]) {
    let Some(id_end) = memchr2(b' ', b'\t', line) else {
        return (line, &[]);
    };

    let mut desc_start = id_end;
    while desc_start < line.len() && matches!(line[desc_start], b' ' | b'\t') {
        desc_start += 1;
    }

    (&line[..id_end], &line[desc_start..])
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::fs::{self, File};
    use std::io::{BufReader, Cursor};
    use std::io::Write;
    use std::path::PathBuf;
    use std::process::{Command, Stdio};
    use std::time::{SystemTime, UNIX_EPOCH};

    type OwnedRecord = (String, String, String, Option<String>);

    // --------------------------------------------------------------------------

    fn read_to_owned_from_reader<R: BufRead>(reader: R) -> Result<Vec<OwnedRecord>, FastxErr> {
        let mut reader: Reader<R> = Reader::from_reader(reader);
        let mut results: Vec<OwnedRecord> = Vec::new();

        while let Some(res) = reader.next() {
            let seq = res?;

            results.push((
                String::from_utf8_lossy(seq.id).to_string(),
                String::from_utf8_lossy(seq.desc).to_string(),
                String::from_utf8_lossy(seq.seq).to_string(),
                seq.qual.map(|q| String::from_utf8_lossy(q).to_string()),
            ));
        }
        Ok(results)
    }

    fn read_to_owned(input: &str) -> Result<Vec<OwnedRecord>, FastxErr> {
        read_to_owned_from_reader(Cursor::new(input.as_bytes()))
    }

    fn temp_path(suffix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("fastx-test-{}-{nanos}{suffix}", std::process::id()))
    }

    fn read_file_to_owned(path: &str) -> Result<Vec<OwnedRecord>, FastxErr> {
        let mut reader = Reader::new(path)?;
        let mut results = Vec::new();
        while let Some(res) = reader.next() {
            let seq = res?;
            results.push((
                String::from_utf8_lossy(seq.id).to_string(),
                String::from_utf8_lossy(seq.desc).to_string(),
                String::from_utf8_lossy(seq.seq).to_string(),
                seq.qual.map(|q| String::from_utf8_lossy(q).to_string()),
            ));
        }
        Ok(results)
    }

    // --------------------------------------------------------------------------
    // empty file or invalid format edge cases

    #[test]
    fn test_fasta_edge_cases_empty_file() {
        let input = "";
        let results = read_to_owned(input);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 0);
    }

    #[test]
    fn test_fasta_edge_cases_blank_file_lf() {
        let input = "\n";
        let results = read_to_owned(input);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 0);
    }

    #[test]
    fn test_fasta_edge_cases_blank_file_crlf() {
        let input = "\r\n";
        let results = read_to_owned(input);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 0);
    }

    #[test]
    fn test_fasta_edge_cases_blank_file_space() {
        let input = " ";
        let results = read_to_owned(input);

        assert!(matches!(results.unwrap_err(), FastxErr::InvalidFormat));
    }

    #[test]
    fn test_fasta_edge_cases_blank_file_space_and_lf() {
        let input = " \n";
        let results = read_to_owned(input);

        assert!(matches!(results.unwrap_err(), FastxErr::InvalidFormat));
    }

    #[test]
    fn test_fasta_edge_cases_blank_file_lf_and_space_and_lf() {
        let input = "\n \n";
        let results = read_to_owned(input);

        assert!(matches!(results.unwrap_err(), FastxErr::InvalidFormat));
    }

    // --------------------------------------------------------------------------
    // valid formats

    #[test]
    fn test_fasta_standard_single_and_multi_line() {
        let input = "\
>seq1 desc
ACTG
>seq2  desc2
AAAA
TTTT
GGGG
CCCC
>\tdesc3
AAAA
>    \t
AAAA
";
        let results = read_to_owned(input);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 4);

        // Record 1: single-line Seq
        assert_eq!(seqs[0].0, "seq1");
        assert_eq!(seqs[0].1, "desc");
        assert_eq!(seqs[0].2, "ACTG");
        assert_eq!(seqs[0].3, None);

        // Record 2: multi-line Seq
        assert_eq!(seqs[1].0, "seq2");
        assert_eq!(seqs[1].1, "desc2");
        assert_eq!(seqs[1].3, None);

        // Record 2: empty id and non-empty description
        assert_eq!(seqs[2].0, "");
        assert_eq!(seqs[2].1, "desc3");
        assert_eq!(seqs[2].2, "AAAA");

        // Record 2: empty id and empty description
        assert_eq!(seqs[3].0, "");
        assert_eq!(seqs[3].1, "");
        assert_eq!(seqs[3].2, "AAAA");
    }

    #[test]
    fn test_fastq_standard_single_and_multi_line() {
        let input = "\
@read1 desc
ACGT
+
IIII
@read2 \tdesc2
AC
GT
+
II
II
";
        let results = read_to_owned(input);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 2);

        // Record 1: 4-line FASTQ
        assert_eq!(seqs[0].0, "read1");
        assert_eq!(seqs[0].1, "desc");
        assert_eq!(seqs[0].2, "ACGT");
        assert_eq!(seqs[0].3, Some("IIII".to_string()));

        // Record 2: multi-line FASTQ
        assert_eq!(seqs[1].0, "read2");
        assert_eq!(seqs[1].1, "desc2");
        assert_eq!(seqs[1].2, "ACGT");
        assert_eq!(seqs[1].3, Some("IIII".to_string()));
    }

    // --------------------------------------------------------------------------
    // valid edge cases

    #[test]
    fn test_fasta_edge_cases_empty_seq() {
        let input = "\
>seq1
";
        let results = read_to_owned(input);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 1);

        assert_eq!(seqs[0].0, "seq1");
        assert_eq!(seqs[0].2, "");
    }

    #[test]
    fn test_fasta_edge_cases_beginning_lf() {
        let input = "
>seq1
actg
";
        let results = read_to_owned(input);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 1);

        assert_eq!(seqs[0].0, "seq1");
        assert_eq!(seqs[0].2, "actg");
    }

    #[test]
    fn test_fasta_edge_cases_empty_seq2() {
        let input = "\
>seq1

";
        let results = read_to_owned(input);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 1);

        assert_eq!(seqs[0].0, "seq1");
        assert_eq!(seqs[0].2, "");
    }

    #[test]
    fn test_fasta_edge_cases_empty_seq_and_no_lf() {
        let input = "\
>seq1";
        let results = read_to_owned(input);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 1);

        assert_eq!(seqs[0].0, "seq1");
        assert_eq!(seqs[0].2, "");
    }

    #[test]
    fn test_fasta_edge_cases_empty_fields() {
        // 1. ID is empty
        // 2. Seq is empty
        // 3. Both ID and Seq are empty
        // 4. No \n in the end
        let input = "\
>
ACGT
>seq_empty
>
>last";
        let results = read_to_owned(input);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 4);

        // Case 1: Empty ID, Normal Seq
        assert_eq!(seqs[0].0, "");
        assert_eq!(seqs[0].2, "ACGT");

        // Case 2: Normal ID, Empty Seq
        assert_eq!(seqs[1].0, "seq_empty");
        assert_eq!(seqs[1].2, "");

        // Case 3: Empty ID, Empty Seq
        assert_eq!(seqs[2].0, "");
        assert_eq!(seqs[2].2, "");

        // Case 4: Last check
        assert_eq!(seqs[3].0, "last");
        assert_eq!(seqs[3].2, "");
    }

    // --------------------------------------------------------------------------
    // invalid formats

    #[test]
    fn test_fastq_unequal_seq_qual_length() {
        let input = "\
@long_seq
AAAA
BBBB
+
JJJJ
";
        let results = read_to_owned(input);

        assert!(matches!(
            results.unwrap_err(),
            FastxErr::UnequalSeqAndQual(_, _),
        ));
    }

    #[test]
    fn test_fastq_unequal_seq_qual_length2() {
        let input = "\
@long_seq
AAAA
BBBB
+
JJJJ
KKKKK
";
        let results = read_to_owned(input);

        assert!(matches!(
            results.unwrap_err(),
            FastxErr::UnequalSeqAndQual(_, _),
        ));
    }

    // --------------------------------------------------------------------------
    // windows line endings and blank lines

    #[test]
    fn test_fasta_crlf_windows_line_endings() {
        let input = ">win\r\nACGT\r\n>win2\r\nTGCA\r\n";
        let results = read_to_owned(input);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].0, "win");
        assert_eq!(seqs[0].2, "ACGT");
        assert_eq!(seqs[1].2, "TGCA");
    }

    #[test]
    fn test_fasta_crlf_blank_lines_are_skipped() {
        let input = "\r\n>win\r\nACGT\r\n\r\n>win2\r\nTGCA\r\n";
        let results = read_to_owned(input);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].0, "win");
        assert_eq!(seqs[0].2, "ACGT");
        assert_eq!(seqs[1].0, "win2");
        assert_eq!(seqs[1].2, "TGCA");
    }

    // --------------------------------------------------------------------------
    // small buffer edge cases

    #[test]
    fn test_fastq_small_buffer_and_no_final_lf() {
        let input = "\
@read1 desc
ACGTACGT
+
IIIIIIII";
        let reader = BufReader::with_capacity(3, Cursor::new(input.as_bytes()));
        let results = read_to_owned_from_reader(reader);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 1);
        assert_eq!(seqs[0].0, "read1");
        assert_eq!(seqs[0].1, "desc");
        assert_eq!(seqs[0].2, "ACGTACGT");
        assert_eq!(seqs[0].3, Some("IIIIIIII".to_string()));
    }

    #[test]
    fn test_fasta_small_buffer_with_lookahead() {
        let input = "\
>seq1
ACGT
>seq2
TGCA";
        let reader = BufReader::with_capacity(3, Cursor::new(input.as_bytes()));
        let results = read_to_owned_from_reader(reader);

        assert!(results.is_ok());
        let seqs = results.unwrap();

        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].0, "seq1");
        assert_eq!(seqs[0].2, "ACGT");
        assert_eq!(seqs[1].0, "seq2");
        assert_eq!(seqs[1].2, "TGCA");
    }

    #[test]
    fn test_fasta_long_single_line_keeps_line_buf_small() {
        let seq = "ACGT".repeat(4096);
        let input = format!(">seq1\n{seq}\n>seq2\nTGCA\n");
        let reader = BufReader::with_capacity(3, Cursor::new(input.as_bytes()));
        let mut reader: Reader<_> = Reader::from_reader(reader);

        {
            let record = reader.next().unwrap().unwrap();
            assert_eq!(record.id, b"seq1");
            assert_eq!(record.seq, seq.as_bytes());
        }

        assert!(reader.line_buf.capacity() < seq.len() / 2);

        let record = reader.next().unwrap().unwrap();
        assert_eq!(record.id, b"seq2");
        assert_eq!(record.seq, b"TGCA");
    }

    #[test]
    fn test_reader_new_reads_last_fasta_record_without_final_lf() {
        let path = temp_path(".fa");
        fs::write(&path, b">chr1\nACGT\n>chrM\nTGCA").unwrap();

        let results = read_file_to_owned(path.to_str().unwrap());
        fs::remove_file(&path).unwrap();

        assert!(results.is_ok());
        let seqs = results.unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].0, "chr1");
        assert_eq!(seqs[0].2, "ACGT");
        assert_eq!(seqs[1].0, "chrM");
        assert_eq!(seqs[1].2, "TGCA");
    }

    #[test]
    fn test_reader_new_reads_last_gzip_fasta_record_without_final_lf() {
        let path = temp_path(".fa.gz");
        let file = File::create(&path).unwrap();
        let mut writer = GzEncoder::new(file, Compression::default());
        writer.write_all(b">chr1\nACGT\n>chrM\nTGCA").unwrap();
        writer.finish().unwrap();

        let results = read_file_to_owned(path.to_str().unwrap());
        fs::remove_file(&path).unwrap();

        assert!(results.is_ok());
        let seqs = results.unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].0, "chr1");
        assert_eq!(seqs[0].2, "ACGT");
        assert_eq!(seqs[1].0, "chrM");
        assert_eq!(seqs[1].2, "TGCA");
    }

    #[test]
    fn test_reader_new_reads_from_stdin() {
        const HELPER_ENV: &str = "FASTX_TEST_STDIN_HELPER";

        if std::env::var_os(HELPER_ENV).is_some() {
            let mut reader = Reader::new("-").unwrap();
            let mut count = 0usize;
            let mut last_id = String::new();
            let mut last_seq = String::new();

            while let Some(record) = reader.next() {
                let record = record.unwrap();
                count += 1;
                last_id = String::from_utf8_lossy(record.id).to_string();
                last_seq = String::from_utf8_lossy(record.seq).to_string();
            }

            println!("{count}\t{last_id}\t{last_seq}");
            return;
        }

        let mut child = Command::new(std::env::current_exe().unwrap())
            .env(HELPER_ENV, "1")
            .arg("reader::tests::test_reader_new_reads_from_stdin")
            .arg("--exact")
            .arg("--nocapture")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();

        child
            .stdin
            .take()
            .unwrap()
            .write_all(b">chr1\nACGT\n>chrM\nTGCA")
            .unwrap();

        let output = child.wait_with_output().unwrap();
        assert!(output.status.success());

        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(
            stdout.contains("2\tchrM\tTGCA"),
            "unexpected helper output: {stdout}"
        );
    }

    #[test]
    fn test_reader_new_reads_gzip_from_stdin() {
        const HELPER_ENV: &str = "FASTX_TEST_GZIP_STDIN_HELPER";

        if std::env::var_os(HELPER_ENV).is_some() {
            let mut reader = Reader::new("-").unwrap();
            let mut count = 0usize;
            let mut last_id = String::new();
            let mut last_seq = String::new();

            while let Some(record) = reader.next() {
                let record = record.unwrap();
                count += 1;
                last_id = String::from_utf8_lossy(record.id).to_string();
                last_seq = String::from_utf8_lossy(record.seq).to_string();
            }

            println!("{count}\t{last_id}\t{last_seq}");
            return;
        }

        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b">chr1\nACGT\n>chrM\nTGCA").unwrap();
        let payload = encoder.finish().unwrap();

        let mut child = Command::new(std::env::current_exe().unwrap())
            .env(HELPER_ENV, "1")
            .arg("reader::tests::test_reader_new_reads_gzip_from_stdin")
            .arg("--exact")
            .arg("--nocapture")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();

        child.stdin.take().unwrap().write_all(&payload).unwrap();

        let output = child.wait_with_output().unwrap();
        assert!(output.status.success());

        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(
            stdout.contains("2\tchrM\tTGCA"),
            "unexpected helper output: {stdout}"
        );
    }
}
