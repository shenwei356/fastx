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

    // skip parsing the ID and description fields in the header, and just return the raw header line as the ID
    pub fn skip_id_parsing(&mut self) {
        self.parse_id = false
    }

    // Read a line into line_buf, stripping any trailing "\r\n" or "\n".
    // Returns the number of raw bytes consumed from the reader (line + line ending).
    // 0 means EOF was reached without consuming anything.
    #[inline(always)]
    fn read_line_fill_buf(&mut self) -> std::io::Result<usize> {
        self.line_buf.clear();

        let mut total = 0;
        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                // EOF — strip any trailing \r that might have been buffered without a following \n
                if self.line_buf.last() == Some(&b'\r') {
                    self.line_buf.pop();
                }
                return Ok(total);
            }

            // Look for a newline in the buffer. If found, we can copy the line directly into line_buf and be done.
            let (consumed, done) = match memchr(b'\n', buf) {
                Some(pos) => {
                    // found a line ending; copy the line (including '\n') into line_buf
                    let end = pos + 1;
                    self.line_buf.extend_from_slice(&buf[..end]);
                    (end, true)
                }
                None => {
                    // no line ending found, consume the entire buffer and continue reading
                    self.line_buf.extend_from_slice(buf);
                    (buf.len(), false)
                }
            };

            self.reader.consume(consumed);
            total += consumed;

            //
            if done {
                // strip trailing '\n' and any preceding '\r' so callers see a clean line
                self.line_buf.pop();
                if self.line_buf.last() == Some(&b'\r') {
                    self.line_buf.pop();
                }
                return Ok(total);
            }
        }
    }

    // only for the first record
    #[inline(always)]
    fn read_next_nonempty_line(&mut self) -> Result<bool, FastxErr> {
        loop {
            match self.read_line_fill_buf() {
                Ok(0) => return Ok(false),                     // EOF
                Ok(_) if self.line_buf.is_empty() => continue, // skip blank lines
                Ok(_) => return Ok(true),                      // non-empty line read successfully
                Err(e) => return Err(FastxErr::IOError(e)),    // I/O error occurred
            }
        }
    }

    // Read the next non-empty line and append it to record_buf.
    //
    // STOP_ON_FASTA_HEADER and STOP_ON_FASTQ_SEP are const generics so each call site is
    // monomorphized: the compiler folds the boolean checks at compile time and removes
    // unreachable branches.
    #[inline(always)]
    fn read_next_nonempty_line_into_record_buf<
        const STOP_ON_FASTA_HEADER: bool,
        const STOP_ON_FASTQ_SEP: bool,
    >(
        &mut self,
    ) -> Result<ReadLineOutcome, FastxErr> {
        loop {
            // Snapshot the buffer once and try to handle the line in-place. This relies on
            // disjoint-field borrows: `self.reader` (held by `buf`) and `self.record_buf` /
            // `self.lookahead_line` are independent fields, so they can be touched while
            // the slice lives. After the last use of `buf`, NLL releases the borrow on
            // `self.reader` so `consume` can run.
            let buf = self.reader.fill_buf().map_err(FastxErr::IOError)?;
            if buf.is_empty() {
                return Ok(ReadLineOutcome::Eof);
            }

            // Look for a newline in the buffer. If found, we can handle the entire line without copying if it doesn't start with a header/sep char.
            if let Some(pos) = memchr(b'\n', buf) {
                // Fast path: the entire line is contained in the current buffer.
                let consumed = pos + 1;
                // `trim_crlf` only adjusts the tail, so `line_len` doubles as the trimmed end.
                let line_len = trim_crlf(&buf[..consumed]).len();

                if line_len == 0 {
                    // blank line — drop the buf borrow and consume
                    self.reader.consume(consumed);
                    continue;
                }

                let first_char = buf[0];

                if STOP_ON_FASTA_HEADER && first_char == b'>' {
                    // stash the header (already trimmed) into lookahead, last use of `buf`
                    self.lookahead_line.clear();
                    self.lookahead_line.extend_from_slice(&buf[..line_len]);
                    self.reader.consume(consumed);
                    self.has_lookahead = true;
                    return Ok(ReadLineOutcome::NextHeader);
                }

                if STOP_ON_FASTQ_SEP && first_char == b'+' {
                    // separator — no copy needed
                    self.reader.consume(consumed);
                    return Ok(ReadLineOutcome::FastqSep);
                }

                // normal line: append the trimmed content directly from the reader's buffer
                self.record_buf.extend_from_slice(&buf[..line_len]);
                self.reader.consume(consumed);
                return Ok(ReadLineOutcome::Appended(line_len));
            }

            // No '\n' in this buffer fill. Either the line is long, or this is the last
            // line of the input without a trailing newline.
            let first_char = buf[0];

            // If the first byte cannot be a header / separator, hand off to the long-line
            // appender — it streams directly into record_buf without lookahead concerns.
            if !((STOP_ON_FASTA_HEADER && first_char == b'>')
                || (STOP_ON_FASTQ_SEP && first_char == b'+'))
            {
                return self.read_long_line_into_record_buf();
            }

            // Slow path: assemble the full line into self.line_buf so we can preserve it as
            // a lookahead header if needed.
            match self.read_line_fill_buf() {
                Ok(0) => return Ok(ReadLineOutcome::Eof),
                Ok(_) => {
                    if self.line_buf.is_empty() {
                        continue;
                    }
                    let first_char = self.line_buf[0];
                    if STOP_ON_FASTA_HEADER && first_char == b'>' {
                        // line_buf is already trimmed; swap it into lookahead_line
                        std::mem::swap(&mut self.line_buf, &mut self.lookahead_line);
                        self.has_lookahead = true;
                        return Ok(ReadLineOutcome::NextHeader);
                    }
                    if STOP_ON_FASTQ_SEP && first_char == b'+' {
                        return Ok(ReadLineOutcome::FastqSep);
                    }
                    let len = self.line_buf.len();
                    self.record_buf.extend_from_slice(&self.line_buf);
                    return Ok(ReadLineOutcome::Appended(len));
                }
                Err(e) => return Err(FastxErr::IOError(e)),
            }
        }
    }

    // Specialized reader for FASTQ quality lines: never produces NextHeader / FastqSep,
    // so it skips all the header / separator decision logic.
    #[inline(always)]
    fn read_qual_line_into_record_buf(&mut self) -> Result<ReadLineOutcome, FastxErr> {
        loop {
            let buf = self.reader.fill_buf().map_err(FastxErr::IOError)?;
            if buf.is_empty() {
                return Ok(ReadLineOutcome::Eof);
            }

            if let Some(pos) = memchr(b'\n', buf) {
                let consumed = pos + 1;
                let line_len = trim_crlf(&buf[..consumed]).len();
                if line_len == 0 {
                    self.reader.consume(consumed);
                    continue;
                }
                self.record_buf.extend_from_slice(&buf[..line_len]);
                self.reader.consume(consumed);
                return Ok(ReadLineOutcome::Appended(line_len));
            }

            // long line that spans buffer boundaries — reuse the existing CR/LF state machine
            return self.read_long_line_into_record_buf();
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
        // line_buf and lookahead_line are kept already-trimmed (no trailing \r\n).

        if !self.has_lookahead {
            // first record
            match self.read_next_nonempty_line() {
                Ok(false) => return None,
                Ok(true) => {}
                Err(e) => return Some(Err(e)),
            }
            self.is_fastq = match self.line_buf[0] {
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
        // (line_buf is already trimmed; just skip the leading '>' or '@')
        let header: &[u8] = &self.line_buf[1..];
        self.record_buf.extend_from_slice(header);
        let header_end = self.record_buf.len();

        // --- Step 2: read Sequence ---
        // The const generics specialize each call site so the dead branches are folded away.

        if self.is_fastq {
            loop {
                match self.read_next_nonempty_line_into_record_buf::<true, true>() {
                    Ok(
                        ReadLineOutcome::Eof
                        | ReadLineOutcome::NextHeader
                        | ReadLineOutcome::FastqSep,
                    ) => break,
                    Ok(ReadLineOutcome::Appended(_)) => {}
                    Err(e) => return Some(Err(e)),
                }
            }
        } else {
            loop {
                match self.read_next_nonempty_line_into_record_buf::<true, false>() {
                    Ok(
                        ReadLineOutcome::Eof
                        | ReadLineOutcome::NextHeader
                        | ReadLineOutcome::FastqSep,
                    ) => break,
                    Ok(ReadLineOutcome::Appended(_)) => {}
                    Err(e) => return Some(Err(e)),
                }
            }
        }

        let seq_end = self.record_buf.len();

        // --- Step 3: read Quality ---
        // Use the dedicated quality reader: it never produces NextHeader / FastqSep, so the
        // header/separator decision logic is eliminated for the entire quality block.

        if self.is_fastq {
            let seq_len = seq_end - header_end;
            let mut qual_read_len = 0;

            while qual_read_len < seq_len {
                match self.read_qual_line_into_record_buf() {
                    Ok(ReadLineOutcome::Eof) => break, // unexpected EOF mid-FASTQ
                    Ok(ReadLineOutcome::Appended(len)) => {
                        qual_read_len += len;

                        if qual_read_len > seq_len {
                            return Some(Err(FastxErr::UnequalSeqAndQual(seq_len, qual_read_len)));
                        }
                    }
                    // read_qual_line_into_record_buf never produces these variants
                    Ok(ReadLineOutcome::NextHeader | ReadLineOutcome::FastqSep) => unreachable!(),
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
    use std::io::Write;
    use std::io::{BufReader, Cursor};
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
        std::env::temp_dir().join(format!("fastseq-test-{}-{nanos}{suffix}", std::process::id()))
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

    #[test]
    fn test_fastq_standard_single_and_multi_lin_with_comments_in_3rd_line() {
        let input = "\
@read1 desc
ACGT
+read1 desc
IIII
@read2 \tdesc2
AC
GT
+ read2 \tdesc2
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

#[cfg(test)]
mod bgzf_tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::fs;
    use std::io::Write;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(suffix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("fastseq-bgzf-{}-{nanos}{suffix}", std::process::id()))
    }

    #[test]
    fn test_multi_member_gzip_fastq() {
        // Simulate a BGZF-like file: multiple concatenated gzip members
        let path = temp_path(".fq.gz");
        let mut file = fs::File::create(&path).unwrap();

        let total_records = 500;
        let records_per_block = 100;
        let lengths = [50, 80, 30, 150, 10, 120, 5, 200, 75, 45];

        for block in 0..(total_records / records_per_block) {
            let mut block_data = Vec::new();
            for i in 0..records_per_block {
                let idx = block * records_per_block + i;
                let len = lengths[idx % lengths.len()];
                let seq: String = std::iter::repeat('A').take(len).collect();
                let qual: String = std::iter::repeat('I').take(len).collect();
                write!(block_data, "@read{}\n{}\n+\n{}\n", idx, seq, qual).unwrap();
            }
            // Each block is a separate gzip member
            let mut gz = GzEncoder::new(Vec::new(), Compression::default());
            gz.write_all(&block_data).unwrap();
            let compressed = gz.finish().unwrap();
            file.write_all(&compressed).unwrap();
        }
        drop(file);

        // Now try to read with our Reader
        let mut reader = Reader::new(path.to_str().unwrap()).unwrap();
        let mut count = 0;
        while let Some(res) = reader.next() {
            let seq = res.unwrap();
            let expected_len = lengths[count % lengths.len()];
            assert_eq!(seq.seq.len(), expected_len, "record {} seq len", count);
            count += 1;
        }
        fs::remove_file(&path).unwrap();
        assert_eq!(count, total_records, "expected {} records, got {}", total_records, count);
    }
}
