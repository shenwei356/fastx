use crate::errors::*;
use crate::seq::*;
use crate::util::*;
use memchr::{memchr, memchr2};
use std::io::BufRead;

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
    pub fn new(file: &str) -> Result<Self, FastxErr> {
        let r: Box<dyn BufRead> = xopen(file, 65536).map_err(FastxErr::IOError)?;
        Ok(Self::from_reader(r))
    }
}

impl<R: BufRead> Reader<R> {
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

    #[inline]
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
                        self.line_buf.extend_from_slice(&buf[..end]);
                        (end, true)
                    }
                    None => {
                        // no line ending found, consume the entire buffer and continue reading
                        self.line_buf.extend_from_slice(buf);
                        (buf.len(), false)
                    }
                }
            };

            self.reader.consume(consumed);
            total += consumed;

            if done {
                return Ok(total);
            }
            // else continue reading until we find a line ending or reach EOF
        }
    }

    #[inline]
    fn read_next_nonempty_line(&mut self) -> Result<bool, FastxErr> {
        loop {
            match self.read_line_fill_buf() {
                Ok(0) => return Ok(false),
                Ok(_) if trim_crlf(&self.line_buf).is_empty() => continue,
                Ok(_) => return Ok(true),
                Err(e) => return Err(FastxErr::IOError(e)),
            }
        }
    }

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
        self.record_buf.extend_from_slice(header);
        let header_end = self.record_buf.len();

        // --- Step 2: read Sequence ---

        loop {
            match self.read_next_nonempty_line() {
                Ok(false) => break,
                Ok(true) => {
                    let line = trim_crlf(&self.line_buf);
                    let first = line[0];

                    // next FASTA record
                    if first == b'>' {
                        std::mem::swap(&mut self.line_buf, &mut self.lookahead_line);
                        self.has_lookahead = true;
                        break;
                    }

                    // next FASTQ record
                    if self.is_fastq && first == b'+' {
                        break;
                    }

                    self.record_buf.extend_from_slice(line);
                }
                Err(e) => return Some(Err(e)),
            }
        }

        let seq_end = self.record_buf.len();

        // --- Step 3: read Quality ---

        if self.is_fastq {
            let seq_len = seq_end - header_end;
            let mut qual_read_len = 0;

            while qual_read_len < seq_len {
                match self.read_next_nonempty_line() {
                    Ok(false) => break, // Unexpected EOF in FASTQ
                    Ok(true) => {
                        let line = trim_crlf(&self.line_buf);
                        self.record_buf.extend_from_slice(line);
                        qual_read_len += line.len();

                        if qual_read_len > seq_len {
                            return Some(Err(FastxErr::UnequalSeqAndQual(seq_len, qual_read_len)));
                        }
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
    use std::io::{BufReader, Cursor};

    // --------------------------------------------------------------------------

    fn read_to_owned_from_reader<R: BufRead>(
        reader: R,
    ) -> Result<Vec<(String, String, String, Option<String>)>, FastxErr> {
        let mut reader: Reader<R> = Reader::from_reader(reader);
        let mut results: Vec<(String, String, String, Option<String>)> = Vec::new();

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

    fn read_to_owned(
        input: &str,
    ) -> Result<Vec<(String, String, String, Option<String>)>, FastxErr> {
        read_to_owned_from_reader(Cursor::new(input.as_bytes()))
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
}
