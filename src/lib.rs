use std::io::BufRead;

pub mod util;
use util::*;

mod seq;
use seq::*;

mod errors;
use errors::*;

pub struct Reader<R: BufRead> {
    reader: R,

    record_buf: Vec<u8>,
    line_buf: Vec<u8>,
    lookahead_line: Option<Vec<u8>>,
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
            reader: reader,
            record_buf: Vec::with_capacity(1 << 20),
            line_buf: Vec::with_capacity(128),
            lookahead_line: None,
        }
    }

    pub fn next(&mut self) -> Option<Result<Seq<'_>, FastxErr>> {
        self.record_buf.clear();

        // --- Step 1: load or read Header into self.line_buf ---

        if self.lookahead_line.is_some() {
            std::mem::swap(&mut self.line_buf, self.lookahead_line.as_mut().unwrap());
            self.lookahead_line = None;
        } else {
            loop {
                self.line_buf.clear();
                match self.reader.read_until(b'\n', &mut self.line_buf) {
                    Ok(0) => return None,                           // EOF
                    Ok(_) if self.line_buf[0] == b'\n' => continue, // empty line with only '\n'
                    Ok(_) => break,                                 // a non-empty line
                    Err(e) => return Some(Err(FastxErr::IOError(e))),
                }
            }
        }

        let is_fastq = match self.line_buf[0] {
            b'>' => false,
            b'@' => true,
            _ => return Some(Err(FastxErr::InvalidFormat)), // not a valid fasta/q record
        };

        // extract header from the header line and store it into record_buf
        let header: &[u8] = &trim_crlf(&self.line_buf)[1..];
        self.record_buf.extend_from_slice(header);
        let header_end = self.record_buf.len();

        // --- Step 2: read Sequence ---

        loop {
            self.line_buf.clear();
            match self.reader.read_until(b'\n', &mut self.line_buf) {
                Ok(0) => break,                                 // EOF,
                Ok(_) if self.line_buf[0] == b'\n' => continue, // empty line with only '\n'
                Ok(_) => {
                    let first = self.line_buf[0];

                    // next FASTA record
                    if first == b'>' {
                        // self.lookahead_line = Some(self.line_buf.clone());
                        self.lookahead_line = Some(std::mem::take(&mut self.line_buf));
                        break;
                    }

                    // next FASTQ record
                    if is_fastq && first == b'+' {
                        break;
                    }

                    self.record_buf.extend_from_slice(trim_crlf(&self.line_buf));
                }
                Err(e) => return Some(Err(FastxErr::IOError(e))),
            }
        }

        let seq_end = self.record_buf.len();

        // --- Step 3: read Quality ---

        if is_fastq {
            let seq_len = seq_end - header_end;
            let mut qual_read_len = 0;

            while qual_read_len < seq_len {
                self.line_buf.clear();
                match self.reader.read_until(b'\n', &mut self.line_buf) {
                    Ok(0) => break,                                 // Unexpected EOF in FASTQ
                    Ok(_) if self.line_buf[0] == b'\n' => continue, // empty line with only '\n'
                    Ok(_) => {
                        let clean_qual: &[u8] = trim_crlf(&self.line_buf);
                        self.record_buf.extend_from_slice(clean_qual);
                        qual_read_len += clean_qual.len();
                    }
                    Err(e) => return Some(Err(FastxErr::IOError(e))),
                }
            }

            if self.record_buf.len() - seq_end > seq_len {
                // self.record_buf.truncate(seq_end + seq_len);
                return Some(Err(FastxErr::UnequalSeqAndQual(
                    seq_len,
                    self.record_buf.len() - seq_end,
                )));
            }
        }

        let buf_slice: &Vec<u8> = &self.record_buf;
        let id_slice: &[u8] = &buf_slice[0..header_end];
        let seq_slice: &[u8] = &buf_slice[header_end..seq_end];

        let qual_slice: Option<&[u8]> = if is_fastq {
            let q_len: usize = buf_slice.len() - seq_end;
            if q_len != seq_slice.len() {
                return Some(Err(FastxErr::UnequalSeqAndQual(seq_slice.len(), q_len)));
            }
            Some(&buf_slice[seq_end..])
        } else {
            None
        };

        let (id, desc) = parse_header(id_slice);
        Some(Ok(Seq {
            id: id,
            desc: desc,
            seq: seq_slice,
            qual: qual_slice,
        }))
    }
}

fn parse_header(line: &[u8]) -> (&[u8], &[u8]) {
    // match line.iter().position(|&b| b == b' ' || b == b'\t') {
    //     Some(id_end) => {
    //         let id_slice = &line[0..id_end]; // id_end might be 0
    //         let remainder = &line[id_end..];
    //         let desc_start_offset = remainder
    //             .iter()
    //             .position(|&b| b != b' ' && b != b'\t')
    //             .unwrap_or(remainder.len());
    //         (id_slice, &line[id_end + desc_start_offset..])
    //     }
    //     None => (line, &[]), // no blank or tab
    // }

    let mut i = 0;
    let n = line.len();

    while i < n && line[i] != b' ' && line[i] != b'\t' {
        i += 1;
    }
    let id = &line[..i];

    while i < n && (line[i] == b' ' || line[i] == b'\t') {
        i += 1;
    }

    (id, &line[i..])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn read_to_owned(
        input: &str,
    ) -> Result<Vec<(String, String, String, Option<String>)>, FastxErr> {
        let cursor: Cursor<&[u8]> = Cursor::new(input.as_bytes());
        let mut reader: Reader<Cursor<&[u8]>> = Reader::from_reader(cursor);
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
    fn test_fasta_edge_cases_blank_file_space() {
        let input = " ";
        let results = read_to_owned(input);

        assert!(matches!(results.unwrap_err(), FastxErr::InvalidFormat));
    }

    #[test]
    fn test_fasta_edge_cases_blank_file_space_and_lf() {
        let input = "\n \n";
        let results = read_to_owned(input);

        assert!(matches!(results.unwrap_err(), FastxErr::InvalidFormat));
    }

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
}
