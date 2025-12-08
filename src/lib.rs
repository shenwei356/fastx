use std::io::{self, BufRead};

mod util;
use util::*;

#[derive(Debug, Clone, Copy)]
pub struct Sequence<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>,
}

impl<'a> Sequence<'a> {
    pub fn is_fastq(&self) -> bool {
        self.qual.is_some()
    }
}

pub struct Reader<R: BufRead> {
    reader: R,

    record_buf: Vec<u8>,
    line_buf: Vec<u8>,
    lookahead_line: Option<Vec<u8>>,
}

impl Reader<Box<dyn BufRead>> {
    pub fn new(file: &str) -> io::Result<Self> {
        let r: Box<dyn BufRead> = xopen(file)?;
        Ok(Self::from_reader(r))
    }
}

impl<R: BufRead> Reader<R> {
    pub fn from_reader(reader: R) -> Self {
        Self {
            reader: reader,
            record_buf: Vec::new(),
            line_buf: Vec::new(),
            lookahead_line: None,
        }
    }

    pub fn next(&mut self) -> Option<Sequence<'_>> {
        self.record_buf.clear();

        // --- Step 1: load or read Header ---

        if self.lookahead_line.is_some() {
            std::mem::swap(&mut self.line_buf, self.lookahead_line.as_mut().unwrap());
            self.lookahead_line = None;
        } else {
            self.line_buf.clear();
            match self.reader.read_until(b'\n', &mut self.line_buf) {
                Ok(0) => return None, // EOF
                Ok(_) => {}           // continue
                Err(_) => return None,
            }
        }

        let line_len: usize = self.line_buf.len();
        if line_len == 0 {
            // shoud not reach here
            return None;
        }

        let is_fastq = match self.line_buf[0] {
            b'>' => false,
            b'@' => true,
            _ => return None, // not a valid fasta/q record
        };

        let header_content: &[u8] = trim_crlf(&self.line_buf);
        let id_content: &[u8] = if header_content.len() > 0 {
            &header_content[1..]
        } else {
            &[]
        };

        self.record_buf.extend_from_slice(id_content);
        let id_end = self.record_buf.len();

        // --- Step 2: read Sequence ---

        loop {
            self.line_buf.clear();
            match self.reader.read_until(b'\n', &mut self.line_buf) {
                Ok(0) => break, // EOF
                Ok(_) => {
                    let first = self.line_buf[0];

                    // next FASTA record
                    if first == b'>' {
                        self.lookahead_line = Some(self.line_buf.clone());
                        break;
                    }

                    // next FASTQ record
                    if is_fastq && first == b'+' {
                        break;
                    }

                    let valid_len = self.line_buf.len();
                    if valid_len > 0 {
                        let clean_seq = trim_crlf(&self.line_buf);
                        self.record_buf.extend_from_slice(clean_seq);
                    }
                }
                Err(_) => break,
            }
        }

        let seq_end = self.record_buf.len();

        // --- Step 3: read Quality ---

        if is_fastq {
            let seq_len = seq_end - id_end;
            let mut qual_read_len = 0;

            while qual_read_len < seq_len {
                self.line_buf.clear();
                match self.reader.read_until(b'\n', &mut self.line_buf) {
                    Ok(0) => break, // Unexpected EOF in FASTQ
                    Ok(_) => {
                        let clean_qual: &[u8] = trim_crlf(&self.line_buf);
                        self.record_buf.extend_from_slice(clean_qual);
                        qual_read_len += clean_qual.len();
                    }
                    Err(_) => break,
                }
            }

            if self.record_buf.len() - seq_end > seq_len {
                self.record_buf.truncate(seq_end + seq_len);
            }
        }

        let buf_slice: &Vec<u8> = &self.record_buf;
        let id_slice: &[u8] = &buf_slice[0..id_end];
        let seq_slice: &[u8] = &buf_slice[id_end..seq_end];

        let qual_slice: Option<&[u8]> = if is_fastq {
            let q_len: usize = buf_slice.len() - seq_end;
            if q_len != seq_slice.len() {
                eprintln!("Warning: Seq len {} != Qual len {}", seq_slice.len(), q_len);
            }
            Some(&buf_slice[seq_end..])
        } else {
            None
        };

        Some(Sequence {
            id: id_slice,
            seq: seq_slice,
            qual: qual_slice,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn read_to_owned(input: &str) -> Vec<(String, String, Option<String>)> {
        let cursor: Cursor<&[u8]> = Cursor::new(input.as_bytes());
        let mut reader: Reader<Cursor<&[u8]>> = Reader::from_reader(cursor);
        let mut results: Vec<(String, String, Option<String>)> = Vec::new();

        while let Some(seq) = reader.next() {
            results.push((
                String::from_utf8_lossy(seq.id).to_string(),
                String::from_utf8_lossy(seq.seq).to_string(),
                seq.qual.map(|q| String::from_utf8_lossy(q).to_string()),
            ));
        }
        results
    }

    #[test]
    fn test_fasta_empty_file() {
        let input = "";
        let results = read_to_owned(input);

        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_fasta_blank_file() {
        let input = " ";
        let results = read_to_owned(input);

        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_fasta_standard_single_and_multi_line() {
        let input = "\
>seq1
ACTG
>seq2 description
AAAA
TTTT
GGGG
CCCC
";
        let results = read_to_owned(input);

        assert_eq!(results.len(), 2);

        // Record 1: single-line Seq
        assert_eq!(results[0].0, "seq1");
        assert_eq!(results[0].1, "ACTG");
        assert_eq!(results[0].2, None);

        // Record 2: multi-line Seq
        assert_eq!(results[1].0, "seq2 description");
        assert_eq!(results[1].1, "AAAATTTTGGGGCCCC");
        assert_eq!(results[1].2, None);
    }

    #[test]
    fn test_fastq_standard_single_and_multi_line() {
        let input = "\
@read1
ACGT
+
IIII
@read2
AC
GT
+
II
II
";
        let results = read_to_owned(input);

        assert_eq!(results.len(), 2);

        // Record 1: 4-line FASTQ
        assert_eq!(results[0].0, "read1");
        assert_eq!(results[0].1, "ACGT");
        assert_eq!(results[0].2, Some("IIII".to_string()));

        // Record 2: multi-line FASTQ
        assert_eq!(results[1].0, "read2");
        assert_eq!(results[1].1, "ACGT");
        assert_eq!(results[1].2, Some("IIII".to_string()));
    }

    #[test]
    fn test_fasta_edge_cases_empty_seq_and_no_lf() {
        let input = "\
>seq1";
        let results = read_to_owned(input);
        assert_eq!(results.len(), 1);

        assert_eq!(results[0].0, "seq1");
        assert_eq!(results[0].1, "");
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
        assert_eq!(results.len(), 4);

        // Case 1: Empty ID, Normal Seq
        assert_eq!(results[0].0, "");
        assert_eq!(results[0].1, "ACGT");

        // Case 2: Normal ID, Empty Seq
        assert_eq!(results[1].0, "seq_empty");
        assert_eq!(results[1].1, "");

        // Case 3: Empty ID, Empty Seq
        assert_eq!(results[2].0, "");
        assert_eq!(results[2].1, "");

        // Case 4: Last check
        assert_eq!(results[3].0, "last");
        assert_eq!(results[3].1, "");
    }

    #[test]
    fn test_fastq_seq_qual_length_consistency() {
        let input = "\
@long_seq
AAAA
BBBB
+
JJJJ
JJJJ
";
        let results = read_to_owned(input);
        assert_eq!(results.len(), 1);

        assert_eq!(results[0].0, "long_seq");
        assert_eq!(results[0].1, "AAAABBBB"); // 8 chars
        assert_eq!(results[0].2, Some("JJJJJJJJ".to_string())); // 8 chars
    }

    #[test]
    fn test_fasta_crlf_windows_line_endings() {
        let input = ">win\r\nACGT\r\n>win2\r\nTGCA\r\n";
        let results = read_to_owned(input);

        assert_eq!(results.len(), 2);
        assert_eq!(results[0].0, "win");
        assert_eq!(results[0].1, "ACGT");
        assert_eq!(results[1].1, "TGCA");
    }
}
