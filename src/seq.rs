#[derive(Debug, Clone, Copy)]
pub struct Seq<'a> {
    pub id: &'a [u8],
    pub desc: &'a [u8],
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>,
}

impl<'a> Seq<'a> {
    pub fn is_fastq(&self) -> bool {
        self.qual.is_some()
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn rc(&self) -> Vec<u8> {
        self.seq
            .iter()
            .rev()
            .map(|&base| RC_TABLE[base as usize])
            .collect()
    }

    pub fn rc_unsafe(&self) -> Vec<u8> {
        let mut result = Vec::with_capacity(self.seq.len());
        let rc_ptr = RC_TABLE.as_ptr();
        for &base in self.seq.iter().rev() {
            result.push(unsafe { *rc_ptr.add(base as usize) });
        }
        result
    }

    pub fn count_base(&self, base: u8) -> usize {
        self.seq.iter().copied().filter(|&b| b == base).count()
    }

    pub fn count_base_fn(&self, f: fn(&u8) -> bool) -> usize {
        self.seq.iter().copied().filter(f).count()
    }

    pub fn count_bases(&self, bases: &[u8]) -> usize {
        let mut table = [0u8; 256];
        for b in bases {
            table[*b as usize] = 1;
        }

        self.seq.iter().map(|&b| table[b as usize] as usize).sum()
    }

    pub fn gc_content(&self) -> f32 {
        if self.seq.len() == 0 {
            return 0.0;
        }
        self.count_base_fn(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) as f32 / self.len() as f32
    }
}

const RC_TABLE: [u8; 256] = make_rc_table();

const fn make_rc_table() -> [u8; 256] {
    let mut table = [b'N'; 256];

    // ACGT
    table[b'A' as usize] = b'T';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'T' as usize] = b'A';
    table[b'N' as usize] = b'N';

    table[b'a' as usize] = b't';
    table[b'c' as usize] = b'g';
    table[b'g' as usize] = b'c';
    table[b't' as usize] = b'a';
    table[b'n' as usize] = b'n';

    // IUPAC
    table[b'M' as usize] = b'K'; // M: A/C
    table[b'R' as usize] = b'Y'; // R: A/G
    table[b'W' as usize] = b'W'; // W: A/T
    table[b'S' as usize] = b'S'; // S: C/G
    table[b'Y' as usize] = b'R'; // Y: C/T
    table[b'K' as usize] = b'M'; // K: G/T

    table[b'V' as usize] = b'B'; // V: A/C/G
    table[b'H' as usize] = b'D'; // H: A/C/T
    table[b'D' as usize] = b'H'; // D: A/G/T
    table[b'B' as usize] = b'V'; // B: C/G/T

    table[b'm' as usize] = b'k';
    table[b'r' as usize] = b'y';
    table[b'w' as usize] = b'w';
    table[b's' as usize] = b's';
    table[b'y' as usize] = b'r';
    table[b'k' as usize] = b'm';

    table[b'v' as usize] = b'v';
    table[b'h' as usize] = b'h';
    table[b'd' as usize] = b'd';
    table[b'b' as usize] = b'b';

    // gaps
    table[b'.' as usize] = b'.';
    table[b'-' as usize] = b'-';
    table[b' ' as usize] = b' ';

    table
}

#[cfg(test)]
mod tests {
    use super::*;

    fn a_seq(seq: &'_ [u8]) -> Seq<'_> {
        Seq {
            id: b"",
            desc: b"",
            seq: seq,
            qual: None,
        }
    }

    #[test]
    fn test_rc() {
        // even bases
        let seq = b"ACGTNacgtn";
        let expected = b"nacgtNACGT";
        assert_eq!(a_seq(seq).rc(), expected);

        // odd bases
        let seq_odd = b"ACGTGa";
        let expected_odd = b"tCACGT";
        assert_eq!(a_seq(seq_odd).rc(), expected_odd);

        // iupac
        let seq_iupac = b"MRSWKY";
        let expected_iupac = b"RMWSYK";
        assert_eq!(a_seq(seq_iupac).rc(), expected_iupac);

        // empty seq
        let seq_empty = b"";
        let expected_empty = b"";
        assert_eq!(a_seq(seq_empty).rc(), expected_empty);

        // gaps
        let seq_gaps = b"a-c";
        let expected_gaps = b"g-t";
        assert_eq!(a_seq(seq_gaps).rc(), expected_gaps);

        // unknown base
        let seq_unknown = b"N?";
        let expected_unknown = b"NN";
        assert_eq!(a_seq(seq_unknown).rc(), expected_unknown);
    }

    #[test]
    fn test_base_counting() {
        let seq = b"ACGTNacgtn";
        assert_eq!(a_seq(seq).count_base(b'a'), 1);

        let seq = b"ACGTNacgtn";
        assert_eq!(a_seq(seq).count_bases(b"Aa"), 2);

        let seq = b"ACGTNacgtn";
        assert_eq!(a_seq(seq).count_base_fn(|&b| b == b'a' || b == b'A'), 2);

        let seq = b"ACGTNacgtn";
        assert_eq!(a_seq(seq).gc_content(), 0.4);

        let seq = b"";
        assert_eq!(a_seq(seq).gc_content(), 0.0);
    }
}
