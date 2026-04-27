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
