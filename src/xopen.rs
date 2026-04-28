use bzip2::read::{BzDecoder, BzEncoder};
use flate2::Compression;
use flate2::read::{GzDecoder, GzEncoder};
use liblzma::read::{XzDecoder, XzEncoder};
use std::alloc::{Layout, alloc_zeroed, dealloc};
use std::fs::File;
use std::io::IsTerminal;
use std::io::{self, BufRead, Read, Write};
use std::ptr::NonNull;
use std::slice;
use zstd::stream::read::Decoder as ZstdDecoder;
use zstd::stream::write::Encoder as ZstdEncoder;

pub const DEFAULT_IO_BUFFER_ALIGNMENT: usize = 4096;

/// AlignedBytes is a helper struct that manages an aligned buffer for I/O operations.
struct AlignedBytes {
    ptr: NonNull<u8>, // a non-null pointer to the allocated buffer memory
    len: usize,       // the length of the buffer in bytes
    align: usize,     // the alignment of the buffer in bytes, which must be a non-zero power of two
}

impl AlignedBytes {
    /// Create a new AlignedBytes with the specified length and alignment.
    fn new(len: usize, align: usize) -> io::Result<Self> {
        let align = normalized_alignment(align)?;
        let len = len.max(1);
        // ensure the layout is valid for the given length and alignment
        let layout = Layout::from_size_align(len, align).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "buffer size and alignment form an invalid layout",
            )
        })?;

        // allocate zeroed memory for the buffer and ensure it's non-null
        let ptr = unsafe { alloc_zeroed(layout) };
        // check if allocation succeeded and create a NonNull pointer
        let ptr = NonNull::new(ptr)
            .ok_or_else(|| io::Error::other("failed to allocate aligned I/O buffer"))?;

        Ok(Self { ptr, len, align })
    }

    // safely create a slice from the raw pointer and length,
    // ensuring the buffer is properly aligned and valid for the specified length
    #[inline]
    fn as_slice(&self) -> &[u8] {
        unsafe { slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }

    // safely create a mutable slice from the raw pointer and length,
    // ensuring the buffer is properly aligned and valid for the specified length
    #[inline]
    fn as_mut_slice(&mut self) -> &mut [u8] {
        unsafe { slice::from_raw_parts_mut(self.ptr.as_ptr(), self.len) }
    }
}

impl Drop for AlignedBytes {
    // when the AlignedBytes is dropped, we need to deallocate the memory using the same layout
    fn drop(&mut self) {
        let layout = Layout::from_size_align(self.len, self.align)
            .expect("aligned buffer layout should stay valid");
        unsafe { dealloc(self.ptr.as_ptr(), layout) };
    }
}

// validate that the provided alignment is a non-zero power of two,
// which is required for proper memory alignment in I/O operations
fn normalized_alignment(align: usize) -> io::Result<usize> {
    if align == 0 || !align.is_power_of_two() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "buffer alignment must be a non-zero power of two",
        ));
    }
    Ok(align)
}

/// AlignedBufReader is a buffered reader that uses an aligned buffer for efficient I/O operations.
struct AlignedBufReader<R> {
    inner: R,          // the underlying reader that provides the data
    buf: AlignedBytes, // the internal aligned buffer used for reading data
    pos: usize,        // the current position in the buffer that has been consumed
    filled: usize,     // the amount of data currently filled in the buffer
}

impl<R> AlignedBufReader<R> {
    /// Create a new AlignedBufReader with the specified capacity, alignment, and inner reader.
    fn with_capacity_and_alignment(
        capacity: usize,
        alignment: usize,
        inner: R,
    ) -> io::Result<Self> {
        Ok(Self {
            inner,
            buf: AlignedBytes::new(capacity.max(1), alignment)?,
            pos: 0,
            filled: 0,
        })
    }
}

impl<R: Read> Read for AlignedBufReader<R> {
    /// Read data into the provided buffer, using the internal aligned buffer for efficient reads.
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        // if the requested read size is larger than the internal buffer and the buffer is currently empty,
        // we can bypass the internal buffer and read directly into the output buffer for better performance
        if self.pos >= self.filled && out.len() >= self.buf.len {
            self.pos = 0;
            self.filled = 0;
            return self.inner.read(out);
        }

        // otherwise, we need to fill the internal buffer and copy data from it to the output buffer
        let available = self.fill_buf()?;
        if available.is_empty() {
            return Ok(0);
        }

        // copy as much data as possible from the internal buffer to the output buffer,
        // and then consume the internal buffer accordingly
        let len = available.len().min(out.len());
        out[..len].copy_from_slice(&available[..len]);
        self.consume(len);
        Ok(len)
    }
}

impl<R: Read> BufRead for AlignedBufReader<R> {
    // fill the internal buffer if it's empty and return a slice of the available data
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        // if the internal buffer has been fully consumed, we need to read more data into it
        if self.pos >= self.filled {
            self.filled = self.inner.read(self.buf.as_mut_slice())?;
            self.pos = 0;
        }
        Ok(&self.buf.as_slice()[self.pos..self.filled])
    }

    // consume a specified amount of data from the internal buffer by advancing the position,
    // ensuring we do not consume more than what is currently filled in the buffer
    fn consume(&mut self, amt: usize) {
        self.pos = (self.pos + amt).min(self.filled);
    }
}

/// AlignedBufWriter is a buffered writer that uses an aligned buffer for efficient I/O operations.
struct AlignedBufWriter<W: Write> {
    inner: W, // the underlying writer that will receive the data when the buffer is flushed
    buf: AlignedBytes, // the internal aligned buffer used for writing data before flushing to the inner writer
    filled: usize, // the amount of data currently filled in the buffer that has not yet been flushed to the inner writer
}

impl<W: Write> AlignedBufWriter<W> {
    // Create a new AlignedBufWriter with the specified capacity, alignment, and inner writer.
    fn with_capacity_and_alignment(
        capacity: usize,
        alignment: usize,
        inner: W,
    ) -> io::Result<Self> {
        Ok(Self {
            inner,
            buf: AlignedBytes::new(capacity.max(1), alignment)?,
            filled: 0,
        })
    }
}

impl<W: Write> AlignedBufWriter<W> {
    // flush the internal buffer to the inner writer, ensuring that all filled data is written out,
    // and then reset the filled count to zero for the next round of buffering
    fn flush_buf(&mut self) -> io::Result<()> {
        let mut written = 0;
        while written < self.filled {
            let n = self
                .inner
                .write(&self.buf.as_slice()[written..self.filled])?;
            if n == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::WriteZero,
                    "failed to flush aligned I/O buffer",
                ));
            }
            written += n;
        }
        self.filled = 0;
        Ok(())
    }
}

impl<W: Write> Drop for AlignedBufWriter<W> {
    fn drop(&mut self) {
        let _ = self.flush();
    }
}

impl<W: Write> Write for AlignedBufWriter<W> {
    // write data to the internal buffer, and if the buffer is full or the incoming data is larger than the buffer,
    // flush the buffer to the inner writer before writing more data, ensuring efficient use of the
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        // if the internal buffer is currently empty and the incoming data is larger than or equal to the buffer size,
        // we can bypass the internal buffer and write directly to the inner writer for better performance
        if self.filled == 0 && data.len() >= self.buf.len {
            return self.inner.write(data);
        }

        // if the incoming data does not fit in the remaining space of the internal buffer,
        // we need to flush the buffer first
        if data.len() > self.buf.len - self.filled {
            self.flush_buf()?;
            if data.len() >= self.buf.len {
                return self.inner.write(data);
            }
        }

        // copy the incoming data to the internal buffer and update the filled count accordingly
        let end = self.filled + data.len();
        self.buf.as_mut_slice()[self.filled..end].copy_from_slice(data);
        self.filled = end;
        Ok(data.len())
    }

    // flush the internal buffer to the inner writer and then flush the inner writer to ensure all data is written out
    fn flush(&mut self) -> io::Result<()> {
        self.flush_buf()?;
        self.inner.flush()
    }
}

/// xopen is a helper function that opens a file for reading and returns a buffered reader that automatically detects compression formats.
pub fn xopen(file: &str, buf_size: usize) -> io::Result<Box<dyn BufRead>> {
    xopen_with_alignment(file, buf_size, DEFAULT_IO_BUFFER_ALIGNMENT)
}

/// xopen_with_alignment is a helper function that opens a file for reading and returns a buffered reader with
/// the specified buffer size and alignment, automatically detecting compression formats.
pub fn xopen_with_alignment(
    file: &str,
    buf_size: usize,
    buf_align: usize,
) -> io::Result<Box<dyn BufRead>> {
    let buf_size = buf_size.max(4096);

    let mut r: Box<dyn BufRead> = if file == "-" {
        if io::stdin().is_terminal() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "no data detected in STDIN",
            ));
        }
        // if reading from STDIN, we need to use an aligned buffer reader to ensure proper alignment for efficient I/O operations,
        // and we can directly wrap the locked STDIN in the aligned buffer reader without an additional
        Box::new(AlignedBufReader::with_capacity_and_alignment(
            buf_size,
            buf_align,
            io::stdin().lock(),
        )?)
    } else {
        // if reading from a file, we can open the file and wrap it in an aligned buffer reader to ensure proper alignment for efficient I/O operations,
        // and we can directly wrap the file in the aligned buffer reader without an additional layer of buffering since the aligned buffer reader already provides buffering functionality
        Box::new(AlignedBufReader::with_capacity_and_alignment(
            buf_size,
            buf_align,
            File::open(file)?,
        )?)
    };

    // check compression formats
    let buf = r.fill_buf()?; // peek without consuming

    let reader: Box<dyn BufRead> = if buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b {
        // gzip
        Box::new(AlignedBufReader::with_capacity_and_alignment(
            buf_size,
            buf_align,
            GzDecoder::new(r),
        )?)
    } else if buf.starts_with(&[0xFD, b'7', b'z', b'X', b'Z', 0x00]) {
        // xz
        Box::new(AlignedBufReader::with_capacity_and_alignment(
            buf_size,
            buf_align,
            XzDecoder::new(r),
        )?)
    } else if buf.starts_with(b"BZh") {
        // bzip2
        Box::new(AlignedBufReader::with_capacity_and_alignment(
            buf_size,
            buf_align,
            BzDecoder::new(r),
        )?)
    } else if buf.starts_with(&[0x28, 0xB5, 0x2F, 0xFD]) {
        // zstd
        Box::new(AlignedBufReader::with_capacity_and_alignment(
            buf_size,
            buf_align,
            ZstdDecoder::new(r)?,
        )?)
    } else {
        // no compression
        r
    };

    Ok(reader)
}

/// xwrite is a helper function that opens a file for writing
/// and returns a buffered writer that automatically detects compression formats based on the file extension.
pub fn xwrite(path: &str, buf_size: usize) -> io::Result<Box<dyn Write>> {
    xwrite_with_alignment(path, buf_size, DEFAULT_IO_BUFFER_ALIGNMENT)
}

/// xwrite_with_alignment is a helper function that opens a file for writing and returns a buffered writer with
/// the specified buffer size and alignment, automatically detecting compression formats based on the file extension.
pub fn xwrite_with_alignment(
    path: &str,
    buf_size: usize,
    buf_align: usize,
) -> io::Result<Box<dyn Write>> {
    let buf_size = buf_size.max(4096);

    if path == "-" {
        return Ok(Box::new(AlignedBufWriter::with_capacity_and_alignment(
            buf_size,
            buf_align,
            io::stdout(),
        )?));
    }

    let file = File::create(path)?;

    let writer: Box<dyn Write> = if path.ends_with(".gz") {
        Box::new(AlignedBufWriter::with_capacity_and_alignment(
            buf_size,
            buf_align,
            GzEncoder::new(file, Compression::default()),
        )?)
    } else if path.ends_with(".xz") {
        Box::new(AlignedBufWriter::with_capacity_and_alignment(
            buf_size,
            buf_align,
            XzEncoder::new(file, 6),
        )?)
    } else if path.ends_with(".bz2") {
        Box::new(AlignedBufWriter::with_capacity_and_alignment(
            buf_size,
            buf_align,
            BzEncoder::new(file, bzip2::Compression::default()),
        )?)
    } else if path.ends_with(".zst") || path.ends_with(".zstd") {
        let encoder = ZstdEncoder::new(file, 0)?; // level 0 = default
        Box::new(AlignedBufWriter::with_capacity_and_alignment(
            buf_size,
            buf_align,
            encoder.auto_finish(),
        )?)
    } else {
        // no compression
        Box::new(AlignedBufWriter::with_capacity_and_alignment(
            buf_size, buf_align, file,
        )?)
    };

    Ok(writer)
}

#[cfg(test)]
mod xwrite_drop_tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(suffix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "fastx-xopen-test-{}-{nanos}{suffix}",
            std::process::id()
        ))
    }

    #[test]
    fn test_xwrite_flushes_on_drop_for_plain_file() {
        let path = temp_path(".fa");
        {
            let mut writer = xwrite(path.to_str().unwrap(), 8192).unwrap();
            writer.write_all(b">chr1\nACGT\n>chrM\nTGCA\n").unwrap();
        }

        let content = fs::read(&path).unwrap();
        fs::remove_file(&path).unwrap();
        assert_eq!(content, b">chr1\nACGT\n>chrM\nTGCA\n");
    }

    #[test]
    fn test_xwrite_flushes_on_drop_for_gzip_file() {
        let path = temp_path(".fa.gz");
        {
            let mut writer = xwrite(path.to_str().unwrap(), 8192).unwrap();
            writer.write_all(b">chr1\nACGT\n>chrM\nTGCA\n").unwrap();
        }

        let mut reader = xopen(path.to_str().unwrap(), 8192).unwrap();
        let mut content = Vec::new();
        reader.read_to_end(&mut content).unwrap();
        fs::remove_file(&path).unwrap();
        assert_eq!(content, b">chr1\nACGT\n>chrM\nTGCA\n");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use std::io::Cursor;
    use std::rc::Rc;

    #[test]
    fn test_aligned_buf_reader_alignment() {
        let reader =
            AlignedBufReader::with_capacity_and_alignment(1024, 256, Cursor::new(Vec::<u8>::new()))
                .unwrap();
        assert_eq!(reader.buf.ptr.as_ptr() as usize % 256, 0);
    }

    #[test]
    fn test_aligned_buf_writer_alignment() {
        let writer =
            AlignedBufWriter::with_capacity_and_alignment(1024, 256, Vec::<u8>::new()).unwrap();
        assert_eq!(writer.buf.ptr.as_ptr() as usize % 256, 0);
    }

    #[test]
    fn test_invalid_alignment_is_rejected() {
        assert!(
            AlignedBufReader::with_capacity_and_alignment(
                1024,
                3000,
                Cursor::new(Vec::<u8>::new())
            )
            .is_err()
        );
    }

    struct SpyReader {
        data: Vec<u8>,
        pos: usize,
        read_sizes: Rc<RefCell<Vec<usize>>>,
    }

    impl Read for SpyReader {
        fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
            self.read_sizes.borrow_mut().push(out.len());
            if self.pos >= self.data.len() {
                return Ok(0);
            }

            let len = out.len().min(self.data.len() - self.pos);
            out[..len].copy_from_slice(&self.data[self.pos..self.pos + len]);
            self.pos += len;
            Ok(len)
        }
    }

    #[test]
    fn test_aligned_buf_reader_large_read_bypasses_internal_buffer() {
        let read_sizes = Rc::new(RefCell::new(Vec::new()));
        let inner = SpyReader {
            data: b"abcdefgh".to_vec(),
            pos: 0,
            read_sizes: Rc::clone(&read_sizes),
        };
        let mut reader = AlignedBufReader::with_capacity_and_alignment(4, 256, inner).unwrap();
        let mut out = [0u8; 8];

        let n = reader.read(&mut out).unwrap();

        assert_eq!(n, 8);
        assert_eq!(&out, b"abcdefgh");
        assert_eq!(read_sizes.borrow().as_slice(), &[8]);
    }
}
