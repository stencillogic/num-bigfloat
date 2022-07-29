//! Utility components.


/// writable buffer.
pub struct WritableBuf<'a> {
    buf: &'a mut [u8],
    offset: usize,
}

impl<'a> WritableBuf<'a> {
    pub fn new(buf: &'a mut [u8]) -> Self {
        WritableBuf {
            buf,
            offset: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.offset
    }
}

impl<'a> core::fmt::Write for WritableBuf<'a> {
    fn write_str(&mut self, s: &str) -> core::fmt::Result {
        let bytes = s.as_bytes();
        let rem = &mut self.buf[self.offset..];
        if rem.len() < bytes.len() {
            Err(core::fmt::Error)
        } else {
            let rem = &mut rem[..bytes.len()];
            rem.copy_from_slice(bytes);
            self.offset += bytes.len();
            Ok(())
        }
    }
}


/// Concatenate several &str instances from `str_list` into one unsing buffer `buf`.
pub fn concat_str<'a>(buf: &'a mut [u8], str_list: &[&str]) -> &'a str {
    let mut p = 0;
    for s in str_list {
        let bytes = s.as_bytes();
        buf[p..p+bytes.len()].copy_from_slice(bytes);
        p += bytes.len();
    }
    if buf.len() > p + 1 {
        buf[p] = 0;
    }
    core::str::from_utf8(&buf[0..p]).unwrap()
}