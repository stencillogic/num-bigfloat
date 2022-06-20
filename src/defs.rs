///!  Definitions.

/// Number of "digits" in BigFloat number.
pub const DECIMAL_PARTS: usize = 10;

/// Number representation.
#[derive(Copy, Clone, Debug)]
pub struct BigFloat {
    pub (crate) sign: i8,                // sign
    pub (crate) e: i8,                   // exponent
    pub (crate) n: i16,                  // the number of decimal positions in the mantissa excluding leading zeroes
    pub (crate) m: [i16; DECIMAL_PARTS], // mantissa
}

/// Possible errors.
#[derive(Eq, PartialEq, Debug)]
pub enum Error {
    /// Exponent value becomes greater than the upper bound or smaller than the lower bound for exponent value.
    ExponentOverflow,   

    /// Divizor is zero.
    DivisionByZero,     

    /// Argument must not be a negative number.
    ArgumentIsNegative,

    /// Invalid argument.
    InvalidArgument,
}


pub const DECIMAL_BASE_LOG10: usize = 4;    // number decimal positions in a digit = log10(DECIMAL_BASE)
pub const DECIMAL_POSITIONS: usize = DECIMAL_PARTS * DECIMAL_BASE_LOG10;
pub const DECIMAL_BASE: usize = 10000;      // 9999 is the maximum of a digit
pub const DECIMAL_SIGN_POS: i8 = 1;         // + sign
pub const DECIMAL_SIGN_NEG: i8 = -1;        // - sign
pub const DECIMAL_MIN_EXPONENT: i8 = -128;  // min exponent value
pub const DECIMAL_MAX_EXPONENT: i8 = 127;   // max exponent value
pub const DECIMAL_MAX_EXPONENT_POSITIONS: i16 = 3;  // max decimal positions in exponent
pub const ZEROED_MANTISSA: [i16; DECIMAL_PARTS] = [0; DECIMAL_PARTS];


/// Eulers number.
pub const E: BigFloat = BigFloat {
    m: [7757, 6249, 3526, 7471, 6028, 2353, 9045, 2845, 2818, 2718],
    n: DECIMAL_POSITIONS as i16, 
    sign: DECIMAL_SIGN_POS, 
    e: 1 - (DECIMAL_POSITIONS as i8),
};

/// Pi number.
pub const PI: BigFloat = BigFloat {
    m: [4197, 0288, 2795, 3383, 6264, 2384, 9793, 5358, 5926, 3141],
    n: DECIMAL_POSITIONS as i16, 
    sign: DECIMAL_SIGN_POS, 
    e: 1 - (DECIMAL_POSITIONS as i8),
};


/// Creation and number manipulation functions.
impl BigFloat {

    /// Return new BigFloat with value zero.
    pub fn new() -> Self {
        return BigFloat {
            sign: DECIMAL_SIGN_POS,
            e: 0,
            n: 0,
            m: ZEROED_MANTISSA,
        };
    }
    
    /// Return BigFloat with the value of 1.
    pub fn one() -> Self {
        let mut val = Self::new();
        val.m[DECIMAL_PARTS-1] = DECIMAL_BASE as i16/10;
        val.n = DECIMAL_POSITIONS as i16;
        val.e = 1 - DECIMAL_POSITIONS as i8;
        return val;
    }

    /// Return new BigFloat with value two.
    pub fn two() -> Self {
        let mut val = Self::new();
        val.m[DECIMAL_PARTS-1] = DECIMAL_BASE as i16/5;
        val.n = DECIMAL_POSITIONS as i16;
        val.e = 1 - DECIMAL_POSITIONS as i8;
        return val;
    }

    /// Create a BigFloat value from a sequence of `bytes`. Each byte must represent a decimal digit.
    /// First byte is the most significant. The length of `bytes` can be any. If the length of
    /// `bytes` is greater than required, then the remaining part is ignored.
    /// If `sign` is negative, then the resulting BigFloat will be
    /// negative.
    pub fn from_bytes(bytes: &[u8], sign: i8, exponent: i8) -> BigFloat {
        let mut mantissa = ZEROED_MANTISSA;
        let mut n: usize = 0;
        let mut p: i16 = 1;
        let d = if bytes.len() > DECIMAL_POSITIONS { DECIMAL_POSITIONS } else { bytes.len() };
        for i in 1..d+1 {
            mantissa[n] += (bytes[d - i] % 10) as i16 * p;
            p *= 10;
            if p == DECIMAL_BASE as i16 {
                n += 1;
                p = 1;
            }
        }

        return BigFloat {
            sign: if sign >= 0 { DECIMAL_SIGN_POS } else { DECIMAL_SIGN_NEG },
            e: exponent,
            n: d as i16,
            m: mantissa,
        };
    }

    /// Construct BigFloat from f64.
    ///
    /// # Errors
    ///
    /// ExponentOverflow - when result is too big or too small.
    pub fn from_f64(mut f: f64) -> Result<Self, Error> {
        let mut e: i32 = 0;
        let mut ret = BigFloat::new();
        if f == 0f64 {
            return Ok(ret);
        }
        if f == f64::INFINITY || f == f64::NAN {
            return Err(Error::ExponentOverflow);
        } 
        if f < 0f64 {
            ret.sign = DECIMAL_SIGN_NEG;
            f = -f;
        }
        // bring to 0.xxxxxxxxx
        while f >= 1.0f64 {
            f /= 10f64;
            e += 1;
        }
        while f < 0.1f64 {
            f *= 10f64;
            e -= 1;
        }
        // fill-in mantissa
        ret.n = DECIMAL_POSITIONS as i16;
        let mut p = DECIMAL_PARTS - 1;
        loop {
            f *= DECIMAL_BASE as f64;
            let d = f as i16;
            f = f - d as f64;
            ret.m[p] = d;
            p -= 1;
            if f == 0f64 || p == 0 {
                break;
            }
        }
        
        e -= DECIMAL_POSITIONS as i32;
        if e < DECIMAL_MIN_EXPONENT as i32 || e > DECIMAL_MAX_EXPONENT as i32 {
            return Err(Error::ExponentOverflow);
        }
        ret.e = e as i8;

        return Ok(ret);
    }

    /// Convert BigFloat to f64.
    pub fn to_f64(&self) -> f64 {
        let mut f: f64 = 0f64;
        for i in 0..DECIMAL_PARTS {
            f += self.m[i] as f64;
            f /= DECIMAL_BASE as f64;
        }
        let mut e = self.n + self.e as i16;
        while e < 0 {
            f /= 10f64;
            e += 1;
        }
        while e > 0 {
            f *= 10f64;
            e -= 1;
        }
        if self.sign == DECIMAL_SIGN_NEG {
            f = -f;
        }
        return f;
    }

    /// Construct BigFloat from f32. Wrapper for from_f64.
    ///
    /// # Errors
    ///
    /// ExponentOverflow - when result is too big or too small.
    pub fn from_f32(f: f32) -> Result<Self, Error> {
        Self::from_f64(f as f64)
    }

    /// Convert BigFloat to f32. Wrapper for to_f64
    pub fn to_f32(&self) -> f32 {
        self.to_f64() as f32
    }

    /// Get BigFloat's mantissa as bytes. Each byte represents a decimal digit.
    /// First byte is the most significant. The length of `bytes` can be any. If the length of
    /// `bytes` is smaller than required, then remaining part of mantissa will be omitted.
    ///
    /// The length of mantissa can be determined using `get_mantissa_len`.
    pub fn get_mantissa_bytes(&self, bytes: &mut [u8]) {
        let mut n: usize = 0;
        let mut p: i16 = 1;
        let d = if bytes.len() < self.n as usize { bytes.len() } else { self.n as usize };
        for i in 1..d+1 {
            bytes[d - i] = ((self.m[n] / p) % 10) as u8;
            p *= 10;
            if p == DECIMAL_BASE as i16 {
                n += 1;
                p = 1;
            }
        }
    }

    /// Return the number of decimal positions filled in the mantissa.
    pub fn get_mantissa_len(&self) -> usize {
        self.n as usize
    }

    /// Return 1 if BigFloat is positive, -1 otherwise.
    pub fn get_sign(&self) -> i8 {
        self.sign
    }

    /// Return exponent part.
    pub fn get_exponent(&self) -> i8 {
        self.e
    }

    /// Return raw parts of BigFloat: mantissa, number of decimal positions in mantissa, sing, and
    /// exponent.
    pub fn to_raw_parts(&self) -> ([i16; DECIMAL_PARTS], i16, i8, i8) {
        (self.m, self.n, self.sign, self.e)
    }

    /// Construct BigFloat from raw parts.
    pub fn from_raw_parts(mantissa: [i16; DECIMAL_PARTS], mantissa_len: i16, sign: i8, exponent: i8) -> Self {
        return BigFloat {
            sign: sign,
            e: exponent,
            n: mantissa_len,
            m: mantissa,
        };
    }
}