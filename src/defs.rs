//!  Definitions.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Number of "digits" in BigFloat number.
pub const DECIMAL_PARTS: usize = 10;

/// i64::MAX.
pub const I64_MAX: BigFloatNum = BigFloatNum {
    m: [0, 0, 0, 0, 0, 8070, 4775, 3685, 3720, 9223],
    n: DECIMAL_POSITIONS as i16,
    sign: DECIMAL_SIGN_POS,
    e: -(DECIMAL_POSITIONS as i8 - 19),
};

/// i64::MIN.
pub const I64_MIN: BigFloatNum = BigFloatNum {
    m: [0, 0, 0, 0, 0, 8080, 4775, 3685, 3720, 9223],
    n: DECIMAL_POSITIONS as i16,
    sign: DECIMAL_SIGN_NEG,
    e: -(DECIMAL_POSITIONS as i8 - 19),
};

/// u64::MAX.
pub const U64_MAX: BigFloatNum = BigFloatNum {
    m: [0, 0, 0, 0, 0, 1615, 955, 737, 6744, 1844],
    n: DECIMAL_POSITIONS as i16,
    sign: DECIMAL_SIGN_POS,
    e: -(DECIMAL_POSITIONS as i8 - 20),
};

// i128::MAX
pub const I128_MAX: BigFloatNum = BigFloatNum {
    m: [7270, 4105, 1588, 3037, 1687, 3173, 4692, 3460, 4118, 1701],
    n: DECIMAL_POSITIONS as i16,
    sign: DECIMAL_SIGN_POS,
    e: -1,
};

/// i128::MIN.
pub const I128_MIN: BigFloatNum = BigFloatNum {
    m: [7280, 4105, 1588, 3037, 1687, 3173, 4692, 3460, 4118, 1701],
    n: DECIMAL_POSITIONS as i16,
    sign: DECIMAL_SIGN_NEG,
    e: -1,
};

/// u128::MAX.
pub const U128_MAX: BigFloatNum = BigFloatNum {
    m: [4550, 8211, 3176, 6074, 3374, 6346, 9384, 6920, 8236, 3402],
    n: DECIMAL_POSITIONS as i16,
    sign: DECIMAL_SIGN_POS,
    e: -1,
};

/// Number representation.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BigFloatNum {
    pub(crate) sign: i8,                // sign
    pub(crate) e: i8,                   // exponent
    pub(crate) n: i16, // the number of decimal positions in the mantissa excluding leading zeroes
    pub(crate) m: [i16; DECIMAL_PARTS], // mantissa
}

/// Possible errors.
#[derive(Eq, PartialEq, Debug, Copy, Clone)]
pub enum Error {
    /// Exponent value becomes greater than the upper bound.
    /// Error stores sign of resulting number.
    ExponentOverflow(i8),

    /// Divizor is zero.
    DivisionByZero,

    /// Argument must not be a negative number.
    ArgumentIsNegative,

    /// Invalid argument.
    InvalidArgument,
}

/// Possible errors.
#[derive(Eq, PartialEq, Debug, Copy, Clone)]
pub enum RoundingMode {
    /// Round half toward positive infinity.
    Up,

    /// Round half toward negative infinity.
    Down,

    /// Round half toward zero.
    ToZero,

    /// Round half away from zero.
    FromZero,

    /// Round half to even.
    ToEven,

    /// Round half to odd.
    ToOdd,
}

pub const DECIMAL_BASE_LOG10: usize = 4; // number decimal positions in a digit = log10(DECIMAL_BASE)
pub const DECIMAL_POSITIONS: usize = DECIMAL_PARTS * DECIMAL_BASE_LOG10;
pub const DECIMAL_BASE: usize = 10000; // 9999 is the maximum of a digit
pub const DECIMAL_SIGN_POS: i8 = 1; // + sign
pub const DECIMAL_SIGN_NEG: i8 = -1; // - sign
pub const DECIMAL_MIN_EXPONENT: i8 = -128; // min exponent value
pub const DECIMAL_MAX_EXPONENT: i8 = 127; // max exponent value
pub const ZEROED_MANTISSA: [i16; DECIMAL_PARTS] = [0; DECIMAL_PARTS];

/// Zero.
pub const ZERO: BigFloatNum = BigFloatNum {
    m: ZEROED_MANTISSA,
    n: 0,
    sign: DECIMAL_SIGN_POS,
    e: 0,
};

/// One.
pub const ONE: BigFloatNum = BigFloatNum {
    m: [0, 0, 0, 0, 0, 0, 0, 0, 0, 1000],
    n: DECIMAL_POSITIONS as i16,
    sign: DECIMAL_SIGN_POS,
    e: 1 - (DECIMAL_POSITIONS as i8),
};

/// Two.
pub const TWO: BigFloatNum = BigFloatNum {
    m: [0, 0, 0, 0, 0, 0, 0, 0, 0, 2000],
    n: DECIMAL_POSITIONS as i16,
    sign: DECIMAL_SIGN_POS,
    e: 1 - (DECIMAL_POSITIONS as i8),
};

/// Eulers number.
pub const E: BigFloatNum = BigFloatNum {
    m: [7757, 6249, 3526, 7471, 6028, 2353, 9045, 2845, 2818, 2718],
    n: DECIMAL_POSITIONS as i16,
    sign: DECIMAL_SIGN_POS,
    e: 1 - (DECIMAL_POSITIONS as i8),
};

/// Pi number.
pub const PI: BigFloatNum = BigFloatNum {
    m: [4197, 288, 2795, 3383, 6264, 2384, 9793, 5358, 5926, 3141],
    n: DECIMAL_POSITIONS as i16,
    sign: DECIMAL_SIGN_POS,
    e: 1 - (DECIMAL_POSITIONS as i8),
};

/// Largest value possible.
pub const MAX: BigFloatNum = BigFloatNum {
    m: [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999],
    n: DECIMAL_POSITIONS as i16,
    sign: DECIMAL_SIGN_POS,
    e: DECIMAL_MAX_EXPONENT,
};

/// Smalles value possible.
pub const MIN: BigFloatNum = BigFloatNum {
    m: [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999],
    n: DECIMAL_POSITIONS as i16,
    sign: DECIMAL_SIGN_NEG,
    e: DECIMAL_MAX_EXPONENT,
};

/// Smalles positive number.
pub const MIN_POSITIVE: BigFloatNum = BigFloatNum {
    m: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    n: 1,
    sign: DECIMAL_SIGN_POS,
    e: DECIMAL_MIN_EXPONENT,
};

/// Smalles positive normal number.
pub const MIN_POSITIVE_NORMAL: BigFloatNum = BigFloatNum {
    m: [0, 0, 0, 0, 0, 0, 0, 0, 0, 1000],
    n: 1,
    sign: DECIMAL_SIGN_POS,
    e: DECIMAL_MIN_EXPONENT,
};

/// Creation and number manipulation functions.
impl BigFloatNum {
    /// Return new BigFloat with value zero.
    pub fn new() -> Self {
        BigFloatNum {
            sign: DECIMAL_SIGN_POS,
            e: 0,
            n: 0,
            m: ZEROED_MANTISSA,
        }
    }

    /// Return BigFloat with the value of 1.
    pub fn one() -> Self {
        let mut val = Self::new();
        val.m[DECIMAL_PARTS - 1] = DECIMAL_BASE as i16 / 10;
        val.n = DECIMAL_POSITIONS as i16;
        val.e = 1 - DECIMAL_POSITIONS as i8;
        val
    }

    /// Create a BigFloat value from a sequence of `bytes`. Each byte must represent a decimal digit.
    /// First byte is the most significant. The length of `bytes` can be any. If the length of
    /// `bytes` is greater than required, then the remaining part is ignored.
    /// If `sign` is negative, then the resulting BigFloat will be
    /// negative.
    pub fn from_bytes(bytes: &[u8], sign: i8, exponent: i8) -> BigFloatNum {
        let mut mantissa = ZEROED_MANTISSA;
        let mut n: usize = 0;
        let mut p: i16 = 1;
        let d = if bytes.len() > DECIMAL_POSITIONS { DECIMAL_POSITIONS } else { bytes.len() };
        for i in 1..d + 1 {
            mantissa[n] += (bytes[d - i] % 10) as i16 * p;
            p *= 10;
            if p == DECIMAL_BASE as i16 {
                n += 1;
                p = 1;
            }
        }

        BigFloatNum {
            sign: if sign >= 0 { DECIMAL_SIGN_POS } else { DECIMAL_SIGN_NEG },
            e: exponent,
            n: Self::num_digits(&mantissa),
            m: mantissa,
        }
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
        for i in 1..d + 1 {
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

    /// Return true if the number is zero.
    pub fn is_zero(&self) -> bool {
        self.n == 0
    }

    /// Return true if integer part of number is zero.
    pub fn is_int_even(&self) -> bool {
        let int = self.int();
        if int.e < 0 {
            let p = int.n + int.e as i16;
            let mut d = int.m[p as usize / DECIMAL_BASE_LOG10];
            let mut i = p % DECIMAL_BASE_LOG10 as i16;
            while i > 0 {
                d /= 10;
                i -= 1;
            }
            d & 1 == 0
        } else if int.e == 0 {
            int.m[0] & 1 == 0
        } else {
            true
        }
    }

    /// Returns true if `self` is subnormal.
    pub fn is_subnormal(&self) -> bool {
        self.n < DECIMAL_POSITIONS as i16 && self.e == DECIMAL_MIN_EXPONENT
    }
}
