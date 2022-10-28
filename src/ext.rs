//! BigFloat number whith support of `NaN`, and `Inf` 
//! values, and implementation of `std::ops` traits.

use crate::defs::{BigFloatNum, Error, DECIMAL_SIGN_POS, DECIMAL_PARTS, DECIMAL_SIGN_NEG, DECIMAL_POSITIONS,
    DECIMAL_MAX_EXPONENT, DECIMAL_MIN_EXPONENT, RoundingMode, DECIMAL_BASE_LOG10, DECIMAL_BASE, I64_MAX, I64_MIN, U64_MAX, I128_MIN, I128_MAX, U128_MAX};
use crate::util::WritableBuf;
use core::num::FpCategory;

#[cfg(feature="std")]
use std::fmt::Write;

#[cfg(not(feature="std"))]
use core::fmt::Write;

#[cfg(feature="rand")]
use rand::random;

#[cfg(feature="serde")]
use serde::{Serialize, Deserialize};

/// Maximum value possible.
pub const MAX: BigFloat = BigFloat {inner: Flavor::Value(crate::defs::MAX)};

/// Maximum possible exponent.
pub const MAX_EXP: i8 = DECIMAL_MAX_EXPONENT;

/// Minumum value possible.
pub const MIN: BigFloat = BigFloat {inner: Flavor::Value(crate::defs::MIN)};

/// Minumum possible exponent.
pub const MIN_EXP: i8 = DECIMAL_MIN_EXPONENT;

/// The smalles positive number.
pub const MIN_POSITIVE: BigFloat = BigFloat {inner: Flavor::Value(crate::defs::MIN_POSITIVE)};

/// The smalles positive normal number.
pub const MIN_POSITIVE_NORMAL: BigFloat = BigFloat {inner: Flavor::Value(crate::defs::MIN_POSITIVE_NORMAL)};

/// Radix of BigFloat
pub const RADIX: u32 = 10;

/// NaN representation.
pub const NAN: BigFloat = BigFloat { inner: Flavor::NaN };

/// Positive infinity.
pub const INF_POS: BigFloat = BigFloat { inner: Flavor::Inf(DECIMAL_SIGN_POS) };

/// Negative infinity.
pub const INF_NEG: BigFloat = BigFloat { inner: Flavor::Inf(DECIMAL_SIGN_NEG) };

/// Value of zero.
pub const ZERO: BigFloat = BigFloat { inner: Flavor::Value(crate::defs::ZERO) };

/// Value of one.
pub const ONE: BigFloat = BigFloat { inner: Flavor::Value(crate::defs::ONE) };

/// Value of two.
pub const TWO: BigFloat = BigFloat { inner: Flavor::Value(crate::defs::TWO) };

/// Euler's number.
pub const E: BigFloat = BigFloat { inner: Flavor::Value(crate::defs::E) };

/// PI number.
pub const PI: BigFloat = BigFloat { inner: Flavor::Value(crate::defs::PI) };

/// PI / 2.
pub const HALF_PI: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    m: [2099, 5144, 6397, 1691, 3132, 6192, 4896, 2679, 7963, 1570],
    n: DECIMAL_POSITIONS as i16, 
    sign: DECIMAL_SIGN_POS, 
    e: -(DECIMAL_POSITIONS as i8 - 1),
})};

/// SQRT(2)
pub const SQRT_2: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: -(DECIMAL_POSITIONS as i8) + 1, 
    n: DECIMAL_POSITIONS as i16, 
    m: [8570, 9807, 2096, 8724, 168, 488, 3095, 6237, 2135, 1414] 
})};

/// 1 / PI
pub const FRAC_1_PI: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: -(DECIMAL_POSITIONS as i8), 
    n: DECIMAL_POSITIONS as i16, 
    m: [689, 8724, 4502, 5267, 7767, 7153, 7906, 6183, 988, 3183]
})};

/// 1 / SQRT(2)
pub const FRAC_1_SQRT_2: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: -(DECIMAL_POSITIONS as i8), 
    n: DECIMAL_POSITIONS as i16, 
    m: [2848, 9039, 484, 3621, 844, 2440, 5475, 1186, 678, 7071]
})};

/// 2 / PI
pub const FRAC_2_PI: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: -(DECIMAL_POSITIONS as i8), 
    n: DECIMAL_POSITIONS as i16, 
    m: [1378, 7448, 9005, 534, 5535, 4307, 5813, 2367, 1977, 6366]
})};

/// 2 / SQRT(PI)
pub const FRAC_2_SQRT_PI: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: -(DECIMAL_POSITIONS as i8), 
    n: DECIMAL_POSITIONS as i16, 
    m: [8440, 2585, 6077, 4515, 8079, 8694, 7562, 3547, 8958, 5641]
})};

/// PI / 3
pub const FRAC_PI_3: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: -(DECIMAL_POSITIONS as i8) + 1, 
    n: DECIMAL_POSITIONS as i16, 
    m: [8066, 6762, 931, 4461, 5421, 7461, 6597, 5119, 1975, 1047]
})};

/// PI / 4
pub const FRAC_PI_4: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: -(DECIMAL_POSITIONS as i8), 
    n: DECIMAL_POSITIONS as i16, 
    m: [493, 5721, 1987, 8458, 5660, 961, 4483, 3397, 9816, 7853]
})};

/// PI / 6
pub const FRAC_PI_6: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: -(DECIMAL_POSITIONS as i8), 
    n: DECIMAL_POSITIONS as i16, 
    m: [329, 3814, 4658, 2305, 7107, 7307, 2988, 5598, 9877, 5235]
})};

/// PI / 8
pub const FRAC_PI_8: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: -(DECIMAL_POSITIONS as i8), 
    n: DECIMAL_POSITIONS as i16, 
    m: [5246, 7860, 993, 4229, 7830, 5480, 7241, 1698, 9908, 3926]
})};

/// ln(10)
pub const LN_10: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: -(DECIMAL_POSITIONS as i8) + 1, 
    n: DECIMAL_POSITIONS as i16, 
    m: [7601, 6420, 6843, 1454, 1799, 6840, 4045, 9299, 5850, 2302]
})};

/// ln(2)
pub const LN_2: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: -(DECIMAL_POSITIONS as i8), 
    n: DECIMAL_POSITIONS as i16, 
    m: [755, 6568, 5817, 1214, 7232, 941, 9453, 559, 4718, 6931]
})};

/// log10(E)
pub const LOG10_E: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: -1,
    e: -(DECIMAL_POSITIONS as i8), 
    n: DECIMAL_POSITIONS as i16, 
    m: [318, 6507, 4129, 4666, 288, 8946, 8323, 1216, 6189, 2386]
})};

/// log2(E)
pub const LOG2_E: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: -1,
    e: -(DECIMAL_POSITIONS as i8), 
    n: DECIMAL_POSITIONS as i16, 
    m: [6474, 6789, 169, 9107, 4620, 3637, 1469, 1612, 1764, 7928]
})};

/// The difference between 1 and the smallest floating point number greater than 1.
pub const EPSILON: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: 0,
    n: 1,
    m: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
})};

/// 180 / PI
#[cfg(not(feature="std"))]
#[cfg(feature="num-traits")]
pub const RAD_TO_DEG_FACTOR: BigFloat = BigFloat { inner: Flavor::Value(BigFloatNum {
    sign: 1,
    e: -(DECIMAL_POSITIONS as i8) + 2,
    n: DECIMAL_POSITIONS as i16,
    m: [3240, 1703, 4105, 5481, 7981, 876, 8232, 5130, 5779, 5729]
})};

/// Number representation.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature="serde", derive(Serialize, Deserialize))]
pub struct BigFloat {
    inner: Flavor,
}

#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature="serde", derive(Serialize, Deserialize))]
enum Flavor {
    Value(BigFloatNum),
    NaN,
    Inf(i8)         // signed Inf
}

impl BigFloat {

    /// Returns a new BigFloat with the value of zero.
    pub fn new() -> Self {
        BigFloat {
            inner: Flavor::Value(BigFloatNum::new())
        }
    }
 

    /// Creates a BigFloat value from a sequence of `bytes`. Each byte must represent a decimal digit.
    /// First byte is the most significant. `bytes` can be of any length. 
    /// If `bytes` is longer than required, then the remaining bytes are ignored.
    /// If the "sign" is negative, then the resulting BigFloat will be negative.
    /// 
    /// ## Examples
    /// 
    /// ```
    /// # use num_bigfloat::BigFloat;
    /// let n1 = BigFloat::from_bytes(&[1,2,3,4,2,0,0,0], 1, -3);
    /// let n2 = BigFloat::from_u32(12342);
    /// assert!(n1.cmp(&n2) == Some(0));
    /// ```
    pub fn from_bytes(bytes: &[u8], sign: i8, exponent: i8) -> Self {
        BigFloat {
            inner: Flavor::Value(BigFloatNum::from_bytes(bytes, sign, exponent))
        }
    }

    /// Creates a BigFloat from f64.
    /// The conversion is not guaranteed to be lossless since BigFloat and f64 have different bases.
    pub fn from_f64(f: f64) -> Self {
        #[cfg(feature = "std")] {
            let strepr = format!("{:e}", f);
            Self::parse(&strepr).unwrap()
        }
        #[cfg(not(feature = "std"))] {
            let mut buf = [0u8; 64];
            let mut strepr = WritableBuf::new(&mut buf); // sign + number + dot + "e" + sign + exponent
            write!(strepr, "{:e}", f).unwrap();
            Self::parse(core::str::from_utf8(&buf).unwrap()).unwrap()
        }
    }

    /// Creates a BigFloat from f32.
    /// The conversion is not guaranteed to be lossless since BigFloat and f64 have different bases.
    pub fn from_f32(f: f32) -> Self {
        Self::from_f64(f as f64)
    }

    /// Converts BigFloat to f64.
    pub fn to_f64(&self) -> f64 {
        match self.inner {
            Flavor::Value(_) => {
                let mut buf = [0; 64];
                let mut w = WritableBuf::new(&mut buf);
                self.write_str(&mut w).unwrap();
                let written_len = w.len();
                let s = core::str::from_utf8(&buf[0..written_len]).unwrap();
                str::parse::<f64>(s).unwrap()
            },
            Flavor::Inf(s) => if s == DECIMAL_SIGN_POS {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            },
            Flavor::NaN => f64::NAN,
        }
    }

    /// Converts BigFloat to f32.
    pub fn to_f32(&self) -> f32 {
        self.to_f64() as f32
    }

    /// Converts BigFloat to i64. 
    /// The function retursn None if `self` is Inf, NaN, or out of range of i64.
    pub fn to_i64(&self) -> Option<i64> {
        self.to_int::<i64>(&I64_MIN, &I64_MAX)
    }

    /// Converts BigFloat to i128. 
    /// The function retursn None if `self` is Inf, NaN, or out of range of i128.
    pub fn to_i128(&self) -> Option<i128> {
        self.to_int::<i128>(&I128_MIN, &I128_MAX)
    }

    /// Converts absolute value of `self` to u64. 
    /// The function retursn None if `self` is Inf, NaN, or out of range of u64.
    pub fn to_u64(&self) -> Option<u64> {
        self.to_uint::<u64>(&U64_MAX)
    }

    /// Converts absolute value of `self` to u128. 
    /// The function retursn None if `self` is Inf, NaN, or out of range of u128.
    pub fn to_u128(&self) -> Option<u128> {
        self.to_uint::<u128>(&U128_MAX)
    }

    fn to_int<T>(&self, min: &BigFloatNum, max: &BigFloatNum) -> Option<T> 
    where T: core::ops::DivAssign<T> + core::ops::AddAssign<T> + core::ops::MulAssign<T> 
            + core::ops::SubAssign<T> + core::convert::From<u32> + core::convert::From<i16> {
        match self.inner {
            Flavor::Value(v) => {

                let int = v.int();

                if int.cmp(min) >= 0 && int.cmp(max) <= 0 {

                    let mut ret: T = 0i16.into();
                    let mut n = int.n as usize + DECIMAL_BASE_LOG10 - 1;
                    n /= DECIMAL_BASE_LOG10;
                    let mut miter = int.m[..n].iter().rev();
                    n *= DECIMAL_BASE_LOG10;
                    n = (n as i16 + int.e as i16) as usize;

                    while n >= DECIMAL_BASE_LOG10 {
                        ret *= (DECIMAL_BASE as u32).into();
                        if v.sign == DECIMAL_SIGN_POS {
                            ret += (*miter.next().unwrap()).into();
                        } else {
                            ret -= (*miter.next().unwrap()).into();
                        }
                        n -= DECIMAL_BASE_LOG10;
                    }

                    if n > 0 {
                        let mut d: T = (DECIMAL_BASE as u32).into();
                        while n > 0 {
                            d /= 10i16.into();
                            n -= 1;
                        }
                        let mut p: T = (*miter.next().unwrap()).into();
                        p /= d;
                        if v.sign == DECIMAL_SIGN_POS {
                            ret += p;
                        } else {
                            ret -= p;
                        }
                    }

                    Some(ret)

                } else {

                    None
                }
            },
            Flavor::Inf(_) => None,
            Flavor::NaN => None,
        }
    }

    fn to_uint<T>(&self, max: &BigFloatNum) -> Option<T> 
    where T: core::ops::DivAssign<T> + core::ops::AddAssign<T> + core::ops::MulAssign<T> 
            + core::convert::From<u32> + core::convert::From<u16> {
        match self.inner {
            Flavor::Value(v) => {

                let int = v.int().abs();

                if int.cmp(max) <= 0 {

                    let mut ret: T = 0u16.into();
                    let mut n = int.n as usize + DECIMAL_BASE_LOG10 - 1;
                    n /= DECIMAL_BASE_LOG10;
                    let mut miter = int.m[..n].iter().rev();
                    n *= DECIMAL_BASE_LOG10;
                    n = (n as i16 + int.e as i16) as usize;

                    while n >= DECIMAL_BASE_LOG10 {
                        ret *= (DECIMAL_BASE as u32).into();
                        ret += (*miter.next().unwrap() as u16).into();
                        n -= DECIMAL_BASE_LOG10;
                    }

                    if n > 0 {
                        let mut d: T = (DECIMAL_BASE as u32).into();
                        while n > 0 {
                            d /= 10u16.into();
                            n -= 1;
                        }
                        let mut p: T = (*miter.next().unwrap() as u16).into();
                        p /= d;
                        ret += p;
                    }

                    Some(ret)

                } else {

                    None
                }
            },
            Flavor::Inf(_) => None,
            Flavor::NaN => None,
        }
    }

    /// Returns the mantissa of the BigFloat in `bytes`. Each byte represents a decimal digit.
    /// The first byte is the most significant. `bytes` can be of any length. 
    /// If the length of `bytes` is less than the number of decimal positions filled in the mantissa, 
    /// then the rest of the mantissa will be omitted.
    ///
    /// The length of the mantissa can be determined using `get_mantissa_len`.
    /// If `self` is Inf or NaN, nothing is returned.
    /// 
    /// ## Examples
    /// 
    /// ```
    /// # use num_bigfloat::BigFloat;
    /// let n = BigFloat::from_f64(123.42);
    /// let mut m = [0; 40];
    /// n.get_mantissa_bytes(&mut m);
    /// // compare m[0..10] to [1,2,3,4,2,0,0,0,0,0]
    /// assert!(m[0..10].iter().zip([1,2,3,4,2,0,0,0,0,0].iter()).filter(|x| { x.0 != x.1 }).count() == 0);
    /// ```
    pub fn get_mantissa_bytes(&self, bytes: &mut [u8]) {
        if let Flavor::Value(v) = self.inner {
            v.get_mantissa_bytes(bytes);
        }
    }

    /// Returns the number of decimal places filled in the mantissa.
    /// If `self` is Inf or NaN, 0 is returned.
    pub fn get_mantissa_len(&self) -> usize {
        match self.inner {
            Flavor::Value(v) => v.get_mantissa_len(),
            _ => 0,
        }
    }

    /// Returns 1 if BigFloat is positive, -1 otherwise.
    /// If `self` is NaN, 0 is returned.
    pub fn get_sign(&self) -> i8 {
        match self.inner {
            Flavor::Value(v) => v.sign,
            Flavor::Inf(s) => s,
            _ => 0,
        }
    }

    /// Returns the exponent part of the number.
    /// If `self` is Inf or NaN, 0 is returned.
    pub fn get_exponent(&self) -> i8 {
        match self.inner {
            Flavor::Value(v) => v.e,
            _ => 0,
        }
    }


    /// Sets the exponent part of the number.
    /// The function has no effect on Inf and NaN values.
    pub fn set_exponent(&mut self, e: i8) {
        if let Flavor::Value(mut v) = self.inner { 
            v.e = e;
            self.inner = Flavor::Value(v);
        }
    }

    /// Returns the raw parts of the number: the mantissa, the number of decimal places in the mantissa, 
    /// the sign, and the exponent.
    /// If `self` is Inf or NaN, None is returned.
    pub fn to_raw_parts(&self) -> Option<([i16; DECIMAL_PARTS], i16, i8, i8)> {
        match self.inner {
            Flavor::Value(v) => Some((v.m, v.n, v.sign, v.e)),
            _ => None,
        }
    }

    /// Creates a BigFloat from the raw parts. `to_raw_parts` can be used to get the raw parts of a number.
    pub fn from_raw_parts(mantissa: [i16; DECIMAL_PARTS], mantissa_len: i16, sign: i8, exponent: i8) -> Self {
        let val = BigFloatNum {
            sign,
            e: exponent,
            n: mantissa_len,
            m: mantissa,
        };
        BigFloat { 
            inner: Flavor::Value(val) 
        }
    }

    /// Returns true if `self` is positive infinity.
    pub fn is_inf_pos(&self) -> bool {
        matches!(self.inner, Flavor::Inf(DECIMAL_SIGN_POS))
    }

    /// Returns true if `self` is negative infinity.
    pub fn is_inf_neg(&self) -> bool {
        matches!(self.inner, Flavor::Inf(DECIMAL_SIGN_NEG))
    }

    /// Returns true if `self` is infinite.
    pub fn is_inf(&self) -> bool {
        self.is_inf_pos() || self.is_inf_neg()
    }

    /// Return true if `self` is not a number.
    pub fn is_nan(&self) -> bool {
        matches!(self.inner, Flavor::NaN)
    }

    /// Adds `d2` to `self` and returns the result of the addition.
    pub fn add(&self, d2: &Self) -> Self {
        match self.inner {
            Flavor::Value(v1) => {
                match d2.inner {
                    Flavor::Value(v2) => {
                        Self::result_to_ext(v1.add(&v2), v1.is_zero(), v1.sign == v2.sign)
                    },
                    Flavor::Inf(s2) => {
                        BigFloat { inner: Flavor::Inf(s2) }
                    },
                    Flavor::NaN => NAN,
                }
            },
            Flavor::Inf(s1) => {
                match d2.inner {
                    Flavor::Value(_) => {
                        BigFloat { inner: Flavor::Inf(s1) }
                    },
                    Flavor::Inf(s2) => {
                        if s1 != s2 {
                            NAN
                        } else {
                            BigFloat { inner: Flavor::Inf(s2) }
                        }
                    },
                    Flavor::NaN => NAN,
                }
            },
            Flavor::NaN => NAN,
        }
    }

    /// Subtracts `d2` from `self` and return the result of the subtraction.
    pub fn sub(&self, d2: &Self) -> Self {
        match self.inner {
            Flavor::Value(v1) => {
                match d2.inner {
                    Flavor::Value(v2) => {
                        Self::result_to_ext(v1.sub(&v2), v1.is_zero(), v1.sign == v2.sign)
                    },
                    Flavor::Inf(s2) => {
                        if s2 == DECIMAL_SIGN_POS {
                            INF_NEG
                        } else {
                            INF_POS
                        }
                    },
                    Flavor::NaN => NAN,
                }
            },
            Flavor::Inf(s1) => {
                match d2.inner {
                    Flavor::Value(_) => {
                        BigFloat { inner: Flavor::Inf(s1) }
                    },
                    Flavor::Inf(s2) => {
                        if s1 == s2 {
                            NAN
                        } else {
                            BigFloat { inner: Flavor::Inf(s1) }
                        }
                    },
                    Flavor::NaN => NAN,
                }
            },
            Flavor::NaN => NAN,
        }
    }

    /// Multiplies `self` by `d2` and returns the result of the multiplication.
    pub fn mul(&self, d2: &Self) -> Self {
        match self.inner {
            Flavor::Value(v1) => {
                match d2.inner {
                    Flavor::Value(v2) => {
                        Self::result_to_ext(v1.mul(&v2), v1.is_zero(), v1.sign == v2.sign)
                    },
                    Flavor::Inf(s2) => {
                        if v1.is_zero() { // 0*inf
                            NAN
                        } else {
                            let s = if v1.sign == s2 {
                                DECIMAL_SIGN_POS
                            } else {
                                DECIMAL_SIGN_NEG
                            };
                            BigFloat { inner: Flavor::Inf(s) }
                        }
                    },
                    Flavor::NaN => NAN,
                }
            },
            Flavor::Inf(s1) => {
                match d2.inner {
                    Flavor::Value(v2) => {
                        if v2.is_zero() { // inf*0
                            NAN
                        } else {
                            let s = if v2.sign == s1 {
                                DECIMAL_SIGN_POS
                            } else {
                                DECIMAL_SIGN_NEG
                            };
                            BigFloat { inner: Flavor::Inf(s) }
                        }
                    },
                    Flavor::Inf(s2) => {
                        let s = if s1 == s2 {
                            DECIMAL_SIGN_POS
                        } else {
                            DECIMAL_SIGN_NEG
                        };
                        BigFloat { inner: Flavor::Inf(s) }
                    },
                    Flavor::NaN => NAN,
                }
            },
            Flavor::NaN => NAN,
        }
    }

    /// Divides `self` by `d2` and returns the result of the division.
    pub fn div(&self, d2: &Self) -> Self {
        match self.inner {
            Flavor::Value(v1) => {
                match d2.inner {
                    Flavor::Value(v2) => {
                        Self::result_to_ext(v1.div(&v2), v1.is_zero(), v1.sign == v2.sign)
                    },
                    Flavor::Inf(_) => {
                        ZERO
                    },
                    Flavor::NaN => NAN,
                }
            },
            Flavor::Inf(s1) => {
                match d2.inner {
                    Flavor::Value(v) => {
                        if s1 == v.sign {
                            INF_POS
                        } else {
                            INF_NEG
                        }
                    },
                    Flavor::Inf(_) => {
                        NAN
                    },
                    Flavor::NaN => NAN,
                }
            },
            Flavor::NaN => NAN,
        }
    }

    /// Compares `self` to `d2`.
    /// Returns positive if `self` > `d2`, negative if `self` < `d2`, zero if `self` == `d2`, None if `self` or `d2` is NaN.
    pub fn cmp(&self, d2: &BigFloat) -> Option<i16> {
        match self.inner {
            Flavor::Value(v1) => {
                match d2.inner {
                    Flavor::Value(v2) => {
                        Some(v1.cmp(&v2))
                    }
                    Flavor::Inf(s2) => {
                        if s2 == DECIMAL_SIGN_POS {
                            Some(-1)
                        } else {
                            Some(1)
                        }
                    },
                    Flavor::NaN => None,
                }
            },
            Flavor::Inf(s1) => {
                match d2.inner {
                    Flavor::Value(_) => {
                        Some(s1 as i16)
                    }
                    Flavor::Inf(s2) => {
                        Some((s1 - s2) as i16)
                    },
                    Flavor::NaN => None,
                }
            },
            Flavor::NaN => None,
        }
    }

    /// Reverses the sign of a number.
    pub fn inv_sign(&self) -> BigFloat {
        match self.inner {
            Flavor::Value(mut v1) => {
                if v1.sign == DECIMAL_SIGN_POS {
                    v1.sign = DECIMAL_SIGN_NEG;
                } else {
                    v1.sign = DECIMAL_SIGN_POS;
                }
                BigFloat {inner: Flavor::Value(v1)}
            },
            Flavor::Inf(s1) => if s1 == DECIMAL_SIGN_POS {INF_NEG} else {INF_POS},
            Flavor::NaN => NAN,
        }
    }

    /// Returns `self` to the power of `d1`.
    pub fn pow(&self, d1: &Self) -> Self {
        match self.inner {
            Flavor::Value(v1) => {
                match d1.inner {
                    Flavor::Value(v2) => {
                        Self::result_to_ext(v1.pow(&v2), v1.is_zero(), v1.sign == v2.sign)
                    },
                    Flavor::Inf(s2) => {
                        // v1^inf
                        let val = v1.cmp(&BigFloatNum::one());
                        if val > 0 {
                            BigFloat { inner: Flavor::Inf(s2) }
                        } else if val < 0 {
                            ZERO
                        } else {
                            ONE
                        }
                    },
                    Flavor::NaN => NAN,
                }
            },
            Flavor::Inf(s1) => {
                match d1.inner {
                    Flavor::Value(v2) => {
                        // inf ^ v2
                        let val = v2.cmp(&BigFloatNum::new());
                        if val > 0 {
                            if s1 == DECIMAL_SIGN_NEG &&
                                v2.frac().is_zero() &&
                                !v2.is_int_even() {
                                    // v2 is odd without fractional part. 
                                    INF_NEG
                            } else {
                                INF_POS
                            }
                        } else if val < 0 {
                            ZERO
                        } else {
                            ONE
                        }
                    },
                    Flavor::Inf(s2) => {
                        // inf^inf
                        if s2 == DECIMAL_SIGN_POS {
                            INF_POS
                        } else {
                            ZERO
                        }
                    },
                    Flavor::NaN => NAN,
                }
            },
            Flavor::NaN => NAN,
        }
    }

    /// Returns the logarithm base `b` of a number.
    pub fn log(&self, b: &Self) -> Self {
        match self.inner {
            Flavor::Value(v1) => {
                match b.inner {
                    Flavor::Value(v2) => {
                        Self::result_to_ext(v1.log(&v2), false, true)
                    },
                    Flavor::Inf(s2) => {
                        // v1.log(inf)
                        if s2 == DECIMAL_SIGN_POS {
                            ZERO
                        } else {
                            NAN
                        }
                    },
                    Flavor::NaN => NAN,
                }
            },
            Flavor::Inf(s1) => {
                if s1 == DECIMAL_SIGN_NEG {
                    // -inf.log(any)
                    NAN
                } else {
                    match b.inner {
                        Flavor::Value(v2) => {
                            // +inf.log(v2)
                            let val = v2.cmp(&BigFloatNum::one());
                            if val < 0 {
                                INF_NEG
                            } else {
                                INF_POS
                            }
                        },
                        Flavor::Inf(_) => NAN, // +inf.log(inf)
                        Flavor::NaN => NAN,
                    }
                }
            },
            Flavor::NaN => NAN,
        }
    }

    fn result_to_ext(res: Result<BigFloatNum, Error>, is_dividend_zero: bool, is_same_sign: bool) -> BigFloat {
        match res {
            Err(e) => match e {
                Error::ExponentOverflow(s) => if s == DECIMAL_SIGN_POS { INF_POS } else { INF_NEG },
                Error::DivisionByZero => {
                    if is_dividend_zero {
                        NAN
                    } else if is_same_sign {
                        INF_POS
                    } else {
                        INF_NEG
                    }
                },
                Error::ArgumentIsNegative => NAN,
                Error::InvalidArgument => NAN,
            },
            Ok(v) => BigFloat {inner: Flavor::Value(v)},
        }
    }

    /// Returns true if `self` is positive. 
    /// The function returns false if `self` is NaN. 
    pub fn is_positive(&self) -> bool {
        match self.inner {
            Flavor::Value(v) => v.sign == DECIMAL_SIGN_POS,
            Flavor::Inf(s) => s == DECIMAL_SIGN_POS,
            Flavor::NaN => false,
        }
    }

    /// Returns true if `self` is negative.
    /// The function returns false if `self` is NaN. 
    pub fn is_negative(&self) -> bool {
        match self.inner {
            Flavor::Value(v) => v.sign == DECIMAL_SIGN_NEG,
            Flavor::Inf(s) => s == DECIMAL_SIGN_NEG,
            Flavor::NaN => false,
        }
    }

    /// Returns true if `self` is subnormal.
    /// A number is considered subnormal if not all digits of the mantissa are used, and the exponent has the minimum possible value.
    pub fn is_subnormal(&self) -> bool {
        if let Flavor::Value(v) = self.inner {
            return v.is_subnormal()
        }
        false
    }

    /// Returns true if `self` is zero.
    pub fn is_zero(&self) -> bool {
        match self.inner {
            Flavor::Value(v) => v.is_zero(),
            Flavor::Inf(_) => false,
            Flavor::NaN => false,
        }
    }

    /// Restricts the value of `self` to an interval determined by the values of `min` and `max`.
    /// The function returns `max` if `self` is greater than `max`, `min` if `self` is less than `min`, and `self` otherwise.
    /// If either argument is NaN or `min` is greater than `max`, the function returns NaN.
    pub fn clamp(&self, min: &Self, max: &Self) -> Self {
        if self.is_nan() || min.is_nan() || max.is_nan() || max.cmp(min).unwrap() < 0 {
            NAN
        } else if self.cmp(min).unwrap() < 0 {
            *min
        } else if self.cmp(max).unwrap() > 0 {
            *max
        } else {
            *self
        }
    }

    /// Returns the value of `d1` if `d1` is greater than `self`, or the value of `self` otherwise.
    /// If either argument is NaN, the function returns NaN.
    pub fn max(&self, d1: &Self) -> Self {
        if self.is_nan() || d1.is_nan() {
            NAN
        } else if self.cmp(d1).unwrap() < 0 {
            *d1
        } else {
            *self
        }
    }

    /// Returns value of `d1` if `d1` is less than `self`, or the value of `self` otherwise.
    /// If either argument is NaN, the function returns NaN.
    pub fn min(&self, d1: &Self) -> Self {
        if self.is_nan() || d1.is_nan() {
            NAN
        } else if self.cmp(d1).unwrap() > 0 {
            *d1
        } else {
            *self
        }
    }

    /// Returns a BigFloat with the value -1 if `self` is negative, 1 if `self` is positive, zero otherwise.
    /// The function returns NaN If `self` is NaN.
    pub fn signum(&self) -> Self {
        if self.is_nan() {
            NAN
        } else if self.is_negative() {
            ONE.inv_sign()
        } else {
            ONE
        }
    }


    /// Parses a number from the string `s`.
    /// The function expects `s` to be a number in scientific format in base 10, or +-Inf, or NaN.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use num_bigfloat::BigFloat;
    /// 
    /// let n = BigFloat::parse("0.0").unwrap();
    /// assert!(n.to_f64() == 0.0);
    /// let n = BigFloat::parse("1.123e-123").unwrap();
    /// assert!((n.to_f64() - 1.123e-123).abs() < f64::EPSILON);
    /// let n = BigFloat::parse("-Inf").unwrap();
    /// assert!(n.to_f64() == f64::NEG_INFINITY);
    /// let n = BigFloat::parse("NaN").unwrap();
    /// assert!(n.to_f64().is_nan());
    /// ```
    pub fn parse(s: &str) -> Option<Self> {
        let ps = crate::parser::parse(s);
        if ps.is_valid() {
            if ps.is_inf() {
                if ps.sign() == DECIMAL_SIGN_POS {
                    Some(INF_POS)
                } else {
                    Some(INF_NEG)
                }
            } else if ps.is_nan() {
                Some(NAN)
            } else {
                let (m, _n, s, e) = ps.raw_parts();
                let mut num = BigFloatNum::from_bytes(&m, s, 0);
                if num.n == 0 {
                    return Some(ZERO);
                }
                if e < DECIMAL_MIN_EXPONENT as i32 {
                    if e > DECIMAL_MIN_EXPONENT as i32 - num.n as i32 {
                        BigFloatNum::shift_right(&mut num.m, (DECIMAL_MIN_EXPONENT as i32 - e) as usize);
                        num.n = (num.n as i32 - DECIMAL_MIN_EXPONENT as i32 + e) as i16;
                        num.e = DECIMAL_MIN_EXPONENT;
                    } else {
                        return Some(ZERO);
                    }
                } else if e > DECIMAL_MAX_EXPONENT as i32 {
                    return Some(BigFloat {inner: Flavor::Inf(s)});
                } else {
                    num.e = e as i8;
                }
                Some(BigFloat {inner: Flavor::Value(num)})
            }
        } else {
            None
        }
    }

    pub(crate) fn write_str<T: Write>(&self, w: &mut T) -> Result<(), core::fmt::Error> {
        match self.inner {
            Flavor::Value(v) => {
                if v.is_zero() {
                    w.write_str("0.0")
                } else {
                    let part1 = if v.sign == DECIMAL_SIGN_NEG {
                        "-"
                    } else {
                        ""
                    };
                    let mut bytes = [0; DECIMAL_POSITIONS];
                    v.get_mantissa_bytes(&mut bytes);
                    let len = v.get_mantissa_len();
                    let mut part2 = [48u8; DECIMAL_POSITIONS + 1];
                    part2[0] = bytes[0] + 48;
                    part2[1] = 46;
                    for i in 1..len {
                        part2[i + 1] = bytes[i] + 48;
                    }
                    let part2: &str = core::str::from_utf8_mut(&mut part2).unwrap();
                    let mut buf = [0u8; 64];
                    let s = crate::util::concat_str(&mut buf, &[part1, part2]);
                    let e = v.n as i32 + v.e as i32 - 1;
                    if e != 0 {
                        let exp_sign = if e > 0 {
                            "+"
                        } else {
                            ""
                        };
                        w.write_fmt(format_args!("{}e{}{}", s, exp_sign, e))
                    } else {
                        w.write_str(s)
                    }
                }
            },
            Flavor::Inf(sign) => {
                let s = if sign == DECIMAL_SIGN_NEG {
                    "-Inf"
                } else {
                    "Inf"
                };
                w.write_str(s)
            },
            crate::ext::Flavor::NaN => {
                w.write_str("NaN")
            },
        }
    }


    /// Returns a random normalized (not subnormal) BigFloat number with exponent in the range
    /// from `exp_from` to `exp_to` inclusive. The sign can be positive and negative. Zero is excluded.
    /// Function does not follow any specific distribution law.
    /// The intended use of this function is for testing.
    /// 
    /// # Errors
    /// 
    /// InvalidArgument - when `exp_from` is greater than `exp_to`.
    #[cfg(feature = "rand")]
    pub fn random_normal(exp_from: i8, exp_to: i8) -> Result<Self, Error> {

        if exp_from > exp_to {
            return Err(Error::InvalidArgument);
        }

        // build mantissa
        let mut mantissa = [0i16; DECIMAL_PARTS];
        for v in mantissa.iter_mut() {
            *v = (random::<u16>() % crate::defs::DECIMAL_BASE as u16) as i16;
        }

        if mantissa[DECIMAL_PARTS - 1] == 0 {
            mantissa[DECIMAL_PARTS - 1] = (crate::defs::DECIMAL_BASE - 1) as i16;
        }

        while mantissa[DECIMAL_PARTS - 1] / 1000 == 0 {
            mantissa[DECIMAL_PARTS - 1] *= 10;
        }

        // sign & exponent
        let sign = if random::<i8>() & 1 == 0 { DECIMAL_SIGN_POS } else { DECIMAL_SIGN_NEG };
        let exp_range = exp_to as i32 - exp_from  as i32;
        let exp = (if exp_range != 0 { random::<i32>().abs() % exp_range } else { 0 } + exp_from as i32) as i8;

        Ok(BigFloat::from_raw_parts(mantissa, DECIMAL_POSITIONS as i16, sign, exp))
    }

    pub fn classify(&self) -> FpCategory {
        match self.inner {
            Flavor::Value(v) => {
                if v.is_subnormal() {
                    FpCategory::Subnormal
                } else if v.is_zero() {
                    FpCategory::Zero
                } else {
                    FpCategory::Normal
                }
            },
            Flavor::Inf(_) => FpCategory::Infinite,
            Flavor::NaN => FpCategory::Nan,
        }
    }

    /// Returns the remainder of division of `self` by `d1`.
    pub fn rem(&self, d1: &Self) -> Self {
        self.sub(&(self.div(&d1)).int().mul(&d1))
    }
}


macro_rules! gen_wrapper2 {
    // unwrap error, function requires self as argument
    ($comment:literal, $fname:ident, $ret:ty, $pos_inf:block, $neg_inf:block, $($arg:ident, $arg_type:ty),*) => {
        #[doc=$comment]
        pub fn $fname(&self$(,$arg: $arg_type)*) -> $ret {
            match self.inner {
                Flavor::Value(v) => Self::result_to_ext(v.$fname($($arg,)*), v.is_zero(), true),
                Flavor::Inf(s) => if s == DECIMAL_SIGN_POS $pos_inf else $neg_inf,
                Flavor::NaN => NAN,
            }
        }
    };
}

macro_rules! gen_wrapper4 {
    // function requires self as argument
    ($comment:literal, $fname:ident, $ret:ty, $pos_inf:block, $neg_inf:block, $($arg:ident, $arg_type:ty),*) => {
        #[doc=$comment]
        pub fn $fname(&self$(,$arg: $arg_type)*) -> $ret {
            let inner = match self.inner {
                Flavor::Value(v) => Flavor::Value(v.$fname($($arg,)*)),
                Flavor::Inf(s) => if s == DECIMAL_SIGN_POS $pos_inf else $neg_inf,
                Flavor::NaN => Flavor::NaN,
            };
            BigFloat {
                inner
            }
        }
    };
}

impl BigFloat {

    gen_wrapper4!("Returns the absolute value of `self`.", abs, Self, {Flavor::Inf(DECIMAL_SIGN_POS)}, {Flavor::Inf(DECIMAL_SIGN_POS)},);
    gen_wrapper4!("Returns the integer part of `self`.", int, Self, {Flavor::NaN}, {Flavor::NaN},);
    gen_wrapper4!("Returns the fractional part of `self`.", frac, Self, {Flavor::NaN}, {Flavor::NaN},);
    gen_wrapper2!("Returns the smallest integer greater than or equal to `self`.", ceil, Self, {INF_POS}, {INF_NEG},);
    gen_wrapper2!("Returns the largest integer less than or equal to `self`.", floor, Self, {INF_POS}, {INF_NEG},);
    gen_wrapper2!("Returns a rounded number with `n` decimal positions in the fractional part of the number using the rounding mode `rm`.", round, Self, {INF_POS}, {INF_NEG}, n, usize, rm, RoundingMode);

    gen_wrapper2!("Returns the square root of `self`.", sqrt, Self, {INF_POS}, {NAN},);
    gen_wrapper2!("Returns the cube root of `self`.", cbrt, Self, {INF_POS}, {INF_NEG},);
    gen_wrapper2!("Returns the natural logarithm of `self`.", ln, Self, {INF_POS}, {NAN},);
    gen_wrapper2!("Returns the logarithm base 2 of `self`.", log2, Self, {INF_POS}, {NAN},);
    gen_wrapper2!("Returns the logarithm base 10 of `self`.", log10, Self, {INF_POS}, {NAN},);
    gen_wrapper2!("Returns `e` to the power of `self`.", exp, Self, {INF_POS}, {INF_NEG},);

    gen_wrapper2!("Returns the sine of `self`. The function takes an angle in radians as an argument.", sin, Self, {NAN}, {NAN},);
    gen_wrapper2!("Returns the cosine of `self`. The function takes an angle in radians as an argument.", cos, Self, {NAN}, {NAN},);
    gen_wrapper2!("Returns the tangent of `self`. The function takes an angle in radians as an argument.", tan, Self, {NAN}, {NAN},);
    gen_wrapper2!("Returns the arcsine of `self`. The result is an angle in radians ranging from -pi/2 to pi/2.", asin, Self, {NAN}, {NAN},);
    gen_wrapper2!("Returns the arccosine of `self`. The result is an angle in radians ranging from 0 to pi.", acos, Self, {NAN}, {NAN},);
    gen_wrapper2!("Returns the arctangent of `self`. The result is an angle in radians ranging from -pi/2 to pi/2.", atan, Self, {HALF_PI}, {HALF_PI.inv_sign()},);

    gen_wrapper2!("Returns the hyperbolic sine of `self`.", sinh, Self, {INF_POS}, {INF_NEG},);
    gen_wrapper2!("Returns the hyperbolic cosine of `self`.", cosh, Self, {INF_POS}, {INF_POS},);
    gen_wrapper2!("Returns the hyperbolic tangent of `self`.", tanh, Self, {ONE}, {ONE.inv_sign()},);
    gen_wrapper2!("Returns the inverse hyperbolic sine of `self`.", asinh, Self, {INF_POS}, {INF_NEG},);
    gen_wrapper2!("Returns the inverse hyperbolic cosine of `self`.", acosh, Self, {ZERO}, {ZERO},);
    gen_wrapper2!("Returns the inverse hyperbolic tangent of `self`.", atanh, Self, {ZERO}, {ZERO},);
}

/// Standard library features
pub mod ops {

    use crate::ONE;
    use crate::ZERO;
    use crate::NAN;
    use crate::BigFloat;

    #[cfg(feature = "std")]
    use std::{
        iter::Product,
        iter::Sum,
        ops::Add,
        ops::AddAssign,
        ops::Div,
        ops::DivAssign,
        ops::Mul,
        ops::MulAssign,
        ops::Neg,
        ops::Sub,
        ops::SubAssign,
        cmp::PartialEq,
        cmp::Eq,
        cmp::PartialOrd,
        cmp::Ordering,
        fmt::Display,
        fmt::Formatter,
        str::FromStr,
        ops::Rem
    };


    #[cfg(not(feature = "std"))]
    use core::{
        iter::Product,
        iter::Sum,
        ops::Add,
        ops::AddAssign,
        ops::Div,
        ops::DivAssign,
        ops::Mul,
        ops::MulAssign,
        ops::Neg,
        ops::Sub,
        ops::SubAssign,
        cmp::PartialEq,
        cmp::Eq,
        cmp::PartialOrd,
        cmp::Ordering,
        fmt::Display,
        fmt::Formatter,
        str::FromStr,
        ops::Rem
    };

    //
    // ops traits
    //

    impl Add for BigFloat {
        type Output = Self;
        fn add(self, rhs: Self) -> Self::Output {
            BigFloat::add(&self, &rhs)
        }
    }
    
    impl AddAssign for BigFloat {
        fn add_assign(&mut self, rhs: Self) {
            *self = BigFloat::add(self, &rhs)
        }
    }

    impl Div for BigFloat {
        type Output = Self;
        fn div(self, rhs: Self) -> Self::Output {
            BigFloat::div(&self, &rhs)
        }
    }
    
    impl DivAssign for BigFloat {
        fn div_assign(&mut self, rhs: Self) {
            *self = BigFloat::div(self, &rhs)
        }
    }

    impl Rem for BigFloat {
        type Output = Self;
        fn rem(self, rhs: Self) -> Self::Output {
            BigFloat::rem(&self, &rhs)
        }
    }
    
    impl Mul for BigFloat {
        type Output = Self;
        fn mul(self, rhs: Self) -> Self::Output {
            BigFloat::mul(&self, &rhs)
        }
    }
    
    impl MulAssign for BigFloat {
        fn mul_assign(&mut self, rhs: Self) {
            *self = BigFloat::mul(self, &rhs)
        }
    }

    impl Neg for BigFloat {
        type Output = Self;
        fn neg(self) -> Self::Output {
            self.inv_sign()
        }
    }

    impl Neg for &BigFloat {
        type Output = BigFloat;
        fn neg(self) -> Self::Output {
            (*self).inv_sign()
        }
    }

    impl Sub for BigFloat {
        type Output = Self;
        fn sub(self, rhs: Self) -> Self::Output {
            BigFloat::sub(&self, &rhs)
        }
    }
    
    impl SubAssign for BigFloat {
        fn sub_assign(&mut self, rhs: Self) {
            *self = BigFloat::sub(self, &rhs)
        }
    }

    impl Add<&BigFloat> for BigFloat {
        type Output = Self;
        fn add(self, rhs: &BigFloat) -> Self::Output {
            BigFloat::add(&self, rhs)
        }
    }
    
    impl AddAssign<&BigFloat> for BigFloat {
        fn add_assign(&mut self, rhs: &BigFloat) {
            *self = BigFloat::add(self, &rhs)
        }
    }
    
    impl Div<&BigFloat> for BigFloat {
        type Output = Self;
        fn div(self, rhs: &BigFloat) -> Self::Output {
            BigFloat::div(&self, &rhs)
        }
    }
    
    impl DivAssign<&BigFloat> for BigFloat {
        fn div_assign(&mut self, rhs: &BigFloat) {
            *self = BigFloat::div(self, &rhs)
        }
    }
    
    impl Mul<&BigFloat> for BigFloat {
        type Output = Self;
        fn mul(self, rhs: &BigFloat) -> Self::Output {
            BigFloat::mul(&self, &rhs)
        }
    }
    
    impl MulAssign<&BigFloat> for BigFloat {
        fn mul_assign(&mut self, rhs: &BigFloat) {
            *self = BigFloat::mul(self, &rhs)
        }
    }

    impl Sub<&BigFloat> for BigFloat {
        type Output = Self;
        fn sub(self, rhs: &BigFloat) -> Self::Output {
            BigFloat::sub(&self, &rhs)
        }
    }
    
    impl SubAssign<&BigFloat> for BigFloat {
        fn sub_assign(&mut self, rhs: &BigFloat) {
            *self = BigFloat::sub(self, &rhs)
        }
    }
    
    //
    // ordering traits
    //
    
    impl PartialEq for BigFloat {
        fn eq(&self, other: &Self) -> bool {
            let cmp_result = BigFloat::cmp(self, other);
            matches!(cmp_result, Some(0))
        }
    }
    
    impl Eq for BigFloat {}
    
    impl PartialOrd for BigFloat {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            let cmp_result = BigFloat::cmp(self, other);
            match cmp_result {
                Some(v) => {
                    if v > 0 {
                        Some(Ordering::Greater)
                    } else if v < 0 {
                        Some(Ordering::Less)
                    } else {
                        Some(Ordering::Equal)
                    }
                },
                None => None,
            }
        }
    }

    impl From<f64> for BigFloat {
        fn from(f: f64) -> Self {
            BigFloat::from_f64(f)
        }
    }

    impl From<f32> for BigFloat {
        fn from(f: f32) -> Self {
            BigFloat::from_f32(f)
        }
    }

    impl Display for BigFloat {

        #[cfg(feature="std")]
        fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
            self.write_str(f)
        }

        #[cfg(not(feature="std"))]
        fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), core::fmt::Error> {
            self.write_str(f)
        }
    }

    impl Default for BigFloat {
        fn default() -> BigFloat {
            BigFloat::new()
        }
    }

    impl FromStr for BigFloat {
        type Err = BigFloat;

        /// Returns parsed number or NAN in case of error.
        fn from_str(src: &str) -> Result<BigFloat, Self::Err> {
            BigFloat::parse(src).ok_or(NAN)
        }
    }

    impl Product for BigFloat {
        fn product<I: Iterator<Item = BigFloat>>(iter: I) -> Self {
            let mut acc = ONE;
            for v in iter {
                acc *= v;
            }
            acc
        }
    }

    impl Sum for BigFloat {
        fn sum<I: Iterator<Item = BigFloat>>(iter: I) -> Self {
            let mut acc = ZERO;
            for v in iter {
                acc += v;
            }
            acc
        }
    }


    impl<'a> Product<&'a BigFloat> for BigFloat {
        fn product<I: Iterator<Item = &'a BigFloat>>(iter: I) -> Self {
            let mut acc = ONE;
            for v in iter {
                acc *= v;
            }
            acc
        }
    }

    impl<'a> Sum<&'a BigFloat> for BigFloat {
        fn sum<I: Iterator<Item = &'a BigFloat>>(iter: I) -> Self {
            let mut acc = ZERO;
            for v in iter {
                acc += v;
            }
            acc
        }
    }
}

macro_rules! impl_int_conv {
    ($s:ty, $u:ty, $from_s:ident, $from_u:ident, $from_int:ident) => {
        impl BigFloat {

            /// Construct BigFloat from integer value.
            pub fn $from_s(i: $s) -> Self {
                let sign = if i < 0 {
                    DECIMAL_SIGN_NEG
                } else {
                    DECIMAL_SIGN_POS
                };
                Self::$from_int(i.abs() as $u, sign)
            }
        
            /// Construct BigFloat from integer value.
            pub fn $from_u(i: $u) -> Self {
                Self::$from_int(i, DECIMAL_SIGN_POS)
            }
        
            fn $from_int(mut v: $u, sign: i8) -> Self {
                let mut d = [0u8; DECIMAL_POSITIONS];
                let mut p = DECIMAL_POSITIONS;
                while v > 0 {
                    p -= 1;
                    d[p] = (v % 10) as u8;
                    v /= 10;
                }
                if p < DECIMAL_POSITIONS {
                    Self::from_bytes(&d[p..], sign, 0)
                } else {
                    Self::new()
                }
            }
        }

        #[cfg(feature = "std")]
        impl From<$s> for BigFloat {
            fn from(i: $s) -> Self {
                BigFloat::$from_s(i)
            }
        }

        #[cfg(feature = "std")]
        impl From<$u> for BigFloat {
            fn from(i: $u) -> Self {
                BigFloat::$from_u(i)
            }
        }



        #[cfg(not(feature = "std"))]
        impl core::convert::From<$s> for BigFloat {
            fn from(i: $s) -> Self {
                BigFloat::$from_s(i)
            }
        }

        #[cfg(not(feature = "std"))]
        impl core::convert::From<$u> for BigFloat {
            fn from(i: $u) -> Self {
                BigFloat::$from_u(i)
            }
        }
    };
}

impl_int_conv!(i8, u8, from_i8, from_u8, from_int_u8);
impl_int_conv!(i16, u16, from_i16, from_u16, from_int_u16);
impl_int_conv!(i32, u32, from_i32, from_u32, from_int_u32);
impl_int_conv!(i64, u64, from_i64, from_u64, from_int_u64);
impl_int_conv!(i128, u128, from_i128, from_u128, from_int_u128);


#[cfg(test)]
mod tests {

    use crate::defs::{DECIMAL_PARTS, RoundingMode, I64_MAX, I64_MIN, U64_MAX, I128_MAX, I128_MIN, U128_MAX};
    use crate::ext::Flavor;
    use super::*;

    #[cfg(feature = "std")]
    use std::str::FromStr;

    #[cfg(not(feature = "std"))]
    use core::str::FromStr;

    #[test]
    fn test_ext() {

        // Inf & NaN
        let d1 = ONE;
        assert!(!d1.is_inf());
        assert!(!d1.is_nan());
        assert!(!d1.is_inf_pos());
        assert!(!d1.is_inf_neg());
        assert!(d1.get_sign() > 0);

        let d1 = ONE.div(&BigFloat::new());
        assert!(d1.is_inf());
        assert!(!d1.is_nan());
        assert!(d1.is_inf_pos());
        assert!(!d1.is_inf_neg());
        assert!(d1.get_sign() > 0);

        let d1 = d1.inv_sign();
        assert!(d1.is_inf());
        assert!(!d1.is_nan());
        assert!(!d1.is_inf_pos());
        assert!(d1.is_inf_neg());
        assert!(d1.get_sign() < 0);

        let d1 = BigFloat::new().div(&BigFloat::new());
        assert!(!d1.is_inf());
        assert!(d1.is_nan());
        assert!(!d1.is_inf_pos());
        assert!(!d1.is_inf_neg());
        assert!(d1.get_sign() == 0);

        // conversions
        let d1 = ONE;
        assert!(d1.to_f64() == 1.0);
        assert!(d1.to_f32() == 1.0);
        assert!(d1.to_i64() == Some(1));
        assert!(d1.to_u64() == Some(1));
        assert!(d1.to_i128() == Some(1));
        assert!(d1.to_u128() == Some(1));
        let d1 = BigFloat::new().div(&BigFloat::new());
        assert!(d1.to_f64().is_nan());
        assert!(d1.to_f32().is_nan());
        assert!(d1.to_i64().is_none());
        assert!(d1.to_u64().is_none());
        assert!(d1.to_i128().is_none());
        assert!(d1.to_u128().is_none());
        let d1 = ONE.div(&BigFloat::new());
        assert!(d1.to_f64().is_infinite());
        assert!(d1.to_f32().is_infinite());
        assert!(d1.to_f64().is_sign_positive());
        assert!(d1.to_f32().is_sign_positive());
        assert!(d1.to_i64().is_none());
        assert!(d1.to_u64().is_none());
        assert!(d1.to_i128().is_none());
        assert!(d1.to_u128().is_none());
        let d1 = d1.inv_sign();
        assert!(d1.to_f64().is_sign_negative());
        assert!(d1.to_f32().is_sign_negative());
        assert!(d1.to_f64().is_infinite());
        assert!(d1.to_f32().is_infinite());
        assert!(d1.to_i64().is_none());
        assert!(d1.to_u64().is_none());
        assert!(d1.to_i128().is_none());
        assert!(d1.to_u128().is_none());

        let d1 = BigFloat { inner: Flavor::Value(I64_MAX) };
        assert!(d1.to_i64() == Some(i64::MAX));
        assert!(d1.add(&ONE).to_i64() == None);
        let d1 = BigFloat { inner: Flavor::Value(I64_MIN) };
        assert!(d1.to_i64() == Some(i64::MIN));
        assert!(d1.sub(&ONE).to_i64() == None);
        let d1 = BigFloat { inner: Flavor::Value(U64_MAX) };
        assert!(d1.to_u64() == Some(u64::MAX));
        assert!(d1.add(&ONE).to_u64() == None);
        let d1 = BigFloat { inner: Flavor::Value(I128_MAX) };
        assert!(d1.to_i128() == Some(i128::MAX));
        assert!(d1.add(&ONE).to_i128() == None);
        let d1 = BigFloat { inner: Flavor::Value(I128_MIN) };
        assert!(d1.to_i128() == Some(i128::MIN));
        assert!(d1.sub(&ONE).to_i128() == None);
        let d1 = BigFloat { inner: Flavor::Value(U128_MAX) };
        assert!(d1.to_u128() == Some(u128::MAX));
        assert!(d1.add(&ONE).to_u128() == None);

        for _ in 0..1000 {
            let i = rand::random::<i64>();
            let d1 = BigFloat::from_i64(i);
            assert!(d1.to_i64() == Some(i));

            let i = rand::random::<u64>();
            let d1 = BigFloat::from_u64(i);
            assert!(d1.to_u64() == Some(i));

            let i = rand::random::<i128>();
            let d1 = BigFloat::from_i128(i);
            assert!(d1.to_i128() == Some(i));

            let i = rand::random::<u128>();
            let d1 = BigFloat::from_u128(i);
            assert!(d1.to_u128() == Some(i));
        }

        let d1 = ONE;
        let mut bytes = [1; DECIMAL_PARTS];
        d1.get_mantissa_bytes(&mut bytes);
        assert!(bytes != [1; DECIMAL_PARTS]);
        assert!(d1.get_mantissa_len() != 0);
        let mut bytes = [1; DECIMAL_PARTS];
        let d1 = INF_POS;
        d1.get_mantissa_bytes(&mut bytes);
        assert!(d1.get_mantissa_len() == 0);
        assert!(bytes == [1; DECIMAL_PARTS]);
        let d1 = INF_NEG;
        d1.get_mantissa_bytes(&mut bytes);
        assert!(d1.get_mantissa_len() == 0);
        assert!(bytes == [1; DECIMAL_PARTS]);
        let d1 = NAN;
        d1.get_mantissa_bytes(&mut bytes);
        assert!(d1.get_mantissa_len() == 0);
        assert!(bytes == [1; DECIMAL_PARTS]);

        assert!(ONE.get_exponent() < 0);
        assert!(INF_POS.get_exponent() == 0);
        assert!(INF_NEG.get_exponent() == 0);
        assert!(NAN.get_exponent() == 0);
    
        assert!(ONE.to_raw_parts().is_some());
        assert!(INF_POS.to_raw_parts().is_none());
        assert!(INF_NEG.to_raw_parts().is_none());
        assert!(NAN.to_raw_parts().is_none());

        assert!(ONE.add(&ONE).cmp(&TWO) == Some(0));
        assert!(ONE.add(&INF_POS).is_inf_pos());
        assert!(INF_POS.add(&ONE).is_inf_pos());
        assert!(ONE.add(&INF_NEG).is_inf_neg());
        assert!(INF_NEG.add(&ONE).is_inf_neg());
        assert!(INF_POS.add(&INF_POS).is_inf_pos());
        assert!(INF_POS.add(&INF_NEG).is_nan());
        assert!(INF_NEG.add(&INF_NEG).is_inf_neg());
        assert!(INF_NEG.add(&INF_POS).is_nan());

        assert!(TWO.sub(&ONE).cmp(&ONE) == Some(0));
        assert!(ONE.sub(&INF_POS).is_inf_neg());
        assert!(INF_POS.sub(&ONE).is_inf_pos());
        assert!(ONE.sub(&INF_NEG).is_inf_pos());
        assert!(INF_NEG.sub(&ONE).is_inf_neg());
        assert!(INF_POS.sub(&INF_POS).is_nan());
        assert!(INF_POS.sub(&INF_NEG).is_inf_pos());
        assert!(INF_NEG.sub(&INF_NEG).is_nan());
        assert!(INF_NEG.sub(&INF_POS).is_inf_neg());

        assert!(TWO.mul(&ONE).cmp(&TWO) == Some(0));
        assert!(ONE.mul(&INF_POS).is_inf_pos());
        assert!(INF_POS.mul(&ONE).is_inf_pos());
        assert!(ONE.mul(&INF_NEG).is_inf_neg());
        assert!(INF_NEG.mul(&ONE).is_inf_neg());
        assert!(ONE.inv_sign().mul(&INF_POS).is_inf_neg());
        assert!(ONE.inv_sign().mul(&INF_NEG).is_inf_pos());
        assert!(INF_POS.mul(&ONE.inv_sign()).is_inf_neg());
        assert!(INF_NEG.mul(&ONE.inv_sign()).is_inf_pos());
        assert!(INF_POS.mul(&INF_POS).is_inf_pos());
        assert!(INF_POS.mul(&INF_NEG).is_inf_neg());
        assert!(INF_NEG.mul(&INF_NEG).is_inf_pos());
        assert!(INF_NEG.mul(&INF_POS).is_inf_neg());
        assert!(INF_POS.mul(&BigFloat::new()).is_nan());
        assert!(INF_NEG.mul(&BigFloat::new()).is_nan());
        assert!(BigFloat::new().mul(&INF_POS).is_nan());
        assert!(BigFloat::new().mul(&INF_NEG).is_nan());

        assert!(TWO.div(&TWO).cmp(&ONE) == Some(0));
        assert!(TWO.div(&INF_POS).is_zero());
        assert!(INF_POS.div(&TWO).is_inf_pos());
        assert!(TWO.div(&INF_NEG).is_zero());
        assert!(INF_NEG.div(&TWO).is_inf_neg());
        assert!(TWO.inv_sign().div(&INF_POS).is_zero());
        assert!(TWO.inv_sign().div(&INF_NEG).is_zero());
        assert!(INF_POS.div(&TWO.inv_sign()).is_inf_neg());
        assert!(INF_NEG.div(&TWO.inv_sign()).is_inf_pos());
        assert!(INF_POS.div(&INF_POS).is_nan());
        assert!(INF_POS.div(&INF_NEG).is_nan());
        assert!(INF_NEG.div(&INF_NEG).is_nan());
        assert!(INF_NEG.div(&INF_POS).is_nan());
        assert!(INF_POS.div(&BigFloat::new()).is_inf_pos());
        assert!(INF_NEG.div(&BigFloat::new()).is_inf_neg());
        assert!(BigFloat::new().div(&INF_POS).is_zero());
        assert!(BigFloat::new().div(&INF_NEG).is_zero());

        for op in [BigFloat::add, 
            BigFloat::sub, 
            BigFloat::mul, 
            BigFloat::div, ] {
            assert!(op(&NAN, &ONE).is_nan());
            assert!(op(&ONE, &NAN).is_nan());
            assert!(op(&NAN, &INF_POS).is_nan());
            assert!(op(&INF_POS, &NAN).is_nan());
            assert!(op(&NAN, &INF_NEG).is_nan());
            assert!(op(&INF_NEG, &NAN).is_nan());
            assert!(op(&NAN, &NAN).is_nan());
        }

        assert!(ONE.cmp(&ONE).unwrap() == 0);
        assert!(ONE.cmp(&INF_POS).unwrap() < 0);
        assert!(INF_POS.cmp(&ONE).unwrap() > 0);
        assert!(INF_POS.cmp(&INF_POS).unwrap() == 0);
        assert!(ONE.cmp(&INF_NEG).unwrap() > 0);
        assert!(INF_NEG.cmp(&ONE).unwrap() < 0);
        assert!(INF_NEG.cmp(&INF_NEG).unwrap() == 0);
        assert!(ONE.cmp(&NAN).is_none());
        assert!(NAN.cmp(&ONE).is_none());
        assert!(INF_POS.cmp(&NAN).is_none());
        assert!(NAN.cmp(&INF_POS).is_none());
        assert!(INF_NEG.cmp(&NAN).is_none());
        assert!(NAN.cmp(&INF_NEG).is_none());
        assert!(NAN.cmp(&NAN).is_none());

        assert!(ONE.is_positive());
        assert!(!ONE.is_negative());

        assert!(ONE.inv_sign().is_negative());
        assert!(!ONE.inv_sign().is_positive());
        assert!(!INF_POS.is_negative());
        assert!(INF_POS.is_positive());
        assert!(INF_NEG.is_negative());
        assert!(!INF_NEG.is_positive());
        assert!(!NAN.is_positive());
        assert!(!NAN.is_negative());


        assert!(ONE.pow(&ONE).cmp(&ONE) == Some(0));
        assert!(BigFloat::new().pow(&INF_POS).is_zero());
        assert!(BigFloat::new().pow(&INF_NEG).is_zero());
        assert!(ONE.pow(&INF_POS).cmp(&ONE) == Some(0));
        assert!(ONE.pow(&INF_NEG).cmp(&ONE) == Some(0));
        assert!(TWO.pow(&INF_POS).is_inf_pos());
        assert!(TWO.pow(&INF_NEG).is_inf_neg());
        assert!(INF_POS.pow(&ONE).is_inf_pos());
        assert!(INF_NEG.pow(&ONE).is_inf_neg());
        assert!(INF_NEG.pow(&TWO).is_inf_pos());
        assert!(INF_NEG.pow(&BigFloat::from_f64(10.2)).is_inf_pos());
        assert!(INF_NEG.pow(&BigFloat::from_f64(3.0)).is_inf_neg());
        assert!(INF_POS.pow(&ONE.inv_sign()).is_zero());
        assert!(INF_NEG.pow(&ONE.inv_sign()).is_zero());
        assert!(INF_POS.pow(&BigFloat::new()).cmp(&ONE) == Some(0));
        assert!(INF_NEG.pow(&BigFloat::new()).cmp(&ONE) == Some(0));
        assert!(INF_POS.pow(&INF_POS).is_inf_pos());
        assert!(INF_NEG.pow(&INF_POS).is_inf_pos());
        assert!(INF_POS.pow(&INF_NEG).is_zero());
        assert!(INF_NEG.pow(&INF_NEG).is_zero());

        let half = ONE.div(&TWO);
        assert!(TWO.log(&TWO).cmp(&ONE) == Some(0));
        assert!(TWO.log(&INF_POS).is_zero());
        assert!(TWO.log(&INF_NEG).is_nan());
        assert!(INF_POS.log(&TWO).is_inf_pos());
        assert!(INF_NEG.log(&TWO).is_nan());
        assert!(half.log(&half).cmp(&ONE) == Some(0));
        assert!(half.log(&INF_POS).is_zero());
        assert!(half.log(&INF_NEG).is_nan());
        assert!(INF_POS.log(&half).is_inf_neg());
        assert!(INF_NEG.log(&half).is_nan());
        assert!(INF_POS.log(&INF_POS).is_nan());
        assert!(INF_POS.log(&INF_NEG).is_nan());
        assert!(INF_NEG.log(&INF_POS).is_nan());
        assert!(INF_NEG.log(&INF_NEG).is_nan());
        assert!(TWO.log(&ONE).is_inf_pos());
        assert!(half.log(&ONE).is_inf_pos());
        assert!(ONE.log(&ONE).is_nan());

        assert!(BigFloat::from_f32(f32::NAN).is_nan());
        assert!(BigFloat::from_f32(f32::INFINITY).is_inf_pos());
        assert!(BigFloat::from_f32(f32::NEG_INFINITY).is_inf_neg());
        assert!(!BigFloat::from_f32(1.0).is_nan());
        assert!(BigFloat::from_f64(f64::NAN).is_nan());
        assert!(BigFloat::from_f64(f64::INFINITY).is_inf_pos());
        assert!(BigFloat::from_f64(f64::NEG_INFINITY).is_inf_neg());
        assert!(!BigFloat::from_f64(1.0).is_nan());

        assert!(ONE.pow(&NAN).is_nan());
        assert!(NAN.pow(&ONE).is_nan());
        assert!(INF_POS.pow(&NAN).is_nan());
        assert!(NAN.pow(&INF_POS).is_nan());
        assert!(INF_NEG.pow(&NAN).is_nan());
        assert!(NAN.pow(&INF_NEG).is_nan());
        assert!(NAN.pow(&NAN).is_nan());

        assert!(TWO.log(&NAN).is_nan());
        assert!(NAN.log(&TWO).is_nan());
        assert!(INF_POS.log(&NAN).is_nan());
        assert!(NAN.log(&INF_POS).is_nan());
        assert!(INF_NEG.log(&NAN).is_nan());
        assert!(NAN.log(&INF_NEG).is_nan());
        assert!(NAN.log(&NAN).is_nan());

        assert!(INF_NEG.abs().is_inf_pos());
        assert!(INF_POS.abs().is_inf_pos());
        assert!(NAN.abs().is_nan());

        assert!(INF_NEG.int().is_nan());
        assert!(INF_POS.int().is_nan());
        assert!(NAN.int().is_nan());

        assert!(INF_NEG.frac().is_nan());
        assert!(INF_POS.frac().is_nan());
        assert!(NAN.frac().is_nan());

        assert!(INF_NEG.ceil().is_inf_neg());
        assert!(INF_POS.ceil().is_inf_pos());
        assert!(NAN.ceil().is_nan());

        assert!(INF_NEG.floor().is_inf_neg());
        assert!(INF_POS.floor().is_inf_pos());
        assert!(NAN.floor().is_nan());

        for rm in [RoundingMode::Up, RoundingMode::Down, RoundingMode::ToZero, 
                RoundingMode::FromZero, RoundingMode::ToEven, RoundingMode::ToOdd] {
            assert!(INF_NEG.round(0, rm).is_inf_neg());
            assert!(INF_POS.round(0, rm).is_inf_pos());
            assert!(NAN.round(0, rm).is_nan());
        }

        assert!(INF_NEG.sqrt().is_nan());
        assert!(INF_POS.sqrt().is_inf_pos());
        assert!(NAN.sqrt().is_nan());

        assert!(INF_NEG.cbrt().is_inf_neg());
        assert!(INF_POS.cbrt().is_inf_pos());
        assert!(NAN.cbrt().is_nan());

        for op in [BigFloat::ln, 
            BigFloat::log2, 
            BigFloat::log10,] {
            assert!(op(&INF_NEG).is_nan());
            assert!(op(&INF_POS).is_inf_pos());
            assert!(op(&NAN).is_nan());
        }

        assert!(INF_NEG.exp().is_inf_neg());
        assert!(INF_POS.exp().is_inf_pos());
        assert!(NAN.exp().is_nan());

        assert!(INF_NEG.sin().is_nan());
        assert!(INF_POS.sin().is_nan());
        assert!(NAN.sin().is_nan());

        assert!(INF_NEG.cos().is_nan());
        assert!(INF_POS.cos().is_nan());
        assert!(NAN.cos().is_nan());

        assert!(INF_NEG.tan().is_nan());
        assert!(INF_POS.tan().is_nan());
        assert!(NAN.tan().is_nan());

        assert!(INF_NEG.asin().is_nan());
        assert!(INF_POS.asin().is_nan());
        assert!(NAN.asin().is_nan());

        assert!(INF_NEG.acos().is_nan());
        assert!(INF_POS.acos().is_nan());
        assert!(NAN.acos().is_nan());

        assert!(INF_NEG.atan().cmp(&HALF_PI.inv_sign()) == Some(0));
        assert!(INF_POS.atan().cmp(&HALF_PI) == Some(0));
        assert!(NAN.atan().is_nan());

        assert!(INF_NEG.sinh().is_inf_neg());
        assert!(INF_POS.sinh().is_inf_pos());
        assert!(NAN.sinh().is_nan());
        
        assert!(INF_NEG.cosh().is_inf_pos());
        assert!(INF_POS.cosh().is_inf_pos());
        assert!(NAN.cosh().is_nan());

        assert!(INF_NEG.tanh().cmp(&ONE.inv_sign()) == Some(0));
        assert!(INF_POS.tanh().cmp(&ONE) == Some(0));
        assert!(NAN.tanh().is_nan());
        
        assert!(INF_NEG.asinh().is_inf_neg());
        assert!(INF_POS.asinh().is_inf_pos());
        assert!(NAN.asinh().is_nan());
        
        assert!(INF_NEG.acosh().is_zero());
        assert!(INF_POS.acosh().is_zero());
        assert!(NAN.acosh().is_nan());
        
        assert!(INF_NEG.atanh().is_zero());
        assert!(INF_POS.atanh().is_zero());
        assert!(NAN.atanh().is_nan());
    }

    #[test]
    pub fn test_ops() {

        let d1 = ONE;
        let d2 = BigFloat::new();
        assert!(d1 - d2 == d1);
        assert!(d1 + d2 == d1);
        let mut d3 = BigFloat::new();
        d3 += d1;
        assert!(d1 == d3);
        d3 -= d1;
        assert!(d1 > d3);
        d3 = TWO;
        d3 *= TWO;
        assert!(d3 == TWO*TWO);
        d3 /= TWO;
        assert!(TWO == d3);
        assert!(ONE < d3);
        assert!(ONE == TWO/TWO);

        let d1 = -TWO;
        assert!(d1.is_negative());

        let d1 = ONE;
        let d2 = BigFloat::new();
        assert!(d1 - &d2 == d1);
        assert!(d1 + &d2 == d1);
        let mut d3 = BigFloat::new();
        d3 += &d1;
        assert!(d1 == d3);
        d3 -= &d1;
        assert!(d1 > d3);
        d3 = TWO;
        d3 *= &TWO;
        assert!(d3 == TWO*&TWO);
        d3 /= &TWO;
        assert!(TWO == d3);
        assert!(ONE < d3);
        assert!(ONE == TWO/&TWO);

        let d1 = -&TWO;
        assert!(d1.is_negative());

        let d1 = BigFloat::from_f64(0.0123456789);

        let mut buf = [0u8; 256];
        let wblen = fmt_to_str(&d1, &mut buf).len();
        let d1str = core::str::from_utf8(&buf[..wblen]).unwrap();
        assert_eq!(d1str, "1.234567890000000000000000000000000000000e-2");
        assert!(BigFloat::from_str(d1str).unwrap() == d1);

        let d1 = BigFloat::from_f64(-123.456789);
        let wblen = fmt_to_str(&d1, &mut buf).len();
        let d1str = core::str::from_utf8(&buf[..wblen]).unwrap();
        assert!(d1str == "-1.234567890000000000000000000000000000000e+2");
        assert!(BigFloat::from_str(d1str).unwrap() == d1);

        let wblen = fmt_to_str(&INF_POS, &mut buf).len();
        let d1str = core::str::from_utf8(&buf[..wblen]).unwrap();
        assert!(d1str == "Inf");

        let wblen = fmt_to_str(&INF_NEG, &mut buf).len();
        let d1str = core::str::from_utf8(&buf[..wblen]).unwrap();
        assert!(d1str == "-Inf");

        let wblen = fmt_to_str(&NAN, &mut buf).len();
        let d1str = core::str::from_utf8(&buf[..wblen]).unwrap();
        assert!(d1str == "NaN");

        assert!(BigFloat::from_str("abc").unwrap_err().is_nan());

        let arr = [TWO, ONE, TWO];
        assert!(arr.into_iter().product::<BigFloat>() == TWO * TWO);
        assert!(arr.into_iter().sum::<BigFloat>() == TWO + ONE + TWO);

        assert!(BigFloat::from_i8(-123) == BigFloat::parse("-1.23e+2").unwrap());
        assert!(BigFloat::from_u8(123) == BigFloat::parse("1.23e+2").unwrap());
        assert!(BigFloat::from_i16(-12312) == BigFloat::parse("-1.2312e+4").unwrap());
        assert!(BigFloat::from_u16(12312) == BigFloat::parse("1.2312e+4").unwrap());
        assert!(BigFloat::from_i32(-123456789) == BigFloat::parse("-1.23456789e+8").unwrap());
        assert!(BigFloat::from_u32(123456789) == BigFloat::parse("1.23456789e+8").unwrap());
        assert!(BigFloat::from_i64(-1234567890123456789) == BigFloat::parse("-1.234567890123456789e+18").unwrap());
        assert!(BigFloat::from_u64(1234567890123456789) == BigFloat::parse("1.234567890123456789e+18").unwrap());
        assert!(BigFloat::from_i128(-123456789012345678901234567890123456789) == BigFloat::parse("-1.23456789012345678901234567890123456789e+38").unwrap());
        assert!(BigFloat::from_u128(123456789012345678901234567890123456789) == BigFloat::parse("1.23456789012345678901234567890123456789e+38").unwrap());
    }

    fn fmt_to_str<'a>(f: &BigFloat, buf: &'a mut [u8]) -> WritableBuf<'a> {
        buf.fill(0);
        let mut strepr = WritableBuf::new(buf);
        write!(strepr, "{}", f).unwrap();
        strepr
    }
}


#[cfg(feature="serde")]
#[cfg(test)]
mod serde_tests {

    use super::*;

    #[test]
    fn test_serde() {

        let d1 = E;

        let json = serde_json::to_string(&d1).unwrap();

        assert_eq!("{\"inner\":{\"Value\":{\"sign\":1,\"e\":-39,\"n\":40,\"m\":[7757,6249,3526,7471,6028,2353,9045,2845,2818,2718]}}}", json);

        let json = "{
            \"inner\": {
                \"Value\": {
                    \"sign\": -1,
                    \"e\": -39,
                    \"n\": 40,
                    \"m\": [7757, 6249, 3526, 7471, 6028, 2353, 9045, 2845, 2818, 2718]
                }
            }
        }";

        let d1 = d1.inv_sign();
        let d2: BigFloat = serde_json::from_str(json).unwrap();

        assert!(d1.cmp(&d2).unwrap() == 0);
    }
}

#[cfg(feature="rand")]
#[cfg(test)]
mod rand_tests {

    use super::*;

    #[test]
    fn test_rand() {

        for _ in 0..1000 {

            let exp_from = rand::random::<i8>();
            let exp_shift = if DECIMAL_MAX_EXPONENT > exp_from {
                rand::random::<u8>() % (DECIMAL_MAX_EXPONENT as i16 - exp_from as i16) as u8
            } else {
                0
            };
            let exp_to = (exp_from as i16 + exp_shift as i16) as i8;

            let n = BigFloat::random_normal(exp_from, exp_to).unwrap();

            assert!(!n.is_subnormal());
            assert!(n.get_exponent() >= exp_from && n.get_exponent() <= exp_to);
        }
    }
}
