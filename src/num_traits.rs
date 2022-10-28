//! `num-traits` implementation.

use crate::BigFloat;
use crate::INF_NEG;
use crate::INF_POS;
use crate::MAX;
use crate::MIN;
use crate::MIN_POSITIVE_NORMAL;
use crate::NAN;
use crate::ONE;
#[cfg(feature="std")] use crate::PI;
#[cfg(feature="std")] use crate::RoundingMode;
#[cfg(feature="std")] use crate::TWO;
use crate::ZERO;
#[cfg(not(feature="std"))] use crate::EPSILON;
use crate::Error;
#[cfg(not(feature="std"))] use crate::ext::RAD_TO_DEG_FACTOR;
use num_traits::Num;
use num_traits::bounds::Bounded;
#[cfg(feature="std")] use num_traits::float::Float;
#[cfg(not(feature="std"))] use num_traits::float::FloatCore;
use num_traits::float::FloatConst;
use num_traits::cast::FromPrimitive;
use num_traits::cast::NumCast;
use num_traits::cast::ToPrimitive;
use num_traits::identities::One;
use num_traits::identities::Zero;
use num_traits::ops::euclid::Euclid;
use num_traits::ops::inv::Inv;
use num_traits::ops::mul_add::MulAdd;
use num_traits::ops::mul_add::MulAddAssign;
use num_traits::pow::Pow;
use num_traits::sign::Signed;
use core::num::FpCategory;

impl Bounded for BigFloat {

    fn min_value() -> Self {
        MIN
    }

    fn max_value() -> Self {
        MAX
    }
}

impl Zero for BigFloat {

    fn zero() -> Self {
        ZERO
    }

    fn is_zero(&self) -> bool {
        self.is_zero()
    }
}

impl One for BigFloat {

    fn one() -> Self {
        ONE
    }
}


impl Num for BigFloat {

    type FromStrRadixErr = Error;

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        if radix != 10 {
            return Err(Error::InvalidArgument);
        }

        BigFloat::parse(str).ok_or(Error::InvalidArgument)
    }
}

impl ToPrimitive for BigFloat {

    fn to_i64(&self) -> Option<i64> {
        BigFloat::to_i64(self)
    }

    fn to_u64(&self) -> Option<u64> {
        BigFloat::to_u64(self)
    }

    fn to_i128(&self) -> Option<i128> {
        BigFloat::to_i128(self)
    }

    fn to_u128(&self) -> Option<u128> {
        BigFloat::to_u128(self)
    }

    fn to_f64(&self) -> Option<f64> {
        Some(BigFloat::to_f64(self))
    }
}

impl NumCast for BigFloat {

    fn from<T: ToPrimitive>(n: T) -> Option<Self> {
        n.to_f64().map(BigFloat::from_f64)
    }
}

#[cfg(feature="std")]
impl Float for BigFloat {

    fn nan() -> Self {
        NAN
    }

    fn infinity() -> Self {
        INF_POS
    }

    fn neg_infinity() -> Self {
        INF_NEG
    }

    /// This function is provided only for compatibility since `-0.0` is not implemented.
    fn neg_zero() -> Self {
        ZERO
    }

    fn min_value() -> Self {
        MIN
    }

    fn min_positive_value() -> Self {
        MIN_POSITIVE_NORMAL
    }

    fn max_value() -> Self {
        MAX
    }

    fn is_nan(self) -> bool {
        BigFloat::is_nan(&self)
    }

    fn is_infinite(self) -> bool {
        BigFloat::is_inf(&self)
    }

    fn is_finite(self) -> bool {
        !(BigFloat::is_inf(&self) || BigFloat::is_nan(&self))
    }

    fn is_normal(self) -> bool {
        !(self.is_subnormal() || self.is_inf() || self.is_nan() || self.is_zero())
    }

    fn classify(self) -> FpCategory {
        BigFloat::classify(&self)
    }

    fn floor(self) -> Self {
        BigFloat::floor(&self)
    }

    fn ceil(self) -> Self {
        BigFloat::ceil(&self)
    }

    fn round(self) -> Self {
        BigFloat::round(&self, 0, RoundingMode::FromZero)
    }

    fn trunc(self) -> Self {
        self.int()
    }

    fn fract(self) -> Self {
        self.frac()
    }

    fn abs(self) -> Self {
        BigFloat::abs(&self)
    }

    fn signum(self) -> Self {
        let ret = BigFloat::signum(&self);
        if ret.is_zero() {
            ONE
        } else {
            ret
        }
    }

    /// Note: BigFloat NaN has no sign.
    fn is_sign_positive(self) -> bool {
        BigFloat::is_positive(&self)
    }

    /// Note: BigFloat NaN has no sign.
    fn is_sign_negative(self) -> bool {
        BigFloat::is_negative(&self)
    }

    /// This function is provided only for compatibility. It is not faster than separate multiplication and addition.
    fn mul_add(self, a: Self, b: Self) -> Self {
        self * a + b
    }

    fn recip(self) -> Self {
        ONE / self
    }

    /// This function is provided only for compatibility. It is not faster than `powf`.
    fn powi(self, n: i32) -> Self {
        let p = BigFloat::from_i32(n);
        BigFloat::pow(&self, &p)
    }

    fn powf(self, n: Self) -> Self {
        BigFloat::pow(&self, &n)
    }

    fn sqrt(self) -> Self {
        BigFloat::sqrt(&self)
    }

    fn exp(self) -> Self {
        BigFloat::exp(&self)
    }

    fn exp2(self) -> Self {
        BigFloat::pow(&TWO, &self)
    }

    fn ln(self) -> Self {
        BigFloat::ln(&self)
    }

    fn log(self, base: Self) -> Self {
        BigFloat::log(&self, &base)
    }

    fn log2(self) -> Self {
        BigFloat::log2(&self)
    }

    fn log10(self) -> Self {
        BigFloat::log10(&self)
    }

    fn max(self, other: Self) -> Self {
        BigFloat::max(&self, &other)
    }

    fn min(self, other: Self) -> Self {
        BigFloat::min(&self, &other)
    }

    fn abs_sub(self, other: Self) -> Self {
        let ret = self.sub(&other);
        if ret.is_negative() {
            ZERO
        } else {
            ret
        }
    }

    fn cbrt(self) -> Self {
        BigFloat::cbrt(&self)
    }

    fn hypot(self, other: Self) -> Self {
        (self * self + other * other).sqrt()
    }

    fn sin(self) -> Self {
        BigFloat::sin(&self)
    }

    fn cos(self) -> Self {
        BigFloat::cos(&self)
    }

    fn tan(self) -> Self {
        BigFloat::tan(&self)
    }

    fn asin(self) -> Self {
        BigFloat::asin(&self)
    }

    fn acos(self) -> Self {
        BigFloat::acos(&self)
    }

    fn atan(self) -> Self {
        BigFloat::atan(&self)
    }

    fn atan2(self, other: Self) -> Self {
        if self.is_zero() && other.is_zero() {
            ZERO
        } else if self.is_positive() || self.is_zero() {
            BigFloat::atan(&other.div(&self))
        } else if other.is_positive() || other.is_zero() {
            BigFloat::atan(&other.div(&self)) + PI
        } else {
            BigFloat::atan(&other.div(&self)) - PI
        }
    }

    /// This function is provided only for compatibility.
    fn sin_cos(self) -> (Self, Self) {
        (BigFloat::sin(&self), BigFloat::cos(&self))
    }

    /// This function is provided only for compatibility.
    fn exp_m1(self) -> Self {
        self.exp().sub(&ONE)
    }

    /// This function is provided only for compatibility.
    fn ln_1p(self) -> Self {
        self.add(&ONE).ln()
    }

    fn sinh(self) -> Self {
        BigFloat::sinh(&self)
    }

    fn cosh(self) -> Self {
        BigFloat::cosh(&self)
    }

    fn tanh(self) -> Self {
        BigFloat::tanh(&self)
    }

    fn asinh(self) -> Self {
        BigFloat::asinh(&self)
    }

    fn acosh(self) -> Self {
        BigFloat::acosh(&self)
    }

    fn atanh(self) -> Self {
        BigFloat::atanh(&self)
    }

    /// This function converts BigFloat to f64 and decomposes it.
    fn integer_decode(self) -> (u64, i16, i8) {

        let f = self.to_f64();

        let bits: u64 = f.to_bits();

        let sign: i8 = if bits >> 63 == 0 { 1 } else { -1 };

        let mut exponent: i16 = ((bits >> 52) & 0x7ff) as i16;

        let mantissa = if exponent == 0 {
            (bits & 0xfffffffffffff) << 1
        } else {
            (bits & 0xfffffffffffff) | 0x10000000000000
        };

        exponent -= 1023 + 52;

        (mantissa, exponent, sign)
    }
}

impl FloatConst for BigFloat {
    fn E() -> Self {
        crate::E
    } 
 
    fn FRAC_1_PI() -> Self {
        crate::PI
    } 
 
    fn FRAC_1_SQRT_2() -> Self {
        crate::FRAC_1_SQRT_2
    } 
 
    fn FRAC_2_PI() -> Self {
        crate::FRAC_2_PI
    } 
 
    fn FRAC_2_SQRT_PI() -> Self {
        crate::FRAC_2_SQRT_PI
    } 
 
    fn FRAC_PI_2() -> Self {
        crate::HALF_PI
    } 
 
    fn FRAC_PI_3() -> Self {
        crate::FRAC_PI_3
    } 
 
    fn FRAC_PI_4() -> Self {
        crate::FRAC_PI_4
    }

    fn FRAC_PI_6() -> Self {
        crate::FRAC_PI_6
    } 
 
    fn FRAC_PI_8() -> Self {
        crate::FRAC_PI_8
    } 
 
    fn LN_10() -> Self {
        crate::LN_10
    } 
 
    fn LN_2() -> Self {
        crate::LN_2
    } 
 
    fn LOG10_E() -> Self {
        crate::LOG10_E
    }

    fn LOG2_E() -> Self {
        crate::LOG2_E
    } 
 
    fn PI() -> Self {
        crate::PI
    } 
 
    fn SQRT_2() -> Self {
        crate::SQRT_2
    }
}

impl FromPrimitive for BigFloat {
    fn from_i64(n: i64) -> Option<Self> {
        Some(BigFloat::from_i64(n))
    }

    fn from_u64(n: u64) -> Option<Self> {
        Some(BigFloat::from_u64(n))
    }

    fn from_i128(n: i128) -> Option<Self> {
        Some(BigFloat::from_i128(n))
    }

    fn from_u128(n: u128) -> Option<Self> {
        Some(BigFloat::from_u128(n))
    }

    fn from_f32(n: f32) -> Option<Self> {
        Some(BigFloat::from_f32(n))
    }

    fn from_f64(n: f64) -> Option<Self> {
        Some(BigFloat::from_f64(n))
    }
}

impl Inv for BigFloat {
    type Output = BigFloat;

    fn inv(self) -> Self::Output {
        ONE.div(&self)
    }
}

/// This trait is provided only for compatibility. It does not provide performance benefits.
impl MulAdd for BigFloat {
    type Output = BigFloat;

    fn mul_add(self, a: Self, b: Self) -> Self::Output {
        self.mul(&a).add(&b)
    }
}

/// This trait is provided only for compatibility. It does not provide performance benefits.
impl MulAddAssign for BigFloat {
    fn mul_add_assign(&mut self, a: Self, b: Self) {
        *self = self.mul(&a).add(&b)
    }
}

impl Pow<BigFloat> for BigFloat {
    type Output = BigFloat;

    fn pow(self, rhs: BigFloat) -> Self::Output {
        BigFloat::pow(&self, &rhs)
    }
}

impl Signed for BigFloat {

    /// Note: BigFloat NaN has no sign.
    fn abs(&self) -> Self {
        BigFloat::abs(self)
    }

    fn abs_sub(&self, other: &Self) -> Self {
        let ret = self.sub(other);
        if ret.is_negative() {
            ZERO
        } else {
            ret
        }
    }

    /// Note: BigFloat NaN has no sign.
    fn signum(&self) -> Self {
        let ret = BigFloat::signum(self);
        if ret.is_zero() {
            ONE
        } else {
            ret
        }
    }

    fn is_positive(&self) -> bool {
        self.is_positive() && !self.is_zero()
    }

    fn is_negative(&self) -> bool {
        self.is_negative() && !self.is_zero()
    }
}

impl Euclid for BigFloat {    
    fn div_euclid(&self, v: &BigFloat) -> BigFloat {
        let q = BigFloat::int(&self.div(&v));
        if BigFloat::rem(self, v).is_negative() {
            return if v.is_negative() { q.add(&ONE) } else { q.sub(&ONE) };
        }
        q
    }

    fn rem_euclid(&self, v: &BigFloat) -> BigFloat {
        let r = BigFloat::rem(self, v);
        if r.is_negative() {
            v.abs().add(&r)
        } else {
            r
        }
    }
}

#[cfg(not(feature="std"))]
impl FloatCore for BigFloat {
    fn infinity() -> Self {
        INF_POS
    }

    fn neg_infinity() -> Self {
        INF_NEG
    }

    fn nan() -> Self {
        NAN
    }

    fn neg_zero() -> Self {
        ZERO
    }

    fn min_value() -> Self {
        MIN
    }

    fn min_positive_value() -> Self {
        MIN_POSITIVE_NORMAL
    }

    fn epsilon() -> Self {
        EPSILON
    }

    fn max_value() -> Self {
        MAX
    }

    fn classify(self) -> FpCategory {
        BigFloat::classify(&self)
    }

    fn to_degrees(self) -> Self {
        self.mul(&RAD_TO_DEG_FACTOR)
    }

    fn to_radians(self) -> Self {
        self.div(&RAD_TO_DEG_FACTOR)
    }

    fn integer_decode(self) -> (u64, i16, i8) {
        let f = self.to_f64();

        let bits: u64 = f.to_bits();

        let sign: i8 = if bits >> 63 == 0 { 1 } else { -1 };

        let mut exponent: i16 = ((bits >> 52) & 0x7ff) as i16;

        let mantissa = if exponent == 0 {
            (bits & 0xfffffffffffff) << 1
        } else {
            (bits & 0xfffffffffffff) | 0x10000000000000
        };

        exponent -= 1023 + 52;

        (mantissa, exponent, sign)
    }
}