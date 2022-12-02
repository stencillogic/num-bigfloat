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

    fn from<T: ToPrimitive>(n: T) -> Option<BigFloat> {
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
        } else if other.is_negative() {
            if self.is_negative() {
                BigFloat::atan(&self.div(&other)) - PI
            } else {
                BigFloat::atan(&self.div(&other)) + PI
            }
        } else {
            BigFloat::atan(&self.div(&other))
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
        crate::FRAC_1_PI
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

#[cfg(test)]
mod tests {

    use core::num::FpCategory;
    use num_traits::{Num, Bounded, Zero, FloatConst, Euclid, MulAdd, MulAddAssign, Signed, Inv, Pow, ToPrimitive, NumCast, FromPrimitive};
    use crate::{BigFloat, MIN_POSITIVE_NORMAL, MIN_POSITIVE, ONE, E, EPSILON, PI, HALF_PI, TWO, INF_NEG, INF_POS, NAN, RoundingMode};

    #[cfg(feature = "std")]
    use num_traits::Float;

    #[cfg(not(feature = "std"))]
    use num_traits::float::FloatCore;
    

    #[test]
    fn test_num_traits() {

        // Num
        let s1 = "1.234";
        let d1 = BigFloat::from_str_radix(s1, 10).unwrap();
        let d2 = BigFloat::parse(s1).unwrap();
        assert!(d1.cmp(&d2) == Some(0));
        assert!(BigFloat::from_str_radix(s1, 123).is_err());

        // Float
        #[cfg(feature = "std")] {
            let nan = <BigFloat as Float>::nan();
            assert!(Float::is_nan(nan));
            assert!(!Float::is_nan(d1));
    
            let infinity = <BigFloat as Float>::infinity();
            assert!(Float::is_infinite(infinity));
            assert!(!Float::is_finite(infinity));
            assert!(infinity > <BigFloat as Bounded>::max_value());
    
            let neg_infinity = <BigFloat as Float>::neg_infinity();
            assert!(neg_infinity.is_infinite());
            assert!(!neg_infinity.is_finite());
            assert!(neg_infinity < <BigFloat as Bounded>::min_value());
    
            let zero = <BigFloat as Zero>::zero();
            let neg_zero = <BigFloat as Float>::neg_zero();
    
            assert_eq!(zero, neg_zero);
            assert_eq!(BigFloat::from_f32(7.0)/infinity, zero);
            assert_eq!(zero * BigFloat::from_f32(10.0), zero);
    
            assert_eq!(<BigFloat as Float>::min_value(), <BigFloat as Bounded>::min_value());
            assert_eq!(<BigFloat as Float>::min_positive_value(), MIN_POSITIVE_NORMAL);
            assert_eq!(<BigFloat as Float>::max_value(), <BigFloat as Bounded>::max_value());
    
            assert!(<BigFloat as Float>::min_value().is_normal());
            assert!(<BigFloat as Float>::max_value().is_normal());
    
            let subnormal = MIN_POSITIVE;
            assert!(!Float::is_normal(zero));
            assert!(!Float::is_normal(nan));
            assert!(!Float::is_normal(infinity));
            assert!(!Float::is_normal(subnormal));
    
            assert_eq!(Float::classify(d1), FpCategory::Normal);
            assert_eq!(Float::classify(infinity), FpCategory::Infinite);
            assert_eq!(Float::classify(nan), FpCategory::Nan);
            assert_eq!(Float::classify(zero), FpCategory::Zero);
            assert_eq!(Float::classify(subnormal), FpCategory::Subnormal);
    
            let d1 = BigFloat::parse("3.3").unwrap();
            let d2 = BigFloat::parse("3.0").unwrap();
            assert_eq!(Float::floor(d1), d2);
            assert_eq!(Float::floor(d2), d2);
    
            let d1 = BigFloat::parse("3.3").unwrap();
            let d2 = BigFloat::parse("4.0").unwrap();
            assert_eq!(Float::ceil(d1), d2);
            assert_eq!(Float::ceil(d2), d2);
    
            let d1 = BigFloat::parse("3.3").unwrap();
            let d2 = BigFloat::parse("3.0").unwrap();
            assert_eq!(Float::round(d1), d2);
    
            let d1 = BigFloat::parse("3.5").unwrap();
            let d2 = BigFloat::parse("4.0").unwrap();
            assert_eq!(Float::round(d1), d2);
    
            let d1 = BigFloat::parse("-3.3").unwrap();
            let d2 = BigFloat::parse("-3.0").unwrap();
            assert_eq!(Float::round(d1), d2);
    
            let d1 = BigFloat::parse("-3.5").unwrap();
            let d2 = BigFloat::parse("-4.0").unwrap();
            assert_eq!(Float::round(d1), d2);
    
            let d1 = BigFloat::parse("3.7").unwrap();
            let d2 = BigFloat::parse("3.0").unwrap();
            assert_eq!(Float::trunc(d1), d2);
    
            let d1 = BigFloat::parse("-3.7").unwrap();
            let d2 = BigFloat::parse("-3.0").unwrap();
            assert_eq!(Float::trunc(d1), d2);
    
            let d1 = BigFloat::parse("-11.234").unwrap();
            let d2 = BigFloat::parse("-0.234").unwrap();
            assert_eq!(Float::fract(d1), d2);
    
            let d1 = BigFloat::parse("-0.234").unwrap();
            let d2 = BigFloat::parse("0.234").unwrap();
            assert_eq!(Float::abs(d1), d2);
            assert_eq!(Float::abs(d2), d2);
    
            assert_eq!(Float::signum(d2), ONE);
            assert_eq!(Float::signum(d1), -ONE);
            assert_eq!(Float::signum(infinity), ONE);
            assert_eq!(Float::signum(neg_infinity), -ONE);
            assert!(Float::signum(nan).is_nan());
    
            assert!(Float::is_sign_positive(d2));
            assert!(Float::is_sign_positive(infinity));
            assert!(!Float::is_sign_positive(d1));
            assert!(!Float::is_sign_positive(neg_infinity));
            assert!(!Float::is_sign_positive(nan));
            assert!(Float::is_sign_negative(d1));
            assert!(Float::is_sign_negative(neg_infinity));
            assert!(!Float::is_sign_negative(d2));
            assert!(!Float::is_sign_negative(infinity));
            assert!(!Float::is_sign_negative(nan));
    
            let d1 = BigFloat::parse("-2.1").unwrap();
            let d2 = BigFloat::parse("3.34").unwrap();
            let d3 = BigFloat::parse("43.657").unwrap();
            assert_eq!(Float::mul_add(d1, d2, d3), d1 * d2 + d3);
    
            assert_eq!(Float::recip(d1), ONE / d1);
    
            let d1 = BigFloat::parse("3.0").unwrap();
            let d2 = BigFloat::parse("81.0").unwrap();
            let d3 = BigFloat::parse("4.0").unwrap();
            assert_eq!(Float::powi(d1, 4), d2);
            assert_eq!(Float::powf(d1, d3), d2);
    
            let d1 = BigFloat::parse("9.0").unwrap();
            let d2 = BigFloat::parse("3.0").unwrap();
            assert_eq!(Float::sqrt(d1), d2);
    
            let d1 = BigFloat::parse("3.0").unwrap();
            assert!(Float::exp(d1).sub(&E.powf(d1)).get_exponent() <= -38);
    
            let d1 = BigFloat::parse("3.0").unwrap();
            let d2 = BigFloat::parse("2.0").unwrap();
            assert!(Float::exp2(d1).sub(&d2.powf(d1)).abs() <= EPSILON);
    
            let d2 = BigFloat::parse("4.0").unwrap();
            assert_eq!(Float::ln(d1), BigFloat::ln(&d1));
            assert_eq!(Float::log10(d1), BigFloat::log10(&d1));
            assert_eq!(Float::log2(d1), BigFloat::log2(&d1));
            assert_eq!(Float::log(d1, d2), BigFloat::log(&d1, &d2));
    
            assert_eq!(Float::max(d1, d2), d2);
            assert_eq!(Float::min(d1, d2), d1);
    
            let d1 = -d1;
            let d2 = -d2;
            assert_eq!(Float::max(d1, d2), d1);
            assert_eq!(Float::min(d1, d2), d2);
    
            let d1 = BigFloat::parse("3.0").unwrap();
            let d2 = BigFloat::parse("4.0").unwrap();
            assert!(Float::abs_sub(d1, d2).is_zero());
            assert_eq!(Float::abs_sub(d2, d1), ONE);
    
            let d2 = BigFloat::parse("27.0").unwrap();
            assert!(Float::cbrt(d2).sub(&d1).abs() <= EPSILON);
    
            let d1 = BigFloat::parse("3.0").unwrap();
            let d2 = BigFloat::parse("4.0").unwrap();
            let d3 = BigFloat::parse("5.0").unwrap();
            assert!(Float::hypot(d2, d1).sub(&d3).abs() <= EPSILON);
    
    
            let d1 = BigFloat::parse("0.5").unwrap();
            assert_eq!(Float::sin(d1), BigFloat::sin(&d1));
            assert_eq!(Float::cos(d1), BigFloat::cos(&d1));
            assert_eq!(Float::sin_cos(d1), (BigFloat::sin(&d1), BigFloat::cos(&d1)));
            assert_eq!(Float::tan(d1), BigFloat::tan(&d1));
            assert_eq!(Float::asin(d1), BigFloat::asin(&d1));
            assert_eq!(Float::acos(d1), BigFloat::acos(&d1));
            assert_eq!(Float::atan(d1), BigFloat::atan(&d1));
    
            let d1 = BigFloat::parse("1.5").unwrap();
            assert_eq!(Float::sinh(d1), BigFloat::sinh(&d1));
            assert_eq!(Float::cosh(d1), BigFloat::cosh(&d1));
            assert_eq!(Float::tanh(d1), BigFloat::tanh(&d1));
            assert_eq!(Float::asinh(d1), BigFloat::asinh(&d1));
            assert_eq!(Float::acosh(d1), BigFloat::acosh(&d1));
            let d1 = BigFloat::parse("0.5").unwrap();
            assert_eq!(Float::atanh(d1), BigFloat::atanh(&d1));
    
            let d1 = BigFloat::parse("0.0").unwrap();
            let d2 = BigFloat::parse("0.0").unwrap();
            assert!(Float::atan2(d1, d2).is_zero());
    
            let d1 = BigFloat::parse("2.0").unwrap();
            let d2 = BigFloat::parse("3.0").unwrap();
            assert_eq!(Float::atan2(d1, d2), d1.div(&d2).atan());
    
            let d1 = BigFloat::parse("2.0").unwrap();
            let d2 = BigFloat::parse("-3.0").unwrap();
            assert_eq!(Float::atan2(d1, d2), d1.div(&d2).atan().add(&PI));
    
            let d1 = BigFloat::parse("-2.0").unwrap();
            let d2 = BigFloat::parse("-3.0").unwrap();
            assert_eq!(Float::atan2(d1, d2), d1.div(&d2).atan().sub(&PI));
    
            let d1 = BigFloat::parse("2.0").unwrap();
            let d2 = BigFloat::parse("0.0").unwrap();
            assert_eq!(Float::atan2(d1, d2), HALF_PI);
    
            let d1 = BigFloat::parse("-2.0").unwrap();
            let d2 = BigFloat::parse("0.0").unwrap();
            assert_eq!(Float::atan2(d1, d2), -HALF_PI);
    
            assert_eq!(Float::exp_m1(d1), BigFloat::exp(&d1).sub(&ONE));
            assert_eq!(Float::ln_1p(E), BigFloat::ln(&E.add(&ONE)));
    
            assert_eq!(Float::integer_decode(d1), Float::integer_decode(-2.0f64));
        }

        // FloatCore
        #[cfg(not(feature="std"))] {
            let nan = <BigFloat as FloatCore>::nan();
            assert!(FloatCore::is_nan(nan));
            assert!(!FloatCore::is_nan(d1));
    
            let infinity = <BigFloat as FloatCore>::infinity();
            assert!(FloatCore::is_infinite(infinity));
            assert!(!FloatCore::is_finite(infinity));
            assert!(infinity > <BigFloat as Bounded>::max_value());
    
            let neg_infinity = <BigFloat as FloatCore>::neg_infinity();
            assert!(neg_infinity.is_infinite());
            assert!(!neg_infinity.is_finite());
            assert!(neg_infinity < <BigFloat as Bounded>::min_value());
    
            let zero = <BigFloat as Zero>::zero();
            let neg_zero = <BigFloat as FloatCore>::neg_zero();
    
            assert_eq!(zero, neg_zero);
            assert_eq!(BigFloat::from_f32(7.0)/infinity, zero);
            assert_eq!(zero * BigFloat::from_f32(10.0), zero);
    
            assert_eq!(<BigFloat as FloatCore>::min_value(), <BigFloat as Bounded>::min_value());
            assert_eq!(<BigFloat as FloatCore>::min_positive_value(), MIN_POSITIVE_NORMAL);
            assert_eq!(<BigFloat as FloatCore>::max_value(), <BigFloat as Bounded>::max_value());
    
            assert!(<BigFloat as FloatCore>::min_value().is_normal());
            assert!(<BigFloat as FloatCore>::max_value().is_normal());
    
            let subnormal = MIN_POSITIVE;
            assert!(!FloatCore::is_normal(zero));
            assert!(!FloatCore::is_normal(nan));
            assert!(!FloatCore::is_normal(infinity));
            assert!(!FloatCore::is_normal(subnormal));
    
            assert_eq!(FloatCore::classify(d1), FpCategory::Normal);
            assert_eq!(FloatCore::classify(infinity), FpCategory::Infinite);
            assert_eq!(FloatCore::classify(nan), FpCategory::Nan);
            assert_eq!(FloatCore::classify(zero), FpCategory::Zero);
            assert_eq!(FloatCore::classify(subnormal), FpCategory::Subnormal);
    
            let d1 = BigFloat::parse("3.3").unwrap();
            let d2 = BigFloat::parse("3.0").unwrap();
            assert_eq!(FloatCore::floor(d1), d2);
            assert_eq!(FloatCore::floor(d2), d2);
    
            let d1 = BigFloat::parse("3.3").unwrap();
            let d2 = BigFloat::parse("4.0").unwrap();
            assert_eq!(FloatCore::ceil(d1), d2);
            assert_eq!(FloatCore::ceil(d2), d2);
    
            let d1 = BigFloat::parse("3.3").unwrap();
            let d2 = BigFloat::parse("3.0").unwrap();
            assert_eq!(FloatCore::round(d1), d2);
    
            let d1 = BigFloat::parse("3.5").unwrap();
            let d2 = BigFloat::parse("4.0").unwrap();
            assert_eq!(FloatCore::round(d1), d2);
    
            let d1 = BigFloat::parse("-3.3").unwrap();
            let d2 = BigFloat::parse("-3.0").unwrap();
            assert_eq!(FloatCore::round(d1), d2);
    
            let d1 = BigFloat::parse("-3.5").unwrap();
            let d2 = BigFloat::parse("-4.0").unwrap();
            assert_eq!(FloatCore::round(d1), d2);
    
            let d1 = BigFloat::parse("3.7").unwrap();
            let d2 = BigFloat::parse("3.0").unwrap();
            assert_eq!(FloatCore::trunc(d1), d2);
    
            let d1 = BigFloat::parse("-3.7").unwrap();
            let d2 = BigFloat::parse("-3.0").unwrap();
            assert_eq!(FloatCore::trunc(d1), d2);
    
            let d1 = BigFloat::parse("-11.234").unwrap();
            let d2 = BigFloat::parse("-0.234").unwrap();
            assert_eq!(FloatCore::fract(d1), d2);
    
            let d1 = BigFloat::parse("-0.234").unwrap();
            let d2 = BigFloat::parse("0.234").unwrap();
            assert_eq!(FloatCore::abs(d1), d2);
            assert_eq!(FloatCore::abs(d2), d2);
    
            assert_eq!(FloatCore::signum(d2), ONE);
            assert_eq!(FloatCore::signum(d1), -ONE);
            assert_eq!(FloatCore::signum(infinity), ONE);
            assert_eq!(FloatCore::signum(neg_infinity), -ONE);
            assert!(FloatCore::signum(nan).is_nan());
    
            assert!(FloatCore::is_sign_positive(d2));
            assert!(FloatCore::is_sign_positive(infinity));
            assert!(!FloatCore::is_sign_positive(d1));
            assert!(!FloatCore::is_sign_positive(neg_infinity));
            assert!(FloatCore::is_sign_positive(nan));
            assert!(FloatCore::is_sign_negative(d1));
            assert!(FloatCore::is_sign_negative(neg_infinity));
            assert!(!FloatCore::is_sign_negative(d2));
            assert!(!FloatCore::is_sign_negative(infinity));
            assert!(!FloatCore::is_sign_negative(nan));

            let d1 = BigFloat::parse("90.0").unwrap();
            assert!(FloatCore::to_radians(d1).sub(&HALF_PI).abs() <= EPSILON);
            assert!(FloatCore::to_degrees(HALF_PI).sub(&d1).abs() <= EPSILON);

            let d1 = BigFloat::parse("-123.123").unwrap();
            assert_eq!(FloatCore::to_radians(d1).to_degrees(), d1);

            let d1 = BigFloat::parse("-2.0").unwrap();
            assert_eq!(FloatCore::integer_decode(d1), FloatCore::integer_decode(-2.0f64));
        }

        // FloatConst
        assert!(<BigFloat as FloatConst>::SQRT_2().mul(&FloatConst::SQRT_2()).sub(&TWO) <= EPSILON);
        assert!(<BigFloat as FloatConst>::FRAC_1_PI().mul(&PI).sub(&ONE) <= EPSILON);
        assert!(<BigFloat as FloatConst>::FRAC_1_SQRT_2().mul(&<BigFloat as FloatConst>::SQRT_2()).sub(&ONE) <= EPSILON);
        assert!(<BigFloat as FloatConst>::FRAC_2_PI().mul(&PI).sub(&TWO) <= EPSILON);
        assert!(<BigFloat as FloatConst>::FRAC_2_SQRT_PI().mul(&PI.sqrt()).sub(&TWO) <= EPSILON);
        assert!(<BigFloat as FloatConst>::FRAC_PI_2().mul(&TWO).sub(&PI) <= EPSILON);
        assert!(<BigFloat as FloatConst>::FRAC_PI_3().mul(&BigFloat::from_i8(3)).sub(&PI) <= EPSILON);
        assert!(<BigFloat as FloatConst>::FRAC_PI_4().mul(&BigFloat::from_i8(4)).sub(&PI) <= EPSILON);
        assert!(<BigFloat as FloatConst>::FRAC_PI_6().mul(&BigFloat::from_i8(6)).sub(&PI) <= EPSILON);
        assert!(<BigFloat as FloatConst>::FRAC_PI_8().mul(&BigFloat::from_i8(8)).sub(&PI) <= EPSILON);
        assert!(<BigFloat as FloatConst>::LN_10().sub(&BigFloat::from_i8(10).ln()) <= EPSILON);
        assert!(<BigFloat as FloatConst>::LN_2().sub(&TWO.ln()) <= EPSILON);
        assert!(<BigFloat as FloatConst>::LOG10_E().sub(&E.log10()) <= EPSILON);
        assert!(<BigFloat as FloatConst>::LOG2_E().sub(&E.log2()) <= EPSILON);

        // Euclid
        let a = BigFloat::from_i8(7);
        let b = BigFloat::from_i8(4);
        assert_eq!(Euclid::div_euclid(&a, &b), ONE); // 7 > 4 * 1
        assert_eq!(Euclid::div_euclid(&-a, &b), -TWO); // -7 >= 4 * -2
        assert_eq!(Euclid::div_euclid(&a, &-b), -ONE); // 7 >= -4 * -1
        assert_eq!(Euclid::div_euclid(&-a, &-b), TWO); // -7 >= -4 * 2

        let c = BigFloat::from_i8(3);
        assert_eq!(Euclid::rem_euclid(&a, &b), c);
        assert_eq!(Euclid::rem_euclid(&-a, &b), ONE);
        assert_eq!(Euclid::rem_euclid(&a, &-b), c);
        assert_eq!(Euclid::rem_euclid(&-a, &-b), ONE);

        // MulAdd
        let d1 = BigFloat::parse("-2.1").unwrap();
        let d2 = BigFloat::parse("3.34").unwrap();
        let d3 = BigFloat::parse("43.657").unwrap();
        assert_eq!(MulAdd::mul_add(d1, d2, d3), d1 * d2 + d3);

        // MulAssign
        let mut d1 = BigFloat::parse("-2.1").unwrap();
        let d2 = BigFloat::parse("3.34").unwrap();
        let d3 = BigFloat::parse("43.657").unwrap();
        let ret = d1 * d2 + d3;
        MulAddAssign::mul_add_assign(&mut d1, d2, d3);
        assert_eq!(d1, ret);

        // Inv
        let d1 = BigFloat::parse("-2.1").unwrap();
        assert_eq!(Inv::inv(d1), ONE / d1);

        // Signed
        let d1 = BigFloat::parse("-0.234").unwrap();
        let d2 = BigFloat::parse("0.234").unwrap();
        assert_eq!(Signed::abs(&d1), d2);
        assert_eq!(Signed::abs(&d2), d2);

        let d1 = BigFloat::parse("3.0").unwrap();
        let d2 = BigFloat::parse("4.0").unwrap();
        assert!(Signed::abs_sub(&d1, &d2).is_zero());
        assert_eq!(Signed::abs_sub(&d2, &d1), ONE);

        let d1 = -d1;
        assert_eq!(Signed::signum(&d2), ONE);
        assert_eq!(Signed::signum(&d1), -ONE);
        assert_eq!(Signed::signum(&INF_POS), ONE);
        assert_eq!(Signed::signum(&INF_NEG), -ONE);
        assert!(Signed::signum(&NAN).is_nan());

        assert!(Signed::is_positive(&BigFloat::from_i8(12)));
        assert!(Signed::is_positive(&INF_POS));
        assert!(!Signed::is_positive(&BigFloat::from_i8(0)));
        assert!(!Signed::is_positive(&INF_NEG));
        assert!(!Signed::is_positive(&NAN));

        assert!(Signed::is_negative(&BigFloat::from_i8(-12)));
        assert!(Signed::is_negative(&INF_NEG));
        assert!(!Signed::is_negative(&BigFloat::from_i8(0)));
        assert!(!Signed::is_negative(&INF_POS));
        assert!(!Signed::is_negative(&NAN));

        // Pow
        let d1 = BigFloat::parse("3.45").unwrap();
        let d2 = BigFloat::parse("-4.567").unwrap();
        assert_eq!(Pow::pow(d1, d2), BigFloat::pow(&d1, &d2));

        // ToPrimitive
        let d1 = BigFloat::parse("-356798765.45678").unwrap();
        assert_eq!(ToPrimitive::to_i64(&d1), Some(-356798765));
        assert_eq!(ToPrimitive::to_u64(&d1), Some(356798765));
        assert_eq!(ToPrimitive::to_i128(&d1), Some(-356798765));
        assert_eq!(ToPrimitive::to_u128(&d1), Some(356798765));
        assert_eq!(ToPrimitive::to_f32(&d1), Some(-356798765.45678));
        assert_eq!(ToPrimitive::to_f64(&d1), Some(-356798765.45678));

        // NumCast
        assert_eq!(NumCast::from(1i8), Some(ONE));
        assert_eq!(NumCast::from(1u8), Some(ONE));
        let d2: BigFloat = NumCast::from(-356798765.45678f64).unwrap();
        assert_eq!(d2, d1);

        let d1 = BigFloat::parse("5.45678").unwrap();
        let d2: BigFloat = NumCast::from(5.45678f32).unwrap();
        assert_eq!(BigFloat::round(&d2, 5, RoundingMode::ToEven), d1);

        // FromPrimitive
        let d1 = BigFloat::parse("-356798765").unwrap();
        assert_eq!(FromPrimitive::from_i64(-356798765), Some(d1));
        assert_eq!(FromPrimitive::from_u64(356798765), Some(-d1));
        assert_eq!(FromPrimitive::from_i128(-356798765), Some(d1));
        assert_eq!(FromPrimitive::from_u128(356798765), Some(-d1));
        let d1 = BigFloat::parse("-5.45678").unwrap();
        let d2 = FromPrimitive::from_f32(-5.45678f32).unwrap();
        let d2 = BigFloat::round(&d2, 5, RoundingMode::ToEven);
        assert_eq!(d2, d1);
        let d1 = BigFloat::parse("-356798765.45678").unwrap();
        assert_eq!(FromPrimitive::from_f64(-356798765.45678f64), Some(d1));
    }
}