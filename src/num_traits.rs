//! `num-traits` implementation.

use crate::BigFloat;
use crate::INF_NEG;
use crate::INF_POS;
use crate::MAX;
use crate::MIN;
use crate::MIN_POSITIVE;
use crate::NAN;
use crate::ONE;
use crate::RoundingMode;
use crate::ZERO;
use crate::Error;
use num_traits::Num;
use num_traits::bounds::Bounded;
use num_traits::float::Float;
use num_traits::float::FloatConst;
use num_traits::cast::cast;
use num_traits::cast::AsPrimitive;
use num_traits::cast::FromPrimitive;
use num_traits::cast::NumCast;
use num_traits::cast::ToPrimitive;
use num_traits::identities::one;
use num_traits::identities::zero;
use num_traits::identities::One;
use num_traits::identities::Zero;
use num_traits::int::PrimInt;
use num_traits::ops::checked::CheckedAdd;
use num_traits::ops::checked::CheckedDiv;
use num_traits::ops::checked::CheckedMul;
use num_traits::ops::checked::CheckedNeg;
use num_traits::ops::checked::CheckedRem;
use num_traits::ops::checked::CheckedShl;
use num_traits::ops::checked::CheckedShr;
use num_traits::ops::checked::CheckedSub;
use num_traits::ops::euclid::CheckedEuclid;
use num_traits::ops::euclid::Euclid;
use num_traits::ops::inv::Inv;
use num_traits::ops::mul_add::MulAdd;
use num_traits::ops::mul_add::MulAddAssign;
use num_traits::ops::saturating::Saturating;
use num_traits::ops::saturating::SaturatingAdd;
use num_traits::ops::saturating::SaturatingMul;
use num_traits::ops::saturating::SaturatingSub;
use num_traits::ops::wrapping::WrappingAdd;
use num_traits::ops::wrapping::WrappingMul;
use num_traits::ops::wrapping::WrappingNeg;
use num_traits::ops::wrapping::WrappingShl;
use num_traits::ops::wrapping::WrappingShr;
use num_traits::ops::wrapping::WrappingSub;
use num_traits::pow::checked_pow;
use num_traits::pow::pow;
use num_traits::pow::Pow;
use num_traits::sign::abs;
use num_traits::sign::abs_sub;
use num_traits::sign::signum;
use num_traits::sign::Signed;
use num_traits::sign::Unsigned;


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

#[cfg(feature = "std")]
impl NumCast for BigFloat {

    fn from<T: ToPrimitive>(n: T) -> Option<Self> {
        todo!()
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

    fn neg_zero() -> Self {
        ZERO
    }

    fn min_value() -> Self {
        MIN
    }

    fn min_positive_value() -> Self {
        MIN_POSITIVE
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
        !BigFloat::is_inf(&self)
    }

    fn is_normal(self) -> bool {
        !self.is_subnormal()
    }

    fn classify(self) -> std::num::FpCategory {
        todo!()
    }

    fn floor(self) -> Self {
        BigFloat::floor(&self)
    }

    fn ceil(self) -> Self {
        BigFloat::ceil(&self)
    }

    fn round(self) -> Self {
        BigFloat::round(&self, 0, RoundingMode::ToEven)
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
        BigFloat::signum(&self)
    }

    fn is_sign_positive(self) -> bool {
        BigFloat::is_positive(&self)
    }

    fn is_sign_negative(self) -> bool {
        BigFloat::is_negative(&self)
    }

    fn mul_add(self, a: Self, b: Self) -> Self {
        self * a + b
    }

    fn recip(self) -> Self {
        ONE / self
    }

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
        todo!()
    }

    fn ln(self) -> Self {
        BigFloat::ln(&self)
    }

    fn log(self, base: Self) -> Self {
        todo!()
    }

    fn log2(self) -> Self {
        todo!()
    }

    fn log10(self) -> Self {
        todo!()
    }

    fn max(self, other: Self) -> Self {
        BigFloat::max(&self, &other)
    }

    fn min(self, other: Self) -> Self {
        BigFloat::min(&self, &other)
    }

    fn abs_sub(self, other: Self) -> Self {
        todo!()
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
        todo!()
    }

    fn sin_cos(self) -> (Self, Self) {
        todo!()
    }

    fn exp_m1(self) -> Self {
        todo!()
    }

    fn ln_1p(self) -> Self {
        todo!()
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

    fn integer_decode(self) -> (u64, i16, i8) {
        todo!()
    }
}