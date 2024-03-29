//! Addition and subtraction.

use crate::defs::BigFloatNum;
use crate::defs::Error;

impl BigFloatNum {
    /// Add d2 and return result of addition.
    ///
    /// # Errors
    ///
    /// ExponentOverflow - when result is too big.
    pub fn add(&self, d2: &Self) -> Result<Self, Error> {
        let n = Self::to_big_float_inc(self);
        let m = Self::to_big_float_inc(d2);
        let ret = n.add(&m)?;
        Self::from_big_float_inc(ret)
    }

    /// Subtract d2 and return result of subtraction.
    ///
    /// # Errors
    ///
    /// ExponentOverflow - when result is too big.
    pub fn sub(&self, d2: &Self) -> Result<Self, Error> {
        let n = Self::to_big_float_inc(self);
        let m = Self::to_big_float_inc(d2);
        let ret = n.sub(&m)?;
        Self::from_big_float_inc(ret)
    }
}

#[cfg(test)]
mod tests {

    use crate::*;
    use defs::*;

    #[test]
    fn test_add() {
        let mut d1 = BigFloatNum::new();
        let mut d2 = BigFloatNum::new();
        let mut d3: BigFloatNum;
        let mut ref_num = BigFloatNum::new();

        //
        // addition
        //

        d2.sign = DECIMAL_SIGN_POS;
        d1.sign = DECIMAL_SIGN_POS;
        for i in 0..DECIMAL_PARTS {
            d1.m[i] = 9999;
            d2.m[i] = 0;
        }
        d1.m[0] = 9990;
        d2.m[0] = 10;
        d1.n = DECIMAL_POSITIONS as i16;
        d2.n = 2;
        d2.e = DECIMAL_MAX_EXPONENT;
        d1.e = DECIMAL_MAX_EXPONENT;

        assert!(d1.add(&d2).unwrap_err() == Error::ExponentOverflow(DECIMAL_SIGN_POS));

        d2.sign = DECIMAL_SIGN_NEG;
        d1.sign = DECIMAL_SIGN_NEG;
        assert!(d1.add(&d2).unwrap_err() == Error::ExponentOverflow(DECIMAL_SIGN_NEG));

        d2.e = 0;
        d1.e = 0;
        ref_num.m.iter_mut().for_each(|x| *x = 0);
        ref_num.m[DECIMAL_PARTS - 1] = DECIMAL_BASE as i16 / 10;
        ref_num.sign = DECIMAL_SIGN_NEG;
        ref_num.n = DECIMAL_POSITIONS as i16;
        ref_num.e = 1;
        d3 = d1.add(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        d2.sign = DECIMAL_SIGN_POS;
        ref_num.sign = DECIMAL_SIGN_NEG;
        for i in 0..DECIMAL_PARTS {
            ref_num.m[i] = 9999;
        }
        ref_num.m[0] = 9980;
        ref_num.n = DECIMAL_POSITIONS as i16;
        ref_num.e = 0;
        d3 = d1.add(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        ref_num.sign = DECIMAL_SIGN_POS;
        d1.sign = DECIMAL_SIGN_POS;
        d2.sign = DECIMAL_SIGN_NEG;
        d3 = d1.add(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        d1 = BigFloatNum::new();
        d2 = BigFloatNum::new();
        ref_num = BigFloatNum::new();
        d1.m[1] = 9999;
        d2.m[1] = 9999;
        ref_num.m[1] = 9998;
        ref_num.m[2] = 1;
        d1.n = 8;
        d2.n = 8;
        ref_num.n = 9;
        d3 = d1.add(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        // 0 + d2
        d1.m[1] = 0;
        d1.n = 0;
        d2.e = 3;
        ref_num = d2;
        d3 = d1.add(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        // d1 + 0
        d2.m[1] = 0;
        d2.n = 0;
        d2.e = 0;
        d1.m[1] = 11;
        d1.n = 6;
        d1.e = 123;
        ref_num = d1;
        d3 = d1.add(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        // different exponents and precisions
        d1 = BigFloatNum::new();
        d2 = BigFloatNum::new();

        for i in 0..DECIMAL_PARTS {
            d1.m[i] = 9999;
        }
        d1.n = DECIMAL_POSITIONS as i16;
        d1.e = 5;
        d2.m[0] = 11;
        d2.n = 2;
        d2.e = 3;
        ref_num = d1;
        d3 = d1.add(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        d1.e = 4;
        d2.m[0] = 11;
        d2.n = 2;
        d2.e = 3;
        ref_num.m.iter_mut().for_each(|x| *x = 0);
        ref_num.e = 5;
        ref_num.n = DECIMAL_POSITIONS as i16;
        ref_num.m[DECIMAL_PARTS - 1] = DECIMAL_BASE as i16 / 10;
        d3 = d1.add(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        d1.e = 5; // 0000 0022 2..2 2211 0000e+5 + 1122 3346e+3
        for i in 0..DECIMAL_PARTS {
            ref_num.m[i] = 2222;
            d1.m[i] = 2222;
        }
        d1.m[0] = 0;
        d1.m[1] = 2211;
        d1.m[DECIMAL_PARTS - 1] = 0;
        d1.m[DECIMAL_PARTS - 2] = 22;
        d1.n = DECIMAL_POSITIONS as i16 - 6;
        d2.m[0] = 3346;
        d2.m[1] = 1122;
        d2.n = 8;
        d2.e = 3;
        ref_num.e = 3;
        ref_num.m[0] = 3346;
        ref_num.m[DECIMAL_PARTS - 1] = 0;
        ref_num.n = DECIMAL_POSITIONS as i16 - 4;
        d3 = d1.add(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        d1.m[DECIMAL_PARTS - 1] = 222; // 0222 2..2 2211 0000e+5 + 1122 3346e+3
        d1.m[DECIMAL_PARTS - 2] = 2222;
        d1.n = DECIMAL_POSITIONS as i16 - 1;
        ref_num.e = 4;
        ref_num.m[0] = 2335;
        ref_num.m[DECIMAL_PARTS - 1] = 2222;
        ref_num.n = DECIMAL_POSITIONS as i16;
        d3 = d1.add(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        //
        // subtracton
        //

        d2.sign = DECIMAL_SIGN_POS;
        d1.sign = DECIMAL_SIGN_NEG;
        for i in 0..DECIMAL_PARTS {
            d1.m[i] = 9999;
            d2.m[i] = 0;
        }
        d1.m[0] = 9990;
        d2.m[0] = 10;
        d1.n = DECIMAL_POSITIONS as i16;
        d2.n = 2;
        d1.e = DECIMAL_MAX_EXPONENT;
        d2.e = DECIMAL_MAX_EXPONENT;
        assert!(d1.sub(&d2).unwrap_err() == Error::ExponentOverflow(DECIMAL_SIGN_NEG));

        d2.sign = DECIMAL_SIGN_NEG;
        d1.sign = DECIMAL_SIGN_POS;
        assert!(d1.sub(&d2).unwrap_err() == Error::ExponentOverflow(DECIMAL_SIGN_POS));

        d1.e -= 1;
        ref_num.m.iter_mut().for_each(|x| *x = 0);
        ref_num.m[DECIMAL_PARTS - 1] = DECIMAL_BASE as i16 / 10;
        ref_num.m[0] = 9;
        ref_num.sign = DECIMAL_SIGN_POS;
        ref_num.n = DECIMAL_POSITIONS as i16;
        ref_num.e = DECIMAL_MAX_EXPONENT;
        d3 = d1.sub(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        d1.sign = DECIMAL_SIGN_POS;
        d2.sign = DECIMAL_SIGN_POS;
        d1.e += 1;
        ref_num.sign = DECIMAL_SIGN_POS;
        for i in 0..DECIMAL_PARTS {
            ref_num.m[i] = 9999;
        }
        ref_num.m[0] = 9980;
        ref_num.n = DECIMAL_POSITIONS as i16;
        d3 = d1.sub(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        ref_num.sign = DECIMAL_SIGN_NEG;
        d1.sign = DECIMAL_SIGN_NEG;
        d2.sign = DECIMAL_SIGN_NEG;
        d3 = d1.sub(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        d1 = BigFloatNum::new();
        d2 = BigFloatNum::new();
        ref_num = BigFloatNum::new();
        d1.m[1] = 9998;
        d1.m[2] = 1;
        d1.n = 9;
        d2.m[1] = 9999;
        d2.n = 8;
        ref_num.m[1] = 9999;
        ref_num.n = 8;
        d3 = d1.sub(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        // d1 - 0
        d2.m[1] = 0;
        d2.n = 0;
        d2.e = -5;
        d1.e = 5;
        ref_num = d1;
        d3 = d1.sub(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);

        // 0 - d2
        d1.m[1] = 0;
        d1.m[2] = 0;
        d1.n = 0;
        d2.m[3] = 345;
        d2.n = 15;
        ref_num = d2;
        ref_num.sign = -ref_num.sign;
        d3 = d1.sub(&d2).unwrap();
        assert!(d3.cmp(&ref_num) == 0);
    }
}
