//! Multiplication and division.

use crate::defs::Error;
use crate::defs::DECIMAL_BASE;
use crate::defs::DECIMAL_BASE_LOG10;
use crate::defs::DECIMAL_MAX_EXPONENT;
use crate::defs::DECIMAL_MIN_EXPONENT;
use crate::defs::DECIMAL_SIGN_NEG;
use crate::defs::DECIMAL_SIGN_POS;
use crate::inc::inc::BigFloatInc;
use crate::inc::inc::DECIMAL_PARTS;
use crate::inc::inc::DECIMAL_POSITIONS;
use crate::RoundingMode;

impl BigFloatInc {
    /// Multiply by d2 and return result of multiplication.
    ///
    /// # Errors
    ///
    /// ExponentOverflow - when result is too big.
    pub fn mul(&self, d2: &Self) -> Result<Self, Error> {
        let mut d3 = Self::new();
        let mut m: i32;
        let mut k: i32;
        let mut n: i32;
        let mut nd: i32;
        let mut e: i32 = self.e as i32 + d2.e as i32;
        let mut d1mi: i32;
        let mut m3 = [0i16; DECIMAL_PARTS * 2 + 1];

        if self.n == 0 || d2.n == 0 {
            return Ok(Self::new());
        }

        for (i, v1) in self.m.iter().enumerate() {
            d1mi = *v1 as i32;
            if d1mi == 0 {
                continue;
            }

            k = 0;
            for (m2j, m3ij) in d2.m.iter().zip(m3[i..].iter_mut()) {
                m = d1mi * (*m2j as i32) + *m3ij as i32 + k;

                *m3ij = (m % DECIMAL_BASE as i32) as i16;
                k = m / DECIMAL_BASE as i32;
            }

            m3[i + d2.m.len()] += k as i16;
        }

        n = Self::num_digits(&m3[DECIMAL_PARTS..]) as i32;
        if n > 0 {
            n += DECIMAL_POSITIONS as i32;
        } else {
            n = Self::num_digits(&m3) as i32;
        }

        // take care if result is not fitting in d3.m without truncating
        if n > DECIMAL_POSITIONS as i32 {
            // save as much digits as we can
            e += n - DECIMAL_POSITIONS as i32;
            nd = n % DECIMAL_BASE_LOG10 as i32;
            n = n / DECIMAL_BASE_LOG10 as i32 - DECIMAL_PARTS as i32 + 1;

            Self::shift_left(&mut m3[n as usize..], DECIMAL_BASE_LOG10 - nd as usize);

            k = 1;
            while nd > 0 {
                k *= 10;
                nd -= 1;
            }
            m3[n as usize] += m3[n as usize - 1] / k as i16;

            d3.n = DECIMAL_POSITIONS as i16;
        } else {
            d3.n = n as i16;
            n = 0;
        }

        // round n digits in the beginning before copying [n..n + DECIMAL_PARTS] to d3
        if Self::round_mantissa(
            &mut m3[0..n as usize + DECIMAL_PARTS],
            n as i16 * DECIMAL_BASE_LOG10 as i16,
            RoundingMode::ToEven,
            true,
        ) {
            e += 1;
        }

        for i in 0..d3.m.len() {
            d3.m[i] = m3[n as usize + i];
        }
        d3.sign = if self.sign == d2.sign || self.n == 0 {
            DECIMAL_SIGN_POS
        } else {
            DECIMAL_SIGN_NEG
        };

        if e < DECIMAL_MIN_EXPONENT as i32 {
            return Ok(d3.process_subnormal(e));
        }

        if e > DECIMAL_MAX_EXPONENT as i32 {
            return Err(Error::ExponentOverflow(d3.sign));
        }

        d3.e = e as i8;

        Ok(d3)
    }

    /// Divide by d2 and return result of division.
    ///
    /// # Errors
    ///
    /// ExponentOverflow - when result is too big.
    /// DivisionByZero - in case of d2 equal to zero.
    /// InvalidArgument - in case both d1 and self are zeroes.
    pub fn div(&self, d2: &Self) -> Result<Self, Error> {
        // Knuth's division
        let mut d3 = Self::new();
        let mut m3 = [0; DECIMAL_PARTS + 1];
        let mut d: i32;
        let mut c: i16;
        let mut i: i32;
        let mut j: i32;
        let n: i32;
        let m: i32;
        let e: i32;
        let mut qh: i32;
        let mut k: i32;
        let mut rh: i32;
        let mut p: i32;
        let mut buf = [0i16; DECIMAL_PARTS * 3 + 4];
        let v1: i32;
        let v2: i32;
        let n1: usize = 2 + DECIMAL_PARTS;
        let n2: usize = DECIMAL_PARTS * 2 + 3;
        let mut n1j: usize;

        if d2.n == 0 {
            if self.n == 0 {
                return Err(Error::InvalidArgument);
            } else {
                return Err(Error::DivisionByZero);
            }
        }
        n = (d2.n as i32 - 1) / DECIMAL_BASE_LOG10 as i32;

        if self.n == 0 {
            return Ok(d3); // d1 / d2 = 0
        }
        m = (self.n as i32 - 1) / DECIMAL_BASE_LOG10 as i32;

        i = DECIMAL_PARTS as i32;
        p = 1;

        if n == 0 {
            // division by single digit
            d = d2.m[0] as i32;
            rh = 0;
            j = m;
            if (self.m[j as usize] as i32) < d {
                rh = self.m[j as usize] as i32;
                j -= 1;
                p -= 1;
            }

            let mut m3i = m3.iter_mut().rev();
            let mut m1i = self.m[..(j + 1) as usize].iter().rev();
            for m3v in m3i.by_ref() {
                qh =
                    rh * DECIMAL_BASE as i32 + if j >= 0 { *m1i.next().unwrap() as i32 } else { 0 };
                rh = qh % d;
                *m3v = (qh / d) as i16;

                if rh == 0 && j <= 0 {
                    break;
                }

                j -= 1;
            }

            for m3v in m3i {
                *m3v = 0;
            }
        } else {
            // normalize: n1 = d1 * d, n2 = d2 * d
            d = DECIMAL_BASE as i32 / (d2.m[n as usize] as i32 + 1); // factor d: d * d2[most significant] is close to DECIMAL_BASE

            if d == 1 {
                buf[n1..(self.m.len() + n1)].copy_from_slice(&self.m[..]);
                buf[n2..(d2.m.len() + n2)].copy_from_slice(&d2.m[..]);
            } else {
                Self::mul_by_digit(&self.m, d, &mut buf[n1..]);
                Self::mul_by_digit(&d2.m, d, &mut buf[n2..]);
            }

            v1 = buf[n2 + n as usize] as i32;
            v2 = buf[n2 + n as usize - 1] as i32;

            j = m - n;
            let mut m3i = m3.iter_mut().rev();
            loop {
                n1j = (n1 as i32 + j) as usize;

                let b2 = buf[n1j + n as usize + 1];
                let b1 = buf[n1j + n as usize];
                let b0 = buf[n1j + n as usize - 1];

                qh = b2 as i32 * DECIMAL_BASE as i32 + b1 as i32;
                rh = qh % v1;
                qh /= v1;

                if qh >= DECIMAL_BASE as i32 || (qh * v2 > DECIMAL_BASE as i32 * rh + b0 as i32) {
                    qh -= 1;
                    rh += v1;
                    if rh < DECIMAL_BASE as i32
                        && (qh >= DECIMAL_BASE as i32
                            || (qh * v2 > DECIMAL_BASE as i32 * rh + b0 as i32))
                    {
                        qh -= 1;
                    }
                }

                // n1_j = n1_j - n2 * qh
                c = 0;
                k = 0;
                let (buf1, buf2) = buf.split_at_mut(n2);
                for (a, b) in buf2[..(n + 2) as usize]
                    .iter()
                    .zip(buf1[n1j..n1j + (n + 2) as usize].iter_mut())
                {
                    k = *a as i32 * qh + k / DECIMAL_BASE as i32;
                    let val = k % DECIMAL_BASE as i32 + c as i32;
                    if (*b as i32) < val {
                        *b += (DECIMAL_BASE as i32 - val) as i16;
                        c = 1;
                    } else {
                        *b -= val as i16;
                        c = 0;
                    }
                }

                if c > 0 {
                    // compensate
                    qh -= 1;
                    c = 0;
                    for (a, b) in buf2[..(n + 2) as usize]
                        .iter()
                        .zip(buf1[n1j..n1j + (n + 2) as usize].iter_mut())
                    {
                        *b += *a + c;
                        if *b >= DECIMAL_BASE as i16 {
                            *b -= DECIMAL_BASE as i16;
                            c = 1;
                        } else {
                            c = 0;
                        }
                    }
                    assert!(c > 0);
                }

                if i < DECIMAL_PARTS as i32 || qh > 0 {
                    *m3i.next().unwrap() = qh as i16;
                    i -= 1;
                } else {
                    p -= 1;
                }

                j -= 1;
                if i < 0 || n1 as i32 + j < 0 {
                    break;
                }
            }
        }

        let mut rnd_e = 0;
        if Self::round_mantissa(
            &mut m3[0..DECIMAL_PARTS + 1],
            DECIMAL_BASE_LOG10 as i16,
            RoundingMode::ToEven,
            true,
        ) {
            rnd_e = 1;
        }
        d3.m.copy_from_slice(&m3[1..DECIMAL_PARTS + 1]);

        // exponent
        j = 0;
        d = d3.m[DECIMAL_PARTS - 1] as i32;
        while d > 0 {
            d /= 10;
            j += 1;
        }

        e = self.e as i32
            - d2.e as i32
            - (DECIMAL_PARTS as i32 - m + n - p) * DECIMAL_BASE_LOG10 as i32
            + rnd_e;

        d3.n = DECIMAL_POSITIONS as i16 - DECIMAL_BASE_LOG10 as i16 + j as i16;
        d3.sign = if self.sign == d2.sign { DECIMAL_SIGN_POS } else { DECIMAL_SIGN_NEG };

        if e < DECIMAL_MIN_EXPONENT as i32 {
            return Ok(d3.process_subnormal(e));
        }

        if e > DECIMAL_MAX_EXPONENT as i32 {
            return Err(Error::ExponentOverflow(d3.sign));
        }

        d3.e = e as i8;

        Ok(d3)
    }
}
