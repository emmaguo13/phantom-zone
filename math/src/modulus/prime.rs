use crate::{
    distribution::{Sampler},
    modulus::{ElemFrom, ElemOps, ElemTo, Modulus, ModulusOps},
};
use primality_test::is_prime;
use core::ops::Deref;
use num_bigint_dig::BigUint;
use num_traits::ToPrimitive;
use prime_factorization::Factorization;

/// A `ModulusOps` implementation that supports small prime modulus (less than
/// `1 << 61`) .
///
/// It panics in [`ModulusOps::new`] if `modulus` is not in range `3..1 << 61`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Prime {
    q: u64,
    q_half: u64,
    log_q: usize,
    barrett_mu: u128,
    barrett_alpha: usize,
    montgomery: u64,
}

fn get_distinct_factors(q: u64) -> Vec<u64> {
    let factors = Factorization::run(q).prime_factor_repr();
    let mut distincts_factors: Vec<u64> = Vec::with_capacity(factors.len());
    for factor in factors.iter(){
        distincts_factors.push(factor.0)
    }
    distincts_factors
}

impl Prime {
    pub fn new(q: u64) -> Self {
        assert!(is_prime(q) && q > 2);
        Self::new_unchecked(q)
    }

    pub fn new_unchecked(q: u64) -> Self {
        assert!(q.next_power_of_two().ilog2() <= 61);
        let log_q = q.next_power_of_two().ilog2() as usize;
        let barrett_mu = (1u128 << (log_q * 2 + 3)) / (q as u128);
        let barrett_alpha = log_q + 3;
        let mut montgomery: u64 = 1;
        let mut q_pow = q;
        for _i in 0..63{
            montgomery = montgomery.wrapping_mul(q_pow);
            q_pow = q_pow.wrapping_mul(q_pow);
        }
        Self {
            q,
            q_half: q >> 1,
            log_q,
            barrett_mu,
            barrett_alpha,
            montgomery,
        }
    }

    #[inline(always)]
    pub fn q(&self) -> u64{
        self.q
    }

    #[inline(always)]
    pub fn bits(&self) -> usize {
        self.log_q
    }

    #[inline(always)]
    pub fn max(&self) -> u64 {
        self.q - 1
    }

    #[inline(always)]
    pub fn half(&self) -> u64 {
        self.q_half
    }

    #[inline(always)]
    pub fn as_f64(&self) -> f64 {
        self.q as _
    }

    pub fn pow(&self, b: u64, e: u64) -> u64 {
        BigUint::from(b)
            .modpow(&BigUint::from(e), &BigUint::from(self.q))
            .to_u64()
            .unwrap()
    }

    pub fn two_adic_generator(&self, two_adicity: usize) -> u64 {
        assert_eq!((self.q - 1) % (1 << two_adicity) as u64, 0);
        self.pow(self.multiplicative_generator(), (self.q - 1) >> two_adicity)
    }

    pub fn multiplicative_generator(&self) -> u64 {
        let order = self.q - 1;
        (1..order)
            .find(|g| self.pow(*g, order >> 1) == order)
            .unwrap()
    }

    #[inline(always)]
    pub fn center(&self, v: u64) -> i64 {
        if v >= self.q_half {
            -((self.q - v) as i64)
        } else {
            v as _
        }
    }

    #[inline(always)]
    pub(crate) fn reduce_i128(&self, c: i128) -> u64 {
        // c / (2^{n + \beta})
        // note: \beta is assumed to -2
        let tmp = c >> (self.log_q - 2);
        // k = ((c / (2^{n + \beta})) * \mu) / 2^{\alpha - (-2)}
        let k = (tmp * self.barrett_mu as i128) >> (self.barrett_alpha + 2);
        // c - k*p
        let tmp = k * (self.q as i128);

        let mut c = (c - tmp) as u64;
        self.reduce_once_assign(&mut c);
        c
    }

    #[inline(always)]
    pub(crate) fn reduce_u128(&self, c: u128) -> u64 {
        // c / (2^{n + \beta})
        // note: \beta is assumed to -2
        let tmp = c >> (self.log_q - 2);
        // k = ((c / (2^{n + \beta})) * \mu) / 2^{\alpha - (-2)}
        let k = (tmp * self.barrett_mu) >> (self.barrett_alpha + 2);
        // c - k*p
        let tmp = k * (self.q as u128);

        let mut c = (c - tmp) as u64;
        self.reduce_once_assign(&mut c);
        c
    }

    #[inline(always)]
    pub(crate) fn reduce_once_assign(&self, a: &mut u64) {
        if *a >= self.q {
            *a -= self.q
        }
    }
}

#[cfg(any(test, feature = "dev"))]
impl Prime {
    pub fn gen(bits: usize, two_adicity: usize) -> Self {
        Self::gen_iter(bits, two_adicity).next().unwrap()
    }

    pub fn gen_iter(bits: usize, two_adicity: usize) -> impl Iterator<Item = Self> {
        assert!(bits <= 61);
        assert!(bits > two_adicity);
        let min = 1u64 << (bits - two_adicity - 1);
        let max = if min.leading_zeros() == 0 {
            u64::MAX
        } else {
            min << 1
        };
        let candidates = (min..max).rev().map(move |hi| (hi << two_adicity) + 1);
        candidates
            .into_iter()
            .filter(|v| is_prime(*v))
            .map(Self::new)
    }
}

impl Deref for Prime {
    type Target = u64;

    #[inline(always)]
    fn deref(&self) -> &Self::Target {
        &self.q
    }
}

impl ElemOps for Prime {
    type Elem = u64;
}



impl ElemFrom<u64> for Prime {
    #[inline(always)]
    fn elem_from(&self, v: u64) -> Self::Elem {
        v % self.q
    }
}

impl ElemFrom<i64> for Prime {
    #[inline(always)]
    fn elem_from(&self, v: i64) -> Self::Elem {
        v.rem_euclid(self.q as _) as _
    }
}

impl ElemFrom<u32> for Prime {
    #[inline(always)]
    fn elem_from(&self, v: u32) -> Self::Elem {
        self.elem_from(v as u64)
    }
}

impl ElemFrom<i32> for Prime {
    #[inline(always)]
    fn elem_from(&self, v: i32) -> Self::Elem {
        self.elem_from(v as i64)
    }
}

impl ElemFrom<f64> for Prime {
    #[inline(always)]
    fn elem_from(&self, v: f64) -> Self::Elem {
        self.reduce_i128(v.round() as i128)
    }
}

impl ElemTo<u64> for Prime {
    #[inline(always)]
    fn elem_to(&self, v: Self::Elem) -> u64 {
        v
    }
}

impl ElemTo<i64> for Prime {
    #[inline(always)]
    fn elem_to(&self, v: Self::Elem) -> i64 {
        self.center(v)
    }
}

impl ElemTo<f64> for Prime {
    #[inline(always)]
    fn elem_to(&self, v: Self::Elem) -> f64 {
        self.center(v) as f64
    }
}

impl Sampler for Prime {}

impl From<Prime> for Modulus {
    fn from(value: Prime) -> Self {
        Self::Prime(value.q)
    }
}

impl TryFrom<Modulus> for Prime {
    type Error = String;

    fn try_from(value: Modulus) -> Result<Self, Self::Error> {
        match value {
            Modulus::Prime(prime) => {
                if (3..1 << 61).contains(&prime) {
                    Ok(Self::new(prime))
                } else {
                    Err(format!(
                        "unsupported prime `{prime}`, expected in range `3..1 << 61`"
                    ))
                }
            }
            _ => Err(format!(
                "invalid modulus `{value:?}`, expected `Modulus::Prime(prime)`"
            )),
        }
    }
}

