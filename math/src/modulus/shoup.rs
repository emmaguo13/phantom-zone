use crate::{
    decomposer::PrimeDecomposer,
    distribution::{DistributionSized},
    modulus::{Prime, Modulus, ModulusOps},
};

use rand::distributions::{Distribution, Uniform};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Shoup(u64, u64);

impl Shoup {
    #[inline(always)]
    pub fn new(v: u64, q: u64) -> Self {
        debug_assert!(v < q);
        let quotient = (((v as u128) << 64) / q as u128) as _;
        Self(v, quotient)
    }

    #[inline(always)]
    pub fn value(&self) -> u64 {
        self.0
    }

    #[inline(always)]
    pub fn quotient(&self) -> u64 {
        self.1
    }

    #[inline(always)]
    pub fn mul(&self, a: u64, q: u64) -> u64 {
        let t = ((self.quotient() as u128 * a as u128) >> 64) as _;
        (a.wrapping_mul(self.value())).wrapping_sub(q.wrapping_mul(t))
    }
}

impl ModulusOps for Prime {
    type ElemPrep = Shoup;
    type Decomposer = PrimeDecomposer;

    fn new(modulus: Modulus) -> Self {
        modulus.try_into().unwrap()
    }

    #[inline(always)]
    fn modulus(&self) -> Modulus {
        (*self).into()
    }

    #[inline(always)]
    fn uniform_distribution(
        &self,
    ) -> impl Distribution<Self::Elem> + DistributionSized<Self::Elem> {
        Uniform::new_inclusive(0, self.max())
    }

    #[inline(always)]
    fn zero(&self) -> Self::Elem {
        0
    }

    #[inline(always)]
    fn one(&self) -> Self::Elem {
        1
    }

    #[inline(always)]
    fn neg_one(&self) -> Self::Elem {
        self.max()
    }

    #[inline(always)]
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        debug_assert!(*a < self.q());
        if *a != 0 {
            self.q() - a
        } else {
            0
        }
    }

    #[inline(always)]
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        debug_assert!(*a < self.q());
        debug_assert!(*b < self.q());
        let mut c = a + b;
        self.reduce_once_assign(&mut c);
        c
    }

    #[inline(always)]
    fn sub(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        debug_assert!(*a < self.q());
        debug_assert!(*b < self.q());
        if a >= b {
            a - b
        } else {
            self.q() + a - b
        }
    }

    #[inline(always)]
    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        debug_assert!(*a < self.q() << 1);
        debug_assert!(*b < self.q() << 1);

        self.reduce_u128(*a as u128 * *b as u128)
    }

    fn inv(&self, a: &Self::Elem) -> Option<Self::Elem> {
        (*a != 0).then(|| self.pow(*a, self.q() - 2))
    }

    #[inline(always)]
    fn prepare(&self, a: &Self::Elem) -> Self::ElemPrep {
        Shoup::new(*a, self.q())
    }

    #[inline(always)]
    fn mul_prep(&self, a: &Self::Elem, b: &Self::ElemPrep) -> Self::Elem {
        b.mul(*a, self.q())
    }
}