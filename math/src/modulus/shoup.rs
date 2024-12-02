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