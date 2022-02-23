use num_traits::{One, Zero};
use std::ops::{Add, Div, Mul, Rem, Sub};
pub trait RingOperation: Sized + Add + Sub + Mul {}
impl<T> RingOperation for T where T: Sized + Add + Sub + Mul {}
pub trait RingOperationFrom:
    Sized
    + for<'x> From<<&'x Self as Add>::Output>
    + for<'x> From<<&'x Self as Sub>::Output>
    + for<'x> From<<&'x Self as Mul>::Output>
where
    for<'x> &'x Self: RingOperation,
{
}
impl<T> RingOperationFrom for T
where
    T: Sized
        + for<'x> From<<&'x T as Add>::Output>
        + for<'x> From<<&'x T as Sub>::Output>
        + for<'x> From<<&'x T as Mul>::Output>,
    for<'x> &'x T: RingOperation,
{
}
pub trait EuclideanRingOperation: RingOperation + Div + Rem {}
impl<T> EuclideanRingOperation for T where T: RingOperation + Div + Rem {}
pub trait EuclideanRingOperationFrom:
    RingOperationFrom
    + for<'x> From<<&'x Self as Div>::Output>
    + for<'x> From<<&'x Self as Rem>::Output>
where
    for<'x> &'x Self: EuclideanRingOperation,
{
}
impl<T> EuclideanRingOperationFrom for T
where
    T: RingOperationFrom
        + for<'x> From<<&'x T as Div>::Output>
        + for<'x> From<<&'x T as Rem>::Output>,
    for<'x> &'x T: EuclideanRingOperation,
{
}
/** Normarize ring element

`abs(a)` in $`\mathbb{Z}`$.
`a/lc(a)` in $`R[x]`$ (`lc(x)` is leading coefficent of x).
*/
pub trait RingNormalize {
    #[must_use]
    fn leading_unit(&self) -> Self;
    fn normalize_mut(&mut self);
    #[must_use]
    fn into_normalize(mut self) -> Self
    where
        Self: Sized,
    {
        self.normalize_mut();
        self
    }
    #[must_use]
    fn normalize(&self) -> Self
    where
        Self: Clone,
    {
        self.clone().into_normalize()
    }
    fn is_similar(&self, other: &Self) -> bool
    where
        Self: Clone + Eq,
    {
        self.normalize() == other.normalize()
    }
}

macro_rules! ring_normalize {
    ($t:ty) => {
        impl RingNormalize for $t {
            fn leading_unit(&self) -> Self {
                if self >= &Self::zero() {
                    Self::one()
                } else {
                    -Self::one()
                }
            }
            fn normalize_mut(&mut self) {
                *self = self.abs();
            }
        }
    };
}
macro_rules! ring_normalize_unsigned {
    ($t:ty) => {
        impl RingNormalize for $t {
            fn leading_unit(&self) -> Self {
                Self::one()
            }
            fn normalize_mut(&mut self) {}
        }
    };
}

ring_normalize!(i8);
ring_normalize!(i16);
ring_normalize!(i32);
ring_normalize!(i64);
ring_normalize!(i128);
ring_normalize!(isize);
ring_normalize_unsigned!(u8);
ring_normalize_unsigned!(u16);
ring_normalize_unsigned!(u32);
ring_normalize_unsigned!(u64);
ring_normalize_unsigned!(u128);
ring_normalize_unsigned!(usize);

#[cfg(feature = "num-bigint")]
impl RingNormalize for num_bigint::BigInt {
    fn leading_unit(&self) -> Self {
        if self.sign() == num_bigint::Sign::Minus {
            -Self::one()
        } else {
            Self::one()
        }
    }
    fn normalize_mut(&mut self) {
        if self.sign() == num_bigint::Sign::Minus {
            *self = -&*self
        }
    }
}
#[cfg(feature = "num-bigint")]
impl RingNormalize for num_bigint::BigUint {
    fn leading_unit(&self) -> Self {
        Self::one()
    }
    fn normalize_mut(&mut self) {}
}
#[cfg(feature = "rug")]
impl RingNormalize for rug::Integer {
    fn leading_unit(&self) -> Self {
        if self < &Self::zero() {
            -Self::one()
        } else {
            Self::one()
        }
    }
    fn normalize_mut(&mut self) {
        self.abs_mut()
    }
}
