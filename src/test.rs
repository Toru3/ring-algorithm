use crate::sealed;
use crate::*;
use is_prime_for_primitive_int::IsPrime;
use num_traits::{One, Zero};
use polynomial_ring::Polynomial;

impl<K> RingNormalize for Polynomial<K>
where
    K: Clone + Zero + One + for<'x> std::ops::AddAssign<&'x K> + for<'x> std::ops::DivAssign<&'x K>,
    for<'x> &'x K: std::ops::Mul<Output = K>,
{
    fn leading_unit(&self) -> Self {
        if let Some(x) = self.lc() {
            Self::from_monomial(x.clone(), 0)
        } else {
            Self::one()
        }
    }
    fn normalize_mut(&mut self) {
        self.monic();
    }
}

macro_rules! poly {
    ($($x:expr),*) => {
        Polynomial::new(vec![$(num::Rational64::from_integer($x)),*])
    }
}
macro_rules! expand_poly {
    ($([$($x:expr),*]),*) => {
        vec![$(poly![$($x),*]),*].into_iter().product::<Polynomial<num::Rational64>>()
    }
}

#[test]
fn test_times() {
    assert_eq!(times::<i32>(42, 756), 42 * 756);
    let a = expand_poly![[2], [1, 1], [2, 1], [3, 1]];
    let c = poly![42];
    assert_eq!(times::<Polynomial<num::Rational64>>(a.clone(), 42), c * a);
}
#[test]
fn test_power() {
    type R = Polynomial<num::Rational64>;
    assert_eq!(power::<i32>(3, 4), 81);
    let a = poly![1, 1];
    let b = poly![1, 4, 6, 4, 1];
    assert_eq!(power::<R>(a, 4), b);
}
#[test]
fn test_gcd() {
    assert_eq!(gcd::<i32>(0, 0), 0);
    assert_eq!(gcd::<i32>(42, 0), 42);
    assert_eq!(gcd::<i32>(0, 42), 42);
    assert_eq!(gcd::<i32>(64, 58), 2);
    assert_eq!(gcd::<i32>(97, 89), 1);
}
#[test]
fn test_gcd2() {
    type R = Polynomial<num::Rational64>;
    let z = R::zero();
    assert_eq!(gcd::<R>(z.clone(), z.clone()), z);
    let a = expand_poly![[2], [1, 1], [2, 1], [3, 1]];
    let b = expand_poly![[3], [1, 1], [4, 1]];
    let c = poly![1, 1];
    let d = expand_poly![[4, 1], [5, 1]];
    assert_eq!(gcd::<R>(a.clone(), z.clone()), a);
    assert_eq!(gcd::<R>(z, a.clone()), a);
    let mut m = gcd::<R>(a.clone(), b);
    m.monic();
    assert_eq!(m, c);
    let mut m = gcd::<R>(a, d);
    m.monic();
    assert!(m.is_one());
}
fn check_eea<T>(a: T, b: T) -> bool
where
    T: sealed::Sized + Zero + One + Clone + Eq + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    let g = gcd::<T>(a.clone(), b.clone());
    let (d, x, y) = extended_euclidian_algorithm::<T>(a.clone(), b.clone());
    g.is_similar(&d) && x * a + y * b == d
}
#[test]
fn test_eea() {
    assert!(check_eea::<i32>(0, 0));
    assert!(check_eea::<i32>(42, 0));
    assert!(check_eea::<i32>(0, 42));
    assert!(check_eea::<i32>(64, 58));
    assert!(check_eea::<i32>(97, 89));
}
#[test]
fn test_eea2() {
    type R = Polynomial<num::Rational64>;
    let z = R::zero();
    check_eea::<R>(z.clone(), z.clone());
    let a = expand_poly![[2], [1, 1], [2, 1], [3, 1]];
    let b = expand_poly![[3], [1, 1], [4, 1]];
    let d = expand_poly![[4, 1], [5, 1]];
    assert!(check_eea::<R>(a.clone(), z.clone()));
    assert!(check_eea::<R>(z, a.clone()));
    assert!(check_eea::<R>(a.clone(), b));
    assert!(check_eea::<R>(a, d));
}
fn check_neea<T>(a: T, b: T) -> bool
where
    T: sealed::Sized + Zero + One + Clone + Eq + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    let g = gcd::<T>(a.clone(), b.clone());
    let (d, x, y) = normalized_extended_euclidian_algorithm::<T>(a.clone(), b.clone());
    g.is_similar(&d) && x * a + y * b == d
}
#[test]
fn test_neea() {
    assert!(check_neea::<i32>(0, 0));
    assert!(check_neea::<i32>(42, 0));
    assert!(check_neea::<i32>(0, 42));
    assert!(check_neea::<i32>(64, 58));
    assert!(check_neea::<i32>(97, 89));
}
#[test]
fn test_neea2() {
    type R = Polynomial<num::Rational64>;
    let z = R::zero();
    check_eea::<R>(z.clone(), z.clone());
    let a = expand_poly![[2], [1, 1], [2, 1], [3, 1]];
    let b = expand_poly![[3], [1, 1], [4, 1]];
    let d = expand_poly![[4, 1], [5, 1]];
    assert!(check_neea::<R>(a.clone(), z.clone()));
    assert!(check_neea::<R>(z, a.clone()));
    assert!(check_neea::<R>(a.clone(), b));
    assert!(check_neea::<R>(a, d));
}
fn check_mod_inv<T>(a: T, m: T) -> Option<T>
where
    T: sealed::Sized + Zero + One + Clone + Eq + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    modulo_inverse::<T>(a.clone(), m.clone()).map(|x| T::from(&T::from(&(a * x) - &T::one()) % &m))
}
#[test]
fn test_mod_inv() {
    // not exists inverse
    assert_eq!(check_mod_inv::<i32>(0, 0), None);
    assert_eq!(check_mod_inv::<i32>(42, 0), None);
    assert_eq!(check_mod_inv::<i32>(0, 42), None);
    assert_eq!(check_mod_inv::<i32>(64, 58), None);
    // exists inverse
    assert_eq!(check_mod_inv::<i32>(97, 89), Some(0));
    assert_eq!(check_mod_inv::<i32>(7, 15), Some(0));
    assert_eq!(check_mod_inv::<i32>(42, 55), Some(0));
    assert_eq!(check_mod_inv::<i32>(15, 64), Some(0));
}
#[test]
fn test_mod_inv2() {
    type R = Polynomial<num::Rational64>;
    // not exists inverse
    let z = R::zero();
    let a = expand_poly![[2], [1, 1], [2, 1], [3, 1]];
    let b = expand_poly![[3], [1, 1], [4, 1]];
    let d = expand_poly![[4, 1], [5, 1]];
    assert_eq!(check_mod_inv::<R>(z.clone(), z.clone()), None);
    assert_eq!(check_mod_inv::<R>(a.clone(), z.clone()), None);
    assert_eq!(check_mod_inv::<R>(z, a.clone()), None);
    assert_eq!(check_mod_inv::<R>(b, d.clone()), None);
    // exists inverse
    let sz = Some(R::zero());
    assert_eq!(check_mod_inv::<R>(a, d), sz);
    let a = poly![7, 1];
    let b = expand_poly![[3, 1], [5, 1]];
    assert_eq!(check_mod_inv::<R>(a, b), sz);
    let a = poly![42, 1];
    let b = expand_poly![[5, 1], [11, 1]];
    assert_eq!(check_mod_inv::<R>(a, b), sz);
    let a = expand_poly![[3, 1], [5, 1]];
    let b = power::<R>(poly![2, 1], 6);
    assert_eq!(check_mod_inv::<R>(a, b), sz);
}

fn check_mod_div<T>(a: T, b: T, m: T) -> Option<T>
where
    T: sealed::Sized + Zero + One + Clone + Eq + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    modulo_division::<T>(a.clone(), b.clone(), m.clone())
        .map(|x| T::from(&T::from(&(b * x) - &a) % &m))
}
#[test]
fn test_mod_div() {
    // not exists inverse
    assert_eq!(check_mod_div::<i32>(0, 0, 0), None);
    assert_eq!(check_mod_div::<i32>(0, 42, 0), None);
    assert_eq!(check_mod_div::<i32>(0, 0, 42), None);
    assert_eq!(check_mod_div::<i32>(6, 4, 8), None);
    assert_eq!(check_mod_div::<i32>(1, 3, 6), None);
    // exists inverse
    assert_eq!(check_mod_div::<i32>(1, 97, 89), Some(0));
    assert_eq!(check_mod_div::<i32>(1, 7, 15), Some(0));
    assert_eq!(check_mod_div::<i32>(1, 42, 55), Some(0));
    assert_eq!(check_mod_div::<i32>(1, 15, 64), Some(0));
    assert_eq!(check_mod_div::<i32>(6, 2, 8), Some(0));
    assert_eq!(check_mod_div::<i32>(6, 9, 12), Some(0));
}
fn check_crt<T>(u: &[T], m: &[T])
where
    T: sealed::Sized + Clone + Eq + Zero + One + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    let a = chinese_remainder_theorem::<T>(u, m).unwrap();
    for (u, m) in u.iter().zip(m.iter()) {
        assert!(T::from(&T::from(&a - u) % m).is_zero());
    }
}
#[test]
fn test_crt() {
    type R = Polynomial<num::Rational64>;
    let u = vec![2, 3, 2];
    let m = vec![3, 5, 7];
    check_crt::<i32>(&u, &m);
    let u = vec![
        3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3, 3,
    ];
    let m = vec![
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97, 101,
    ];
    check_crt::<i128>(&u, &m);
    let u = vec![
        poly![3, 1],
        poly![5],
        poly![7, 1],
        poly![1, 1],
        poly![2],
        poly![1, 1],
        poly![3, 3],
        poly![1, 7],
    ];
    let m = vec![
        poly![1, 1, 1],
        poly![1, 2, 1],
        poly![1, 3, 1],
        poly![1, 4, 1],
        poly![1, 5, 1],
        poly![1, 6, 1],
        poly![1, 7, 1],
        poly![1, 8, 1],
    ];
    check_crt::<R>(&u, &m);
}
#[cfg(feature = "num-bigint")]
fn check_fcrt<T>(u: &[T], m: &[T])
where
    T: sealed::Sized + Clone + Eq + Zero + One + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    let a = fast_chinese_remainder_theorem::<T>(u, m);
    for (u, m) in u.iter().zip(m.iter()) {
        assert!(T::from(&T::from(&a - u) % m).is_zero());
    }
}
#[cfg(feature = "num-bigint")]
#[test]
fn test_fcrt() {
    type R = Polynomial<num::Rational64>;
    let u = vec![2, 3, 2, 6];
    let m = vec![3, 5, 7, 11];
    check_fcrt::<i32>(&u, &m);
    let u = vec![3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3]
        .into_iter()
        .map(num::BigInt::from)
        .collect::<Vec<_>>();
    let m = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
        .into_iter()
        .map(num::BigInt::from)
        .collect::<Vec<_>>();
    check_fcrt::<num::BigInt>(&u, &m);
    let u = vec![
        poly![3, 1],
        poly![5],
        poly![7, 1],
        poly![1, 1],
        poly![2],
        poly![1, 1],
        poly![3, 3],
        poly![1, 7],
    ];
    let m = vec![
        poly![1, 1, 1],
        poly![1, 2, 1],
        poly![1, 3, 1],
        poly![1, 4, 1],
        poly![1, 5, 1],
        poly![1, 6, 1],
        poly![1, 7, 1],
        poly![1, 8, 1],
    ];
    check_crt::<R>(&u, &m);
}
#[cfg(feature = "num-bigint")]
fn make_prime_list(n: usize) -> Vec<u64> {
    let mut v = Vec::with_capacity(n);
    for p in 2.. {
        if p.is_prime() {
            v.push(p);
            if v.len() >= n {
                return v;
            }
        }
    }
    unreachable!()
}

#[cfg(feature = "num-bigint")]
#[test]
fn test_fcrt2() {
    type Z = num::BigInt;
    let m = make_prime_list(8192)
        .into_iter()
        .map(Z::from)
        .collect::<Vec<_>>();
    let u = m
        .iter()
        .map(|m| Z::from(rand::random::<u128>()) % m)
        .collect::<Vec<_>>();
    check_fcrt::<Z>(&u, &m);
}

#[cfg(feature = "rug")]
mod rug_test {
    use super::*;
    type Z = rug::Integer;
    type Q = rug::Rational;
    #[test]
    fn test_times() {
        let a = Z::from(42i32);
        let b = 756;
        assert_eq!(times::<Z>(a.clone(), b), a * 756i32);
        let a = Q::from(3i32);
        let b = Q::from(15i32);
        assert_eq!(times::<Q>(a, 5), b);
        //let a = poly![1i32, 2i32, 3i32];
        //let c = poly![42i32];
        //assert_eq!(times::<Polynomial<Q>>(a.clone(), 42), c * a);
    }
    #[test]
    fn test_power() {
        //type R = Polynomial<Q>;
        let a = Z::from(3i32);
        let b = Z::from(81i32);
        assert_eq!(power::<Z>(a, 4), b);
        let a = Q::from(3i32) / Q::from(2i32);
        let b = Q::from(27i32) / Q::from(8i32);
        assert_eq!(power::<Q>(a, 3), b);
        //let a = poly![1, 1];
        //let b = poly![1, 4, 6, 4, 1];
        //assert_eq!(power::<R>(a, 4), b);
    }
    #[test]
    fn test_gcd() {
        let a = Z::from(42i32);
        assert_eq!(gcd::<Z>(Z::zero(), Z::zero()), Z::zero());
        assert_eq!(gcd::<Z>(a.clone(), Z::zero()), a.clone());
        assert_eq!(gcd::<Z>(Z::zero(), a.clone()), a);
        let b = Z::from(64i32);
        let c = Z::from(58i32);
        assert_eq!(gcd::<Z>(b, c), Z::from(2i32));
        let d = Z::from(97i32);
        let e = Z::from(89i32);
        assert_eq!(gcd::<Z>(d, e), Z::one());
    }
    #[test]
    fn test_coprime() {
        let a = Z::from(-57i32);
        let b = Z::from(-49i32);
        assert!(is_coprime::<Z>(a, b));
    }
    fn check_eea(a: Z, b: Z) -> bool {
        let g = gcd::<Z>(a.clone(), b.clone());
        let (d, x, y) = extended_euclidian_algorithm::<Z>(a.clone(), b.clone());
        g.is_similar(&d) && x * a + y * b == d
    }
    #[test]
    fn test_eea() {
        assert!(check_eea(Z::zero(), Z::zero()));
        assert!(check_eea(Z::from(42i32), Z::zero()));
        assert!(check_eea(Z::zero(), Z::from(42i32)));
        assert!(check_eea(Z::from(64i32), Z::from(58i32)));
        assert!(check_eea(Z::from(97i32), Z::from(89i32)));
    }
    fn check_neea(a: Z, b: Z) -> bool {
        let g = gcd::<Z>(a.clone(), b.clone());
        let (d, x, y) = normalized_extended_euclidian_algorithm::<Z>(a.clone(), b.clone());
        g.is_similar(&d) && x * a + y * b == d
    }
    #[test]
    fn test_neea() {
        assert!(check_neea(Z::zero(), Z::zero()));
        assert!(check_neea(Z::from(42i32), Z::zero()));
        assert!(check_neea(Z::zero(), Z::from(42i32)));
        assert!(check_neea(Z::from(64i32), Z::from(58i32)));
        assert!(check_neea(Z::from(97i32), Z::from(89i32)));
    }
    fn check_mod_inv(a: i32, m: i32) -> Option<i32> {
        use std::convert::TryFrom;
        let a = Z::from(a);
        let m = Z::from(m);
        modulo_inverse::<Z>(a.clone(), m.clone())
            .map(|x| (a * x - Z::one()) % m)
            .map(|x| i32::try_from(x).unwrap())
    }
    #[test]
    fn test_mod_inv() {
        // not exists inverse
        assert_eq!(check_mod_inv(0, 0), None);
        assert_eq!(check_mod_inv(42, 0), None);
        assert_eq!(check_mod_inv(0, 42), None);
        assert_eq!(check_mod_inv(64, 58), None);
        // exists inverse
        assert_eq!(check_mod_inv(97, 89), Some(0));
        assert_eq!(check_mod_inv(7, 15), Some(0));
        assert_eq!(check_mod_inv(42, 55), Some(0));
        assert_eq!(check_mod_inv(15, 64), Some(0));
    }
    fn check_mod_div(a: i32, b: i32, m: i32) -> Option<i32> {
        use std::convert::TryFrom;
        let a = Z::from(a);
        let b = Z::from(b);
        let m = Z::from(m);
        modulo_division::<Z>(a.clone(), b.clone(), m.clone())
            .map(|x| i32::try_from((b * x - a) % m).unwrap())
    }
    #[test]
    fn test_mod_div() {
        // not exists inverse
        assert_eq!(check_mod_div(0, 0, 0), None);
        assert_eq!(check_mod_div(0, 42, 0), None);
        assert_eq!(check_mod_div(0, 0, 42), None);
        assert_eq!(check_mod_div(6, 4, 8), None);
        assert_eq!(check_mod_div(1, 3, 6), None);
        // exists inverse
        assert_eq!(check_mod_div(1, 97, 89), Some(0));
        assert_eq!(check_mod_div(1, 7, 15), Some(0));
        assert_eq!(check_mod_div(1, 42, 55), Some(0));
        assert_eq!(check_mod_div(1, 15, 64), Some(0));
        assert_eq!(check_mod_div(6, 2, 8), Some(0));
        assert_eq!(check_mod_div(6, 9, 12), Some(0));
    }
    fn check_fcrt(u: &[Z], m: &[Z]) {
        let a = fast_chinese_remainder_theorem::<Z>(u, m);
        for (u, m) in u.iter().zip(m.iter()) {
            assert!(Z::from(&Z::from(&a - u) % m).is_zero());
        }
    }
    #[test]
    fn test_fcrt() {
        let u = vec![2i32, 3, 2, 6]
            .into_iter()
            .map(Z::from)
            .collect::<Vec<_>>();
        let m = vec![3i32, 5, 7, 11]
            .into_iter()
            .map(Z::from)
            .collect::<Vec<_>>();
        check_fcrt(&u, &m);
        let u = vec![3i32, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3]
            .into_iter()
            .map(Z::from)
            .collect::<Vec<_>>();
        let m = vec![
            2i32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
        ]
        .into_iter()
        .map(Z::from)
        .collect::<Vec<_>>();
        check_fcrt(&u, &m);
    }
    fn make_prime_list(n: usize) -> Vec<u64> {
        let mut v = Vec::with_capacity(n);
        for p in 2.. {
            if p.is_prime() {
                v.push(p);
                if v.len() >= n {
                    return v;
                }
            }
        }
        unreachable!()
    }
    #[test]
    fn test_fcrt2() {
        let m = make_prime_list(8192)
            .into_iter()
            .map(Z::from)
            .collect::<Vec<_>>();
        let u = m
            .iter()
            .map(|m| Z::from(rand::random::<u128>()) % m)
            .collect::<Vec<_>>();
        check_fcrt(&u, &m);
    }
}
