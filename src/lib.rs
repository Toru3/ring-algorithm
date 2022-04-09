#![cfg_attr(feature = "__internal_inject_debug", recursion_limit = "8")]
use num_traits::{One, Zero};
use std::ops::{Add, AddAssign, BitAnd, Mul, MulAssign, Rem, ShrAssign};
mod sealed {
    pub trait SizedExt: std::marker::Sized + std::fmt::Debug + std::fmt::Display {}
    impl<T> SizedExt for T where T: std::marker::Sized + std::fmt::Debug + std::fmt::Display {}
    #[cfg(not(feature = "__internal_inject_debug"))]
    pub use std::marker::Sized;
    #[cfg(feature = "__internal_inject_debug")]
    pub use SizedExt as Sized;
}
mod ring_traits;
#[cfg(test)]
mod test;
pub use ring_traits::{
    EuclideanRingOperation, EuclideanRingOperationFrom, RingNormalize, RingOperation,
    RingOperationFrom,
};

/** calcurate $`pa`$ with mutliprecation by doubling
```
use ring_algorithm::times;
assert_eq!(times::<i32, u64>(2, 16), 32);
```
*/
pub fn times<T, U>(a: T, mut p: U) -> T
where
    T: sealed::Sized + Zero + for<'x> AddAssign<&'x T> + for<'x> From<<&'x T as Add>::Output>,
    for<'x> &'x T: Add,
    U: sealed::Sized
        + Zero
        + One
        + Eq
        + for<'x> ShrAssign<usize>
        + for<'x> From<<&'x U as BitAnd>::Output>,
    for<'x> &'x U: BitAnd,
{
    let mut x = T::zero();
    let mut y = a;
    loop {
        if U::from(&p & &U::one()) == U::one() {
            x += &y;
        }
        p >>= 1;
        if p == U::zero() {
            break;
        }
        y = T::from(&y + &y);
    }
    x
}

/** calcurate $`a^p`$ with exponentiation by squaring
```
use ring_algorithm::power;
assert_eq!(power::<i32, u64>(2, 16), 65536);
```
*/
pub fn power<T, U>(a: T, mut p: U) -> T
where
    T: sealed::Sized + One + for<'x> MulAssign<&'x T> + for<'x> From<<&'x T as Mul>::Output>,
    for<'x> &'x T: Mul,
    U: sealed::Sized
        + Zero
        + One
        + Eq
        + for<'x> ShrAssign<usize>
        + for<'x> From<<&'x U as BitAnd>::Output>,
    for<'x> &'x U: BitAnd,
{
    let mut x = T::one();
    let mut y = a;
    loop {
        if U::from(&p & &U::one()) == U::one() {
            x *= &y;
        }
        p >>= 1;
        if p == U::zero() {
            break;
        }
        y = T::from(&y * &y);
    }
    x
}

/** calcurate greatest common divisor
```
use ring_algorithm::gcd;
assert_eq!(gcd::<i32>(15, 21), 3);
assert_eq!(gcd::<i32>(14, 15), 1);
assert_eq!(gcd::<i32>(0, 42), 42);
assert_eq!(gcd::<i32>(0, 0), 0);
```
*/
pub fn gcd<T>(mut x: T, mut y: T) -> T
where
    T: sealed::Sized + Zero + for<'x> From<<&'x T as Rem>::Output>,
    for<'x> &'x T: Rem,
{
    while !y.is_zero() {
        let r = T::from(&x % &y);
        x = y;
        y = r;
    }
    x
}

/** test $`\gcd(x, y) = 1`$
*/
pub fn is_coprime<T>(x: T, y: T) -> bool
where
    T: sealed::Sized + Eq + Zero + One + RingNormalize + for<'x> From<<&'x T as Rem>::Output>,
    for<'x> &'x T: Rem,
{
    gcd::<T>(x, y).into_normalize().is_one()
}

/** extended euclidian algorithm

calcurate g (`gcd(a, b)`) and x, y ( $`g = ax + by`$ )
```
use ring_algorithm::{gcd, extended_euclidian_algorithm};
let a = 314;
let b = 271;
let (d, x, y) = extended_euclidian_algorithm::<i32>(a, b);
assert_eq!(d, gcd::<i32>(a, b));
assert_eq!(d, x * a + y * b);
```
 */
pub fn extended_euclidian_algorithm<T>(x: T, y: T) -> (T, T, T)
where
    T: sealed::Sized + Zero + One + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    let mut old = (x, T::one(), T::zero());
    let mut now = (y, T::zero(), T::one());
    while !now.0.is_zero() {
        let q = T::from(&old.0 / &now.0);
        let new = (
            T::from(&old.0 - &T::from(&q * &now.0)),
            T::from(&old.1 - &T::from(&q * &now.1)),
            T::from(&old.2 - &T::from(&q * &now.2)),
        );
        old = now;
        now = new;
    }
    old
}

/** extended euclidian algorithm with normalize
```
use ring_algorithm::{gcd, normalized_extended_euclidian_algorithm, RingNormalize};
let a = 314;
let b = 271;
let (d, x, y) = normalized_extended_euclidian_algorithm::<i32>(a, b);
assert_eq!(d, gcd::<i32>(a, b));
assert_eq!(d, x * a + y * b);
```
*/
pub fn normalized_extended_euclidian_algorithm<T>(x: T, y: T) -> (T, T, T)
where
    T: sealed::Sized + Zero + One + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    let lc_x = x.leading_unit();
    let lc_y = y.leading_unit();
    let mut old = (x.into_normalize(), T::from(&T::one() / &lc_x), T::zero());
    let mut now = (y.into_normalize(), T::zero(), T::from(&T::one() / &lc_y));
    while !now.0.is_zero() {
        let q = T::from(&old.0 / &now.0);
        let r = T::from(&old.0 % &now.0);
        let lc_r = r.leading_unit();
        let new = (
            r.into_normalize(),
            T::from(&T::from(&old.1 - &T::from(&q * &now.1)) / &lc_r),
            T::from(&T::from(&old.2 - &T::from(&q * &now.2)) / &lc_r),
        );
        old = now;
        now = new;
    }
    old
}

/** calc inverse in modulo

calc x ($`ax \equiv 1 \pmod{m}`$)
```
use ring_algorithm::modulo_inverse;
let a = 42;
let m = 55;
let b = modulo_inverse::<i32>(a, m).unwrap();
assert_eq!((a * b - 1) % m, 0);
```
*/
pub fn modulo_inverse<T>(a: T, m: T) -> Option<T>
where
    T: sealed::Sized + Eq + Zero + One + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    if m.is_zero() || a.is_zero() {
        return None;
    }
    let (gcd, inv_a, _) = normalized_extended_euclidian_algorithm::<T>(a, m);
    if gcd.is_one() {
        Some(inv_a)
    } else {
        None
    }
}

/** division in modulo

calc x ($`bx \equiv a \pmod{m}`$)
```
use ring_algorithm::modulo_division;
let a = 42;
let b = 32;
let m = 98;
let x = modulo_division::<i32>(a, b, m).unwrap();
assert_eq!((b * x - a) % m, 0);
```
*/
pub fn modulo_division<T>(a: T, b: T, m: T) -> Option<T>
where
    T: sealed::Sized + Clone + Eq + Zero + One + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    if m.is_zero() || b.is_zero() {
        return None;
    }
    let (gcd, inv_b, _) = normalized_extended_euclidian_algorithm::<T>(b, m.clone());
    if T::from(&a % &gcd).is_zero() {
        let q = T::from(&a / &gcd);
        let t = q * inv_b;
        Some(T::from(&t % &m))
    } else {
        None
    }
}

/** Chinese remainder theorem

```
use ring_algorithm::chinese_remainder_theorem;
let u = vec![2, 3, 2];
let m = vec![3, 5, 7];
let a = chinese_remainder_theorem::<i32>(&u, &m).unwrap();
for (u, m) in u.iter().zip(m.iter()) {
    assert_eq!((a - u) % m, 0);
}
```
*/
pub fn chinese_remainder_theorem<T>(u: &[T], m: &[T]) -> Option<T>
where
    T: sealed::Sized + Clone + Eq + Zero + One + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    if u.len() != m.len() {
        return None;
    }
    let len = u.len();
    let mut v = Vec::with_capacity(len);
    v.push(u[0].clone());
    for (i, (u_i, m_i)) in u.iter().zip(m.iter()).enumerate().skip(1) {
        let p = (v[i - 1].clone(), m[i - 1].clone());
        let (t, n) = m
            .iter()
            .zip(v.iter())
            .rev()
            .skip(1)
            .fold(p, |(t, n), (m_j, v_j)| {
                (
                    T::from(&T::from(v_j + &T::from(m_j * &t)) % m_i),
                    T::from(&T::from(&n * m_j) % m_i),
                )
            });
        v.push(modulo_division::<T>(
            T::from(u_i + &T::from(m_i - &t)),
            n,
            m_i.clone(),
        )?);
    }
    let mut ret = v.pop().unwrap();
    for (v_i, m_i) in v.iter().zip(m.iter()).rev() {
        ret = T::from(&T::from(&ret * m_i) + v_i);
    }
    Some(ret)
}

/** calc subproduct tree

```text
 120
  /\
 6  20
/\  /\
2 3 4 5
```

```rust
use ring_algorithm::build_subproduct_tree;
let a = vec![2, 3, 4, 5];
let b = vec![120, 6, 20, 2, 3, 4, 5];
let r = build_subproduct_tree::<i32>(a);
assert_eq!(b, r);
let a = vec![2, 3, 4, 5, 6];
let b = vec![720, 120, 6, 6, 20, 6, 1, 2, 3, 4, 5, 6, 1, 1, 1];
let r = build_subproduct_tree::<i32>(a);
assert_eq!(b, r);
```
*/
pub fn build_subproduct_tree<T>(a: Vec<T>) -> Vec<T>
where
    T: sealed::Sized + Clone + One + for<'x> From<<&'x T as Mul>::Output>,
    for<'x> &'x T: Mul,
{
    let mut len = a.len().next_power_of_two();
    let mut r = vec![T::one(); 2 * len - 1];
    a.into_iter()
        .zip(r.iter_mut().skip(len - 1))
        .for_each(|(a, r)| *r = a);
    while len >= 2 {
        let (ro, ri) = r.split_at_mut(len - 1);
        let (_, ro) = ro.split_at_mut(len / 2 - 1);
        ro.iter_mut()
            .zip(ri.chunks_exact(2))
            .for_each(|(ro, ri)| *ro = T::from(&ri[0] * &ri[1]));
        len >>= 1;
    }
    r
}

/** calc subsum tree

```text
  14
  /\
 5  9
/\  /\
2 3 4 5
```

```rust
use ring_algorithm::build_subsum_tree;
let a = vec![2, 3, 4, 5];
let b = vec![14, 5, 9, 2, 3, 4, 5];
let r = build_subsum_tree::<i32>(a);
assert_eq!(b, r);
let a = vec![2, 3, 4, 5, 6];
let b = vec![20, 14, 6, 5, 9, 6, 0, 2, 3, 4, 5, 6, 0, 0, 0];
let r = build_subsum_tree::<i32>(a);
assert_eq!(b, r);
```
*/
pub fn build_subsum_tree<T>(a: Vec<T>) -> Vec<T>
where
    T: sealed::Sized + Clone + Zero + for<'x> From<<&'x T as Add>::Output>,
    for<'x> &'x T: Add,
{
    let mut len = a.len().next_power_of_two();
    let mut r = vec![T::zero(); 2 * len - 1];
    a.into_iter()
        .zip(r.iter_mut().skip(len - 1))
        .for_each(|(a, r)| *r = a);
    while len >= 2 {
        let (ro, ri) = r.split_at_mut(len - 1);
        let (_, ro) = ro.split_at_mut(len / 2 - 1);
        ro.iter_mut()
            .zip(ri.chunks_exact(2))
            .for_each(|(ro, ri)| *ro = T::from(&ri[0] + &ri[1]));
        len >>= 1;
    }
    r
}

fn modular_reduction_aux<T, U>(f: &U, v: &[T], sv: usize, start: usize, end: usize) -> Vec<U>
where
    T: sealed::Sized,
    U: sealed::Sized + for<'x> From<<&'x U as Rem<&'x T>>::Output>,
    for<'x> &'x U: Rem<&'x T>,
{
    let s = (end - start) / 2;
    let mid = start + s;
    let ev = sv * 2 + 1;
    if s == 1 {
        let f0 = U::from(f % &v[sv + start]);
        let f1 = U::from(f % &v[sv + start + 1]);
        vec![f0, f1]
    } else {
        let mut v0 = {
            let f0 = U::from(f % &v[sv + start / s]);
            modular_reduction_aux::<T, U>(&f0, v, ev, start, mid)
        };
        let mut v1 = {
            let f1 = U::from(f % &v[sv + start / s + 1]);
            modular_reduction_aux::<T, U>(&f1, v, ev, mid, end)
        };
        v0.append(&mut v1);
        v0
    }
}

/** modular reduction

`m.len()` must be $`2^k~(k=1,2,3,\ldots)`$

```
use ring_algorithm::modular_reduction;
let t = 997u128;
let m = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
let v = modular_reduction::<u128, u128>(&t, &m);
let w = m.iter().map(|m| t % m).collect::<Vec<_>>();
assert_eq!(v, w);
```
*/
pub fn modular_reduction<T, U>(f: &U, m: &[T]) -> Vec<U>
where
    T: sealed::Sized + Clone + One + for<'x> From<<&'x T as Mul>::Output>,
    for<'x> &'x T: Mul,
    U: sealed::Sized + for<'x> From<<&'x U as Rem<&'x T>>::Output>,
    for<'x> &'x U: Rem<&'x T>,
{
    let len = m.len();
    assert!(len >= 2 && len.is_power_of_two());
    let v = build_subproduct_tree::<T>(m.to_vec());
    modular_reduction_aux::<T, U>(f, &v, 1, 0, len)
}

fn crt_inverses<T>(m: &[T], big_m: &T) -> Vec<T>
where
    T: sealed::Sized + Clone + Eq + Zero + One + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    let m2 = m.iter().map(|m| T::from(m * m)).collect::<Vec<_>>();
    modular_reduction::<T, T>(big_m, &m2)
        .into_iter()
        .zip(m.iter())
        .map(|(v, m)| modulo_inverse::<T>(T::from(&v / m), m.clone()).unwrap())
        .collect()
}

fn crt_combination<T>(c: &[T], m: &[T], v: &[T], sv: usize, start: usize, end: usize) -> T
where
    T: sealed::Sized + Clone + Eq + Zero + One + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    if end - start == 1 {
        c[0].clone()
    } else {
        let s = (end - start) / 2;
        let mid = start + s;
        let ev = sv * 2 + 1;
        let (c0, c1) = c.split_at(s);
        let (m0, m1) = m.split_at(s);
        let r0 = crt_combination::<T>(c0, m0, v, ev, start, mid);
        let r1 = crt_combination::<T>(c1, m1, v, ev, mid, end);
        T::from(&T::from(&r0 * &v[sv + start / s + 1]) + &T::from(&r1 * &v[sv + start / s]))
    }
}

/** Chinese remainder theorem

`m.len()` must be $`2^k~(k=1,2,3,\ldots)`$ and `u.len() == m.len()`
```
use ring_algorithm::fast_chinese_remainder_theorem;
let u = vec![2, 3, 2, 6];
let m = vec![3, 5, 7, 11];
let a = fast_chinese_remainder_theorem::<i32>(&u, &m);
for (u, m) in u.iter().zip(m.iter()) {
    assert_eq!((a - u) % m, 0);
}
```
*/
pub fn fast_chinese_remainder_theorem<T>(u: &[T], m: &[T]) -> T
where
    T: sealed::Sized + Clone + Eq + Zero + One + RingNormalize + EuclideanRingOperationFrom,
    for<'x> &'x T: EuclideanRingOperation,
{
    assert_eq!(u.len(), m.len());
    let len = m.len();
    assert!(len >= 2 && len.is_power_of_two());
    let v = build_subproduct_tree::<T>(m.to_vec());
    let s = crt_inverses::<T>(m, &v[0]);
    let c = s
        .iter()
        .zip(u.iter())
        .map(|(s, u)| T::from(s * u))
        .collect::<Vec<_>>();
    crt_combination::<T>(&c, m, &v, 1, 0, len)
}

/** power in modulo

```rust
use ring_algorithm::modulo_power;
let a = 314i32;
let p = 271i32;
let m = 42;
assert_eq!(modulo_power(a, p, &m), 20);
```
*/
pub fn modulo_power<T, U>(a: T, mut p: U, m: &T) -> T
where
    T: sealed::Sized
        + One
        + for<'x> MulAssign<&'x T>
        + for<'x> From<<&'x T as Mul>::Output>
        + for<'x> From<<&'x T as Rem>::Output>,
    for<'x> &'x T: Mul + Rem,
    U: sealed::Sized
        + Zero
        + One
        + Eq
        + for<'x> ShrAssign<usize>
        + for<'x> From<<&'x U as BitAnd>::Output>,
    for<'x> &'x U: BitAnd,
{
    let mut x = T::one();
    let mut y = a;
    loop {
        if U::from(&p & &U::one()) == U::one() {
            x = T::from(&T::from(&x * &y) % m);
        }
        p >>= 1;
        if p == U::zero() {
            break;
        }
        y = T::from(&T::from(&y * &y) % m);
    }
    x
}
