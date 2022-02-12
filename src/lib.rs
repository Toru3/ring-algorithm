#![cfg_attr(feature = "__internal_inject_debug", recursion_limit = "8")]
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
pub use ring_traits::{EuclideanRingOperation, RingNormalize, RingOperation};

/** calcurate $`pa`$ with mutliprecation by doubling
```
use ring_algorithm::times;
assert_eq!(times::<i32>(2, 16), 32);
```
*/
pub fn times<T>(a: T, mut p: u64) -> T
where
    T: sealed::Sized + num_traits::Zero + for<'x> std::ops::AddAssign<&'x T>,
    for<'x> &'x T: std::ops::Add<Output = T>,
{
    let mut x = T::zero();
    let mut y = a;
    loop {
        if p % 2 == 1 {
            x += &y;
        }
        p /= 2;
        if p == 0 {
            break;
        }
        y = &y + &y;
    }
    x
}

/** calcurate $`a^p`$ with exponentiation by squaring
```
use ring_algorithm::power;
assert_eq!(power::<i32>(2, 16), 65536);
```
*/
pub fn power<T>(a: T, mut p: u64) -> T
where
    T: sealed::Sized + num_traits::One + for<'x> std::ops::MulAssign<&'x T>,
    for<'x> &'x T: std::ops::Mul<Output = T>,
{
    let mut x = T::one();
    let mut y = a;
    loop {
        if p % 2 == 1 {
            x *= &y;
        }
        p /= 2;
        if p == 0 {
            break;
        }
        y = &y * &y;
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
    T: sealed::Sized + num_traits::Zero,
    for<'x> &'x T: std::ops::Rem<Output = T>,
{
    while !y.is_zero() {
        let r = &x % &y;
        x = y;
        y = r;
    }
    x
}

/** test $`\gcd(x, y) = 1`$
*/
pub fn is_coprime<T>(x: T, y: T) -> bool
where
    T: sealed::Sized + Eq + num_traits::Zero + num_traits::One + RingNormalize,
    for<'x> &'x T: std::ops::Rem<Output = T>,
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
    T: sealed::Sized + num_traits::Zero + num_traits::One,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    let mut old = (x, T::one(), T::zero());
    let mut now = (y, T::zero(), T::one());
    while !now.0.is_zero() {
        let q = &old.0 / &now.0;
        let new = (
            &old.0 - &(&q * &now.0),
            &old.1 - &(&q * &now.1),
            &old.2 - &(&q * &now.2),
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
    T: sealed::Sized + num_traits::Zero + num_traits::One + RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    let lc_x = x.leading_unit();
    let lc_y = y.leading_unit();
    let mut old = (x.into_normalize(), &T::one() / &lc_x, T::zero());
    let mut now = (y.into_normalize(), T::zero(), &T::one() / &lc_y);
    while !now.0.is_zero() {
        let q = &old.0 / &now.0;
        let r = &old.0 % &now.0;
        let lc_r = r.leading_unit();
        let new = (
            r.into_normalize(),
            &(&old.1 - &(&q * &now.1)) / &lc_r,
            &(&old.2 - &(&q * &now.2)) / &lc_r,
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
    T: sealed::Sized + Eq + num_traits::Zero + num_traits::One + RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
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
    T: sealed::Sized + Clone + Eq + num_traits::Zero + num_traits::One + RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    if m.is_zero() || b.is_zero() {
        return None;
    }
    let (gcd, inv_b, _) = normalized_extended_euclidian_algorithm::<T>(b, m.clone());
    if (&a % &gcd).is_zero() {
        Some(&(&a / &gcd * inv_b) % &m)
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
    T: sealed::Sized + Clone + Eq + num_traits::Zero + num_traits::One + RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
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
                (&(v_j + &(m_j * &t)) % m_i, &(&n * m_j) % m_i)
            });
        v.push(modulo_division::<T>(u_i + &(m_i - &t), n, m_i.clone())?);
    }
    let mut ret = v.pop().unwrap();
    for (v_i, m_i) in v.iter().zip(m.iter()).rev() {
        ret = &(&ret * m_i) + v_i;
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
let r = build_subproduct_tree(a);
assert_eq!(b, r);
let a = vec![2, 3, 4, 5, 6];
let b = vec![720, 120, 6, 6, 20, 6, 1, 2, 3, 4, 5, 6, 1, 1, 1];
let r = build_subproduct_tree(a);
assert_eq!(b, r);
```
*/
pub fn build_subproduct_tree<T>(a: Vec<T>) -> Vec<T>
where
    T: sealed::Sized + Clone + num_traits::One,
    for<'x> &'x T: std::ops::Mul<Output = T>,
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
            .for_each(|(ro, ri)| *ro = &ri[0] * &ri[1]);
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
let r = build_subsum_tree(a);
assert_eq!(b, r);
let a = vec![2, 3, 4, 5, 6];
let b = vec![20, 14, 6, 5, 9, 6, 0, 2, 3, 4, 5, 6, 0, 0, 0];
let r = build_subsum_tree(a);
assert_eq!(b, r);
```
*/
pub fn build_subsum_tree<T>(a: Vec<T>) -> Vec<T>
where
    T: sealed::Sized + Clone + num_traits::Zero,
    for<'x> &'x T: std::ops::Add<Output = T>,
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
            .for_each(|(ro, ri)| *ro = &ri[0] + &ri[1]);
        len >>= 1;
    }
    r
}

fn modular_reduction_aux<T, U>(f: &U, v: &[T], sv: usize, start: usize, end: usize) -> Vec<U>
where
    T: sealed::Sized,
    U: sealed::Sized,
    for<'x> &'x U: std::ops::Rem<&'x T, Output = U>,
{
    let s = (end - start) / 2;
    let mid = start + s;
    let ev = sv * 2 + 1;
    if s == 1 {
        let f0 = f % &v[sv + start];
        let f1 = f % &v[sv + start + 1];
        vec![f0, f1]
    } else {
        let mut v0 = {
            let f0 = f % &v[sv + start / s];
            modular_reduction_aux::<T, U>(&f0, v, ev, start, mid)
        };
        let mut v1 = {
            let f1 = f % &v[sv + start / s + 1];
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
let v = modular_reduction(&t, &m);
let w = m.iter().map(|m| t % m).collect::<Vec<_>>();
assert_eq!(v, w);
```
*/
pub fn modular_reduction<T, U>(f: &U, m: &[T]) -> Vec<U>
where
    T: sealed::Sized + Clone + num_traits::One,
    for<'x> &'x T: std::ops::Mul<Output = T>,
    U: sealed::Sized,
    for<'x> &'x U: std::ops::Rem<&'x T, Output = U>,
{
    let len = m.len();
    assert!(len >= 2 && len.is_power_of_two());
    let v = build_subproduct_tree::<T>(m.to_vec());
    modular_reduction_aux::<T, U>(f, &v, 1, 0, len)
}

fn crt_inverses<T>(m: &[T], big_m: &T) -> Vec<T>
where
    T: sealed::Sized + Clone + Eq + num_traits::Zero + num_traits::One + ring_traits::RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    let m2 = m.iter().map(|m| m * m).collect::<Vec<_>>();
    modular_reduction::<T, T>(big_m, &m2)
        .into_iter()
        .zip(m.iter())
        .map(|(v, m)| modulo_inverse::<T>(&v / m, m.clone()).unwrap())
        .collect()
}

fn crt_combination<T>(c: &[T], m: &[T], v: &[T], sv: usize, start: usize, end: usize) -> T
where
    T: sealed::Sized + Clone + Eq + num_traits::Zero + num_traits::One + ring_traits::RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
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
        &(&r0 * &v[sv + start / s + 1]) + &(&r1 * &v[sv + start / s])
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
    T: sealed::Sized + Clone + Eq + num_traits::Zero + num_traits::One + ring_traits::RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    assert_eq!(u.len(), m.len());
    let len = m.len();
    assert!(len >= 2 && len.is_power_of_two());
    let v = build_subproduct_tree::<T>(m.to_vec());
    let s = crt_inverses::<T>(m, &v[0]);
    let c = s
        .iter()
        .zip(u.iter())
        .map(|(s, u)| s * u)
        .collect::<Vec<_>>();
    crt_combination::<T>(&c, m, &v, 1, 0, len)
}
