![Tests](https://github.com/stencillogic/num-bigfloat/workflows/Rust/badge.svg)
![Minimum rustc version](https://img.shields.io/badge/rustc-1.62.1+-blue.svg)


Increased precision floating point numbers implemented purely in Rust. 

## Arbitrary precision numbers

**A fork of the library is available that implements floating-point numbers of arbitrary precision: [astro-float](https://github.com/stencillogic/astro-float).
The new library uses more elaborate algorithms, has higher accuracy, and better test coverage.**

## Rationale

There are several notable implementations of numbers with increased precision for Rust. Among the libraries, one can recall [num-bigint](https://crates.io/crates/num-bigint), [rust_decimal](https://crates.io/crates/rust_decimal).

While these libraries are great in many ways, they don't allow you to perform operations on numbers while still providing fairly high precision.

There are also wrapper libraries, like [rug](https://crates.io/crates/rug), that depend on MPFR for implementing arbitrary precision floating point numbers.

This library is written in pure Rust, provides more precision than f32, f64, and some other data types with increased precision.

## Number characteristics

Number has fixed-size mantissa and exponent.

| Name                          | Value  |
|:------------------------------|-------:|
| Base                          |     10 |
| Decimal positions in mantissa |     40 |
| Exponent minimum value        |   -128 |
| Exponent maximum value        |    127 |


## no_std

Library can be used without the standard Rust library. This can be achieved by turning off `std` feature.


## Performance

Benchmark for `num-bigfloat`, `rug`(`MPFR`), and `astro-float` can be found here: [bigfloat-bench](https://github.com/stencillogic/bigfloat-bench).

## Other features

The library depends on [rand](https://crates.io/crates/rand) which is used by `BigFloat::random_normal`. This dependecy can be excluded by turning off the `rand` feature. `rand` feature requires `std` feature.

The library depends on [serde](https://crates.io/crates/serde) which is also optional and can be eliminated by turning off the `serde` feature.

In addition, the library implements [num-traits](https://crates.io/crates/num_traits). Dependency on `num-traits` can be excluded by turning off `num-traits` feature.