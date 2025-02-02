# ChangeLog
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
## [Unreleased]
### Changed
- Changed MSRV 1.66.1
### Added
- `mul_power`
## [0.7.0] - 2023-03-13
### Changed
- Changed MSRV 1.65.0
- Removed `num-bigint`, `rug` from default features
## [0.6.1] - 2022-11-27
### Added
- Add minimum supported Rust version
## [0.6.0] - 2022-04-10
### Added
- `modulo_power`
### Changed
- second argument type of `times`, `power` is changed from `u64` to `U` (generic type)
- changed edition from 2018 to 2021
## [0.5.0] - 2022-02-26
### Added
- `rug::Integer` support
### Changed
- `RingOperation` and `EuclideanRingOperation` changed
