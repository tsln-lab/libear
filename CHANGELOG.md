# Changelog

## unreleased changes

### Added

- `GainCalculatorObjects` now implements all Objects parameters, ported from the reference implementation: Cartesian positions (allocentric panner and Cartesian extent), `objectDivergence`, `channelLock`, `zoneExclusion`, `screenRef` and `screenEdgeLock` (new `ObjectsTypeMetadata::screenEdgeLock` member).
- `GainCalculatorDirectSpeakers` supports Cartesian positions and `screenEdgeLock`.
- `M+SC` and `M-SC` loudspeakers with azimuths between 35 and 60 degrees.
- Warning codes `DIVERGENCE_POSITIONRANGE_IGNORED` and `DIVERGENCE_AZIMUTHRANGE_IGNORED` for mismatched divergence types.
- `objects_reference_tests` and `direct_speakers_reference_tests`: comparison against gains generated from the reference implementation (`tools/reference/`).

### Changed

- `Layout::screen` defaults to `getDefaultScreen()` to match the EAR. Call `layout.screen(boost::none)` to get the old behaviour.

- added xsimd submodule and updated eigen to 3.4.0; this required changing the eigen remote, so you may need to run `git submodule sync` as well as the usual `git submodule update --init --recursive`

## 0.9.0

Initial release.
