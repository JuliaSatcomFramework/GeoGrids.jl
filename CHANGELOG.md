# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Added a functionality to deal with polygons crossing the antimeridian (180° longitude), especially in the context of ofsetting regions. The implementation is based on the implementation in the python library at https://github.com/gadomski/antimeridian, though it's simplified to only work on a subset of cases for the moment.

### Changed
- The `borders` is now the one defined in `CountriesBorders.jl` and is simply extended in GeoGrids (instead of replacing it with a local function of the same name)