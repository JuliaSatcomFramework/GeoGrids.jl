# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Added a functionality to deal with polygons crossing the antimeridian (180° longitude), especially in the context of ofsetting regions. The implementation is based on the implementation in the python library at https://github.com/gadomski/antimeridian, though it's simplified to only work on a subset of cases for the moment.
- Added methods for `PlotlyBase.scattergeo` and `CountriesBorders.extract_plot_coords` to work with `MultiBorder`, `PolyBorder` and `AbstractRegion`
- Added Distances.jl as dependency and implemented the `GreatCircleMetric` as `Metric` to simplify computing haversine distances between points already in `LatLon` or `Point{🌐, <:LatLon}`.
- Added the `LatLonTree` type, which is a dummy type which only gains functionality (and meaningulf constructors) when loading `NearestNeighbors.jl` via package extensions. This allows to simplify querying _closeness_ and distances between points in `LatLon` efficiently via the functions provided by `NearestNeighbors.jl`.

### Fixed
- Offsetting polygons close to the antimeridian by enough to cause a cross of the antimeridian does not result in degenerate polygons anymore.

### Changed
- The `borders` is now the one defined in `CountriesBorders.jl` and is simply extended in GeoGrids (instead of replacing it with a local function of the same name)