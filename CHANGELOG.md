# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

## [0.5.9] - 2025-03-22

### Added
- The `plot_geo_poly` function now has a `single_trace` keyword argument to allow for plotting multiple polyareas provided as a single trace (default to `false`).
- The `plot_geo_poly` function now also supports generic iterables of PolyAreas and has specific methods to be used with `AbtractRegion` inputs.
- The `plot_geo_cells` function now has a `color_contours` keyword argument to control whether the cell contours are also colored when providing a `color` vector(default to `true`).

### Fixed
- Fixed some issues with Plotting introduced with CountriesBorders.jl 0.4.8 that are now fixed with 0.4.9.
- Removed some code duplication and overspecified arguments in function signatures for internal code.
- The internal `gen_circle_pattern` function now returns a number of points which respects the `n` input keyword argument.


## [0.5.8] - 2025-03-20

### Added
- Added a new `color_greedy` function to color a set of points with a greedy algorithm, useful for assigning colors to cells which are not located on a regular grid.
- Three new subtypes of `AbstractRegion` have been added (check docstrings for more details):
  - `HotspotRegion` which is a region representing all points within a certain distance (radius) of a point.
  - `MultiRegion` which is a region encompassing multiple polyareas, and that can be constructed with a vector of `AbstractRegion`s as input.
  - `ClippedRegion` which is a region obtained by clipping another region with a mask.
- Added a new `BoxBorder` type, which expands the concept of `PolyBorder` and `MultiBorder` to the `Box` type from Meshes.jl
- Added a new abstract type `BorderGeometry{P} <: Geometry{🌐,LATLON{P}}` which is a supertype of `PolyBorder`, `MultiBorder` and `BoxBorder`. This is mostly to simplify dispatch of various methods.

### Changed
- Changed the implementation of `offset_region` to be more generic and simplified, without multiple methods for GeoRegion and PolyRegion as they were already doing basically the same thing.
- The `LatBeltRegion` now has the lim field which is a `NTuple{2, Deg{Float64}}` to be concrete and has an implementation of `bboxes` and `polyareas` to simpify plotting and inclusion within a `MultiRegion`.

### Fixed
- Fixed the tessellation function returning vectors with abstract eltype

## [0.5.7] - 2025-03-08

### Changed
- Sped up check for point inclusion in `AbstractRegion` by saving and pre-filtering on bounding boxes of polyareas, like in CountriesBorders.jl, resulting in a speed up factor of >20x
- Moved back the functionality from the `NearestNeighborsExt.jl` to GeoGrids directly as NearestNeighbors.jl was anyhow an indirect dependency coming from Meshes.jl.

## [0.5.6] - 2025-03-05

### Added
- Added a functionality to deal with polygons crossing the antimeridian (180° longitude), especially in the context of ofsetting regions. The implementation is based on the implementation in the python library at https://github.com/gadomski/antimeridian, though it's simplified to only work on a subset of cases for the moment.
- Added methods for `PlotlyBase.scattergeo` and `CountriesBorders.extract_plot_coords` to work with `MultiBorder`, `PolyBorder` and `AbstractRegion`
- Added Distances.jl as dependency and implemented the `GreatCircleMetric` as `Metric` to simplify computing haversine distances between points already in `LatLon` or `Point{🌐, <:LatLon}`.
- Added the `LatLonTree` type, which is a dummy type which only gains functionality (and meaningulf constructors) when loading `NearestNeighbors.jl` via package extensions. This allows to simplify querying _closeness_ and distances between points in `LatLon` efficiently via the functions provided by `NearestNeighbors.jl`.

### Fixed
- Offsetting polygons close to the antimeridian by enough to cause a cross of the antimeridian does not result in degenerate polygons anymore.

### Changed
- The `borders` is now the one defined in `CountriesBorders.jl` and is simply extended in GeoGrids (instead of replacing it with a local function of the same name)