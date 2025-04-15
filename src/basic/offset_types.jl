"""
    GeoRegionOffset{D,P} <: AbstractRegion

Type representing an enlarged geographical region based on a GeoRegion.

Fields:
- `original::GeoRegion{D,P}`: The original GeoRegion
- `name::String`: Name of the enlarged region
- `domain::MultiBorder{P}`: Domain of the enlarged region
- `convexhull::PolyBorder{P}`: Convex hull of the enlarged region

# Constructors

    GeoRegionOffset(; name="offset_georegion", continent="", subregion="", admin="", delta::Number, resolution=110, refRadius=constants.Re_mean, magnitude=3, precision=7)
    GeoRegionOffset(gr::GeoRegion, delta::Number; name="offset_georegion", refRadius=constants.Re_mean, magnitude=3, precision=7)

Create an enlarged GeoRegion either from scratch or from an existing GeoRegion.

## Arguments
- `delta`: Distance to enlarge the region by, in meters
- `gr::GeoRegion`: The original GeoRegion to enlarge (for the second constructor)

## Keyword Arguments
- `name::String="enlarged_region"`: Name of the enlarged region
- `continent::String=""`: Continent of the region (only for the first constructor)
- `subregion::String=""`: Subregion within the continent (only for the first constructor)
- `admin::String=""`: Administrative area (only for the first constructor)
- `resolution::Int=110`: Resolution of the geographical data (only for the first constructor)
- `refRadius::Float64=constants.Re_mean`: Reference radius of the Earth
- `magnitude::Int=3`: Magnitude for polygon offsetting
- `precision::Int=7`: Precision for polygon offsetting

## Returns
- `GeoRegionOffset`: The enlarged geographical region
"""
mutable struct GeoRegionOffset{D,P} <: AbstractRegion
    original::GeoRegion{D,P}
    name::String
    domain::MultiBorder{P}
    convexhull::PolyBorder{P}
end

function GeoRegionOffset(; name="offset_georegion", continent="", subregion="", admin="", skip_areas=nothing, delta::Number, resolution=110, refRadius=constants.Re_mean, magnitude=3, precision=7)
    gr = GeoRegion(; name, continent, subregion, admin, skip_areas, resolution)

    GeoRegionOffset(gr, delta; name, refRadius, magnitude, precision)
end

function GeoRegionOffset(gr::GeoRegion, delta::Number; name="offset_georegion", refRadius=constants.Re_mean, magnitude=3, precision=7)
    # Resolution and Skip Areas must be set on the original GeoRegion definition, cannot be set a posteriori on the offset.
    orLatLon = offset_region(gr, delta; refRadius, magnitude, precision)
    orCart = cartesian_geometry(orLatLon)
    or = MultiBorder(orLatLon, orCart)

    chCart = convexhull(orCart)
    chLatlon = latlon_geometry(chCart)
    ch = PolyBorder(chLatlon, chCart)

    GeoRegionOffset(gr, name, or, ch)
end

"""
    PolyRegionOffset{P} <: AbstractRegion

Struct representing an enlarged polygonal region.

Fields:
- `original::PolyRegion{P}`: The original PolyRegion
- `name::String`: Name of the enlarged region
- `domain::MultiBorder{P}`: Domain of the enlarged region as a MultiBorder

# Constructors

    PolyRegionOffset(delta::Number; kwargs...)
    PolyRegionOffset(pr::PolyRegion, delta::Number; kwargs...)

Create an enlarged PolyRegion either from scratch or from an existing PolyRegion.

## Arguments
- `deltaDist`: Distance to enlarge the region by, in meters
- `pr::PolyRegion`: The original PolyRegion to enlarge (for the second constructor)

## Keyword Arguments
- `name::String="offset_polyregion"`: Name of the enlarged region
- `domain`: Domain of the region (only for the first constructor)
- `refRadius::Float64=constants.Re_mean`: Reference radius of the Earth
- `magnitude::Int=3`: Magnitude for polygon offsetting
- `precision::Int=7`: Precision for polygon offsetting

## Returns
- `PolyRegionOffset`: The enlarged polygonal region
"""
mutable struct PolyRegionOffset{P} <: AbstractRegion
    original::PolyRegion{P}
    name::String
    domain::MultiBorder{P}
end

function PolyRegionOffset(; name="offset_polyregion", domain, delta::Number, refRadius=constants.Re_mean, magnitude=3, precision=7)
    pr = PolyRegion(name, domain)

    PolyRegionOffset(pr, delta; name, refRadius, magnitude, precision)
end

function PolyRegionOffset(pr::PolyRegion, delta::Number; name="offset_polyregion", refRadius=constants.Re_mean, magnitude=3, precision=7)
    orLatLon = offset_region(pr, delta; refRadius, magnitude, precision)
    orCart = cartesian_geometry(orLatLon)
    or = MultiBorder(orLatLon, orCart)

    PolyRegionOffset(pr, name, or)
end