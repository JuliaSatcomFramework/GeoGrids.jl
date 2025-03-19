## Angle Types
const UnitfulAngleType = Union{typeof(°),typeof(rad)}
const UnitfulAngleQuantity = Quantity{<:Real,<:Any,<:UnitfulAngleType}
const ValidAngle = Union{UnitfulAngleQuantity,Real}
const ValidDistance = Union{Length, Real}

const constants = (
    Re_mean = 6371e3, # Mean Earth Radius [m]
    a = 6378137, # [m] WGS84 semi-major axis
    b = 6356752.315 # [m] WGS84 semi-minor axis
)

"""
    PolyBorder{P} <: Geometry{🌐,LATLON{P}}

Struct representing a PolyArea in both LatLon and Cartesian coordinates.

Fields:
- `latlon::POLY_LATLON{P}`: The borders in LatLon CRS
- `cart::POLY_CART{P}`: The borders in Cartesian2D CRS

Where `P` is the precision type (e.g., Float32, Float64) for the coordinates.
"""
struct PolyBorder{P} <: Geometry{🌐,LATLON{P}}
    latlon::POLY_LATLON{P}
    cart::POLY_CART{P}
    bbox::BOX_CART{P}
end

PolyBorder(latlon::POLY_LATLON) = PolyBorder(latlon, cartesian_geometry(latlon))
function PolyBorder(latlon::POLY_LATLON, cart::POLY_CART) 
    bbox = boundingbox(cart)
    PolyBorder(latlon, cart, bbox)
end

"""
    MultiBorder{P} <: Geometry{🌐,LATLON{P}}

Struct representing a Multi in both LatLon and Cartesian coordinates.

Fields:
- `latlon::MULTI_LATLON{P}`: The borders in LatLon CRS
- `cart::MULTI_CART{P}`: The borders in Cartesian2D CRS

Where `P` is the precision type (e.g., Float32, Float64) for the coordinates.
"""
struct MultiBorder{P} <: Geometry{🌐,LATLON{P}}
    latlon::MULTI_LATLON{P}
    cart::MULTI_CART{P}
    bboxes::Vector{BOX_CART{P}}
end

MultiBorder(latlon::MULTI_LATLON) = 
    MultiBorder(latlon, cartesian_geometry(latlon))
MultiBorder(latlon::MULTI_LATLON, cart::MULTI_CART) =
    MultiBorder(latlon, cart, bboxes(cart))

abstract type AbstractRegion end

"""
    GlobalRegion <: AbstractRegion

Type representing a global region.
"""
struct GlobalRegion <: AbstractRegion end

"""
    GeoRegion{D,P} <: AbstractRegion

Type representing a geographical region based on CountriesBorders.

Fields:
- `name::String`: Name of the region
- `continent::String`: Continent of the region
- `subregion::String`: Subregion within the continent
- `admin::String`: Administrative area
- `domain::D`: Domain of the region
- `convexhull::PolyBorder{P}`: Convex hull of the region

Where `D` is the domain type and `P` is the precision type for coordinates.
"""
mutable struct GeoRegion{D,P} <: AbstractRegion
    name::String
    continent::String
    subregion::String
    admin::String
    domain::D
    convexhull::PolyBorder{P}
end

function GeoRegion(; name="", continent="", subregion="", admin="", skip_areas=nothing, resolution=110)
    all(isempty(v) for v in (continent, subregion, admin)) && error("Input at least one argument between continent, subregion and admin...")

    nt = (; continent, subregion, admin)
    kwargs = (k => v for (k, v) in pairs(nt) if !isempty(v))
    d = CountriesBorders.extract_countries(CountriesBorders.get_geotable(; resolution); skip_areas, kwargs...)
    cart = convexhull(d) # Using convexhull() method from CountriesBorders. Convexhull always give a PolyArea.
    latlon = latlon_geometry(cart)
    ch = PolyBorder(latlon, cart)
    name = isempty(name) ? (!isempty(admin) ? admin : "region_name") : name
    GeoRegion(name, continent, subregion, admin, d, ch)
end

"""
    PolyRegion{P} <: AbstractRegion

Type representing a polygonal region based on PolyArea.

Fields:
- `name::String`: Name of the region
- `domain::PolyBorder{P}`: Domain of the region as a PolyBorder

Where `P` is the precision type for coordinates.
"""
mutable struct PolyRegion{P} <: AbstractRegion
    name::String
    domain::PolyBorder{P}
end

PolyRegion(name, domain::Vector{<:LatLon}) = PolyRegion(name, PolyBorder(PolyArea(map(Point, domain))))
PolyRegion(; name::String="region_name", domain) = PolyRegion(name, domain)

"""
    LatBeltRegion <: AbstractRegion

Type representing a latitude belt region.

Fields:
- `name::String`: Name of the region
- `lim::Tuple{ValidAngle,ValidAngle}`: Latitude limits of the belt in degrees
"""
mutable struct LatBeltRegion <: AbstractRegion
    name::String
    lim::Tuple{Deg{Float64},Deg{Float64}} # [°]

    function LatBeltRegion(name::String, lim::Tuple{ValidAngle,ValidAngle})
        # Inputs validation    
        _lim = map(lim) do l
            l isa Real ? l * ° : l |> u"°" # Convert to Unitful °
        end

        for x in _lim
            abs(x) ≤ 90° || error(
#! format: off
"LAT provided as numbers must be expressed in degrees and satisfy -90 ≤ x ≤ 90. 
Consider using `°` (or `rad`) from `Unitful` if you want to pass numbers in degrees (or rad), by doing `x * °` (or `x * rad`)."
#! format: on   
            )
        end
        _lim[1] > _lim[2] && error("The first LAT limit must be lower than the second one...")
        _lim[1] == _lim[2] && error("The first LAT limit must be different than the second one...")

        new(name, _lim)
    end
end

LatBeltRegion(; name::String="region_name", lim) = LatBeltRegion(name, lim)


# HotSpotRegion <: AbstractRegion
"""
    HotSpotRegion{P} <: AbstractRegion

Type representing a hot spot region, which is defined as the set of points which are less than `radius` distance from a given `center` in LatLon CRS.

Fields:
- `name::String`: Name identifying the hot spot
- `center::POINT_LATLON{P}`: Center of the hot spot
- `radius::Float64`: Radius of the hot spot [m]
- `domain::PolyBorder{P}`: Polygon identifying the polygon in latlon which represents a _circle_ of specified radius from the center

# Constructor
    HotSpotRegion(; name::String, center::Union{LATLON, POINT_LATLON}, radius::Number)

Create a HotSpotRegion from a name, a center and a radius.
"""
struct HotSpotRegion{P} <: AbstractRegion
    "Name identifying the hot spot"
    name::String
    "Center of the hot spot"
    center::POINT_LATLON{P}
    "Radius of the hot spot [m]"
    radius::Float64
    "Polygon identifying the polygon in latlon which represents a _circle_ of specified radius from the center"
    domain::PolyBorder{P}
end
function HotSpotRegion(; name::String, center::Union{LATLON, POINT_LATLON}, radius::Number)
    c = center isa LATLON ? Point(center) : center
    circ_poly = gen_circle_pattern(c, radius; n = 51) |> only
    domain = PolyBorder(circ_poly[1:end-1] |> PolyArea)
    return HotSpotRegion{floattype(c)}(name, c, radius, domain)
end

"""
    MultiRegion{P} <: AbstractRegion

Type representing a region which is defined as the union of multiple PolyAreas.

Fields:
- `name::String`: Name identifying the multi region
- `domain::MultiBorder{P}`: MultiPolygon identifying the various regions included in the multi region

# Constructor
    MultiRegion(areas::Vector; name::String)

The constructor takes a vector of polyareas, Multi, or other AbstractRegions and creates a single MultiBorder which encompasses all the polyareas of the provided `areas`.
"""
@kwdef struct MultiRegion{P} <: AbstractRegion
    "Name identifying the multi region"
    name::String
    "List of Polygons identifying the various regions included in the multi region"
    domain::MultiBorder{P}
end
function MultiRegion(areas::Vector; name::String)
    # We convert all polygons to Float32 machine precision
    f(reg) = Iterators.map(change_floattype(Float32), polyareas(reg))
    f(reg::Union{POLY_LATLON, MULTI_LATLON}) = cartesian_geometry(reg) |> f

    multi_cart = Iterators.flatten((f(area) for area in areas)) |> Multi
    multi_latlon = latlon_geometry(multi_cart)
    domain = MultiBorder(multi_latlon, multi_cart)
    return MultiRegion(name, domain)
end

## Define Tessellation Types
abstract type AbstractTiling end

"""
    ICO <: AbstractTiling

Struct representing an icosahedral tiling.

Fields:
- `correction::Number`: Default correction factor for the icosahedral cell grid partial overlap
- `pattern::Symbol`: Default pattern shape to be used with this type of tiling (:circ or :hex)
"""
struct ICO <: AbstractTiling 
    correction::Number
    pattern::Symbol

    function ICO(correction::Number, pattern::Symbol)
        pattern in (:circ, :hex) || error("Pattern must be :circ or :hex")
        new(correction, pattern)
    end
end

ICO(; correction::Number=3/2, pattern::Symbol=:circ) = ICO(correction, pattern)

"""
    HEX <: AbstractTiling

Struct representing a hexagonal tiling.

Fields:
- `direction::Symbol`: Default direction of hexagons in the tiling (:pointy or :flat)
- `pattern::Symbol`: Default pattern shape to be used with this type of tiling (:circ or :hex)
"""
struct HEX <: AbstractTiling 
    direction::Symbol
    pattern::Symbol

    function HEX(direction::Symbol, pattern::Symbol)
        direction in (:pointy, :flat) || error("Direction must be :pointy or :flat")
        pattern in (:circ, :hex) || error("Pattern must be :circ or :hex")
        new(direction, pattern)
    end
end

HEX(; direction::Symbol=:pointy, pattern::Symbol=:hex) = HEX(direction, pattern)

"""
    H3 <: AbstractTiling

Struct representing an H3 tiling.
"""
struct H3 <: AbstractTiling end

"""
    EO

Struct used to create function methods that return more than one output.
Used within multiple methods of the GeoGrids API, usually given as last optional argument.
"""
struct EO end