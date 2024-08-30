## Angle Types
const UnitfulAngleType = Union{typeof(°),typeof(rad)}
const UnitfulAngleQuantity = Quantity{<:Real,<:Any,<:UnitfulAngleType}
const ValidAngle = Union{UnitfulAngleQuantity,Real}

const constants = (
    Re_mean = 6371e3, # Mean Earth Radius [m]
    a = 6378137, # [m] WGS84 semi-major axis
    b = 6356752.315 # [m] WGS84 semi-minor axis
)

"""
    PolyBorder <: Geometry{🌐,LATLON}

Struct representing a PolyArea in both LatLon and Cartesian coordinates.

Fields:
- `latlon::POLY_LATLON`: The borders in LatLon CRS
- `cart::POLY_CART`: The borders in Cartesian2D CRS
"""
struct PolyBorder{P} <: Geometry{🌐,LATLON} # Use parameetric precision (e.g., Float32, Float64) for the coordinates.
    latlon::POLY_LATLON{P}
    cart::POLY_CART{P}
end
function PolyBorder(latlon::POLY_LATLON) 
    cart = cartesian_geometry(latlon)
    PolyBorder(latlon, cart)
end

"""
    MultiBorder <: Geometry{🌐,LATLON}

Struct representing a Multi in both LatLon and Cartesian coordinates.

Fields:
- `latlon::MULTI_LATLON`: The borders in LatLon CRS
- `cart::MULTI_CART`: The borders in Cartesian2D CRS
"""
struct MultiBorder{P} <: Geometry{🌐,LATLON} # Use parameetric precision (e.g., Float32, Float64) for the coordinates.
    latlon::MULTI_LATLON{P}
    cart::MULTI_CART{P}
end
function MultiBorder(latlon::MULTI_LATLON) 
    cart = cartesian_geometry(latlon)
    MultiBorder(latlon, cart)
end

abstract type AbstractRegion end

"""
    GlobalRegion

Type representing a global region.
"""
struct GlobalRegion <: AbstractRegion end

"""
    GeoRegion{D} <: AbstractRegion

Type representing a geographical region based on CountriesBorders.

Fields:
- `name::String`: Name of the region
- `continent::String`: Continent of the region
- `subregion::String`: Subregion within the continent
- `admin::String`: Administrative area
- `domain::D`: Domain of the region
- `convexhull::PolyArea`: Convex hull of the region
"""
mutable struct GeoRegion{D,P} <: AbstractRegion
    name::String
    continent::String
    subregion::String
    admin::String
    domain::D
    convexhull::PolyBorder{P}
end
function GeoRegion(; name="region_name", continent="", subregion="", admin="")
    all(isempty(v) for v in (continent, subregion, admin)) && error("Input at least one argument between continent, subregion and admin...")

    nt = (; continent, subregion, admin)
    kwargs = (k => v for (k, v) in pairs(nt) if !isempty(v))
    d = CountriesBorders.extract_countries(; kwargs...)
    cart = convexhull(d) # Using convexhull() method from CountriesBorders. Convexhull always give a PolyArea.
    latlon = latlon_geometry(cart)
    ch = PolyBorder(latlon, cart)
    GeoRegion(name, continent, subregion, admin, d, ch)
end

"""
    PolyRegion <: AbstractRegion

Type representing a polygonal region based on PolyArea.

Fields:
- `name::String`: Name of the region
- `domain::PolyBorder`: Domain of the region as a PolyBorder
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
- `lim::Tuple{ValidAngle,ValidAngle}`: Latitude limits of the belt in radians
"""
mutable struct LatBeltRegion <: AbstractRegion
    name::String
    lim::Tuple{ValidAngle,ValidAngle} # [rad] 

    function LatBeltRegion(name::String, lim::Tuple{ValidAngle,ValidAngle})
        # Inputs validation    
        _lim = map(lim) do l
            l isa Real ? l * ° : l |> u"°" # Convert to Uniful °
        end

        for x in _lim
            abs(x) ≤ 90° || error(
#! format: off
"LAT provided as numbers must be expressed in radians and satisfy -90 ≤ x ≤ 90. 
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

## Define Tessellation Types
abstract type AbstractTiling end

"""
    ICO <: AbstractTiling

Struct representing an icosahedral tiling.

Fields:
- `correction::Number`: Default correction factor for the icosahedral cell grid partial overlap
- `pattern::Symbol`: Default pattern shape to be used with this type of tiling
"""
struct ICO <: AbstractTiling 
    "Default correction factor for the icosahedral cell grid partial overlap"
    correction::Number
    "Default pattern shape to be used with this type of tiling"
    pattern::Symbol

    function ICO(correction::Number, pattern::Symbol)
        # Input validation
        pattern in (:circ, :hex) || error("Pattern must be :circ or :hex...")

        new(correction, pattern)
    end
end
ICO(; correction::Number=3/2, pattern::Symbol=:circ) = ICO(correction, pattern)

"""
    HEX <: AbstractTiling

Struct representing a hexagonal tiling.

Fields:
- `direction::Symbol`: Default direction of hexagons in the tiling
- `pattern::Symbol`: Default pattern shape to be used with this type of tiling
"""
struct HEX <: AbstractTiling 
    "Default direction of hexagons in the tiling"
    direction::Symbol
    "Default pattern shape to be used with this type of tiling"
    pattern::Symbol

    function HEX(direction::Symbol, pattern::Symbol)
        # Input validation
        direction in (:pointy, :flat) || error("Direction must be :pointy or :flat...")
        pattern in (:circ, :hex) || error("Pattern must be :circ or :hex...")

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