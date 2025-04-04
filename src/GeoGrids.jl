module GeoGrids

using AngleBetweenVectors
using Clipper
using CoordRefSystems
using CoordRefSystems: Deg, RevolutionEllipsoid, Cartesian
using CountriesBorders
using CountriesBorders: borders, cartesian_geometry, latlon_geometry, change_geometry, bboxes, polyareas, in_exit_early, to_cart_point, VALID_CRS, floattype, to_latlon_point
using CountriesBorders: LATLON, CART, POLY_LATLON, POLY_CART, MULTI_LATLON, MULTI_CART, RING_LATLON, RING_CART, POINT_LATLON, BOX_CART, DOMAIN, BOX_LATLON
using Dictionaries
using Distances: Distances, Metric, result_type
using GeoPlottingHelpers: GeoPlottingHelpers, extract_latlon_coords, to_raw_lonlat, geom_iterable
using LinearAlgebra
using Meshes
using Meshes: 🌐, WGS84Latest, coords
using NearestNeighbors
using NearestNeighbors: always_false, check_k, knn_point!, inrange_point!, check_radius
using PlotlyExtensionsHelper
using StaticArrays
using Unitful: °, rad, Quantity, Length, @u_str, ustrip
import Graphs
using SimpleWeightedGraphs

include("basic_types.jl")
include("distances.jl")
include("offset_types.jl")
include("interface_func.jl")
include("helper_func.jl")
include("offsetting_func.jl")
include("ico_func.jl")
include("rect_func.jl")
include("filtering_func.jl")
include("tessellation_func.jl")
include("plot_func.jl")

export AbstractRegion, GlobalRegion, GeoRegion, PolyRegion, LatBeltRegion,
    GeoRegionOffset, PolyRegionOffset,
    MultiBorder, PolyBorder, BoxBorder,
    HotSpotRegion, MultiRegion, ClippedRegion,
    AbstractTiling, ICO, HEX, H3,
    EO
    # UnitfulAngleType, UnitfulAngleQuantity, ValidAngle,

export LatLonTree, GreatCircleMetric, great_circle_distance
    
export icogrid, rectgrid, vecgrid,
    extract_countries, SKIP_NONCONTINENTAL_EU, filter_points, group_by_domain,
    gen_hex_lattice, generate_tesselation, _tesselate, gen_circle_pattern, gen_hex_pattern,
    borders, centroid, in, get_lat, get_lon, latlon_geometry, cartesian_geometry,
    offset_region

export °, rad, ustrip,
    LatLon, Cartesian, WGS84Latest, coords, PolyArea, Point,
    SVector

include("coloring.jl")
export color_greedy

end # module GeoGrids 

# //TODO:
# [x] Add tests for enlarged types and offsetting functions.
# [x] Add tests for plotting functions for enlarged regions.
# [x] Add documetnation on README for region enlargement.
# [x] Add interface functions for in for offset regions
# [x] Add interface functions for centroid for offset regions
# [x] Add interface functions for borders for offset regions (?)
# [x] Add filtering and grouping functions for enlarged regions 
# [x] Merge PR.
# [x] Version bump to 0.5.0
# [x] Add tessellation for offset regions
# [] Add ValidDistance support (Issue)
# [] Add tests for all the functionalities added for offset regions