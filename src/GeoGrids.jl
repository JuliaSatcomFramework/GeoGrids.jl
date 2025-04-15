module GeoGrids

using AngleBetweenVectors: AngleBetweenVectors
using Clipper: Clipper, ClipperOffset, EndTypeClosedPolygon, IntPoint,
    JoinTypeMiter, add_path!, execute, tofloat
using CoordRefSystems: CoordRefSystems, Deg, RevolutionEllipsoid, Cartesian, 
    CoordRefSystems, Datum
using CountriesBorders: CountriesBorders, borders, cartesian_geometry, 
    latlon_geometry, change_geometry, bboxes, polyareas, in_exit_early, 
    to_cart_point, VALID_CRS, floattype, to_latlon_point, LATLON, 
    POLY_LATLON, POLY_CART, MULTI_LATLON, MULTI_CART, RING_LATLON, 
    POINT_LATLON, BOX_CART, BOX_LATLON, LatLon, SKIP_NONCONTINENTAL_EU, 
    extract_countries
using Dictionaries: Dictionaries, Dictionary
using Distances: Distances, Metric, result_type
using GeoPlottingHelpers: GeoPlottingHelpers, 
    to_raw_lonlat, geom_iterable
using LinearAlgebra: LinearAlgebra, convert, norm
using Meshes: Meshes, 🌐, WGS84Latest, coords, Meshes, Box, CCW, 
    Geometry, Multi, Point, PointSet, PolyArea, Ring, SimpleMesh, 
    TesselationMethod, VoronoiTesselation, boundingbox, centroid, 
    convexhull, crs, direction, radius, rings, segments, tesselate, 
    vertices
using NearestNeighbors: NearestNeighbors, always_false, check_k, 
    knn_point!, inrange_point!, check_radius, BallTree, Euclidean, 
    KDTree, NNTree, inrange!, knn!
using PlotlyExtensionsHelper: PlotlyExtensionsHelper
using SimpleWeightedGraphs: SimpleWeightedGraphs, SimpleWeightedGraph
using StaticArrays: StaticArrays, SA_F64, SVector
using Unitful: °, rad, Quantity, Length, @u_str, ustrip
using Graphs: reverse

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