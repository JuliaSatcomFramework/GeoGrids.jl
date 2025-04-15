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
using Graphs: Graphs, reverse

include("basic/basic_types.jl")
export AbstractRegion, GlobalRegion, GeoRegion, PolyRegion, LatBeltRegion,
    MultiBorder, PolyBorder, BoxBorder, HotSpotRegion, MultiRegion, 
    ClippedRegion, AbstractTiling, ICO, HEX, H3, EO

include("basic/offset_types.jl")
export GeoRegionOffset, PolyRegionOffset

include("misc/distances.jl")
export LatLonTree, GreatCircleMetric, great_circle_distance

include("misc/interface_func.jl")
export extract_countries, borders, centroid, in, get_lat, get_lon

include("basic/offsetting_func.jl")
export offset_region

include("grids/ico_func.jl")
export icogrid

include("grids/rect_func.jl")
export rectgrid, vecgrid

include("grids/filtering_func.jl")
export filter_points, group_by_domain

include("grids/tessellation_func.jl")
export gen_hex_lattice, generate_tesselation, gen_circle_pattern, gen_hex_pattern

include("misc/coloring.jl")
export color_greedy

# Export from other modules
export °, rad, ustrip, LatLon, Cartesian, WGS84Latest, coords, 
    PolyArea, Point, SVector, SKIP_NONCONTINENTAL_EU, latlon_geometry, 
    cartesian_geometry

include("misc/helper_func.jl")
include("misc/plot_func.jl")

end # module GeoGrids 