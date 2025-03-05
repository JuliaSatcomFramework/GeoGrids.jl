module NearestNeighborsExt

using NearestNeighbors
using NearestNeighbors: always_false, check_k, knn_point!, inrange_point!, check_radius
using GeoGrids
using GeoGrids.Distances
using GeoGrids.StaticArrays
using GeoGrids.CoordRefSystems
using GeoGrids.Meshes
using GeoGrids: LatLonTree, VALID_COORD, authalic_radius, raw_svector, GreatCircleMetric, mactype

function GeoGrids.LatLonTree{BallTree}(points::AbstractVector{<:VALID_COORD}, metric::Metric = GreatCircleMetric(first(points)); kwargs...) 
    normalized_points = map(raw_svector, points)
    tree = BallTree(normalized_points, metric; kwargs...)
    return LatLonTree{typeof(tree), eltype(points)}(tree, points)
end
function GeoGrids.LatLonTree{KDTree}(points::AbstractVector{<:VALID_COORD}, metric::Metric = Euclidean(); kwargs...) 
    normalized_points = map(raw_svector, points)
    tree = KDTree(normalized_points, metric; kwargs...)
    return LatLonTree{typeof(tree), eltype(points)}(tree, points)
end

GeoGrids.LatLonTree(points::AbstractVector{<:VALID_COORD}; kwargs...) = LatLonTree{BallTree}(points; kwargs...)

# Single point
"""
    knn(lltree::LatLonTree, points::Vector, k::Int [, sortres=false]) -> indices, distances

Performs a lookup of the `k` nearest neigbours to the `points` from the data
in the `lltree`. `points` can either be a single coordinate expressed in `LatLon` or `Point{🌐, <:LatLon}` or a vector of such coordinates. `skip` is an optional predicate
to determine if a point that would be returned should be skipped based on its
index.

See also: `knn!`, `nn`.
"""
function NearestNeighbors.knn(lltree::LatLonTree, point::VALID_COORD, k::Int, sortres=false, skip::F=always_false) where {F<:Function}
    idx = Vector{Int}(undef, k)
    (; tree) = lltree
    _T = mactype(point)
    T = result_type(tree.metric, _T, _T)
    dist = Vector{T}(undef, k)
    point = raw_svector(point)
    return knn!(idx, dist, tree, point, k, sortres, skip)
end
# Vector of points
function NearestNeighbors.knn(lltree::LatLonTree, points::AbstractVector{<:VALID_COORD}, k::Int, sortres=false, skip::F=always_false) where {F<:Function}
    (; tree) = lltree
    check_k(tree, k)
    n_points = length(points)
    _T = mactype(eltype(points))
    T = result_type(tree.metric, _T, _T)
    dists = [Vector{T}(undef, k) for _ in 1:n_points]
    idxs = [Vector{Int}(undef, k) for _ in 1:n_points]
    for i in 1:n_points
        knn_point!(tree, raw_svector(points[i]), sortres, dists[i], idxs[i], skip)
    end
    return idxs, dists
end
# Version with provided vectors to store indices and distances for a point
function NearestNeighbors.knn!(idxs::AbstractVector{<:Integer}, dists::AbstractVector, lltree::LatLonTree, point::VALID_COORD, k::Int, sortres=false, skip::F=always_false) where {F<:Function}
    (; tree) = lltree
    check_k(tree, k)
    length(idxs) == k || throw(ArgumentError("idxs must be of length k"))
    length(dists) == k || throw(ArgumentError("dists must be of length k"))
    knn_point!(tree, raw_svector(point), sortres, dists, idxs, skip)
    return idxs, dists
end

# Check for range
function NearestNeighbors.inrange!(idxs::AbstractVector, lltree::LatLonTree{<:BallTree}, point::VALID_COORD, radius::Number, sortres=false)
    (; tree) = lltree
    check_radius(radius)
    length(idxs) == 0 || throw(ArgumentError("idxs must be empty"))
    inrange_point!(tree, raw_svector(point), radius, sortres, idxs)
    return idxs
end
function NearestNeighbors.inrange(lltree::LatLonTree{<:BallTree}, point::VALID_COORD, radius::Number, sortres=false)
    return inrange!(Int[], lltree, point, radius, sortres)
end

# Version with multiple points
function NearestNeighbors.inrange(lltree::LatLonTree{<:BallTree},
                 points::AbstractVector{<:VALID_COORD},
                 radius::Number,
                 sortres=false)
    (; tree) = lltree
    check_radius(radius)

    idxs = [Vector{Int}() for _ in 1:length(points)]

    for i in 1:length(points)
        inrange_point!(tree, raw_svector(points[i]), radius, sortres, idxs[i])
    end
    return idxs
end

function NearestNeighbors.inrangecount(lltree::LatLonTree{<:BallTree}, point::VALID_COORD, radius::Number)
    (; tree) = lltree
    check_radius(radius)
    return inrange_point!(tree, raw_svector(point), radius, false, nothing)
end

function NearestNeighbors.inrangecount(lltree::LatLonTree{<:BallTree},
        points::AbstractVector{<:VALID_COORD},
        radius::Number)
    (; tree) = lltree
    check_radius(radius)
    out= Vector{Int}(undef, length(points))
    for i in eachindex(out, points)
        out[i] = inrange_point!(tree, raw_svector(points[i]), radius, false, nothing)
    end
    return out
end

end
