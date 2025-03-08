# This is not exactly just POINT_LATLON as it allows also different datums
const VALID_COORD = Union{Point{🌐, <:LatLon}, LatLon}


# Authalic radius computation based on https://en.wikipedia.org/wiki/Earth_radius#Authalic_radius
function authalic_radius(a, b, e)
    lg = log((1 + e) / (b/ a))
    inner = a^2 + b^2 / e * lg
    return sqrt(inner/2)
end
function authalic_radius(::Type{🌎}) where 🌎 <: RevolutionEllipsoid
    (; a, b, e) = CoordRefSystems.ellipsoidparams(🌎)
    return authalic_radius(ustrip(a), ustrip(b), e)
end

for f in (:authalic_radius,)
    @eval $f(D::Type{<:Datum}) = $f(CoordRefSystems.ellipsoid(D))
    @eval $f(L::Type{<:LatLon}) = $f(CoordRefSystems.datum(L))
    @eval $f(P::Type{<:Point{🌐, <:LatLon}}) = $f(Meshes.crs(P))
    @eval $f(p::VALID_COORD) = $f(typeof(p))
end

# Create a SVector with lon and lat as first and second elements translated in radians
function raw_svector(p::LatLon)
    (; lat, lon) = p
    f(x) = ustrip(x) |> deg2rad
    SVector(f(lon), f(lat))
end
raw_svector(p::Point{🌐, <:LatLon}) = raw_svector(coords(p))

# Now we create our custom distance based on Haversine but working on radians or directly on LatLon or Point{🌐, <:LatLon} values
"""
    GreatCircleMetric(radius=6_371_000)

The great circle distance between two points on a sphere of given `radius`. This is basically equivalent to `Haversine` from Distances.jl but allows computing distance between points already expressed in `LatLon` or `Point{🌐, <:LatLon}` types.
"""
struct GreatCircleMetric{T<:Number} <: Metric
    radius::T
end
GreatCircleMetric{T}() where {T<:Number} = GreatCircleMetric(T(6_371_000))
GreatCircleMetric() = GreatCircleMetric{Float64}()
GreatCircleMetric(p::VALID_COORD) = GreatCircleMetric{floattype(p)}(authalic_radius(p))

function (dist::GreatCircleMetric)(x, y)
    length(x) == length(y) == 2 || throw(ArgumentError("x and y must have length 2 for GreatCircleMetric distance computation"))
    @inbounds λ₁, φ₁ = x
    @inbounds λ₂, φ₂ = y

    Δλ = λ₂ - λ₁  # longitudes
    Δφ = φ₂ - φ₁  # latitudes

    # haversine formula
    a = sin(Δφ/2)^2 + cos(φ₁)*cos(φ₂)*sin(Δλ/2)^2

    # distance on the sphere
    2 * (dist.radius * asin( min(√a, one(a)) )) # take care of floating point errors
end
(dist::GreatCircleMetric)(x::VALID_COORD, y::VALID_COORD) = dist(raw_svector(x), raw_svector(y))

Distances.result_type(::GreatCircleMetric{T1}, ::Type{T2}, ::Type{T3}) where {T1<:Number,T2<:Number,T3<:Number} =
    float(promote_type(T1, T2, T3))

great_circle_distance(p1::VALID_COORD, p2::VALID_COORD, r = authalic_radius(p1)) = GreatCircleMetric(r)(p1, p2)

"""
    LatLonTree{T <: NNTree, P <: Union{LatLon, Point{🌐, <:LatLon}}}

Structure providing a modified `NNTree` from NearestNeighbors.jl that allows querying with Points in LatLon coordinates.

The LatLonTree internally stores an `NNTree` where the provided points are converted to SVector{2, Float64} and store the longitude and latitude values of the provided points in radians.

# Constructors
    LatLonTree{BallTree}(points::AbstractVector, metric::Metric = GreatCircleMetric(first(points)); kwargs...)
    LatLonTree{KDTree}(points::AbstractVector, metric::Metric = Euclidean(); kwargs...)
    LatLonTree(points::AbstractVector; kwargs...)

The first two constructors explicitly specify the type of `NNTree` to use internally, and allow customizing the metric used while the latter method simply fall backs to the `BallTree` with default metric.

The default preferred tree structure is the `BallTree` one as it computes distance between points using the Great Circle distance assuming spherical earth, which is a good approximation for actual distances between points on the surface of the earth (and much faster than the actual inverse geodesic computation)

Instances of LatLonTree support the following functions from NearestNeighbors.jl:
- `knn`
- `knn!`
- `inrangecount`
- `inrange`
- `inrange!`

All these methods allow specifing the query points as one of the following:
- A `LatLon` coordinate
- A `Point{🌐, <:LatLon}` point
- A vector with either of the two previous types as `eltype`

The `inrangecount`, `inrange` and `inrange!` methods are only callable with instances of `LatLonTree{<:BallTree}` as the distance with `KDTree` is not really meaningful.
"""
struct LatLonTree{T <: NNTree, P <: VALID_COORD}
    tree::T
    points::Vector{P}
end

function LatLonTree{BallTree}(points::AbstractVector{<:VALID_COORD}, metric::Metric = GreatCircleMetric(first(points)); kwargs...) 
    normalized_points = map(raw_svector, points)
    tree = BallTree(normalized_points, metric; kwargs...)
    return LatLonTree{typeof(tree), eltype(points)}(tree, points)
end
function LatLonTree{KDTree}(points::AbstractVector{<:VALID_COORD}, metric::Metric = Euclidean(); kwargs...) 
    normalized_points = map(raw_svector, points)
    tree = KDTree(normalized_points, metric; kwargs...)
    return LatLonTree{typeof(tree), eltype(points)}(tree, points)
end

LatLonTree(points::AbstractVector{<:VALID_COORD}; kwargs...) = LatLonTree{BallTree}(points; kwargs...)

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
    _T = floattype(point)
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
    _T = floattype(eltype(points))
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