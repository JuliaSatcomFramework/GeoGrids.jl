# This is not exactly just POINT_LATLON as it allows also different datums
const VALID_COORD = Union{Point{🌐, <:LatLon}, LatLon}

"""
    LatLonTree{T <: NNTree, P <: Union{LatLon, Point{🌐, <:LatLon}}}

Structure providing a modified `NNTree` from NearestNeighbors.jl that allows querying with Points in LatLon coordinates.

The LatLonTree internally stores an `NNTree` where the provided points are converted to SVector{2, Float64} and store the longitude and latitude values of the provided points in radians.

While the structure is defined in GeoGrids.jl, its constructor and its usefulness is only realized through an extension when loading the NearestNeighbors.jl package.

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
struct LatLonTree{T, P <: POINT_LATLON}
    tree::T
    points::Vector{P}
end

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
    @eval $f(v::Vector{<:VALID_COORD}) = $f(eltype(v))
end

mactype(x) = CoordRefSystems.mactype(x)
mactype(P::Type{<:Point{🌐, <:LatLon}}) = mactype(P |> crs)
mactype(p::VALID_COORD) = mactype(typeof(p))
mactype(v::Vector{<:VALID_COORD}) = mactype(eltype(v))

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
GreatCircleMetric(p::VALID_COORD) = GreatCircleMetric{mactype(p)}(authalic_radius(p))

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