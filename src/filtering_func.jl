"""
    in_domain(p::LLA, domain::Union{GeometrySet,PolyArea})
    in_domain(p::LLA, domain::Union{GeoRegion, PolyRegion}) = in_domain(p, domain.domain)
    in_domain(p::Union{Tuple{Float64, Float64},SVector{2,Float64}}, domain::Union{GeometrySet,PolyArea})
    in_domain(p::Union{Tuple{Float64, Float64},SVector{2,Float64}}, domain::Union{GeoRegion, PolyRegion}) = in_domain(p, domain.domain)
    in_domain(p::Point2, domain::Union{GeometrySet,PolyArea})
    in_domain(p::Point2, domain::Union{GeoRegion, PolyRegion}) = in_domain(p, domain.domain)

This function determines if a given point belongs to a 2-dimensional `Meshes.Domain` object. The `Meshes.Domain` object represents a geometric domain, which is essentially a 2D region in space, specified by its bounds and discretization. 

The function first converts the input tuple into a `Meshes.Point` object, which is then checked if it falls inside the given `Meshes.Domain` object.
The `Meshes.Domain` can be either a `GeometrySet` or a `PolyArea` object.

## Arguments
* `p`: A tuple of two numbers `(x,y)` representing a point in 2D space. 
* `domain::Union{GeometrySet,PolyArea}`: A `Meshes.Domain` object representing a 2D region in space. 

### Output
The function returns a boolean value: `true` if the point represented by the input tuple falls inside the `Meshes.Domain` object and `false` otherwise. 
"""
function in_domain(p::LLA, domain::Union{GeometrySet,PolyArea})
    _p = (rad2deg(p.lon), rad2deg(p.lat))
    Meshes.Point2(_p) in domain # Meshes.Point in Meshes.Geometry
end
in_domain(p::LLA, domain::Union{GeoRegion, PolyRegion}) = in_domain(p, domain.domain)

function in_domain(p::Union{Tuple{Float64, Float64},SVector{2,Float64}}, domain::Union{GeometrySet,PolyArea})
    _p = _point_check(p) # Input check
	Meshes.Point2(_p) in domain
end
in_domain(p::Union{Tuple{Float64, Float64},SVector{2,Float64}}, domain::Union{GeoRegion, PolyRegion}) = in_domain(p, domain.domain)

function in_domain(p::Point2, domain::Union{GeometrySet,PolyArea})
    _p = _point_check(p.coords) # Input check
	return Meshes.Point2(_p) in domain
end
in_domain(p::Point2, domain::Union{GeoRegion, PolyRegion}) = in_domain(p, domain.domain)

"""
    filter_points(points::Union{Vector{LLA}, Vector{SVector{2,Float64}}, Vector{Point2}, Vector{Tuple{Float64,Float64}}}, domain::Union{GeoRegion, PolyRegion})
    filter_points(points::Union{Vector{LLA}, Vector{SVector{2,Float64}}, Vector{Point2}, Vector{Tuple{Float64,Float64}}}, domain::Union{GeometrySet, PolyArea})
    
Filters a list of points based on whether they fall within a specified geographical domain.

## Arguments
- `points`: A vector of points. The points can be of type `LLA`, `SVector{2,Float64}`, `Point2`, or `Tuple{Float64,Float64}`.
- `domain`: A geographical domain which can be of type `GeoRegion` or `PolyRegion`, in alternative a `Meshes.Domain` of type `GeometrySet` or `PolyArea`.

## Returns
- A vector of points that fall within the specified domain.
"""
function filter_points(points::Union{Vector{LLA}, Vector{SVector{2,Float64}}, Vector{Point2}, Vector{Tuple{Float64,Float64}}}, domain::Union{GeoRegion, PolyRegion})
    mask = map(x -> in_domain(x, domain.domain), points)
    return points[mask]
end

function filter_points(points::Union{Vector{LLA}, Vector{SVector{2,Float64}}, Vector{Point2}, Vector{Tuple{Float64,Float64}}}, domain::Union{GeometrySet, PolyArea})
    mask = map(x -> in_domain(x, domain), points)
    return points[mask]
end