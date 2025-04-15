"""
    offset_region(originalRegion::GeoRegion, deltaDist; refRadius=constants.Re_mean, magnitude=3, precision=7)
    offset_region(originalRegion::PolyRegion, deltaDist; refRadius=constants.Re_mean, magnitude=3, precision=7)

Offset a GeoRegion or PolyRegion by a given distance. This function offsets each polygon in
the region separately and combines the results into a Multi geometry.

## Arguments
- `originalRegion::Union{GeoRegion,PolyRegion}`: The original region to be offset.
- `deltaDist`: The distance to offset the region by, in meters. Positive for \
enlargement, negative for shrinking.
- `refRadius::Float64=constants.Re_mean`: The reference radius to use for the \
Earth.
- `magnitude::Int=3`: The number of integer digits for IntPoint conversion.
- `precision::Int=7`: The total number of digits to be considered for each \
coordinate in IntPoint conversion.

## Returns
- `Multi`: A Multi geometry containing the offset polygons.

## Notes
- For GeoRegion, only the outer ring of each polygon is considered for offsetting.
- For PolyRegion, if multiple outer rings are produced, inner rings are ignored \
and separate PolyAreas are created for each outer ring.
"""
function offset_region(originalRegion, deltaDist; refRadius=constants.Re_mean, magnitude=3, precision=7)
    # `magnitude` represents the number of integer digits while `precision` the
    # total number of digits that will be considered for each of the coordinates
    # for the `IntPoint` conversion. Look at Clipper documentation for more
    # details. `delta` is the distance to offset the polygon by, it is a
    # positive value for enlargement or negative number for shrinking. The value
    # should be expressed in m.

    # Compute the delta value to be used in the offsetting process (in deg)
    delta = rad2deg(deltaDist / refRadius)
    intDelta = Float64(IntPoint(delta, delta, magnitude, precision).X) # We use IntPoint to exploit the conversion to IntPoint in Clipping, then we can use either X or Y as delta value. Clipper wants Float64

    # We directly offset every polyarea in the region
    polys = map(polyareas(originalRegion)) do polyarea
        # Only outer ring (i.e., the first is considered for the enlargement of the
        # GeoRegion. We avoid considering the inner rings since they are not
        # relevant for our application.
        outerRing = rings(latlon_geometry(polyarea)) |> first
        offsetRings = _offset_ring(outerRing, intDelta; magnitude, precision)
        # Create a separate PolyArea for each of the offset rings.
        map(offsetRings) do ring
            PolyArea(ring)
        end
    end |> splat(vcat)

    return Multi(polys)
end

"""
    _offset_ring(ring::Ring{🌐,<:LatLon{WGS84Latest}}, delta; magnitude=3, precision=7)

Offset a ring by a given delta value. This function uses the Clipper library
for polygon offsetting. It may return multiple rings even when starting from
a single Ring.

## Arguments
- `ring::Ring{🌐,<:LatLon{WGS84Latest}}`: The ring to be offset.
- `delta`: The distance to offset the ring by. Positive for enlargement, \
negative for shrinking.
- `magnitude::Int=3`: The number of integer digits for IntPoint conversion.
- `precision::Int=7`: The total number of digits to be considered for each \
coordinate in IntPoint conversion.

## Returns
- Vector of `Ring{🌐,<:LatLon{WGS84Latest}}`: The resulting offset rings.
"""
function _offset_ring(ring::RING_LATLON{T}, delta; magnitude=3, precision=7) where {T}
    # delta translated in deg wrt the Earrth radius
    intPoly = map([vertices(ring)...]) do vertex # Use splat to avoid CircularVector as output from map
        y = get_lat(vertex) |> ustrip
        x = get_lon(vertex) |> ustrip
        IntPoint(x, y, magnitude, precision) # Consider LON as X and LAT as Y
    end
    co = ClipperOffset()
    add_path!(co, intPoly, JoinTypeMiter, EndTypeClosedPolygon) # We fix JoinTypeMiter, EndTypeClosedPolygon because it works well with complex polygons, look at Clipper documentation for details.
    offset_polygons = execute(co, delta) # Clipper polygon
    # We can end up with multiple polygons even when starting from a single
    # Ring. So we will return all the polygons generate for this Ring as A#
    # a vector of Rings. Afterwards they will be used as outer or inner rings.
    outRings = map(eachindex(offset_polygons)) do i
        map(offset_polygons[i]) do vertex
            lonlat = tofloat(vertex, magnitude, precision)
            # We force latitude to be in the range [-90°, 90°]
            lat = min(max(lonlat[2] |> T, -90), 90)
            lon = lonlat[1] |> T
            LatLon{WGS84Latest}(lat, lon) |> Point # Use consistent type precision for the coordinates.
        end |> Ring
    end

    # We now eventually split the rings if they contain antimeridian crossings
    split_antimeridian!(outRings)

    # Return a vector of Ring.
    return outRings
end

# Check if a ring contains a crossing of the antimeridian (180° longitude)
function has_antimeridian(ring)
    any(segments(ring)) do segment
        p1, p2 = extrema(segment)
        abs(get_lon(p1) - get_lon(p2)) > 180°
    end
end

function split_antimeridian(poly::Union{PolyArea, Multi})
	rs = rings(poly)
    split_antimeridian!(rs)
	return Multi(map(PolyArea, rs))
end

function split_antimeridian!(rs::Vector{<:Ring})
	for i in reverse(eachindex(rs))
		r = rs[i]
		has_antimeridian(r) || continue
		new_rings = split_antimeridian(r)
		splice!(rs, i, new_rings)
	end
end

# Split a ring crossing the antimeridian into multiple rings. Code is based on the implementation in the python library at https://github.com/gadomski/antimeridian, though it's simplified to only work on a subset of cases for the moment
function split_antimeridian(ring::Ring)
    ptype = eltype(vertices(ring))
    segs = Vector{ptype}[]
    seg = ptype[]
    δlon(p1, p2) = get_lon(p2) - get_lon(p1)
    for s in segments(ring)
        p1, p2 = extrema(s)
        push!(seg, p1)
        Δlon = δlon(p1, p2)
        if 180° < Δlon < 360° # Left
            lat = crossing_latitude_flat(p1, p2)
            push!(seg, ptype(LatLon(lat, -180°)))
            push!(segs, seg)
            seg = [ptype(LatLon(lat, 180°))]
        elseif -360° < Δlon < -180° # Right
            lat = crossing_latitude_flat(p2, p1)
            push!(seg, ptype(LatLon(lat, 180°)))
            push!(segs, seg)
            seg = [ptype(LatLon(lat, -180°))]
        end
    end
    # This is slightly different from the python library implementation, but I am not sure what is the python library doing there
	prepend!(first(segs), seg)
    # We always force rings to be CCW
    return map(segs) do seg
		r = Ring(seg)
		Meshes.orientation(r) == CCW ? r : reverse(r)
	end
end

# We only assume flat computation for the lat crossing instead of also geodetic as in the python library
function crossing_latitude_flat(p1, p2)
    Δlat = get_lat(p2) - get_lat(p1)
    coeff = 180° - get_lon(p1)
    den = get_lon(p2) + 360° - get_lon(p1)
    if get_lon(p1) <= 0°
        coeff = 180° + get_lon(p1)
        den = get_lon(p1) + 360° - get_lon(p2)
    end
    return get_lat(p1) + coeff * Δlat / den
end