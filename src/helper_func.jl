"""
    _wrap_latlon(lat::Number, lon::Number)

The `_wrap_latlon` function normalizes and wraps geographic coordinates,
latitude (`lat`) and longitude (`lon`). It ensures that the latitude is within
the range [-90, 90] degrees and the longitude is within the range [-180, 180)
degrees. This function is useful for handling geographic data where coordinates
might exceed their typical bounds.

## Arguments
- `lat::Number`: The latitude value to be normalized and wrapped, expressed in \
degrees.
- `lon::Number`: The longitude value to be normalized and wrapped, expressed \
in degrees.

## Returns
- `Tuple{Number, Number}`: A tuple `(lat, lon)` in degrees where `lat` is in the \
range [-90, 90] and `lon` is in the range [-180, 180).
"""
function _wrap_latlon(lat::Number, lon::Number)
    # Normalize lat to the range [-180, 180)
    lat = rem(lat, 360, RoundNearest)
    lon = rem(lon, 360, RoundNearest)

    # Wrap to the range [-90, 90] and make the longitude "jump"
    if lat > 90
        lat = 180 - lat
        lon = lon + 180
        lon = rem(lon, 360, RoundNearest)
    elseif lat < -90
        lat = -180 - lat
        lon = lon + 180
        lon = rem(lon, 360, RoundNearest)
    end

    return lat, lon
end

"""
    _add_angular_offset(inputθϕ, offsetθϕ) -> NamedTuple{(:θ, :ϕ), Tuple{Float64, Float64}}

Add an angular offset to given spherical coordinates.

## Arguments
- `inputθϕ::NamedTuple{(:θ, :ϕ), Tuple{Float64, Float64}}`: The input spherical \
coordinates with components `θ` (polar angle) and `ϕ` (azimuthal angle) in \
radians.
- `offsetθϕ::NamedTuple{(:θ, :ϕ), Tuple{Float64, Float64}}`: The offset \
spherical coordinates with components `θ` (polar angle) and `ϕ` (azimuthal \
angle) in radians.

## Returns
- `NamedTuple{(:θ, :ϕ), Tuple{Float64, Float64}}`: The new spherical coordinates \
after applying the angular offset, with components `θ` and `ϕ` in radians.
"""
function _add_angular_offset(inputθϕ, offsetθϕ)
    sϕ, cϕ = sincos(inputθϕ.ϕ)
    sθ, cθ = sincos(inputθϕ.θ)
    sΔϕ, cΔϕ = sincos(offsetθϕ.ϕ)
    sΔθ, cΔθ = sincos(offsetθϕ.θ)

    # Define the orthonormal basis describing a new reference system [r̂, θ̂, ϕ̂] 
    # and it's relation with ECEF cartesian [x̂, ŷ, ẑ] as per
    # [Wikipedia](https://en.wikipedia.org/wiki/Spherical_coordinate_system#Integration_and_differentiation_in_spherical_coordinates)
    r̂ = SA_F64[sθ*cϕ, sθ*sϕ, cθ]
    θ̂ = SA_F64[cθ*cϕ, cθ*sϕ, -sθ]
    ϕ̂ = SA_F64[-sϕ, cϕ, 0]

    # Define the rotation matrix to go from the new local reference system to
    # the ECEF cartesian system. We want the z axis of the new system to be
    # aligned with the direction of the input vector [centre], so the order of
    # axis such to have a right-handed reference system will be [θ̂, ϕ̂, r̂]. At
    # this point we can describe the offset rotation as a vector represented ni
    # this new reference system, as described in:
    # https://math.stackexchange.com/questions/4343044/rotate-vector-by-a-random-little-amount
    # Even if we inverted the axes such to get the z axis aligned with the
    # input, the system of equations still hold to get the ordered [x̂, ŷ, ẑ]
    # axes.
    R = hcat(θ̂, ϕ̂, r̂) # [θ̂, ϕ̂, r̂] -> [x̂, ŷ, ẑ]  

    # Write the offset vector wrt the new referene system.
    vᴵ = [sΔθ * cΔϕ, sΔθ * sΔϕ, cΔθ]

    # Transform the offset vector to the ECEF cartesian system.
    v = R * vᴵ

    # Convert back to spherical coordinates
    r = norm(v)
    θ = acos(v[3] / r)
    ϕ = atan(v[2], v[1])

    return (θ=θ, ϕ=ϕ) # [deg] ALBERTO: ?? Is it deg though? as the acos and atan return values in radians
end

# This function takes a box and constructs the corresponding PolyArea by oversampling the lon to avoid artifacts while displaying the PolyArea on scattergeo
function box_to_poly_oversample(b::Union{BOX_CART, BOX_LATLON}, lon_dist = 5)
	lo, hi = extrema(b) .|> coords .|> CoordRefSystems.raw
    lo_lon, lo_lat = lo
    hi_lon, hi_lat = hi
    f = to_cart_point
	Δlon = hi_lon - lo_lon
    # We have points in lon distanced by around 5° mostly to avoid problems in plotting on scattergeo
    nlon = max(2, ceil(Int, Δlon / lon_dist))
	hi_gen = [f(LatLon(hi_lat, lon)) for lon in range(lo_lon, hi_lon, nlon)]
	lo_gen = [f(LatLon(lo_lat, lon)) for lon in range(hi_lon, lo_lon, nlon)]
    p = Ring(vcat(hi_gen, lo_gen)) |> PolyArea
end

# This function is used to take a MULTI_CART or Region and a clipping mask and create a MultiBorder obtained by clipping with the mask using the Sutherland-Hodgman algorithm
function clipped_multiborder(original::MULTI_CART{P}, mask::Union{BoxBorder{P}, PolyBorder{P}}) where P
    # We iterate over all polygons, converting everything to Float32 precision
    clipped_polys = POLY_CART{P}[]
    for poly in polyareas(original)
        clipped = Meshes.clip(poly, mask.cart, Meshes.SutherlandHodgmanClipping())
        if !isnothing(clipped)
            push!(clipped_polys, clipped)
        end
    end
    isempty(clipped_polys) && error("No polygons left after clipping...")
    # We create a Multi from the clipped polygons
    multi_cart = Multi(clipped_polys)
    multi_latlon = latlon_geometry(multi_cart)
    domain = MultiBorder(multi_latlon, multi_cart)
end
function clipped_multiborder(original::MULTI_CART, mask::Union{BoxBorder{P}, PolyBorder{P}}) where P
    clipped_multiborder(cartesian_geometry(P, original), mask)
end
function clipped_multiborder(original, mask::Union{BoxBorder{P}, PolyBorder{P}}) where P
    multi = polyareas(original) |> Multi
    clipped_multiborder(multi, mask)
end
