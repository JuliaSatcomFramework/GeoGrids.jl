## Define Getters
get_lat(p::Point{𝔼{2},<:Cartesian2D{WGS84Latest}}) = coords(p).y |> ustrip |> Deg # LAT is Y
get_lat(p::Point{🌐,<:LatLon{WGS84Latest}}) = coords(p).lat
get_lat(p::LatLon) = p.lat

get_lon(p::Point{𝔼{2},<:Cartesian2D{WGS84Latest}}) = coords(p).x |> ustrip |> Deg # LON is X
get_lon(p::Point{🌐,<:LatLon{WGS84Latest}}) = coords(p).lon
get_lon(p::LatLon) = p.lon

## borders()
# Define borders for PolyBorder
CountriesBorders.borders(::Type{LatLon}, b::Union{PolyBorder, MultiBorder}) = b.latlon
CountriesBorders.borders(::Type{Cartesian}, b::Union{PolyBorder, MultiBorder}) = b.cart

# Define borders for AbstractRegion
CountriesBorders.borders(crs::VALID_CRS, r::AbstractRegion) = borders(crs, r.domain)
CountriesBorders.borders(crs::VALID_CRS, gr::GeoRegion) = 
    map(Base.Fix1(borders, crs), gr.domain)

# polyareas
CountriesBorders.polyareas(r::AbstractRegion) = polyareas(r.domain)
function CountriesBorders.polyareas(r::LatBeltRegion) 
    lo, hi = r.lim
    # We have to reduce the box to exclude the bounds of the region to keep the behavior of the current LatBeltRegion
    lo = ustrip(lo) |> Float32 |> nextfloat
    hi = ustrip(hi) |> Float32 |> prevfloat
    f(ll) = to_cart_point(ll, Float32)
    # We make this not just 4 points to have plotting with scattergeo not result in very ugly thing. We anyhow expect this not to be a problem for the speed of `in`
    pts = reduce(vcat, (
        [LatLon(hi, lon) |> f for lon in range(-180, 180; step = 5)],
        [LatLon(lat, 180) |> f for lat in range(hi, lo; step = -1)],
        [LatLon(lo, lon) |> f for lon in range(180, -180; step = -5)],
        [LatLon(lat, -180) |> f for lat in range(lo, hi; step = 1)],
    ))
    p = PolyArea(pts[1:end-1])
    return (p, )
end

# bboxes
CountriesBorders.bboxes(r::AbstractRegion) = bboxes(r.domain)
function CountriesBorders.bboxes(r::LatBeltRegion)
    lo, hi = r.lim
    # We have to reduce the box to exclude the bounds of the region to keep the behavior of the current LatBeltRegion
    lo = ustrip(lo) |> Float32 |> nextfloat
    hi = ustrip(hi) |> Float32 |> prevfloat
    return (
        Box(
            to_cart_point(LatLon(lo, -180), Float32),
            to_cart_point(LatLon(hi, 180), Float32),
        ),
    )
end
CountriesBorders.bboxes(b::PolyBorder) = (b.bbox,)
CountriesBorders.bboxes(b::MultiBorder) = b.bboxes


## Base.in()
# //NOTE: Interface choice: no possbility to call Base.in on GeoRegion, PolyRegion, or LatBeltRegion with a Cartesian2D point. This is a safe choice of interface of users.
for P in (POINT_LATLON, LATLON)
    @eval Base.in(p::$P, r::AbstractRegion) = in_exit_early(p, r)
    @eval Base.in(p::$P, b::Union{PolyBorder, MultiBorder}) = in_exit_early(p, b)
end

# LatBeltRegion()
CountriesBorders.in_exit_early(p, llr::LatBeltRegion) = to_cart_point(p, Float32) in only(bboxes(llr))

## centroid()
# Define ad-hoc methods for GeoRegion - using centroid definition of CountriesBorders.jl
Meshes.centroid(crs::VALID_CRS, d::GeoRegion) = centroid(crs, d.domain) # Fallback on all the definitions in CountriesBorders.jl for CountryBorder
Meshes.centroid(d::GeoRegion) = centroid(Cartesian, d)

# Define ad-hoc methods for GeoRegionOffset - using centroid definition of CountriesBorders.jl
Meshes.centroid(crs::Type{<:Union{LatLon,Cartesian}}, d::GeoRegionOffset) = centroid(crs, d.domain) # Fallback on all the definitions in CountriesBorders.jl for CountryBorder
Meshes.centroid(d::GeoRegionOffset) = centroid(Cartesian, d)

# Define ad-hoc methods for PolyRegion - using centroid definition of Meshes.jl
Meshes.centroid(::Type{Cartesian}, d::PolyBorder) = centroid(d.cart)
function Meshes.centroid(::Type{LatLon}, d::PolyBorder)
    c = centroid(d.cart)
    LatLon{WGS84Latest}(get_lat(c), get_lon(c)) |> Point
end
Meshes.centroid(crs::Type{<:Union{LatLon,Cartesian}}, d::PolyRegion) = centroid(crs, d.domain)
Meshes.centroid(d::PolyRegion) = centroid(Cartesian, d.domain)

# Define ad-hoc methods for PolyRegionOffset - using centroid definition of Meshes.jl
Meshes.centroid(::Type{Cartesian}, d::MultiBorder) = centroid(d.cart)
function Meshes.centroid(::Type{LatLon}, d::MultiBorder)
    c = centroid(d.cart)
    LatLon{WGS84Latest}(get_lat(c), get_lon(c)) |> Point
end
Meshes.centroid(crs::Type{<:Union{LatLon,Cartesian}}, d::PolyRegionOffset) = centroid(crs, d.domain)
Meshes.centroid(d::PolyRegionOffset) = centroid(Cartesian, d.domain)

## CountriesBorders.extract_countries()
CountriesBorders.extract_countries(r::GeoRegion) = r.domain

## geom_iterable, used by extract_plot_coords
CountriesBorders.geom_iterable(b::Union{PolyBorder, MultiBorder}) = geom_iterable(borders(LatLon, b))

CountriesBorders.geom_iterable(r::AbstractRegion) = geom_iterable(r.domain)
CountriesBorders.geom_iterable(r::LatBeltRegion) = polyareas(r)