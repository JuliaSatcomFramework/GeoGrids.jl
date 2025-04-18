## Define Getters
get_lat(p) = to_raw_lonlat(p)[2] |> Deg 
get_lon(p) = to_raw_lonlat(p)[1] |> Deg 

CountriesBorders.floattype(r::AbstractRegion) = return floattype(r.domain)
CountriesBorders.floattype(::LatBeltRegion) = return Float64

CountriesBorders.latlon_geometry(T::Type{<:Real}, b::BorderGeometry) = latlon_geometry(T, b.latlon)
CountriesBorders.cartesian_geometry(T::Type{<:Real}, b::BorderGeometry) = cartesian_geometry(T, b.cart)

## borders()
# Define borders for PolyBorder
CountriesBorders.borders(::Type{LatLon}, b::BorderGeometry) = b.latlon
CountriesBorders.borders(::Type{Cartesian}, b::BorderGeometry) = b.cart

# Define borders for AbstractRegion
CountriesBorders.borders(crs::VALID_CRS, r::AbstractRegion) = borders(crs, r.domain)
CountriesBorders.borders(crs::VALID_CRS, gr::GeoRegion) = 
    map(Base.Fix1(borders, crs), gr.domain)

# This extracts a Box as the borders of the LatBeltRegion
function CountriesBorders.borders(crs::VALID_CRS, r::LatBeltRegion) 
    lo, hi = r.lim
    # We have to reduce the box to exclude the bounds of the region to keep the behavior of the current LatBeltRegion
    f = crs == LatLon ? to_latlon_point : to_cart_point
    lo = ustrip(lo) |> Float32 |> nextfloat
    hi = ustrip(hi) |> Float32 |> prevfloat
    b = Box(
        f(LatLon(lo, -180)),
        f(LatLon(hi, 180)),
    )
end

# polyareas
CountriesBorders.polyareas(b::BoxBorder) = polyareas(b.cart)
CountriesBorders.polyareas(r::AbstractRegion) = polyareas(r.domain)
function CountriesBorders.polyareas(r::LatBeltRegion) 
    b = borders(Cartesian, r)
    return polyareas(b)
end

# bboxes
CountriesBorders.bboxes(b::BoxBorder) = (b.cart, )
CountriesBorders.bboxes(r::AbstractRegion) = bboxes(r.domain)
CountriesBorders.bboxes(r::LatBeltRegion) = (borders(Cartesian, r), )
CountriesBorders.bboxes(b::PolyBorder) = (b.bbox,)
CountriesBorders.bboxes(b::MultiBorder) = b.bboxes


## Base.in()
# //NOTE: Interface choice: no possbility to call Base.in on GeoRegion, PolyRegion, or LatBeltRegion with a Cartesian2D point. This is a safe choice of interface of users.
for P in (POINT_LATLON, LATLON)
    @eval Base.in(p::$P, r::AbstractRegion) = in_exit_early(p, r)
    @eval Base.in(p::$P, b::BorderGeometry) = in_exit_early(p, b)
end

# LatBeltRegion()
function CountriesBorders.in_exit_early(p, llr::LatBeltRegion)
    lo, hi = llr.lim
    (; lat) = to_latlon_point(p) |> coords
    return lo < lat < hi
end
CountriesBorders.in_exit_early(p, b::BoxBorder{P}) where P = to_cart_point(P, p) in borders(Cartesian, b)

## centroid()
# Methods for BorderGeometry
function Meshes.centroid(crs::VALID_CRS, b::BorderGeometry)
    cart = centroid(borders(Cartesian, b))
    return crs === Cartesian ? cart : to_latlon_point(cart)
end
Meshes.centroid(b::BorderGeometry) = centroid(Cartesian, b)

# Define ad-hoc methods for AbstractRegions - using centroid definition of CountriesBorders.jl
Meshes.centroid(crs::VALID_CRS, d::AbstractRegion) = centroid(crs, d.domain) # Fallback on all the definitions in CountriesBorders.jl for CountryBorder
Meshes.centroid(d::AbstractRegion) = centroid(Cartesian, d)

## CountriesBorders.extract_countries()
CountriesBorders.extract_countries(r::GeoRegion) = r.domain

## geom_iterable, used by extract_latlon_coords
GeoPlottingHelpers.geom_iterable(b::BorderGeometry) = geom_iterable(borders(LatLon, b))
GeoPlottingHelpers.geom_iterable(b::BoxBorder) = polyareas(b) # The box is not directly supported

GeoPlottingHelpers.geom_iterable(r::AbstractRegion) = geom_iterable(r.domain)
GeoPlottingHelpers.geom_iterable(r::LatBeltRegion) = polyareas(r)

# Default plot kwargs to lines for regions
GeoPlottingHelpers.geo_plotly_trace_default_kwargs(::AbstractRegion, tracefunc) = (; mode = "lines")