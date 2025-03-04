## Define Getters
get_lat(p::Point{𝔼{2},<:Cartesian2D{WGS84Latest}}) = coords(p).y |> ustrip |> Deg # LAT is Y
get_lat(p::Point{🌐,<:LatLon{WGS84Latest}}) = coords(p).lat
get_lat(p::LatLon) = p.lat

get_lon(p::Point{𝔼{2},<:Cartesian2D{WGS84Latest}}) = coords(p).x |> ustrip |> Deg # LON is X
get_lon(p::Point{🌐,<:LatLon{WGS84Latest}}) = coords(p).lon
get_lon(p::LatLon) = p.lon

## boders()
# Define borders for PolyBorder
CountriesBorders.borders(::Type{LatLon}, pb::PolyBorder) = pb.latlon
CountriesBorders.borders(::Type{Cartesian}, pb::PolyBorder) = pb.cart

# Define borders for MultiBorder
CountriesBorders.borders(::Type{LatLon}, mb::MultiBorder) = mb.latlon
CountriesBorders.borders(::Type{Cartesian}, mb::MultiBorder) = mb.cart

# Define borders for PolyRegion
function CountriesBorders.borders(::Type{T}, pr::PolyRegion) where {T}
    borders(T, pr.domain)
end

# Define borders for PolyRegionOffset
function CountriesBorders.borders(::Type{T}, pr::PolyRegionOffset) where {T}
    borders(T, pr.domain)
end

# Define borders for GeoRegion
function CountriesBorders.borders(::Type{T}, gr::GeoRegion) where {T}
    map(Base.Fix1(borders, T), gr.domain)
end

# Define borders for GeoRegionOffset
function CountriesBorders.borders(::Type{T}, gr::GeoRegionOffset) where {T}
    borders(T, gr.domain)
end

## Base.in()
# //NOTE: Interface choice: no possbility to call Base.in on GeoRegion, PolyRegion, or LatBeltRegion with a Cartesian2D point. This is a safe choice of interface of users.
# GeoRegion()
# Fallback on Base.in for CountryBorder defined in CountriesBorders
Base.in(p::Point{🌐,<:LatLon{WGS84Latest}}, gr::GeoRegion) = in(p, gr.domain)
Base.in(p::LatLon, gr::GeoRegion) = in(p, gr.domain)

# GeoRegionOffset()
Base.in(p::Point{🌐,<:LatLon{WGS84Latest}}, gr::GeoRegionOffset) = in(p, gr.domain)
Base.in(p::LatLon, gr::GeoRegionOffset) = in(p, gr.domain)

# PolyRegion()
Base.in(p::Point{𝔼{2},<:Cartesian2D{WGS84Latest}}, pb::PolyBorder) = in(p, borders(Cartesian, pb))
Base.in(p::Point{🌐,<:LatLon{WGS84Latest}}, pb::PolyBorder) = in(Meshes.flat(p), pb) # Flatten the point in Cartesian2D and call the method above
Base.in(p::LatLon, pb::PolyBorder) = in(Point(LatLon{WGS84Latest,Deg{Float32}}(p.lat, p.lon)), pb)

Base.in(p::Point{🌐,<:LatLon{WGS84Latest}}, pr::PolyRegion) = in(p, pr.domain)
Base.in(p::LatLon, pr::PolyRegion) = in(p, pr.domain)

# PolyRegionOffset()
Base.in(p::Point{𝔼{2},<:Cartesian2D{WGS84Latest}}, mb::MultiBorder) = in(p, borders(Cartesian, mb))
Base.in(p::Point{🌐,<:LatLon{WGS84Latest}}, mb::MultiBorder) = in(Meshes.flat(p), mb) # Flatten the point in Cartesian2D and call the method above
Base.in(p::LatLon, mb::MultiBorder) = in(Point(LatLon{WGS84Latest,Deg{Float32}}(p.lat, p.lon)), mb)

Base.in(p::Point{🌐,<:LatLon{WGS84Latest}}, pr::PolyRegionOffset) = in(p, pr.domain)
Base.in(p::LatLon, pr::PolyRegionOffset) = in(p, pr.domain)

# LatBeltRegion()
Base.in(p::Point{🌐,<:LatLon{WGS84Latest}}, lbr::LatBeltRegion) = lbr.lim[1] < get_lat(p) < lbr.lim[2]
Base.in(p::LatLon, lbr::LatBeltRegion) = lbr.lim[1] < p.lat < lbr.lim[2]

## centroid()
# Define ad-hoc methods for GeoRegion - using centroid definition of CountriesBorders.jl
Meshes.centroid(crs::Type{<:Union{LatLon,Cartesian}}, d::GeoRegion) = centroid(crs, d.domain) # Fallback on all the definitions in CountriesBorders.jl for CountryBorder
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

## extract_plot_coords
CountriesBorders.extract_plot_coords(b::Union{PolyBorder, MultiBorder}) = extract_plot_coords(borders(LatLon, b))

CountriesBorders.extract_plot_coords(r::AbstractRegion) = extract_plot_coords(r.domain)
