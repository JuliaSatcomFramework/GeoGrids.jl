module PlotlyBaseExt

using PlotlyExtensionsHelper
using PlotlyBase
using Unitful: ustrip
using Meshes: vertices, rings, Multi, Ngon, 🌐, WGS84Latest

using GeoGrids
using GeoGrids: MultiBorder, PolyBorder, AbstractRegion, BorderGeometry
using GeoGrids.CountriesBorders: LATLON, POINT_LATLON, extract_plot_coords, extract_plot_coords!

const VALID_COORD = Union{LATLON, POINT_LATLON}

plot_coords(x) = extract_plot_coords(x)
function plot_coords(ps::AbstractVector{<:VALID_COORD})
    lat = Float32[]
    lon = Float32[]
    extract_plot_coords!(lat, lon, ps)
    return (; lat, lon)
end
function plot_coords(ps::AbstractVector{<:AbstractVector{<:VALID_COORD}})
    lat = Float32[]
    lon = Float32[]
    for pv in ps
        !isempty(lat) && !isempty(lon) && extract_plot_coords!(lat, lon, LatLon(NaN, NaN))
        extract_plot_coords!(lat, lon, pv)
    end
    return (; lat, lon)
end


const DEFAULT_CELL_CONTOUR = (;
    mode="lines",
    marker=attr(;
        size=1,
        color="rgb(92,97,102)",
    ),
    name="Cell Contour",
    showlegend=false
)

const DEFAULT_CELL_CENTER = (;
    mode="text",
    textfont=attr(; size=10),
    name="Cell Number",
    showlegend=false
)

const DEFAULT_POINT = (;
    mode="markers",
    marker_size=5,
)

const DEFAULT_GEOLAYOUT = (;
    geo=attr(
        projection=attr(
            type="robinson",
        ),
        showocean=true,
        # oceancolor =  "rgb(0, 255, 255)",
        oceancolor="rgb(255, 255, 255)",
        showland=true,
        # landcolor =  "rgb(230, 145, 56)",
        landcolor="rgb(217, 217, 217)",
        showlakes=true,
        # lakecolor =  "rgb(0, 255, 255)",
        lakecolor="rgb(255, 255, 255)",
        showcountries=true,
        lonaxis=attr(
            showgrid=true,
            gridcolor="rgb(102, 102, 102)"
        ),
        lataxis=attr(
            showgrid=true,
            gridcolor="rgb(102, 102, 102)"
        )
    )
)

const DEFAULT_COLORS = [
    "rgb(230, 25, 75)",   # Red
    "rgb(60, 180, 75)",   # Green
    "rgb(255, 225, 25)",  # Yellow
    "rgb(0, 130, 200)",   # Blue
    "rgb(245, 130, 48)",  # Orange
    "rgb(145, 30, 180)",  # Purple
    "rgb(70, 240, 240)",  # Cyan
    "rgb(240, 50, 230)",  # Magenta
    "rgb(210, 245, 60)",  # Lime
    "rgb(250, 190, 190)", # Pink
    "rgb(0, 128, 128)",   # Teal
    "rgb(230, 190, 255)", # Lavender
    "rgb(170, 110, 40)",  # Brown
    "rgb(255, 250, 200)", # Beige
    "rgb(128, 128, 128)"  # Grey
]

## Auxiliary Functions
# Internal functions for creating the scatter plots.
"""
    GeoGrids._get_scatter_points(points::AbstractVector{<:Union{LatLon, Point{🌐,<:LatLon{WGS84Latest}}}}; kwargs...)

This function takes an array of geographic points and generates a scatter plot
using the `scattergeo` function from the PlotlyJS package. The points are
converted to a vector of `LatLon` objects if they are not already. The
latitude and longitude of each point are extracted and used to create the
scatter plot.

## Arguments
- `points::AbstractVector{<:Union{LatLon, Point{🌐,<:LatLon{WGS84Latest}}}}`: An array of \
points that can be of type `LatLon`, `AbstractVector`, or `Tuple`. Each \
point represents a geographic location.

## Keyword Arguments
- `kwargs`: Additional keyword arguments to customize the `scattergeo` trace. \
These arguments are passed directly to the `scattergeo` function from \
PlotlyBase.

## Returns
- A `scattergeo` object representing the scatter plot of the provided \
geographic points.
"""
function GeoGrids._get_scatter_points(points::AbstractVector{<:VALID_COORD}; kwargs...)
    # Markers for the points
    return scattergeo(;
        plot_coords(points)...,
        DEFAULT_POINT...,
        kwargs...
    )
end

"""
    GeoGrids._get_scatter_cellcontour(polygons::AbstractVector{<:Union{LatLon,Point{🌐,<:LatLon{WGS84Latest}}}}; kwargs...)

This function creates a geographic scatter plot of cell contours based on the
input polygons. Each polygon is processed to extract its vertices' latitude and
longitude, which are then used to plot the contours on a geographic map.

## Arguments
- `polygons::AbstractVector{<:Union{LatLon,Point{🌐,<:LatLon{WGS84Latest}}}}`: A vector of \
polygons, where each polygon is represented by a vector of `LatLon` \
objects. Each `LatLon` object holds latitude and longitude information.

## Keyword Arguments
- `colors::Union{AbstractVector{<:Integer},Nothing}`: A vector of integers representing the color index of each cell. If `nothing`, all will be colored in the same color.
- `kwargs...`: Additional keyword arguments to customize the scatter plot. These \
are passed directly to the `scattergeo` function, allowing customization of \
the plot's appearance (e.g., color, line style, marker options).

## Returns
- A vector of `scattergeo` plot objects: The scatter plot visualization of the cell \
contours, ready for rendering in a geographic plot.
"""
function GeoGrids._get_scatter_cellcontour(polygons::AbstractVector{<:AbstractVector{<:VALID_COORD}}; colors, kwargs...)
    function f(polys; color = nothing) 
        cc = isnothing(color) ? DEFAULT_CELL_CONTOUR : attr(; DEFAULT_CELL_CONTOUR..., marker_color=color)
        scattergeo(;
            plot_coords(polys)...,
            cc...,
            kwargs...
        )
    end
    isnothing(colors) && return [f(polygons)]
    # We have more colors
    @assert length(colors) == length(polygons) "You have a different number of colors and polygons."
    map(unique(colors)) do c
        polys = polygons[findall(==(c), colors)]
        f(polys; color=DEFAULT_COLORS[c])
    end
end

"""
    GeoGrids._get_scatter_poly(poly::PolyArea{🌐,<:LatLon{WGS84Latest}}; kwargs...)

This function creates a geographic scatter plot of a polygon's boundary.

## Arguments
- `poly::PolyArea{🌐,<:LatLon{WGS84Latest}}`: A polygon object representing the \
area to be plotted.

## Keyword Arguments
- `kwargs...`: Additional keyword arguments to customize the scatter plot. These \
are passed directly to the `scattergeo` function.

## Returns
- A vector of `scattergeo` plot objects: Each element represents a ring of the \
polygon, ready for rendering in a geographic plot.

## Notes
- The function processes each ring of the polygon separately, creating a trace \
for each.
- The first and last points of each ring are connected to close the polygon.
- By default, the lines are colored red and have no legend entry.
"""
function GeoGrids._get_scatter_poly(poly::PolyArea{🌐,<:LatLon{WGS84Latest}}; kwargs...)
    # scatter line
    map(rings(poly)) do ring
        # Each ring will be a separate trace.
        scattergeo(;
            plot_coords(ring)...,
            mode="lines",
            line_color="red",
            showlegend=false,
            kwargs...
        )
    end
end

"""
    _default_geolayout(; title::String="Point Position GEO Map", camera::Symbol=:twodim, kwargs...) -> PlotlyJS.Layout

This function generates a default geographic layout for a PlotlyJS plot. It sets
up the layout with a geographic projection, ocean and land colors, and grid
lines for latitude and longitude. The projection type is chosen based on the
camera perspective: "orthographic" for 3D and "natural earth" for 2D.

## Keyword Arguments
- `title::String="Point Position GEO Map"`: The title of the plot. Defaults to \
"Point Position GEO Map".
- `camera::Symbol=:twodim`: The camera perspective of the plot. Can be `:twodim` \
for a 2D view or `:threedim` for a 3D view.
- `kwargs`: Additional keyword arguments to customize the layout further. These \
arguments are passed directly to the `Layout` constructor.

## Returns
- A `PlotlyJS.Layout` object representing the default geographic layout with the \
specified options.
"""
function _default_geolayout(; title="Point Position GEO Map", camera::Symbol=:twodim, kwargs...)
    projection = camera == :threedim ? "orthographic" : "natural earth"

    # Create the geo layout
    return Layout(;
        DEFAULT_GEOLAYOUT...,
        title=title,
        geo_projection_type=projection,
        kwargs...
    )
end

## Core plotting functions.
"""
    plot_geo_points(points::AbstractVector{<:Union{LatLon, Point{🌐,<:LatLon{WGS84Latest}}}}; title::String="Point Position GEO Map", camera::Symbol=:twodim, kwargs_scatter::NamedTuple=(); kwargs_layout::NamedTuple=())
    plot_geo_points(p::Union{LatLon, Point{🌐,<:LatLon{WGS84Latest}}}; title::String="Point Position GEO Map", camera::Symbol=:twodim, kwargs_scatter::NamedTuple=(); kwargs_layout::NamedTuple=())

This function generates a geographic plot for a given array of points. It
creates a scatter plot of the points using `GeoGrids._get_scatter_points` and sets up the
layout with `_default_geolayout`. The plot is created using the `plotly_plot`
function from PlotlyExtensionsHelper.

## Arguments
- `points::AbstractVector{<:Union{LatLon, Point{🌐,<:LatLon{WGS84Latest}}}}`: An \
array of points to be plotted. Each point can be of type `LatLon`, \
`AbstractVector`, or `Tuple` containing latitude and longitude.

## Keyword Arguments
- `title::String="Point Position GEO Map"`: The title of the plot. Defaults to \
"Point Position GEO Map".
- `camera::Symbol=:twodim`: The camera perspective of the plot. Can be `:twodim` \
for a 2D view or other supported camera views.
- `kwargs_scatter::NamedTuple=()` : Additional keyword arguments for customizing \
the scatter plot. These arguments are passed directly to the \
`GeoGrids._get_scatter_points` function.
- `kwargs_layout::NamedTuple=()` : Additional keyword arguments for customizing \
the plot layout. These arguments are passed directly to the \
`_default_geolayout` function.

## Returns
- A object representing the geographic plot of the provided \
points (the specific type of the object depends on the loaded Plotly package as per `plotly_plot` output).

See also: [`GeoGrids._get_scatter_points`](@ref), [`_default_geolayout`](@ref),
[`plot_geo_cells`](@ref)
"""
function GeoGrids.plot_geo_points(points::AbstractVector{<:VALID_COORD}; title="Point Position GEO Map", camera::Symbol=:twodim, kwargs_scatter=(;), kwargs_layout=(;))
    # Markers for the points
    scatterpoints = GeoGrids._get_scatter_points(points; kwargs_scatter...)
    layout = _default_geolayout(; title, camera, kwargs_layout...)

    plotly_plot([scatterpoints], layout)
end
GeoGrids.plot_geo_points(p::VALID_COORD; kwargs...) = GeoGrids.plot_geo_points([p]; kwargs...)

"""
    plot_geo_cells(cellCenters::AbstractVector{<:Union{LatLon, Point{🌐,<:LatLon{WGS84Latest}}}}; title::String="Cell Layout GEO Map", camera::Symbol=:twodim, kwargs_centers::NamedTuple=(;), kwargs_layout::NamedTuple=(;))
    plot_geo_cells(cc::Union{LatLon, Point{🌐,<:LatLon{WGS84Latest}}}; kwargs...)
    plot_geo_cells(cellCenters::AbstractVector{<:Union{LatLon, Point{🌐,<:LatLon{WGS84Latest}}}}, cellContours::AbstractVector{<:AbstractVector{<:Union{LatLon,Point{🌐,<:LatLon{WGS84Latest}}}}}; title::String="Cell Layout GEO Map", camera::Symbol=:twodim, kwargs_centers::NamedTuple=(;), kwargs_contours::NamedTuple=(;),kwargs_layout::NamedTuple=(;))
    plot_geo_cells(cellCenter::Union{LatLon, Point{🌐,<:LatLon{WGS84Latest}}}, cellContour::AbstractVector{<:Union{LatLon,Point{🌐,<:LatLon{WGS84Latest}}}}; kwargs...)

Plot geographical cell centers and/or contours on a map using various input
configurations.

## Arguments
- `cellCenters::AbstractVector{<:Union{LatLon, \
Point{🌐,<:LatLon{WGS84Latest}}}}`: A vector of geographical coordinates \
representing the centers of the cells.
- `cc::Union{LatLon, Point{🌐,<:LatLon{WGS84Latest}}}`: A single geographical \
coordinate representing the center of a cell.
- `cellContours::AbstractVector{<:Union{LatLon,Point{🌐,<:LatLon{WGS84Latest}}}}`: \
A vector of vectors containing the contours for each cell, represented as \
`LatLon` coordinates.
- `cellContour::AbstractVector{<:Union{LatLon,Point{🌐,<:LatLon{WGS84Latest}}}}`: \
A vector of `LatLon` objects representing the contour of a single cell.
- `title::String`: Title of the plot. Default is `"Cell Layout GEO Map"`.
- `camera::Symbol`: Camera view for the plot. Default is `:twodim`.
- `kwargs_centers::NamedTuple`: Additional options for customizing the \
appearance of cell centers.
- `kwargs_contours::NamedTuple`: Additional options for customizing the \
appearance of cell contours.
- `kwargs_layout::NamedTuple`: Additional options for customizing the overall \
layout of the plot.
- `kwargs...`: Additional keyword arguments forwarded to the underlying plotting \
functions.

## Return Value
- Returns a plot object displaying the specified geographical cell centers \
and/or contours on a map.

See also: [`GeoGrids._get_scatter_points`](@ref), [`GeoGrids._get_scatter_cellcontour`](@ref),
[`_default_geolayout`](@ref), [`plot_geo_points`](@ref)
"""
function GeoGrids.plot_geo_cells(cellCenters::AbstractVector{<:VALID_COORD}; title="Cell Layout GEO Map", camera::Symbol=:twodim, kwargs_centers=(;), kwargs_layout=(;))
    # Fallback method to plot only cell centers
    k = (; DEFAULT_CELL_CENTER..., text=map(x -> string(x), 1:length(cellCenters)), kwargs_centers...) # Default for text mode for cellCenters
    GeoGrids.plot_geo_points(cellCenters; title, camera, kwargs_scatter=k, kwargs_layout)
end
GeoGrids.plot_geo_cells(cc::VALID_COORD; kwargs...) = GeoGrids.plot_geo_points([cc]; kwargs...)

function GeoGrids.plot_geo_cells(cellCenters::AbstractVector{<:VALID_COORD}, cellContours::AbstractVector{<:AbstractVector{<:VALID_COORD}}; title="Cell Layout GEO Map", colors::Union{AbstractVector{<:Integer},Nothing}=nothing, camera::Symbol=:twodim, kwargs_centers=(;), kwargs_contours=(;), kwargs_layout=(;))
    # Create scatter plot for the cells contours.
    scatterContours = GeoGrids._get_scatter_cellcontour(cellContours; colors, kwargs_contours...)

    # Create scatter plot for the cell centers.
    
    fontColors = colors === nothing ? (;) : (;textfont_color=map(x -> DEFAULT_COLORS[x], colors))
    k = (; DEFAULT_CELL_CENTER..., text=map(x -> string(x), 1:length(cellCenters)), fontColors..., kwargs_centers...) # Default for text mode for cellCenters
    scatterCenters = GeoGrids._get_scatter_points(cellCenters; k...)

    # Create layout
    layout = _default_geolayout(; title, camera, kwargs_layout...)

    plotly_plot([scatterContours..., scatterCenters], layout)
end
GeoGrids.plot_geo_cells(cellCenter::VALID_COORD, cellContour::AbstractVector{<:VALID_COORD}; kwargs...) = GeoGrids.plot_geo_cells([cellCenter], [cellContour]; kwargs...)

"""
    plot_geo_poly(polys::AbstractVector{<:PolyArea{🌐,<:LatLon{WGS84Latest}}}; title="Polygon GEO Map", camera::Symbol=:twodim, kwargs_scatter=(;), kwargs_layout=(;))
    plot_geo_poly(poly::PolyArea{🌐,<:LatLon{WGS84Latest}}; kwargs...)
    plot_geo_poly(multi::Multi{🌐,<:LatLon{WGS84Latest}}; kwargs...)
    plot_geo_poly(b::Union{<:GeoGrids.PolyBorder, <:GeoGrids.MultiBorder}; kwargs...)

Plot geographical polygons on a map.

## Arguments
- `polys::AbstractVector{<:PolyArea{🌐,<:LatLon{WGS84Latest}}}`: A vector of \
polygon areas to be plotted.
- `poly::PolyArea{🌐,<:LatLon{WGS84Latest}}`: A single polygon area to be \
plotted.
- `multi::Multi{🌐,<:LatLon{WGS84Latest}}`: A multi-polygon object to be \
plotted.
- `b::Union{<:GeoGrids.PolyBorder, <:GeoGrids.MultiBorder}`: A polygon border or \
multi-border object to be plotted.

- `title::String`: Title of the plot. Default is `"Polygon GEO Map"`.
- `camera::Symbol`: Camera view for the plot. Default is `:twodim`.
- `kwargs_scatter::NamedTuple`: Additional options for customizing the \
appearance of the polygons.
- `kwargs_layout::NamedTuple`: Additional options for customizing the overall \
layout of the plot.
- `kwargs...`: Additional keyword arguments forwarded to the underlying plotting \
functions.

## Return Value
- Returns a plot object displaying the specified geographical polygons on a map.

See also: [`GeoGrids._get_scatter_poly`](@ref), [`_default_geolayout`](@ref),
[`plot_geo_cells`](@ref)
"""
function GeoGrids.plot_geo_poly(polys::AbstractVector{<:PolyArea{🌐,<:LatLon{WGS84Latest}}}; title="Polygon GEO Map", camera::Symbol=:twodim, kwargs_scatter=(;), kwargs_layout=(;))
    # Extract the vertices of the polygon
    scattersPoly = map(x -> GeoGrids._get_scatter_poly(x; kwargs_scatter...) |> splat(vcat), polys) |> splat(vcat)

    layout = _default_geolayout(; title, camera, kwargs_layout...)

    plotly_plot([scattersPoly...], layout)
end
GeoGrids.plot_geo_poly(poly::PolyArea{🌐,<:LatLon{WGS84Latest}}; kwargs...) = GeoGrids.plot_geo_poly([poly]; kwargs...)
GeoGrids.plot_geo_poly(multi::Multi{🌐,<:LatLon{WGS84Latest}}; kwargs...) = GeoGrids.plot_geo_poly(multi.geoms; kwargs...)
GeoGrids.plot_geo_poly(b::Union{<:GeoGrids.PolyBorder,<:GeoGrids.MultiBorder}; kwargs...) = GeoGrids.plot_geo_poly(b.latlon; kwargs...)

## Additional Functions
"""
    plot_unitarysphere(points_cart)

This function takes an SVector{3} of Cartesian coordinates and plots the
corresponding points on a unitary sphere. The sphere is defined by a range of
angles that are discretized into a grid of n_sphere points.

## Arguments:
- `points_cart`: an array of Cartesian coordinates of the points to be plotted \
on the unitary sphere.
## Output:
- Plot of the unitary sphere with the input points represented as markers.
"""
function GeoGrids.plot_unitarysphere(points_cart; kwargs_scatter=(;), kwargs_layout=(;))
    # Reference Sphere
    n_sphere = 100
    u = range(-π, π; length=n_sphere)
    v = range(0, π; length=n_sphere)
    x_sphere = cos.(u) * sin.(v)'
    y_sphere = sin.(u) * sin.(v)'
    z_sphere = ones(n_sphere) * cos.(v)'
    color = ones(size(z_sphere))
    sphere = surface(
        z=z_sphere,
        x=x_sphere,
        y=y_sphere,
        colorscale=[[0, "rgb(2,204,150)"], [1, "rgb(2,204,150)"]],
        showscale=false
    )

    # Take an array of SVector
    markers = scatter3d(;
        x=map(x -> x[1], points_cart),
        y=map(x -> x[2], points_cart),
        z=map(x -> x[3], points_cart),
        mode="markers",
        marker_size=4,
        marker_color="rgb(0,0,0)",
        kwargs_scatter...
    )

    layout = Layout(;
        scene=attr(
            xaxis=attr(
                visible=false,
            ),
            yaxis=attr(
                visible=false,
            ),
            zaxis=attr(
                visible=false,
            ),
        ),
        width=700,
        kwargs_layout...
    )

    # Plot([sphere,markers],layout)
    plotly_plot([sphere, markers], layout)
end

PlotlyBase.scattergeo(b::BorderGeometry; kwargs...) = scattergeo(; extract_plot_coords(b)..., mode = "lines")
PlotlyBase.scattergeo(b::AbstractRegion; kwargs...) = scattergeo(; extract_plot_coords(b)..., mode = "lines")

end # module PlotlyBaseExt