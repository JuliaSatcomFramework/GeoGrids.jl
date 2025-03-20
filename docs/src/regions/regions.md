# Regions

**GeoGrids.jl** provides several region types for defining geographical areas of interest. These regions can be used for filtering points, generating tessellations, and other geospatial operations.

```@setup plot
using PlotlyBase
using PlotlyDocumenter
using GeoGrids
```

## Region Types

GeoGrids.jl defines several types of regions through a type hierarchy:

- `AbstractRegion`: Base abstract type for all regions
  - `GlobalRegion`: Represents the entire globe
  - `LatBeltRegion`: A latitude belt defined by min/max latitudes
  - `GeoRegion`: A region defined by country/countries
  - `PolyRegion`: A region defined by a polygon of coordinates
  - `MultiRegion`: A region represented by a single Multi encompassing polyareas from multiple regions
  - `HotspotRegion`: A region defined by a `center` and a radius in `m`, identifying all points within distance `radius` from the `center`
  - `ClippedRegion`: A region obtained by clipping another region with a mask

### Global Region

The simplest region type is `GlobalRegion`, which represents the entire globe. This is useful when you need to perform operations on a global scale.

```@example plot
region = GlobalRegion()
points = icogrid(sepAng=10)
plot = plot_geo_points(points; title="Global Region with Points")
to_documenter(plot) # hide
```

### Latitude Belt Region

`LatBeltRegion` defines a region between two latitudes. This is useful for studying phenomena that vary with latitude or focusing on specific latitude bands.

```@example plot
region = LatBeltRegion(lim=(-10,10))
points = icogrid(sepAng=4)
filtered = filter_points(points, region)
plot = plot_geo_points(filtered; title="Latitude Belt Region (-10° to 10°)")
to_documenter(plot) # hide
```

### Geographical Region

`GeoRegion` allows you to define regions based on country boundaries. You can specify one or multiple countries.

```@example plot
region = GeoRegion(admin="Italy;Spain")
points = icogrid(sepAng=2)
filtered = filter_points(points, region)
plot = plot_geo_points(filtered; title="Geographical Region (Italy and Spain)", kwargs_layout=(;geo_scope="europe"))
to_documenter(plot) # hide
```

### Polygon Region

`PolyRegion` lets you define custom regions using a list of latitude-longitude coordinates that form a polygon.

```@example plot
region = PolyRegion(domain=[
    LatLon(10°, -5°), 
    LatLon(10°, 15°), 
    LatLon(27°, 15°), 
    LatLon(27°, -5°)
])
points = rectgrid(2)[:]
filtered = filter_points(points, region)
plot = plot_geo_points(filtered; title="Custom Polygon Region")
to_documenter(plot) # hide
```

### HotSpot Region

`HotSpotRegion` defines a region as a circle centered at a specific point with a given radius. It internally still stores a PolyArea for representing the `circle` as the set of points within a given radius from the center is not a circle anymore in the PlateCarree projection used internally by GeoGrids.jl.

!!! note
    The points of the `circle` are computed by assuming a spherical Earth for calculating the distance and use the great circle formulation.

```@example plot
region = HotSpotRegion(center=LatLon(30, 40), radius=1000e3)
plot = scattergeo(region) |> Plot
to_documenter(plot) # hide
```

### MultiRegion

`MultiRegion` is a region that contains multiple regions. It is useful when you need to represent multiple independent polygons as a single region.

```@example plot
hr = HotSpotRegion(center=LatLon(30, 40), radius=1000e3)
italy = GeoRegion(;admin="Italy")
mr = MultiRegion([hr, italy])
plot = scattergeo(mr) |> Plot
to_documenter(plot) # hide
```

### ClippedRegion

`ClippedRegion` is a region that is obtained by clipping another region with a mask, which can be provided either as a `PolyBorder` or as a `BoxBorder`.

```@example plot
box = BoxBorder(LatLon(30, 10), LatLon(45, 35)) # We only keep the portions of the regions that are between 30-45 latitude and 10-35 longitude
reg = GeoRegion(;continent = "Europe")
cr = ClippedRegion(reg, box)
plot = scattergeo(cr) |> Plot
to_documenter(plot) # hide
```

## Enlarged Regions

GeoGrids.jl supports creating enlarged versions of regions through the `GeoRegionOffset` and `PolyRegionOffset` types. These are useful when you need to create buffer zones around existing regions.

### Enlarged Geographical Region

```@example plot
ereg = GeoRegionOffset(delta=50e3, admin="Spain; Italy", resolution=110)
plot = plot_geo_poly(ereg.domain; title="Geo Region Enlarged", kwargs_layout=(;geo_fitbounds="locations"))
to_documenter(plot) # hide
```

### Enlarged Polygon Region

```@example plot
region = PolyRegionOffset(
    delta=100e3, 
    domain=[
        LatLon(10°, -5°), 
        LatLon(10°, 15°), 
        LatLon(27°, 15°), 
        LatLon(27°, -5°)
    ]
)
plot = plot_geo_poly([
    region.domain.latlon.geoms..., 
    region.original.domain.latlon
]; title="Original and Enlarged Polygon Region")
to_documenter(plot) # hide
```


Under the hood, both of these types use the `offset_region` function.
