# Filtering & Grouping

**GeoGrids.jl** provides powerful functionality for filtering points based on geographical regions. This includes both basic point-in-region testing and efficient filtering of large point sets.

```@setup plot
using PlotlyBase
using PlotlyDocumenter
using GeoGrids
```

## Basic Filtering

The `filter_points` function provides a straightforward way to filter points based on regions:

```@example plot
# Create a latitude belt region
region = LatBeltRegion(lim=(-10,10))
points = icogrid(sepAng=4)

# Filter points to only those in the region
filtered = filter_points(points, region)
plot = plot_geo_points(filtered; title="Points Filtered to Equatorial Belt")
to_documenter(plot) # hide
```

## Grouping Points by Region

Sometimes you need to categorize points based on which regions they fall into. The `group_by_domain` function handles this:

```@example plot
# Create multiple regions
regions = [
    GeoRegion(admin="Spain"),
    GeoRegion(admin="France"),
    GeoRegion(admin="Italy")
]
points = icogrid(sepAng=2)

# Group points by region
grouped = group_by_domain(points, regions)

# Plot points in Spain as an example
spain_points = grouped["Spain"]
plot = plot_geo_points(spain_points; title="Points Grouped to Spain")
to_documenter(plot) # hide
```
