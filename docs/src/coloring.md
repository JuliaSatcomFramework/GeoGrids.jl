# Coloring

```@setup plot
using PlotlyBase
using PlotlyDocumenter
using GeoGrids
```

```@docs
color_greedy
```

## Example
```@example plot
region = GeoRegion(admin="Italy")
centers, tiles = generate_tesselation(region, 50e3, HEX(), EO())
colors = color_greedy(centers, 50e3, 4)
plot = plot_geo_cells(centers, tiles; colors, title="Hexagonal Tessellation with Colors of Italy", kwargs_layout=(;geo_fitbounds="locations"))
to_documenter(plot) # hide
```