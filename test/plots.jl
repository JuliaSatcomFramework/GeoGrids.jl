@testitem "plot_geo_points" tags = [:general] begin
    using PlotlyBase

    @test plot_unitarysphere(GeoGrids._icogrid(100; coord=:cart)) isa Plot

    # Multi Point
    @test plot_geo_points(icogrid(sepAng=5°)) isa Plot
    @test plot_geo_points(icogrid(sepAng=deg2rad(4) * rad); camera=:threedim) isa Plot
    @test plot_geo_points([rectgrid(5)...]) isa Plot
    @test plot_geo_points([rectgrid(5°)...]; camera=:threedim) isa Plot
    
    # Single Point
    @test plot_geo_points(LatLon(10,10); camera=:threedim) isa Plot
end

@testitem "plot_geo_cells" tags = [:general] begin
    using PlotlyBase

    reg = GeoRegion(; name="Tassellation", admin="Spain")
	
    # Multi Point
    centers, ngon = generate_tesselation(reg, 40000, HEX(;pattern=:circ), EO())
    @test plot_geo_cells(centers, ngon) isa Plot
    @test plot_geo_cells(centers[1], ngon[1]) isa Plot
    # Single Point
    @test plot_geo_cells(centers) isa Plot
    @test plot_geo_cells(centers[1]) isa Plot
	centers, ngon = generate_tesselation(reg, 40000, HEX(), EO())
    # Multi Point
    @test plot_geo_cells(centers, ngon) isa Plot
    @test plot_geo_cells(centers[1], ngon[1]) isa Plot
    # Single Point
    @test plot_geo_cells(centers) isa Plot
    @test plot_geo_cells(centers[1]) isa Plot

    # Check colors
    colors = color_greedy(centers, 40000, 4)
    plt = plot_geo_cells(centers, ngon; colors)
    @test length(plt.data) == 5
end

@testitem "plot_geo_poly" tags = [:general] begin
    using PlotlyBase
    using Meshes
    
    # Test with a single PolyArea
    poly = PolyArea([Point(LatLon{WGS84Latest}(0, 0)), Point(LatLon{WGS84Latest}(0, 1)), Point(LatLon{WGS84Latest}(1, 1)), Point(LatLon{WGS84Latest}(1, 0))])
    @test plot_geo_poly(poly) isa Plot

    # Test with a vector of PolyAreas
    polys = [
    PolyArea([Point(LatLon{WGS84Latest}(0, 0)), Point(LatLon{WGS84Latest}(0, 1)), Point(LatLon{WGS84Latest}(1, 1)), Point(LatLon{WGS84Latest}(1, 0))]),
    PolyArea([Point(LatLon{WGS84Latest}(2, 2)), Point(LatLon{WGS84Latest}(2, 3)), Point(LatLon{WGS84Latest}(3, 3)), Point(LatLon{WGS84Latest}(3, 2))])
    ]
    @test plot_geo_poly(polys) isa Plot

    # Test with a Multi object
    multi = Multi([
    PolyArea([Point(LatLon{WGS84Latest}(0, 0)), Point(LatLon{WGS84Latest}(0, 1)), Point(LatLon{WGS84Latest}(1, 1)), Point(LatLon{WGS84Latest}(1, 0))]),
    PolyArea([Point(LatLon{WGS84Latest}(2, 2)), Point(LatLon{WGS84Latest}(2, 3)), Point(LatLon{WGS84Latest}(3, 3)), Point(LatLon{WGS84Latest}(3, 2))])
    ])
    @test plot_geo_poly(multi) isa Plot

    # Test with a PolyBorder
    poly_border = PolyBorder(poly)
    @test plot_geo_poly(poly_border) isa Plot

    # Test with a MultiBorder
    multi_border = MultiBorder(multi)
    @test plot_geo_poly(multi_border) isa Plot

    # Test with different camera views
    @test plot_geo_poly(poly; camera=:threedim) isa Plot
    @test plot_geo_poly(polys; camera=:threedim) isa Plot

    # Test with custom title and layout options
    @test plot_geo_poly(poly; title="Custom Polygon Plot", kwargs_layout=(width=800, height=600)) isa Plot
end

@testitem "scattergeo" tags = [:plotting] begin
    using PlotlyBase

    offset_region = GeoRegionOffset(; admin = "Spain", delta = 100e3)

    # This is just a coverage test and sanity check to spot errors. It is not really testing functionality
    @test scattergeo(offset_region) == scattergeo(offset_region.domain)

    # Let's do some stupid coverage tests for some plotting functions
    box = BoxBorder(LatLon(0, 0), LatLon(1, 1))
    llr = LatBeltRegion("test", (0, 1))

    @test_nowarn scattergeo(box)
    @test_nowarn scattergeo(llr)
end