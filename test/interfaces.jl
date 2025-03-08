## Base.in()
@testitem "Test Base.in GeoRegion" tags = [:interface] begin
    sample_ita = [LatLon{WGS84Latest}(43.727878°, 12.843441°), LatLon{WGS84Latest}(43.714933°, 10.399326°), LatLon{WGS84Latest}(37.485829°, 14.328285°), LatLon{WGS84Latest}(39.330460°, 8.430780°), LatLon{WGS84Latest}(45.918388°, 10.886654°)]
    sample_eu = [LatLon{WGS84Latest}(52.218550°, 4.420621°), LatLon{WGS84Latest}(41.353144°, 2.167639°), LatLon{WGS84Latest}(42.670341°, 23.322592°)]

    ita = GeoRegion(name="ITA", admin="Italy")
    eu = GeoRegion(; continent="Europe")

    @test all(map(x -> in(x, ita), sample_ita)) # Test LatLon()
    @test all(map(x -> in(Point(x), ita), sample_ita)) # Test Point(LatLon())

    @test all(map(x -> in(x, eu), sample_eu)) # Test LatLon()
    @test all(map(x -> in(Point(x), eu), sample_eu)) # Test Point(LatLon())

    @test in(LatLon(0.7631954460103929rad, 0.22416033273563304rad), ita)
    @test in(LatLon(0.7631954460103929rad, 0.22416033273563304rad), eu)
    @test !in(LatLon(0.7085271959818754rad, -0.2072522112608427rad), eu)
    @test !in(LatLon(52.218550°, 4.420621°), ita)
end

@testitem "Test Base.in PolyRegion" tags = [:interface] begin
    sample_in = [LatLon(14°, 1°), LatLon(26.9°, -4.9°), LatLon(10.1°, 14.9°)]
    sample_out = [LatLon(0°, 0°), LatLon(10°, -5.2°), LatLon(27°, 15.3°)]
    sample_border = [LatLon(10°, -5°), LatLon(10.1°, 10°), LatLon(27°, 15°)] # Due to the Predicates of Meshes the countour is not exact (acceptable)
    poly = PolyRegion("POLY", [LatLon(10°, -5°), LatLon(10°, 15°), LatLon(27°, 15°), LatLon(27°, -5°)])

    @test all(map(x -> in(x, poly), sample_in))
    @test all(map(x -> in(Point(x), poly), sample_in))

    @test all(map(x -> in(x, poly), sample_border))
    @test all(map(x -> in(Point(x), poly), sample_border))

    @test !all(map(x -> in(x, poly), sample_out))
    @test !all(map(x -> in(Point(x), poly), sample_out))

    @test in(LatLon(0.24434609527920614rad, 0.017453292519943295rad), poly)
end

@testitem "Test Base.in LatBeltRegion" tags = [:interface] begin
    belt = LatBeltRegion(; name="test", lim=(-60°, 60°))
    sample_in = [LatLon(14°, 1°), LatLon(26.9°, -65°), LatLon(10.1°, 70°)]
    sample_out = [LatLon(90°, 1°), LatLon(60.1°, 1°), LatLon(-62°, -4.9°), LatLon(-60.1°, 14.9°)]

    @test all(map(x -> in(x, belt), sample_in))
    @test all(map(x -> in(Point(x), belt), sample_in))

    @test !all(map(x -> in(x, belt), sample_out))
    @test !all(map(x -> in(Point(x), belt), sample_out))

    @test in(LatLon(0.24434609527920614, 0.017453292519943295), belt)
end

## borders()
@testitem "Test borders for GeoRegion" tags = [:interface] begin
    using GeoGrids: MULTI_LATLON, MULTI_CART
    reg = GeoRegion(name="ITA", admin="Italy;Spain")

    @test borders(reg) == map(x -> x.latlon, reg.domain)
    @test borders(LatLon, reg) == map(x -> x.latlon, reg.domain)
    @test borders(Cartesian, reg) == map(x -> x.cart, reg.domain)

    # Georegion Offset
    gro = GeoRegionOffset(reg, 50e3)
    @test borders(gro) isa MULTI_LATLON
    @test borders(Cartesian, gro) isa MULTI_CART

    # We test coverage of a point in the multi-border
    @test LatLon(42,12) in gro.domain
end

@testitem "Test borders for PolyRegion" tags = [:interface] begin
    using GeoGrids: MULTI_LATLON, MULTI_CART

    poly = PolyRegion("POLY", [LatLon(10°, -5°), LatLon(10°, 15°), LatLon(27°, 15°), LatLon(27°, -5°)])

    @test borders(poly) == poly.domain.latlon
    @test borders(LatLon, poly) == poly.domain.latlon
    @test borders(Cartesian, poly) == poly.domain.cart

    # This tests both the PolyRegionOffset borders and the one for MultiBorder
    polyoffset = PolyRegionOffset(poly, 50e3)
    @test borders(polyoffset) isa MULTI_LATLON
    @test borders(Cartesian, polyoffset) isa MULTI_CART
end

## centroid()
@testitem "Test centroid for GeoRegion" tags = [:interface] begin
    reg = GeoRegion(; name="Tassellation", admin="Spain;Italy")
    testPoint_cart = Point(Cartesian{WGS84Latest}(2.2774699f0, 41.078136f0))
    testPoint_latlon = Point(LatLon{WGS84Latest}(41.078136f0, 2.2774699f0))
    let
        c = centroid(reg)
        @test abs(get_lat(c) - get_lat(testPoint_cart)) < 1e-2
        @test abs(get_lon(c) - get_lon(testPoint_cart)) < 1e-2
    end
    let 
        c = centroid( reg.domain)
        @test abs(get_lat(c) - get_lat(testPoint_cart)) < 1e-2
        @test abs(get_lon(c) - get_lon(testPoint_cart)) < 1e-2
    end
    let 
        c = centroid(Cartesian, reg)
        @test abs(get_lat(c) - get_lat(testPoint_cart)) < 1e-2
        @test abs(get_lon(c) - get_lon(testPoint_cart)) < 1e-2
    end
    let
        c = centroid(LatLon, reg)
        @test abs(get_lat(c) - get_lat(testPoint_latlon)) < 1e-2
        @test abs(get_lon(c) - get_lon(testPoint_latlon)) < 1e-2
    end
    let
        c = centroid(Cartesian, reg.domain)
        @test abs(get_lat(c) - get_lat(testPoint_cart)) < 1e-2
        @test abs(get_lon(c) - get_lon(testPoint_cart)) < 1e-2
    end
    let
        c = centroid(LatLon, reg.domain)
        @test abs(get_lat(c) - get_lat(testPoint_latlon)) < 1e-2
        @test abs(get_lon(c) - get_lon(testPoint_latlon)) < 1e-2
    end
end

@testitem "Test centroid for PolyRegion" tags = [:interface] begin
    poly = PolyRegion("POLY", [LatLon(10°, -5°), LatLon(10°, 15°), LatLon(27°, 15°), LatLon(27°, -5°)])
    testPoint_cart = Point(Cartesian{WGS84Latest}(5.0, 18.5))
    testPoint_latlon = Point(LatLon{WGS84Latest}(18.5, 5.0))

    @test centroid(poly) == testPoint_cart
    @test centroid(Cartesian, poly) == testPoint_cart
    @test centroid(Cartesian, poly.domain) == testPoint_cart
    @test centroid(Cartesian, poly.domain) == testPoint_cart
    @test centroid(poly.domain.cart) == testPoint_cart

    @test centroid(LatLon, poly) == testPoint_latlon
    @test centroid(LatLon, poly.domain) == testPoint_latlon
end

## Helpers
@testitem "Helper Functions" tags = [:general] begin
    r = GeoRegion(name="ITA", admin="Italy")
    @test extract_countries(r) == r.domain
end