@testitem "GeoRegion Test" tags = [:filtering] begin
    sample_ita = [LatLon(43.727878°, 12.843441°), LatLon(43.714933°, 10.399326°), LatLon(37.485829°, 14.328285°), LatLon(39.330460°, 8.430780°), LatLon(45.918388°, 10.886654°)]
    sample_eu = [LatLon(52.218550°, 4.420621°), LatLon(41.353144°, 2.167639°), LatLon(42.670341°, 23.322592°)]

    ita = GeoRegion(name="ITA", admin="Italy")
    eu = GeoRegion(; continent="Europe")

    @test ita isa GeoRegion
    @test eu isa GeoRegion
    @test ita isa AbstractRegion
    @test eu isa AbstractRegion
    
    @test_throws "Input at least one argument between continent, subregion and admin..." GeoRegion()

    @test filter_points(sample_ita, ita) == sample_ita
    @test filter_points(sample_eu, eu) == sample_eu
    @test filter_points([LatLon(52.218550°, 4.420621°), LatLon(43.727878°, 12.843441°), LatLon(41.353144°, 2.167639°), LatLon(43.714933°, 10.399326°)], ita) == [LatLon(43.727878°, 12.843441°), LatLon(43.714933°, 10.399326°)]
    
    @test filter_points(map(x -> Point(x), sample_ita), ita) == map(x -> Point(x), sample_ita) # Additional test for type Point(LatLon())
end

@testitem "PolyRegion Test" tags = [:filtering] begin
    sample_in = [LatLon(14°, 1°), LatLon(26.9°, -4.9°), LatLon(10.1°, 14.9°)]
    sample_out = [LatLon(0°, 0°), LatLon(10°, -5.2°), LatLon(27°, 15.3°)]

    poly = PolyRegion("POLY", [LatLon(10°, -5°), LatLon(10°, 15°), LatLon(27°, 15°), LatLon(27°, -5°)])
    vertex = [LatLon(10°, -5°), LatLon(10°, 15°), LatLon(27°, 15°), LatLon(27°, -5°)]

    @test poly isa PolyRegion
    @test PolyRegion(;domain=vertex) isa PolyRegion
    @test PolyRegion(;name="Test",domain=vertex) isa PolyRegion
    @test PolyRegion("Test", vertex) isa PolyRegion
    @test_throws "UndefKeywordError: keyword argument `domain` not assigned" PolyRegion()
    
    @test filter_points(vcat(sample_in, sample_out), poly) == sample_in
    @test filter_points(map(x -> Point(x), vcat(sample_in, sample_out)), poly) == map(x -> Point(x), sample_in) # Additional test for type Point(LatLon())
end

@testitem "LatBeltRegion Test" tags = [:filtering] begin
    belt = LatBeltRegion(; name="test", lim=(-60°, 60°))
    sample_in = [LatLon(14°, 1°), LatLon(26.9°, -65°), LatLon(10.1°, 70°)]
    sample_out = [LatLon(90°, 1°), LatLon(60.1°, 1°), LatLon(-62°, -4.9°), LatLon(-60.1°, 14.9°)]

    @test belt isa LatBeltRegion
    @test LatBeltRegion(; name="test", lim=(0°,90°)) isa LatBeltRegion
    @test LatBeltRegion(; lim=(0°,90°)) isa LatBeltRegion
    @test LatBeltRegion("test", (0°,90°)) isa LatBeltRegion
    
    a = LatBeltRegion("test", (0°,90°))
    b = LatBeltRegion("test", (0,90))
    c = LatBeltRegion("test", (0rad, (π/2)rad))
    @test a isa LatBeltRegion
    @test b isa LatBeltRegion
    @test c isa LatBeltRegion
    @test a.lim == b.lim == c.lim == (0°,90°)

    @test_throws "UndefKeywordError: keyword argument `lim` not assigned" LatBeltRegion()
    @test_throws "LAT provided as numbers must be expressed in degrees and satisfy -90 ≤ x ≤ 90. 
Consider using `°` (or `rad`) from `Unitful` if you want to pass numbers in degrees (or rad), by doing `x * °` (or `x * rad`)." LatBeltRegion("test", (0°,91°))
    @test_throws "LAT provided as numbers must be expressed in degrees and satisfy -90 ≤ x ≤ 90. 
Consider using `°` (or `rad`) from `Unitful` if you want to pass numbers in degrees (or rad), by doing `x * °` (or `x * rad`)." LatBeltRegion("test", (-91°,91°))
    @test_throws "LAT provided as numbers must be expressed in degrees and satisfy -90 ≤ x ≤ 90. 
Consider using `°` (or `rad`) from `Unitful` if you want to pass numbers in degrees (or rad), by doing `x * °` (or `x * rad`)." LatBeltRegion("test", (-91°,0°))
    @test_throws "The first LAT limit must be lower than the second one..." LatBeltRegion(; lim=((π/2)rad, 0rad))
    @test_throws "The first LAT limit must be different than the second one..." LatBeltRegion(; lim=(90, 90))
    
    points = [LatLon(14°, 1°), LatLon(90°, 1°), LatLon(60.1°, 1°), LatLon(26.9°, -65°), LatLon(-62°, -4.9°), LatLon(-60.1°, 14.9°), LatLon(10.1°, 70°)]
    @test filter_points(points, belt) == sample_in

    @test filter_points(map(Point, points), belt) == map(Point, sample_in) # Additional test for type Point(LatLon())
end

@testitem "Group By Test" tags = [:filtering] begin
    using Meshes: 🌐, WGS84Latest

    ita = GeoRegion(; name="ITA", admin="Italy")
    eu = GeoRegion(; name="EU", continent="Europe")
    poly = PolyRegion("POLY", [LatLon(10°, -5°), LatLon(10°, 15°), LatLon(27°, 15°), LatLon(27°, -5°)])
    belt = LatBeltRegion(; name="BELT", lim=(0°, 5°))
    
    sample_in_ita = [LatLon(43.727878°, 12.843441°), LatLon(43.714933°, 10.399326°), LatLon(37.485829°, 14.328285°), LatLon(39.330460°, 8.430780°), LatLon(45.918388°, 10.886654°)]
    sample_in_poly = [LatLon(14°, 1°), LatLon(26.9°, -4.9°), LatLon(10.1°, 14.9°)]
    sample_out_poly = [LatLon(0°, 0°), LatLon(10°, -5.2°), LatLon(27°, 15.3°)]
    sample_in_belt = [LatLon(1°, 1°), LatLon(2.5°, -65°), LatLon(4.9°, 70°)]
    sample_out_belt = [LatLon(90°, 1°), LatLon(60.1°, 1°), LatLon(-62°, -4.9°), LatLon(-60.1°, 14.9°)]

    big_vec = [sample_in_ita..., sample_in_poly..., sample_out_poly..., sample_in_belt..., sample_out_belt...]

    # Test with LatLon
    # Unique test
    groups_unique = group_by_domain(big_vec, [ita, eu, poly, belt])
    @test groups_unique["ITA"] == sample_in_ita
    @test groups_unique["POLY"] == sample_in_poly
    @test groups_unique["BELT"] == sample_in_belt
    @test isempty(groups_unique["EU"])
    @test groups_unique["ITA"] isa Vector{<:LatLon}
    @test groups_unique["POLY"] isa Vector{<:LatLon}
    @test groups_unique["BELT"] isa Vector{<:LatLon}
    # Repeated elements test
    groups = group_by_domain(big_vec, [ita, eu, poly, belt]; flagUnique=false)
    @test groups["ITA"] == sample_in_ita
    @test groups["EU"] == sample_in_ita
    @test groups["POLY"] == sample_in_poly
    @test groups["BELT"] == sample_in_belt
    @test groups["ITA"] isa Vector{<:LatLon}
    @test groups["POLY"] isa Vector{<:LatLon}
    @test groups["BELT"] isa Vector{<:LatLon}

    # Test with Point(LatLon)
    # Unique test
    groups_unique = group_by_domain(map(x -> Point(x), big_vec), [ita, eu, poly, belt])
    @test groups_unique["ITA"] == map(x -> Point(x), sample_in_ita)
    @test groups_unique["POLY"] == map(x -> Point(x), sample_in_poly)
    @test groups_unique["BELT"] == map(x -> Point(x), sample_in_belt)
    @test isempty(groups_unique["EU"])
    @test groups_unique["ITA"] isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test groups_unique["POLY"] isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test groups_unique["BELT"] isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}
    # Repeated elements test
    groups = group_by_domain(map(x -> Point(x), big_vec), [ita, eu, poly, belt]; flagUnique=false)
    @test groups["ITA"] == map(x -> Point(x), sample_in_ita)
    @test groups["EU"] == map(x -> Point(x), sample_in_ita)
    @test groups["POLY"] == map(x -> Point(x), sample_in_poly)
    @test groups["BELT"] == map(x -> Point(x), sample_in_belt)
    @test groups["ITA"] isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test groups["POLY"] isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test groups["BELT"] isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}
end