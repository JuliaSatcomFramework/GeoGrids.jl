@testitem "Theory Notebook Functions" tags = [:general] begin
    @test GeoGrids.fibonaccisphere_optimization1(100) isa Vector{SVector{3,Float64}}
    @test GeoGrids.fibonaccisphere_alternative1(100) isa Vector{SVector{3,Float64}}
end

@testitem "Icogrid Functions" tags = [:general] begin
    using Meshes: 🌐, WGS84Latest
    
    @test GeoGrids._icogrid(100; coord=:cart) isa Vector{SVector{3,Float64}}
    @test GeoGrids._icogrid(100; coord=:sphe) isa Vector{SVector{2,Float64}}

    @test icogrid(N=100) isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}

    a = icogrid(sepAng=5)
    b = icogrid(sepAng=5°)
    c = icogrid(sepAng=deg2rad(5) * rad)

    @test a isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test b isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test c isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test length(a) == 1260
    @test length(b) == 1260
    @test length(c) == 1260

    @test abs(get_lat(a[1]) - 87.7171°) < 1e-4
    @test abs(get_lon(a[1]) - 0.0°) < 1e-4
    @test abs(get_lat(a[end]) - -87.7171°) < 1e-4
    @test abs(get_lon(a[end]) - 37.7251°) < 1e-4
    @test abs(get_lat(b[1]) - 87.7171°) < 1e-4
    @test abs(get_lon(b[1]) - 0.0°) < 1e-4
    @test abs(get_lat(b[end]) - -87.7171°) < 1e-4
    @test abs(get_lon(b[end]) - 37.7251°) < 1e-4
    @test abs(get_lat(c[1]) - 87.7171°) < 1e-4
    @test abs(get_lon(c[1]) - 0.0°) < 1e-4
    @test abs(get_lat(c[end]) - -87.7171°) < 1e-4
    @test abs(get_lon(c[end]) - 37.7251°) < 1e-4

    @test_logs (:warn, "Input sepAng is negative, it will be converted to positive...") icogrid(sepAng=-5°)
    @test_throws "The sepAng provided as numbers must be expressed in radians and satisfy -360° ≤ x ≤ 360°. 
Consider using `°` (or `rad`) from `Unitful` if you want to pass numbers in degrees (or rad), by doing `x * °` (or `x * rad`)." icogrid(sepAng=361°)
end

@testitem "Mesh Grid Functions" tags = [:general] begin
    using Meshes: 🌐, WGS84Latest
    
    @test rectgrid(5) isa Matrix{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test rectgrid(5°) isa Matrix{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test rectgrid(deg2rad(5) * rad) isa Matrix{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test rectgrid(5; yRes=3) isa Matrix{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test rectgrid(5°; yRes=3°) isa Matrix{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test rectgrid(deg2rad(5) * rad; yRes=deg2rad(3) * rad) isa Matrix{<:Point{🌐,<:LatLon{WGS84Latest}}}

    @test_throws "Resolution of x is too large, it must be smaller than 180°..." rectgrid(181°)
    @test_throws "Resolution of y is too large, it must be smaller than 180°..." rectgrid(5°; yRes=181°)

    @test_logs (:warn, "Input xRes is negative, it will be converted to positive...") rectgrid(-5°; yRes=3°)
    @test_logs (:warn, "Input yRes is negative, it will be converted to positive...") rectgrid(5°; yRes=-3°)

    a = rectgrid(5)
    b = rectgrid(5°)
    c = rectgrid(deg2rad(5) * rad)
    @test length(a) == 2664
    @test length(b) == 2664
    @test length(c) == 2664

    @test get_lat(a[1, 1]) ≈ -90°
    @test get_lon(a[1, 1]) ≈ -175°
    @test get_lat(a[end, end]) ≈ 90°
    @test get_lon(a[end, end]) ≈ 180°
    @test abs(get_lon(a[1, 2]) - get_lon(a[1, 1])) ≈ 5°
    @test abs(get_lat(a[1, 2]) - get_lat(a[1, 1])) ≈ 0°
    @test abs(get_lon(a[2, 1]) - get_lon(a[1, 1])) ≈ 0°
    @test abs(get_lat(a[2, 1]) - get_lat(a[1, 1])) ≈ 5°
    @test get_lat(b[1, 1]) ≈ -90°
    @test get_lon(b[1, 1]) ≈ -175°
    @test get_lat(b[end, end]) ≈ 90°
    @test get_lon(b[end, end]) ≈ 180°
    @test abs(get_lon(b[1, 2]) - get_lon(b[1, 1])) ≈ 5°
    @test abs(get_lat(b[1, 2]) - get_lat(b[1, 1])) ≈ 0°
    @test abs(get_lon(b[2, 1]) - get_lon(b[1, 1])) ≈ 0°
    @test abs(get_lat(b[2, 1]) - get_lat(b[1, 1])) ≈ 5°
    @test get_lat(c[1, 1]) ≈ -90°
    @test get_lon(c[1, 1]) ≈ -175°
    @test get_lat(c[end, end]) ≈ 90°
    @test get_lon(c[end, end]) ≈ 180°
    @test abs(get_lon(c[1, 2]) - get_lon(c[1, 1])) ≈ 5°
    @test abs(get_lat(c[1, 2]) - get_lat(c[1, 1])) ≈ 0°
    @test abs(get_lon(c[2, 1]) - get_lon(c[1, 1])) ≈ 0°
    @test abs(get_lat(c[2, 1]) - get_lat(c[1, 1])) ≈ 5°
end

@testitem "Vec Grid Functions" tags = [:general] begin
    using Meshes: 🌐, WGS84Latest

    @test vecgrid(5) isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test vecgrid(5°) isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}
    @test vecgrid(deg2rad(5) * rad) isa Vector{<:Point{🌐,<:LatLon{WGS84Latest}}}

    @test_throws "Resolution of grid is too large, it must be smaller than 90°..." vecgrid(91)
    @test_logs (:warn, "Input gridRes is negative, it will be converted to positive...") vecgrid(-deg2rad(5) * rad)

    a = vecgrid(5)
    @test length(a) == 19
    @test abs(get_lat(a[1]) - get_lat(a[2])) ≈ 5°
end