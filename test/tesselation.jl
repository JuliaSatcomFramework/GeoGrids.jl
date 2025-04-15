@testsnippet setup_tessellation begin
    using GeoGrids
    using GeoGrids: constants, to_raw_lonlat

    function haversine(p1, p2; refRadius = constants.Re_mean)
        # Convert latitude and longitude to radians
        lon1, lat1 = to_raw_lonlat(p1) .|> deg2rad
        lon2, lat2 = to_raw_lonlat(p2) .|> deg2rad

        # Calculate differences in longitude and latitude
        Δlon = lon2 - lon1
        Δlat = lat2 - lat1

        # Calculate the great-circle distance
        a = sin(Δlat / 2)^2 + cos(lat1) * cos(lat2) * sin(Δlon / 2)^2
        c = 2 * atan(sqrt(a), sqrt(1 - a))

        # Return the distance in kilometers
        return refRadius * c
    end
end

## HEX
@testitem "GeoRegion HEX Layout Tesselation and :hex Pattern (no EO)" tags = [:tesselation] begin
    # NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test with predefined points like in the cases with :circ.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.
    samples = [LatLon(43.2693, -8.83691), LatLon(43.6057, -8.11563), LatLon(36.808, -5.06421), LatLon(36.792, -2.3703)] # LatLon{WGS84Latest} coordinates
    corresponding_idxs = [1, 2, 121, 122]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers = generate_tesselation(reg, 40000, HEX())

    @test length(centers) == 122
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1
    end
end

@testitem "GeoRegion HEX Layout Tesselation and :circ Pattern (no EO)" tags = [:tesselation] begin
    samples = [LatLon(43.2693, -8.83691), LatLon(43.6057, -8.11563), LatLon(36.808, -5.06421), LatLon(36.792, -2.3703)] # LatLon{WGS84Latest} coordinates
    corresponding_idxs = [1, 2, 121, 122]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers = generate_tesselation(reg, 40000, HEX(; pattern=:circ))

    @test length(centers) == 122
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1e-4
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1e-4
    end
end

@testitem "GeoRegion HEX Layout Tesselation and :hex Pattern (with EO)" tags = [:tesselation] begin
    # //NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test like the previous cases using a circular tesselation.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.
    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers, ngon = generate_tesselation(reg, 40000, HEX(), EO())

    @test length(centers) == 122
    @test length(ngon) == 122

    samplePoints = [LatLon(43.2693, -8.83691), LatLon(43.6057, -8.11563), LatLon(36.808, -5.06421), LatLon(36.792, -2.3703)] # LatLon{WGS84Latest} coordinates
    corresponding_idxs_points = [1, 2, 121, 122]
    sampleNgons = [[LatLon(43.558,-9.1653), LatLon(43.2648,-9.27284), LatLon(42.937,-9.11636), LatLon(42.9797,-8.51291), LatLon(43.2764,-8.40112), LatLon(43.6025,-8.55324), LatLon(43.558,-9.1653)],
    [LatLon(36.4876,-2.64153), LatLon(36.4735,-2.11589), LatLon(36.7875,-1.96157), LatLon(37.0961,-2.09562), LatLon(37.1107,-2.62807), LatLon(36.7972,-2.77902), LatLon(36.4876,-2.64153)]] 
    corresponding_idxs_ngon = [1, 122]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
	centers, ngon = generate_tesselation(reg, 40000, HEX(), EO())

    @test length(centers) == 122
    for i in eachindex(samplePoints) 
        @test abs(get_lat(centers[corresponding_idxs_points[i]]) - samplePoints[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs_points[i]]) - samplePoints[i].lon) < 1    
    end

    @test length(ngon) == 122
    for i in eachindex(sampleNgons)
        for v in eachindex(sampleNgons[i])
            # @test get_lat(ngon[corresponding_idxs_ngon[i]][v]) == sampleNgons[i][v].lat
            @test abs(get_lat(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lat) < 1
            @test abs(get_lon(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lon) < 1
        end
    end
end

@testitem "GeoRegion HEX Layout Tesselation and :circ Pattern (with EO)" tags = [:tesselation] setup=[setup_tessellation] begin
    samplePoints = [LatLon(43.2693, -8.83691), LatLon(43.6057, -8.11563), LatLon(36.808, -5.06421), LatLon(36.792, -2.3703)] # LatLon{WGS84Latest} coordinates
    corresponding_idxs = [1, 2, 121, 122]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers, ngon = generate_tesselation(reg, 40000, HEX(; pattern=:circ), EO())

    @test length(centers) == 122
    for i in eachindex(samplePoints)
        @test haversine(centers[corresponding_idxs[i]], samplePoints[i]) < 10 # 10m is less than 1e-4 radians
    end

    @test length(ngon) == 122
    for (i, idx) in enumerate(corresponding_idxs)
        p1 = samplePoints[i]
        for p in ngon[idx]
            # We test that the great circle distance to the resulting ngon points is roughly the target radius
            @test haversine(p, p1) - 40000 < 10 # 10m is less than 1e-4 radians
        end
    end
end

@testitem "GeoRegion HEX Layout Tesselation and :hex Pattern :flat (no EO)" tags = [:tesselation] begin
    # NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test with predefined points like in the cases with :circ.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.
    samples = [LatLon(42.8428, -9.06506), LatLon(43.4251, -7.82139), LatLon(37.5057, -1.64231), LatLon(36.459, -4.77773),]
    corresponding_idxs = [1, 2, 117, 118]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers = generate_tesselation(reg, 40000, HEX(; direction=:flat))

    @test length(centers) == 118
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1
    end
end

@testitem "PolyRegion HEX Layout Tesselation and :hex Pattern (no EO)" tags = [:tesselation] begin
    # NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test with predefined points like in the cases with :circ.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.    
    samples = [LatLon(79.9691, 5.0), LatLon(79.6441, 7.98881), LatLon(60.0123, 7.14867), LatLon(60.3001, 8.25311)]
    corresponding_idxs = [1, 2, 318, 319]

    reg = PolyRegion(; domain=[LatLon(60, -5), LatLon(80, 0), LatLon(80, 10), LatLon(60, 15)])
    centers = generate_tesselation(reg, 40000, HEX())

    @test length(centers) == 319
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1
    end
end

## ICO
@testitem "GlobalRegion ICO Layout Tesselation and :circ Pattern (no EO)" tags = [:tesselation] begin
    samples = [LatLon(83.991, 0.0), LatLon(79.5826, -137.508), LatLon(-79.5826, 88.6025), LatLon(-83.991, -48.9053)]
    corresponding_idxs = [1, 2, 181, 182]

    reg = GlobalRegion()
    centers = generate_tesselation(reg, 1000000, ICO())

    @test length(centers) == 182
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1e-4
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1e-4
    end
end

@testitem "GlobalRegion ICO Layout Tesselation and :hex Pattern (with EO)" tags = [:tesselation] begin
     # NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test with predefined points like in the cases with :circ.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.
    samplePoints = [LatLon(83.991, 0.0), LatLon(79.5826, -137.508), LatLon(-79.5826, 88.6025), LatLon(-83.991, -48.9053)]
    corresponding_idxs_points = [1, 2, 181, 182]
    sampleNgons = [[LatLon(66.8411, -2.78182), LatLon(75.8026, 15.7204), LatLon(82.6627, 18.6819), LatLon(83.991, 0.0), LatLon(83.4003, -18.4249), LatLon(66.8411, -2.78182)],
        [LatLon(-82.3641, -67.4582), LatLon(-83.991, -48.9053), LatLon(-83.4003, -30.4803), LatLon(-66.8411, -46.1235), LatLon(-75.8026, -64.6257), LatLon(-82.3641, -67.4582)]]
    corresponding_idxs_ngon = [1, 182]

    reg = GlobalRegion()
    centers, ngon = generate_tesselation(reg, 1000000, ICO(; pattern=:hex), EO())

    @test length(centers) == 182
    for i in eachindex(samplePoints)
        @test abs(get_lat(centers[corresponding_idxs_points[i]]) - samplePoints[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs_points[i]]) - samplePoints[i].lon) < 1
    end

    @test length(ngon) == 182
    for i in eachindex(sampleNgons)
        for v in eachindex(sampleNgons[i])
            @test abs(get_lat(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lat) < 1
            @test abs(get_lon(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lon) < 1
        end
    end
end

@testitem "LatBeltRegion ICO Layout Tesselation and :circ Pattern (no EO)" tags = [:tesselation] begin
    samples = [LatLon(9.89685, 41.5062), LatLon(9.78996, -96.0016), LatLon(-9.78996, 25.0621), LatLon(-9.89685, -112.446)]
    corresponding_idxs = [1, 2, 187, 188]

    reg = LatBeltRegion(lim=(-10, 10))
    centers = generate_tesselation(reg, 400000, ICO())

    @test length(centers) == 188
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1e-4
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1e-4
    end

    @test_throws "10000000 points is not sufficient for the requested separation angle" generate_tesselation(GlobalRegion(), 4e3, ICO())
    
    centers = generate_tesselation(GlobalRegion(), 4e3, ICO(); maxPrec=10^8)
    @test length(centers) == 10779808
end

@testitem "LatBeltRegion ICO Layout Tesselation and :hex Pattern (with EO)" tags = [:tesselation] begin
    # //NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test like the previous cases using a circular tesselation.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.
    samplePoints = [LatLon(9.89685,41.5062), LatLon(9.78996,-96.0016), LatLon(-9.78996,25.0621), LatLon(-9.89685,-112.446)]
    corresponding_idxs_points = [1, 2, 187, 188]

    sampleNgons = [[LatLon(11.8183,45.104),LatLon(13.921,40.8336),LatLon(9.54199,37.4515),LatLon(8.05959,37.8844),LatLon(5.89612,42.1999),LatLon(10.2875,45.5542),LatLon(11.8183,45.104)],
    [LatLon(-11.8183,-116.044),LatLon(-13.921,-111.773),LatLon(-9.54199,-108.391),LatLon(-8.05959,-108.824),LatLon(-5.89611,-113.139),LatLon(-10.2875,-116.494),LatLon(-11.8183,-116.044)]]
    corresponding_idxs_ngon = [1, 188]

    reg = LatBeltRegion(lim=(-10,10))
	centers,ngon = generate_tesselation(reg, 400000, ICO(;pattern=:hex), EO())

    @test length(centers) == 188
    for i in eachindex(samplePoints) 
        @test abs(get_lat(centers[corresponding_idxs_points[i]]) - samplePoints[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs_points[i]]) - samplePoints[i].lon) < 1  
    end

    @test length(ngon) == 188
    for i in eachindex(sampleNgons)
        for v in eachindex(sampleNgons[i])
            @test abs(get_lat(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lat) < 1
            @test abs(get_lon(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lon) < 1
        end
    end
end

@testitem "PolyRegion ICO Layout Tesselation and :circ Pattern (no EO)" tags = [:tesselation] begin
    samples = [LatLon(79.9456, 4.99995), LatLon(79.6831, 9.73597), LatLon(60.0598, -3.47868), LatLon(60.0243, 8.92039)]
    corresponding_idxs = [1, 2, 218, 219]

    reg = PolyRegion(; domain=[LatLon(60, -5), LatLon(80, 0), LatLon(80, 10), LatLon(60, 15)])
    centers = generate_tesselation(reg, 40000, ICO(; correction=1.7))

    @test length(centers) == 219
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1e-4
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1e-4
    end
end

@testitem "PolyRegion ICO Layout Tesselation and :hex Pattern (with EO)" tags = [:tesselation] begin
  # //NOTE: 
    #The order of points entering `DelaunayTriangulation.triangulate`
    # is randomised as default so we end up with a slighlty different
    # tessellation each time. This is not an issue for us but we cannot perform
    # an hardcoded test like the previous cases using a circular tesselation.
    # For the moment being we decreased the tolerance to 1 (degrees of latitude 
    # and longitude) to allow for the randomness in the tessellation.
    samplePoints = [LatLon(79.9456,4.99995), LatLon(79.6831,9.73597), LatLon(60.0598,-3.47868), LatLon(60.0243,8.92039)]
    corresponding_idxs_points = [1, 2, 218, 219]
    sampleNgons = [[LatLon(79.2839,5.7835),LatLon(80.0631,6.07365),LatLon(81.0731,4.99257),LatLon(80.8394,4.30809),LatLon(79.929,3.94491),LatLon(78.9239,4.90667),LatLon(79.2839,5.7835)],
        [LatLon(60.0209,9.55066),LatLon(60.4379,9.40303),LatLon(60.6574,8.86491),LatLon(60.0377,8.28825),LatLon(59.6336,8.42962),LatLon(59.4002,8.98413),LatLon(60.0209,9.55066)]]
    corresponding_idxs_ngon = [1, 219]
    
    reg = PolyRegion(;domain=[LatLon(60,-5), LatLon(80,0), LatLon(80,10), LatLon(60,15)])
	centers, ngon = generate_tesselation(reg, 40000, ICO(;correction=1.7, pattern=:hex), EO())

    @test length(centers) == 219
    for i in eachindex(samplePoints) 
        @test abs(get_lat(centers[corresponding_idxs_points[i]]) - samplePoints[i].lat) < 1
        @test abs(get_lon(centers[corresponding_idxs_points[i]]) - samplePoints[i].lon) < 1  
    end

    @test length(ngon) == 219
    for i in eachindex(sampleNgons)
        for v in eachindex(sampleNgons[i])
            @test abs(get_lat(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lat) < 1
            @test abs(get_lon(ngon[corresponding_idxs_ngon[i]][v]) - sampleNgons[i][v].lon) < 1
        end
    end
end

@testitem "GeoRegion ICO Layout Tesselation and :circ Pattern (no EO)" tags = [:tesselation] begin
    samples = [LatLon(43.5772, -7.50902), LatLon(43.4468, -5.70002), LatLon(36.6657, -3.58858), LatLon(36.285, -5.82463)]
    corresponding_idxs = [1, 2, 103, 104]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers = generate_tesselation(reg, 40000, ICO())

    @test length(centers) == 104
    for i in eachindex(samples)
        @test abs(get_lat(centers[corresponding_idxs[i]]) - samples[i].lat) < 1e-4
        @test abs(get_lon(centers[corresponding_idxs[i]]) - samples[i].lon) < 1e-4
    end
end

@testitem "GeoRegion ICO Layout Tesselation and :circ Pattern (with EO)" tags = [:tesselation] setup=[setup_tessellation] begin
    samplePoints = [LatLon(43.5772, -7.50902), LatLon(43.4468, -5.70002), LatLon(36.6657, -3.58858), LatLon(36.285, -5.82463)]
    corresponding_idxs = [1, 2, 103, 104]

    reg = GeoRegion(; name="Tassellation", admin="Spain")
    centers, ngon = generate_tesselation(reg, 40000, ICO(), EO())

    @test length(centers) == 104
    for i in eachindex(samplePoints)
        @test haversine(centers[corresponding_idxs[i]], samplePoints[i]) < 10 # 10m is less than 1e-4 radians
    end

    @test length(ngon) == 104
    for (i, idx) in enumerate(corresponding_idxs)
        p1 = samplePoints[i]
        for p in ngon[idx]
            # We test that the great circle distance to the resulting ngon points is roughly the target radius
            @test haversine(p, p1) - 40000 < 10 # 10m is less than 1e-4 radians
        end
    end
end

@testitem "Misc coverage" begin
    using GeoGrids: gen_circle_pattern, get_lat, _wrap_latlon

    # We test wrapping around poles
    pts = gen_circle_pattern(LatLon(90, 0), 1000e3) |> first
    lat = get_lat(first(pts))
    @test all(p -> get_lat(p) ≈ lat, pts)

    pts = gen_circle_pattern(LatLon(-90, 0), 1000e3) |> first
    lat = get_lat(first(pts))
    @test all(p -> get_lat(p) ≈ lat, pts)

    @test _wrap_latlon(100, 0) == (80, 180)
    @test _wrap_latlon(-100, 0) == (-80, 180)
end

@testitem "Test missing Functions" tags = [:tesselation] begin
    reg = GlobalRegion()
    @test_throws "H3 tassellation is not yet implemented in this version..." generate_tesselation(reg, 10.0, H3())

    @test_throws ArgumentError generate_tesselation(GlobalRegion(), 10.0, HEX())
    @test_throws ArgumentError generate_tesselation(GlobalRegion(), 10.0, HEX(), EO())
end